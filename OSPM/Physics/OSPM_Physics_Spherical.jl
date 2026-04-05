# ############################################################
# ############################################################
# HOT PATH — A-matrix builder & batch evaluator.
# These are the only two functions called in production.
# ############################################################
# ############################################################
# DO NOT MODIFY without consulting Nate.


module OSPMPhysicsSpherical
@info "OSPMPhysicsSpherical loaded from" @__FILE__
using LinearAlgebra, StaticArrays, Statistics, Random, Base.Threads, Optim
export build_R_halo_physical, rho_interp, halo_from_theta, tables_spherical, make_potential_force_funcs, integrate_orbit_rk4, orbit_to_sigma2_profile, build_A_matrix_julia, ospm_runcheck, build_A_matrix_stellar, build_A_matrix_hybrid, mass_enclosed_two_radii, reset_orbit_cache, evaluate_batch_theta, NTHREADS, stellar_log_likelihood_jl, solve_weights_stellar_jl, ShellRefinementState, effective_rank, shell_support, shell_score, is_saturated, update_shell_state!

# ============================================================
# Support code: constants, types, caches, helpers, halo
# physics, orbit integration, likelihood, weight solver,
# and dead/back-compat code.
# ============================================================
include("OSPM_Physics_Support.jl")

# Main A-matrix builder: maps orbital weights → observables (vlos likelihoods + occupancy). Parallel, numerically hardened.
function build_A_matrix_hybrid(Norbit::Int, R_star_m::Vector{Float64}, has_vlos::AbstractVector{Bool}, v_star_mps::Vector{Float64}, verr_star_mps::Vector{Float64}, sini::Float64, rho_s::Float64, r_s::Float64, MBH::Float64, halo_type::String; nsteps::Int=DEFAULT_NSTEPS, Lfrac::NTuple{5,Float64}=DEFAULT_LFRAC, dt_frac_orbit::Float64=DEFAULT_DT_FRAC, dR_frac::Float64=DEFAULT_DR_FRAC, Nbins_occ::Int=DEFAULT_NBINS_OCC, return_occ::Bool=true, max_attempts_factor::Int=DEFAULT_MAX_ATTEMPTS, diag::Bool=false, threaded::Bool=true, fill_pct::Float64=0.95)

    Nstar = length(R_star_m)
    R_sorted = sort(R_star_m)
    ΔR_med = median(diff(R_sorted))
    @assert length(has_vlos)       == Nstar
    @assert length(v_star_mps)     == Nstar
    @assert length(verr_star_mps)  == Nstar
    Nstar == 0 && return zeros(Float64, 0, Norbit)
    vlos_idx = Int[]
    sizehint!(vlos_idx, Nstar)
    @inbounds for i in 1:Nstar
        has_vlos[i] && push!(vlos_idx, i)
    end

    Nvlos = length(vlos_idx)
    ctx   = get_halo_context(rho_s, r_s, MBH, halo_type)
    sini  = clamp01(f64(sini))
    Rmin  = minimum(R_star_m)
    Rmax  = maximum(R_star_m)
    Nout  = Nvlos + (return_occ ? Nbins_occ : 0)

    if !(isfinite(Rmin) && isfinite(Rmax) && Rmax > Rmin)
        return zeros(Float64, Nout, Norbit)
    end

    shells     = sort(copy(R_star_m))
    Nshells    = length(shells)
    occ_edges  = collect(range(Rmin, Rmax; length=Nbins_occ + 1))
    theta0     = f64(pi / 2)
    orbit_ctx  = (frc=ctx.frc, R_pos=ctx.R, halo=ctx.halo)
    A_vlos = zeros(Float64, Nvlos,    Norbit)
    A_occ  = zeros(Float64, Nbins_occ, Norbit)
    success_flags = falses(Norbit)
    min_r_reached = fill(Inf, Norbit)
    rapo_list     = fill(NaN, Norbit)
    star_hit_flags = falses(Nstar)
    nworkers = threaded ? Threads.nthreads() : 1
    rngs = [MersenneTwister(0x5eed1234 + UInt(t)) for t in 1:nworkers]
    # --- Cost-sorted dispatch: expensive orbits (low Lfrac, small rapo) go first ---
    _orbit_cost = Vector{Float64}(undef, Norbit)
    @inbounds for c in 1:Norbit
        lf   = Lfrac[1 + ((c - 1) % length(Lfrac))]
        rapo = shells[mod1(c, Nshells)]
        _orbit_cost[c] = lf * rapo      # low = expensive (radial + deep)
    end
    cost_order = sortperm(_orbit_cost)
    next_orbit = Threads.Atomic{Int}(1)
    filled_atomic = Threads.Atomic{Int}(0)
    fill_target = round(Int, max(fill_pct, 0.99) * Norbit)

    function _worker!(rng)
        col_occ  = zeros(Float64, Nbins_occ)
        col_vlos = zeros(Float64, Nvlos)
        s_arr    = Vector{Float64}(undef, nsteps)
        vlos_buf = Vector{Float64}(undef, nsteps)

        while true
            filled_atomic[] >= fill_target && break
            c_seq = Threads.atomic_add!(next_orbit, 1)
            c_seq > Norbit && break
            c_claim = cost_order[c_seq]
            idx_local = mod1(c_claim, Nshells)
            rapo = f64(shells[idx_local])
            rapo_list[c_claim] = rapo
            !(isfinite(rapo) && rapo > 0.0) && continue
            lf      = Lfrac[1 + ((c_claim - 1) % length(Lfrac))]
            r0_frac = 0.95 + 0.04 * rand(rng)
            ic, Lz0, E0, vc, st = launch_orbit_apocenter(rapo=rapo, theta0=theta0, Lz_frac=f64(lf), pot=ctx.pot, frc=ctx.frc, r0_frac=r0_frac, dt_frac=dt_frac_orbit)
            st != :ok && continue
            r, vr, theta, vtheta = integrate_orbit_rk4(ic=ic, xLz=Lz0, orbit_ctx=orbit_ctx, nsteps=nsteps)
            isempty(r) && continue
            success_flags[c_claim] = true
            rapo_list[c_claim]     = rapo
            min_r_reached[c_claim] = minimum(r)
            Nhits    = length(r)
            dt_orb   = f64(ic[3])
            resize!(s_arr,Nhits)
            resize!(vlos_buf,Nhits)
            phi = 0.0
            @inbounds for i in 1:Nhits
                si          = _ssin(f64(theta[i]))
                ri          = f64(r[i])
                s_arr[i]    = ri * si
                vphi_i      = f64(Lz0) / (ri * si)
                cp, sp      = cos(phi), sin(phi)
                vlos_buf[i] = sini * (f64(vr[i]) * cp - vphi_i * sp)
                phi        += f64(Lz0) / (ri * ri * si * si) * dt_orb
            end

            fill!(col_occ,0.0)
            fill!(col_vlos,0.0)
            if return_occ && Nhits > 0
                @inbounds for j in 1:Nbins_occ
                    lo  = occ_edges[j]
                    hi  = occ_edges[j + 1]
                    cnt = 0
                    @inbounds for x in s_arr
                        cnt += (x >= lo && x < hi) ? 1 : 0
                    end
                    col_occ[j] = cnt / Nhits
                end
            end

            if Nvlos > 0
                pidx        = sortperm(s_arr)
                s_sorted    = s_arr[pidx]
                vlos_sorted = vlos_buf[pidx]
                inv_sqrt2pi = inv(sqrt(2 * pi))
                @inbounds for row in 1:Nvlos
                    istar = vlos_idx[row]
                    Ri    = f64(R_star_m[istar])
                    vi    = f64(v_star_mps[istar])
                    si    = f64(verr_star_mps[istar])
                    s2    = si * si

                    if !(isfinite(s2) && s2 > 0.0)
                        col_vlos[row] = 0.0
                        continue
                    end

                    dR_local = local_spacing(R_sorted, Ri)
                    dR = max( dR_frac * Ri, 0.5 * dR_local, 0.5 * ΔR_med)
                    lo = Ri - dR
                    hi = Ri + dR
                    j0 = searchsortedfirst(s_sorted, lo)
                    j1 = searchsortedlast(s_sorted,  hi)
                    if j1 >= j0
                        star_hit_flags[istar] = true
                    else
                        col_vlos[row] = 0.0
                        continue
                    end
                    p        = 0.0
                    nhit     = j1 - j0 + 1
                    inv_norm = inv_sqrt2pi / si
                    @inbounds for j in j0:j1
                        dv = vi - vlos_sorted[j]
                        p += exp(-0.5 * (dv * dv) / s2)
                    end
                    col_vlos[row] = (p * inv_norm) / nhit
                end
            end
            @inbounds A_vlos[:, c_claim] .= col_vlos
            @inbounds A_occ[:,  c_claim] .= col_occ
            Threads.atomic_add!(filled_atomic, 1)
        end
    end  # _worker!
    if threaded && nworkers > 1
        Threads.@threads for t in 1:nworkers
            _worker!(rngs[t])
        end
    else
        _worker!(rngs[1])
    end
    filled = filled_atomic[]
    if filled < fill_target
        missing = fill_target - filled
        Rmed = median(R_star_m)
        inner_fail = 0
        outer_fail = 0
        for i in 1:Norbit
            if !success_flags[i]
                rmin_i = min_r_reached[i]
                if !isfinite(rmin_i)
                    rmin_i = rapo_list[i]
                end
                if isfinite(rmin_i)
                    if rmin_i < Rmed
                        inner_fail += 1
                    else
                        outer_fail += 1
                    end
                end
            end
        end
    println("WARNING: build_A_matrix_hybrid filled ", filled, " / ", fill_target, " | missing ", round(100*(fill_target-filled)/fill_target, digits=1), "%",
        " | inner miss ", round(100*inner_fail/Norbit, digits=1), "%", " | outer miss ", round(100*outer_fail/Norbit, digits=1), "%")
    end
    if Nvlos > 0
        hits = sum(star_hit_flags[vlos_idx])
        total = length(vlos_idx)
        println("STAR COVERAGE: ", hits, " / ", total, " (", round(100*hits/total, digits=1), "%)")
    end
    if Nvlos > 0
        Rv = R_star_m[vlos_idx]
        Rmin = minimum(Rv)
        dR_vals = similar(Rv)
        for i in eachindex(Rv)
            dR_vals[i] = max( dR_frac * Rv[i], 0.5 * local_spacing(R_sorted, Rv[i]), 0.5 * ΔR_med )
        end
        Rlo = Rv .- dR_vals
        inner_effective = minimum(Rlo)
        Rmed = median(Rv)
        println( "INNER SCALE: Rmin=", round(Rmin, digits=3), " | Reff=", round(inner_effective, digits=3), " | Rmed=", round(Rmed, digits=3) )
    end
    A = return_occ ? vcat(A_vlos, A_occ) : A_vlos
    return diag ? (A, Dict("filled" => filled, "attempts" => Norbit)) : A
end

# Batch evaluator (@threads): build A + solve weights + chi2_red per theta. status: 0=pass 1=A_fail 2=weights_fail 3=exception
function evaluate_batch_theta(thetas::AbstractMatrix{<:Real}, R_star_m::Vector{Float64}, valid_vlos::AbstractVector{Bool}, v_star_mps::Vector{Float64}, verr_star_mps::Vector{Float64}, sini::Float64, Norbit::Int, halo_type::String; Nocc::Int=0, lambda_occ::Float64=0.0, alpha::Float64=DEFAULT_ALPHA, maxiter::Int=DEFAULT_MAXITER, max_refine::Int=0)

    nrow, nbatch = size(thetas)
    R_valid  = R_star_m[valid_vlos]
    v_valid  = v_star_mps[valid_vlos]
    ve_valid = verr_star_mps[valid_vlos]
    Nvalid   = length(R_valid)
    nu       = max(count(valid_vlos) - 3, 1)
    rv_mask = convert(Vector{Bool}, trues(Nvalid))
    status = Vector{Int}(undef, nbatch)
    refine_used = zeros(Int, nbatch)
    chi2   = Vector{Float64}(undef, nbatch)
    Threads.@threads :dynamic for i in 1:nbatch
        try
            rho_s = Float64(thetas[1, i])
            r_s   = Float64(thetas[2, i])
            MBH   = nrow >= 3 ? Float64(thetas[3, i]) : 0.0
            # --- build full hybrid matrix (velocity + light) ---
            A = build_A_matrix_hybrid(Norbit, R_star_m, valid_vlos, v_star_mps, verr_star_mps, sini, rho_s, r_s, MBH, halo_type; return_occ=true, Nbins_occ=Nocc, threaded=false)
            m, n = size(A)
            if m == 0 || n == 0 || !all(isfinite, A)
                status[i] = 1
                chi2[i]   = Inf
                continue
            end
            # --- solve weights ---
            w, ok = solve_weights_stellar_jl(A, ve_valid; rv_mask=rv_mask, Nocc=Nocc, lambda_occ=lambda_occ, alpha=alpha, maxiter=maxiter, seed=UInt(i))
            if !ok
                status[i] = 2
                chi2[i]   = Inf
                continue
            end
            # --- likelihood ---
            ll = stellar_log_likelihood_jl(A, w, ve_valid; rv_mask=rv_mask, lambda_occ=lambda_occ, Nocc=Nocc)
            val = -2.0 * ll / nu

            # --- shell refinement (optional, off by default) ---
            if max_refine > 0 && isfinite(val)
                shells_sorted = sort(R_star_m)
                Ns = length(shells_sorted)
                shell_states = Vector{ShellRefinementState}(undef, Ns)
                n_cols = size(A, 2)

                for k in 1:Ns
                    cols_k = collect(k:Ns:min(Norbit, n_cols))
                    Ak = isempty(cols_k) ? zeros(Float64, size(A,1), 1) : Matrix(A[:, cols_k])
                    support_k = isempty(cols_k) ? 0.0 : shell_support(w, cols_k)
                    reff_k = effective_rank(Ak)
                    shell_states[k] = ShellRefinementState(val, support_k, reff_k, Float64[], :active)
                end

                passes_done = 0
                for pass in 1:max_refine
                    passes_done = pass
                    A_r = build_A_matrix_hybrid(Norbit, R_star_m, valid_vlos, v_star_mps, verr_star_mps, sini, rho_s, r_s, MBH, halo_type; return_occ=true, Nbins_occ=Nocc, threaded=false)
                    m_r, n_r = size(A_r)
                    (m_r == 0 || n_r == 0 || !all(isfinite, A_r)) && break

                    w_r, ok_r = solve_weights_stellar_jl(A_r, ve_valid; rv_mask=rv_mask, Nocc=Nocc, lambda_occ=lambda_occ, alpha=alpha, maxiter=maxiter, seed=UInt(i + pass * nbatch))
                    ok_r || break

                    ll_r = stellar_log_likelihood_jl(A_r, w_r, ve_valid; rv_mask=rv_mask, lambda_occ=lambda_occ, Nocc=Nocc)
                    val_r = -2.0 * ll_r / nu

                    any_active = false
                    for k in 1:Ns
                        shell_states[k].status == :frozen && continue
                        cols_k = collect(k:Ns:min(Norbit, n_r))
                        isempty(cols_k) && continue
                        Ak = Matrix(A_r[:, cols_k])
                        support_k = shell_support(w_r, cols_k)
                        update_shell_state!(shell_states[k], Ak, val_r, support_k)
                        shell_states[k].status == :active && (any_active = true)
                    end

                    isfinite(val_r) && val_r < val && (val = val_r)
                    any_active || break
                end
                refine_used[i] = passes_done
            end

            if isfinite(val)
                status[i] = 0
                chi2[i]   = val
            else
                status[i] = 2
                chi2[i]   = Inf
            end
        catch e
            status[i] = 3
            chi2[i]   = Inf
            @warn "evaluate_batch_theta exception on i=$i" exception=(e, catch_backtrace())
        end
    end
    return status, chi2, refine_used
end
end # module
