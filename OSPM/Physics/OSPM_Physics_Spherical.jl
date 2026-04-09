# HOT PATH — A-matrix builder & batch evaluator.
#
# WHAT THIS DOES
# Orbit-superposition physics engine for spherical dark-matter halo modelling.
# Given halo parameters θ = (ρ_s, r_s, M_BH, M/L), builds a gravitational potential,
# launches a library of stellar orbits at various apocenter radii and angular-
# momentum fractions, then projects each orbit onto observables (line-of-sight
# velocity likelihoods + radial occupancy histograms) to form the A-matrix.
# A weight solver (L-BFGS, non-negative) finds the orbital mix that best fits
# kinematic data, and the resulting log-likelihood becomes the χ² score the
# Python daemon minimises. Everything is multithreaded with work-stealing so
# idle cores join whichever θ still has orbits in flight.
#
# LAYOUT  (this file — see OSPM_Physics_Support.jl for halo, orbit, likelihood code)
# ──────────────────────────────────────────────────────────────────────────────
# OrbitWorkState             mutable struct          — per-θ scratch: A-matrices, flags, atomics
# _init_orbit_work(…)       args → OrbitWorkState   — allocate buffers, sort cost order, set fill targets
# _orbit_worker!(st,rng)    state,rng → ()          — inner loop: claim orbit slot, integrate, project, fill A columns
# build_A_matrix_hybrid(…)  Norbit,stars,θ → A      — public entry: spawn workers, assemble vlos+occ A-matrix
# evaluate_batch_theta(…)   θ-matrix → status,χ²    — batch evaluator: build A, solve weights, score, work-steal
#
# OSPM_Physics_Support.jl   (included)
#   constants                G, pc, Msun, DEFAULT_NSTEPS, DEFAULT_LFRAC, …
#   HaloContext              struct: potential/force tables + closures
#   halo_from_theta          (ρ_s,r_s,MBH) → halo dict
#   tables_spherical         halo,R → potential/force/Menc tables
#   launch_orbit_apocenter   rapo,Lz_frac → initial conditions
#   integrate_orbit_rk4      ic → (r,vr,θ,vθ) trajectory
#   solve_weights_stellar_jl A,verr → weights (L-BFGS, ≥0)
#   stellar_log_likelihood_jl A,w,verr → log-likelihood
#   ShellRefinementState     struct + update!/score/saturated helpers
# ──────────────────────────────────────────────────────────────────────────────

module OSPMPhysicsSpherical
@info "OSPMPhysicsSpherical loaded from" @__FILE__
using LinearAlgebra, StaticArrays, Statistics, Random, Base.Threads, Optim
export build_R_halo_physical, rho_interp, halo_from_theta, tables_spherical, make_potential_force_funcs, integrate_orbit_rk4, orbit_to_sigma2_profile, build_A_matrix_julia, ospm_runcheck, build_A_matrix_stellar, build_A_matrix_hybrid, mass_enclosed_two_radii, reset_orbit_cache, evaluate_batch_theta, NTHREADS, stellar_log_likelihood_jl, solve_weights_stellar_jl, ShellRefinementState, effective_rank, shell_support, shell_score, is_saturated, update_shell_state!
include("OSPM_Physics_Support.jl") # Support code: constants, types, caches, helpers, halo physics, orbit integration,
# likelihood, weight solver, and dead/back-compat code.

mutable struct OrbitWorkState
    Norbit::Int;  Nstar::Int;  Nvlos::Int;  Nbins_occ::Int;  Nshells::Int;  nsteps::Int
    max_attempts_factor::Int;  fill_target::Int;  orbit_budget::Int;  return_occ::Bool
    theta0::Float64;  sini::Float64;  ΔR_med::Float64;  Ra_match::Float64;  Rb_match::Float64
    R_star_m::Vector{Float64};  v_star_mps::Vector{Float64};  verr_star_mps::Vector{Float64}
    R_sorted::Vector{Float64};  shells::Vector{Float64};  vlos_idx::Vector{Int}
    occ_edges::Vector{Float64};  cost_order::Vector{Int}
    orbit_ctx;  pot;  frc;  Lfrac
    dt_frac_orbit::Float64;  t_deadline::UInt64
    A_vlos::Matrix{Float64};  A_occ::Matrix{Float64}
    success_flags::BitVector;  min_r_reached::Vector{Float64};  rapo_list::Vector{Float64}
    star_hit_flags::BitVector
    next_orbit::Threads.Atomic{Int};  filled_atomic::Threads.Atomic{Int}
    phase::Threads.Atomic{Int}  # lifecycle: 0=not started, 1=orbit phase, 2=post-orbit, 3=done
end

function _init_orbit_work(Norbit::Int, R_star_m::Vector{Float64}, has_vlos::AbstractVector{Bool}, v_star_mps::Vector{Float64}, verr_star_mps::Vector{Float64}, sini::Float64, ctx;
        nsteps::Int, Lfrac, dt_frac_orbit::Float64, Nbins_occ::Int, return_occ::Bool, max_attempts_factor::Int, fill_pct::Float64, t_deadline::UInt64)
    Nstar = length(R_star_m);  R_sorted = sort(R_star_m);  ΔR_med = median(diff(R_sorted))
    vlos_idx = Int[];  sizehint!(vlos_idx, Nstar)
    @inbounds for i in 1:Nstar
        has_vlos[i] && push!(vlos_idx, i)
    end
    Nvlos = length(vlos_idx);  shells = sort(copy(R_star_m));  Nshells = length(shells)
    occ_edges = build_occ_edges_adaptive(R_star_m, vlos_idx);  Nbins_occ = length(occ_edges) - 1
    A_vlos = zeros(Float64, Nvlos, Norbit);  A_occ = zeros(Float64, Nbins_occ, Norbit)
    success_flags = falses(Norbit);  min_r_reached = fill(Inf, Norbit);  rapo_list = fill(NaN, Norbit);  star_hit_flags = falses(Nstar)
    _orbit_cost = Vector{Float64}(undef, Norbit)
    @inbounds for c in 1:Norbit
        lf   = Lfrac[1 + ((c - 1) % length(Lfrac))]
        rapo = shells[mod1(c, Nshells)]
        _orbit_cost[c] = lf * rapo
    end
    cost_order = sortperm(_orbit_cost);  fill_target = max(1, round(Int, clamp(fill_pct, 0.0, 1.0) * Norbit));  orbit_budget = max_attempts_factor * Norbit
    Ra_match = quantile(R_sorted, 0.25);  Rb_match = quantile(R_sorted, 0.60);  orbit_ctx = (frc=ctx.frc, R_pos=ctx.R, halo=ctx.halo)
    return OrbitWorkState(
        Norbit, Nstar, Nvlos, Nbins_occ, Nshells, nsteps, max_attempts_factor, fill_target, orbit_budget, return_occ,
        f64(pi / 2), clamp01(f64(sini)), ΔR_med, Ra_match, Rb_match, R_star_m, v_star_mps, verr_star_mps, R_sorted, shells, vlos_idx,
        occ_edges, cost_order, orbit_ctx, ctx.pot, ctx.frc, Lfrac, dt_frac_orbit, t_deadline, A_vlos, A_occ, success_flags,
        min_r_reached, rapo_list, star_hit_flags, Threads.Atomic{Int}(1), Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))
end

function _orbit_worker!(st::OrbitWorkState, rng)
    col_occ = zeros(Float64, st.Nbins_occ);  col_vlos = zeros(Float64, st.Nvlos)
    s_arr = Vector{Float64}(undef, st.nsteps);  vlos_buf = Vector{Float64}(undef, st.nsteps)
    while true
        time_ns() > st.t_deadline && break
        st.filled_atomic[] >= st.fill_target && break
        st.phase[] != 1 && break
        c_seq = Threads.atomic_add!(st.next_orbit, 1)
        c_seq > st.orbit_budget && break
        attempt = (c_seq - 1) ÷ st.Norbit
        c_claim = st.cost_order[mod1(c_seq, st.Norbit)]
        attempt > 0 && st.success_flags[c_claim] && continue
        idx_local = mod1(c_claim, st.Nshells)
        rapo = f64(st.shells[idx_local])
        st.rapo_list[c_claim] = rapo
        !(isfinite(rapo) && rapo > 0.0) && continue
        lf = st.Lfrac[1 + ((c_claim - 1) % length(st.Lfrac))]
        r0_frac = 0.95 - 0.05 * f64(attempt) / st.max_attempts_factor + 0.04 * rand(rng)
        ic, Lz0, E0, vc, launch_state = launch_orbit_apocenter(rapo=rapo, theta0=st.theta0, Lz_frac=f64(lf), pot=st.pot, frc=st.frc, r0_frac=r0_frac, dt_frac=st.dt_frac_orbit)
        launch_state != :ok && continue
        r, vr, theta, vtheta = integrate_orbit_rk4( ic=ic, xLz=Lz0, orbit_ctx=st.orbit_ctx, nsteps=st.nsteps )
        isempty(r) && continue
        st.success_flags[c_claim] = true;  st.min_r_reached[c_claim] = minimum(r)
        Nhits = length(r);  dt_orb = f64(ic[3]);  resize!(s_arr, Nhits);  resize!(vlos_buf, Nhits)
        phi = 0.0
        @inbounds for i in 1:Nhits
            si = _ssin(f64(theta[i]))
            ri = f64(r[i])
            s_arr[i] = ri * si
            vphi_i = f64(Lz0) / (ri * si)
            cp, sp = cos(phi), sin(phi)
            vlos_buf[i] = st.sini * (f64(vr[i]) * cp - vphi_i * sp)
            phi += f64(Lz0) / (ri * ri * si * si) * dt_orb
        end
        fill!(col_occ, 0.0);  fill!(col_vlos, 0.0)
        if st.return_occ && Nhits > 0
            @inbounds for j in 1:st.Nbins_occ
                lo = st.occ_edges[j]
                hi = st.occ_edges[j + 1]
                cnt = 0
                @inbounds for x in s_arr
                    cnt += (x >= lo && x < hi) ? 1 : 0
                end
                col_occ[j] = cnt / Nhits
            end
        end
        if st.Nvlos > 0
            pidx = sortperm(s_arr);  s_sorted = s_arr[pidx];  vlos_sorted = vlos_buf[pidx];  inv_sqrt2pi = inv(sqrt(2 * pi))
            @inbounds for row in 1:st.Nvlos
                istar = st.vlos_idx[row];  Ri = f64(st.R_star_m[istar]);  vi = f64(st.v_star_mps[istar])
                si = f64(st.verr_star_mps[istar]);  s2 = si * si
                if !(isfinite(s2) && s2 > 0.0)
                    col_vlos[row] = 0.0
                    continue
                end
                dR_local = local_spacing(st.R_sorted, Ri)
                dR_floor = 0.35 * st.ΔR_med
                dR = Ri < st.Ra_match ? max(0.60 * dR_local, dR_floor) :
                    Ri < st.Rb_match ? max(0.80 * dR_local, dR_floor) :
                                        max(1.00 * dR_local, dR_floor)
                lo = Ri - dR;  hi = Ri + dR;  j0 = searchsortedfirst(s_sorted, lo);  j1 = searchsortedlast(s_sorted, hi)
                if j1 >= j0
                    st.star_hit_flags[istar] = true
                else
                    col_vlos[row] = 0.0
                    continue
                end
                p = 0.0;  nhit = j1 - j0 + 1;  inv_norm = inv_sqrt2pi / si
                @inbounds for j in j0:j1
                    dv = vi - vlos_sorted[j]
                    p += exp(-0.5 * (dv * dv) / s2)
                end
                col_vlos[row] = (p * inv_norm) / nhit
            end
        end
        @inbounds st.A_vlos[:, c_claim] .= col_vlos
        @inbounds st.A_occ[:,  c_claim] .= col_occ
        Threads.atomic_add!(st.filled_atomic, 1)
    end
end

# Main A-matrix builder: maps orbital weights → observables (vlos likelihoods + occupancy). Parallel, numerically hardened.
function build_A_matrix_hybrid(Norbit::Int, R_star_m::Vector{Float64}, has_vlos::AbstractVector{Bool}, v_star_mps::Vector{Float64}, verr_star_mps::Vector{Float64}, sini::Float64, rho_s::Float64, r_s::Float64, MBH::Float64, ML::Float64, halo_type::String; stellar_model=nothing, nsteps::Int=DEFAULT_NSTEPS, Lfrac::NTuple{5,Float64}=DEFAULT_LFRAC, dt_frac_orbit::Float64=DEFAULT_DT_FRAC, dR_frac::Float64=DEFAULT_DR_FRAC, Nbins_occ::Int=DEFAULT_NBINS_OCC, return_occ::Bool=true, max_attempts_factor::Int=DEFAULT_MAX_ATTEMPTS, diag::Bool=false, threaded::Bool=true, fill_pct::Float64=0.80, t_deadline::UInt64=typemax(UInt64))
    Nstar = length(R_star_m)
    @assert length(has_vlos) == Nstar; @assert length(v_star_mps) == Nstar; @assert length(verr_star_mps) == Nstar
    Nstar == 0 && return zeros(Float64, 0, Norbit)
    stellar_model_jl = normalize_stellar_model(stellar_model)
    ctx = get_halo_context(rho_s, r_s, MBH, ML, halo_type; stellar_model=stellar_model_jl)
    sini = clamp01(f64(sini)); Rmin = minimum(R_star_m); Rmax = maximum(R_star_m)

    if !(isfinite(Rmin) && isfinite(Rmax) && Rmax > Rmin)
        Nvlos = count(has_vlos);  Nout = Nvlos + (return_occ ? Nbins_occ : 0)
        return zeros(Float64, Nout, Norbit)
    end
    st = _init_orbit_work( Norbit, R_star_m, has_vlos, v_star_mps, verr_star_mps, sini, ctx;
        nsteps=nsteps, Lfrac=Lfrac, dt_frac_orbit=dt_frac_orbit, Nbins_occ=Nbins_occ, return_occ=return_occ,
        max_attempts_factor=max_attempts_factor, fill_pct=fill_pct, t_deadline=t_deadline )
    # standalone path: phase is always 1 (no helpers to coordinate with)
    Threads.atomic_xchg!(st.phase, 1)
    nworkers = threaded ? Threads.nthreads() : 1;  rngs = [MersenneTwister(0x5eed1234 + UInt(t)) for t in 1:nworkers]
    if threaded && nworkers > 1
        Threads.@threads for t in 1:nworkers
            _orbit_worker!(st, rngs[t])
        end
    else
        _orbit_worker!(st, rngs[1])
    end
    filled = st.filled_atomic[]
    if filled < st.fill_target
        Rmed = median(R_star_m)
        inner_fail = 0;  outer_fail = 0
        for i in 1:Norbit
            if !st.success_flags[i]
                rmin_i = st.min_r_reached[i]
                if !isfinite(rmin_i)
                    rmin_i = st.rapo_list[i]
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
        println("WARNING: build_A_matrix_hybrid filled ", filled, " / ", st.fill_target, " | missing ", round(100*(st.fill_target-filled)/st.fill_target, digits=1), "%",
            " | inner miss ", round(100*inner_fail/Norbit, digits=1), "%", " | outer miss ", round(100*outer_fail/Norbit, digits=1), "%")
    end
    if st.Nvlos > 0
        hits = sum(st.star_hit_flags[st.vlos_idx]);  total = length(st.vlos_idx)
        println("STAR COVERAGE: ", hits, " / ", total, " (", round(100*hits/total, digits=1), "%)")
    end
    if st.Nvlos > 0
        Rv = st.R_star_m[st.vlos_idx]
        Rmin_v = minimum(Rv)
        dR_vals = similar(Rv)
        for i in eachindex(Rv)
            dR_local = local_spacing(st.R_sorted, Rv[i])
            if Rv[i] < st.Ra_match
                dR_vals[i] = max(0.35 * dR_local, 0.15 * st.ΔR_med)
            elseif Rv[i] < st.Rb_match
                dR_vals[i] = max(0.45 * dR_local, 0.15 * st.ΔR_med)
            else
                dR_vals[i] = max(0.60 * dR_local, 0.15 * st.ΔR_med)
            end
        end
        Rlo = Rv .- dR_vals
        inner_effective = minimum(Rlo)
        Rmed_v = median(Rv)
        println("INNER SCALE: Rmin=", round(Rmin_v, digits=3), " | Reff=", round(inner_effective, digits=3), " | Rmed=", round(Rmed_v, digits=3))
    end
    A = return_occ ? vcat(st.A_vlos, st.A_occ) : st.A_vlos
    return diag ? (A, Dict("filled" => filled, "attempts" => Norbit)) : A
end

# Batch evaluator w/ work-stealing.
# Phase 1: each thread claims thetas, inits OrbitWorkState, runs orbits, solves weights+likelihood.
# Phase 2: idle threads join any theta still in orbit phase. Keeps all cores saturated.
function evaluate_batch_theta(thetas::AbstractMatrix{<:Real}, R_star_m::Vector{Float64}, valid_vlos::AbstractVector{Bool}, v_star_mps::Vector{Float64}, verr_star_mps::Vector{Float64}, sini::Float64, Norbit::Int, halo_type::String; stellar_model=nothing, Nocc::Int=0, lambda_occ::Float64=0.0, alpha::Float64=DEFAULT_ALPHA, maxiter::Int=DEFAULT_MAXITER, max_refine::Int=0, timeout_s::Float64=120.0)
    stellar_model_jl = normalize_stellar_model(stellar_model)
    nrow, nbatch = size(thetas); ve_valid = verr_star_mps[valid_vlos]; Nvalid = length(ve_valid)
    nu = max(count(valid_vlos) - 4, 1); rv_mask = convert(Vector{Bool}, trues(Nvalid))
   
    status = fill(4, nbatch);  refine_used = zeros(Int, nbatch);  chi2 = fill(Inf, nbatch)  # default = timeout/unprocessed
    # Shared orbit state per theta — visible to all threads for helping
    work_states = Vector{Union{Nothing, OrbitWorkState}}(undef, nbatch);  fill!(work_states, nothing)
    next_theta = Threads.Atomic{Int}(1);  nthreads = Threads.nthreads()
    helper_rngs = [MersenneTwister(0x0E100000 + UInt(t) + UInt(nbatch)) for t in 1:nthreads]

    function _batch_worker!(tid::Int)
        rng_own = MersenneTwister(0x5eed1234 + UInt(tid));  rng_help = helper_rngs[tid]

        while true
            # ---- Phase 1: claim a theta ----
            i = Threads.atomic_add!(next_theta, 1)
            if i <= nbatch
                theta_deadline = time_ns() + UInt64(round(timeout_s * 1e9))
                if time_ns() > theta_deadline
                    status[i] = 4
                    chi2[i] = Inf
                    continue
                end
                try
                    rho_s = Float64(thetas[1, i]); r_s = Float64(thetas[2, i]); MBH = nrow >= 3 ? Float64(thetas[3, i]) : 0.0; ML = nrow >= 4 ? Float64(thetas[4, i]) : 1.0
                    Rmin_v = minimum(R_star_m);  Rmax_v = maximum(R_star_m)
                    if !(isfinite(Rmin_v) && isfinite(Rmax_v) && Rmax_v > Rmin_v) || length(R_star_m) == 0
                        status[i] = 1
                        chi2[i] = Inf
                        continue
                    end
                    ctx = get_halo_context(rho_s, r_s, MBH, ML, halo_type; stellar_model=stellar_model_jl)
                    ws = _init_orbit_work( Norbit, R_star_m, valid_vlos, v_star_mps, verr_star_mps, sini, ctx; nsteps=DEFAULT_NSTEPS, Lfrac=DEFAULT_LFRAC,
                        dt_frac_orbit=DEFAULT_DT_FRAC, Nbins_occ=Nocc, return_occ=true, max_attempts_factor=DEFAULT_MAX_ATTEMPTS, fill_pct=0.80, t_deadline=theta_deadline )
                    work_states[i] = ws  # publish and mark joinable
                    Threads.atomic_xchg!(ws.phase, 1)
                    _orbit_worker!(ws, rng_own)  # owner does bulk of orbit work
                    Threads.atomic_xchg!(ws.phase, 2)  # close orbit phase
                    # Assemble A-matrix
                    A = ws.return_occ ? vcat(ws.A_vlos, ws.A_occ) : ws.A_vlos
                    m, n = size(A)
                    if m == 0 || n == 0 || !all(isfinite, A)
                        status[i] = 1
                        chi2[i] = Inf
                        Threads.atomic_xchg!(ws.phase, 3)
                        continue
                    end
                    # Solve weights
                    Nocc_eff = max(size(A, 1) - length(ve_valid), 0)
                    w, ok = solve_weights_stellar_jl(A, ve_valid; rv_mask=rv_mask, Nocc=Nocc_eff, lambda_occ=lambda_occ, alpha=alpha, maxiter=maxiter, seed=UInt(i))
                    if !ok
                        status[i] = 2
                        chi2[i] = Inf
                        Threads.atomic_xchg!(ws.phase, 3)
                        continue
                    end
                    # Likelihood
                    ll = stellar_log_likelihood_jl(A, w, ve_valid; rv_mask=rv_mask, lambda_occ=lambda_occ, Nocc=Nocc_eff)
                    val = -2.0 * ll / nu
                    # Shell refinement (optional, off by default)
                    if max_refine > 0 && isfinite(val)
                        shells_sorted = sort(R_star_m);  Ns = length(shells_sorted)
                        shell_states = Vector{ShellRefinementState}(undef, Ns);  n_cols = size(A, 2)
                        for k in 1:Ns
                            cols_k = collect(k:Ns:min(Norbit, n_cols))
                            Ak = isempty(cols_k) ? zeros(Float64, size(A,1), 1) : Matrix(A[:, cols_k])
                            support_k = isempty(cols_k) ? 0.0 : shell_support(w, cols_k)
                            reff_k = effective_rank(Ak)
                            shell_states[k] = ShellRefinementState(val, support_k, reff_k, Float64[], :active)
                        end
                        passes_done = 0
                        for pass_i in 1:max_refine
                            passes_done = pass_i
                            time_ns() > theta_deadline && break
                            A_r = build_A_matrix_hybrid(Norbit, R_star_m, valid_vlos, v_star_mps, verr_star_mps, sini, rho_s, r_s, MBH, ML, halo_type; stellar_model=stellar_model_jl, return_occ=true, Nbins_occ=Nocc, threaded=false, fill_pct=0.80, t_deadline=theta_deadline)
                            m_r, n_r = size(A_r);  (m_r == 0 || n_r == 0 || !all(isfinite, A_r)) && break
                            Nocc_eff_r = max(size(A_r, 1) - length(ve_valid), 0)
                            w_r, ok_r = solve_weights_stellar_jl(A_r, ve_valid; rv_mask=rv_mask, Nocc=Nocc_eff_r, lambda_occ=lambda_occ, alpha=alpha, maxiter=maxiter, seed=UInt(i + pass_i * nbatch))
                            ok_r || break
                            ll_r = stellar_log_likelihood_jl(A_r, w_r, ve_valid; rv_mask=rv_mask, lambda_occ=lambda_occ, Nocc=Nocc_eff_r)
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
                        chi2[i] = val
                    else
                        status[i] = 2
                        chi2[i] = Inf
                    end
                    Threads.atomic_xchg!(ws.phase, 3)
                catch e
                    status[i] = 3
                    chi2[i] = Inf
                    ws_i = work_states[i]
                    ws_i !== nothing && Threads.atomic_xchg!(ws_i.phase, 3)
                    @warn "evaluate_batch_theta exception on i=$i" exception=(e, catch_backtrace())
                end
                continue  # try to claim another theta
            end
            # ---- Phase 2: no unclaimed thetas — help stragglers ----
            helped = false
            for scan in 1:nbatch
                ws_scan = work_states[scan]
                ws_scan === nothing && continue
                ws_scan.phase[] != 1 && continue
                time_ns() > ws_scan.t_deadline && continue
                _orbit_worker!(ws_scan, rng_help)
                helped = true
                break  # re-scan (a different theta may now be the bottleneck)
            end
            helped || break  # nothing left to help — this thread is done
        end
    end
    Threads.@threads :static for t in 1:nthreads
        _batch_worker!(t)
    end
    return status, chi2, refine_used
end
end # module