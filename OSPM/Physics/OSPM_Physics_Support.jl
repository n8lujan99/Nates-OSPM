# ============================================================
# OSPM_Physics_Support.jl
# Included by OSPM_Physics_Spherical.jl — do NOT load directly.
# Contains: constants, types, caches, helpers, halo physics, orbit integration, likelihood, weight solver.
# Dead code (not on the evaluate_batch_theta hot path) is
# collected at the bottom and clearly labelled.
# ============================================================
# §1  CONSTANTS
const NTHREADS = Threads.nthreads()
const G    = 6.67430e-11
const c    = 2.99792458e8
const pc   = 3.0856775814913673e16
const Msun = 1.98847e30
# machine floors
const EPS_FORCE = 1e-14
const EPS_VEL   = 1e-14
const EPS_ARG   = 1e-14
# physical geometry gate
const EPS_SIN = 1e-6
# scale-aware force gate
const REL_FORCE    = 1e-10   # loosen to 1e-9 if needed
const BRACKET_FRAC = 1e-6    # MUST be >> eps(Float64)

# TUNABLE KNOBS — adjust these to control resolution, accuracy, and parallelism.
# -- Halo potential grid --
const DEFAULT_NR              = 256       # radial grid points for potential table
const DEFAULT_RMAX_FACTOR     = 300.0     # max radius in units of r_s
# -- Orbit integration --
const DEFAULT_NSTEPS          = 4000      # RK4 steps per orbit
const DEFAULT_STOP_RMIN_FACTOR = 1.001    # orbit stops when r < factor * rmin
const DEFAULT_DT_FRAC         = 0.01      # timestep = dt_frac / orbital_frequency
const DEFAULT_DT_FLOOR        = 1e-30     # floor on orbital-frequency denominator
const DEFAULT_R0_FRAC         = 0.98      # starting radius as fraction of apocenter
# -- A-matrix / orbit library --
const DEFAULT_LFRAC           = (0.05, 0.2, 0.4, 0.7, 1.0)  # angular momentum fractions
const DEFAULT_DR_FRAC         = 0.05      # radial matching tolerance (fraction of R)
const DEFAULT_NBINS_OCC       = 6         # occupancy histogram bins
const DEFAULT_MAX_ATTEMPTS    = 60        # max orbit-launch attempts multiplier
const DEFAULT_DR_FLOOR_FRAC   = 0.01      # floor on dR (fraction)
const DEFAULT_DR_FLOOR_PC     = 0.0       # floor on dR (parsecs)
# -- Weight solver --
const DEFAULT_ALPHA           = 1e-2      # L2 regularization strength
const DEFAULT_MAXITER         = 150       # LBFGS iteration cap
const DEFAULT_LOG_EPS         = 1e-300    # log-likelihood floor to avoid log(0)
const DEFAULT_SIGMA_FLOOR_MPS = 2e3       # velocity-error floor [m/s]

# ============================================================
# §2  TYPES, CACHES, INLINE HELPERS
# ============================================================
@inline f64(x)=Float64(x)
@inline safe_sign(x)=x>0 ? 1.0 : (x<0 ? -1.0 : 0.0)
@inline _ssin(theta::Float64)=begin s=sin(theta); abs(s)>1e-12 ? s : safe_sign(s)*1e-12 end
@inline function _sincos_safe(theta::Float64); s,cc=sincos(theta); abs(s)>1e-12 ? (s,cc) : (safe_sign(s)*1e-12,cc) end
@inline clamp01(x::Float64)=x<0 ? 0.0 : (x>1 ? 1.0 : x)
struct HaloContext
    halo::Dict{Symbol,Any}
    R::Vector{Float64}
    tabv::Vector{Float64}
    tabfr::Vector{Float64}
    Menc::Vector{Float64}
    pot::Function
    frc::Function
end
const _HALO_CTX_CACHE = Dict{Tuple{Float64,Float64,Float64,Float64,UInt64,Symbol,Int,Float64},HaloContext}()
const _HALO_LOCK = ReentrantLock()

# ============================================================
# §3  SMALL UTILITIES
# ============================================================
function build_occ_edges_adaptive(R_star_m::Vector{Float64}, vlos_idx::Vector{Int};
    zone_quantiles::Tuple{Float64,Float64} = (0.25, 0.60),
    base_bins::NTuple{3,Int} = (3, 2, 1),          # inner, mid, outer
    max_extra::NTuple{3,Int} = (2, 1, 1),
    min_stars_per_bin::Int = 3,
    res_factor::Float64 = 2.0)

    # use only the stars that actually speak in the velocity block
    R_use = isempty(vlos_idx) ? copy(R_star_m) : R_star_m[vlos_idx]
    R_use = sort(R_use[isfinite.(R_use)])

    n = length(R_use)
    if n == 0
        return [0.0, 1.0]
    elseif n == 1
        r0 = R_use[1]
        return [max(0.0, 0.9r0), 1.1r0]
    end

    Rmin = R_use[1]
    Rmax = R_use[end]
    Rmax <= Rmin && return [Rmin, Rmax + max(abs(Rmax), 1.0)]

    q_in, q_mid = zone_quantiles
    Ra = quantile(R_use, q_in)
    Rb = quantile(R_use, q_mid)

    # keep the zones ordered and non-degenerate
    ΔR_med = median(diff(R_use))
    floor_sep = max(1e-12, 0.5 * ΔR_med)

    Ra = clamp(Ra, Rmin + floor_sep, Rmax - 2floor_sep)
    Rb = clamp(Rb, Ra   + floor_sep, Rmax - floor_sep)

    function zone_refine(Rz::Vector{Float64}, zlo::Float64, zhi::Float64, base::Int, extra_max::Int)
        nz = length(Rz)
        if nz < 2 || !(zhi > zlo)
            return max(base, 1)
        end
        dz = diff(sort(Rz))
        dz = dz[isfinite.(dz) .& (dz .> 0.0)]
        dmed = isempty(dz) ? (zhi - zlo) : median(dz)
        zwidth = zhi - zlo
        # count-based refinement
        n_by_count = fld(nz, min_stars_per_bin)
        # resolution-based refinement
        n_by_res = floor(Int, zwidth / max(res_factor * dmed, 1e-12))
        n_target = min(n_by_count, n_by_res)
        n_bins = base + clamp(n_target - base, 0, extra_max)
        return max(n_bins, 1)
    end

    Rin  = R_use[R_use .<  Ra]
    Rmid = R_use[(R_use .>= Ra) .& (R_use .<  Rb)]
    Rout = R_use[R_use .>= Rb]
    Nin  = zone_refine(Rin,  Rmin, Ra,  base_bins[1], max_extra[1])
    Nmid = zone_refine(Rmid, Ra,   Rb,  base_bins[2], max_extra[2])
    Nout = zone_refine(Rout, Rb,   Rmax, base_bins[3], max_extra[3])
    # do not let the outer cloud dominate the occupancy language
    Nout = min(Nout, max(1, Nmid))
    ein  = collect(range(Rmin, Ra;   length=Nin  + 1))
    emid = collect(range(Ra,   Rb;   length=Nmid + 1))
    eout = collect(range(Rb,   Rmax; length=Nout + 1))
    occ_edges = vcat(ein, emid[2:end], eout[2:end])
    # final monotonic cleanup
    occ_edges = sort(unique(occ_edges))
    if length(occ_edges) < 2
        occ_edges = [Rmin, Rmax]
    end
    return occ_edges
end

@inline function normalize_halo(halo)
    h=Dict{Symbol,Any}()
    for (k,v) in halo
        h[k isa Symbol ? k : Symbol(String(k))]=v
    end
    if haskey(h,:type) && !(h[:type] isa Symbol)
        h[:type]=Symbol(lowercase(String(h[:type])))
    end
    h
end
logspace10(a,b,n)=n==1 ? [10.0^a] : (da=(b-a)/(n-1); [10.0^(a+(i-1)*da) for i in 1:n])
build_R_halo_physical(n; rmin=1e-3, rmax=300.0)=logspace10(log10(rmin), log10(rmax), n)
@inline function _quant(x::Float64; digits::Int=10)
    return round(x, digits=digits)
end
function local_spacing(R_sorted, Ri)
    idx = searchsortedfirst(R_sorted, Ri)
    if idx <= 1
        return R_sorted[2] - R_sorted[1]
    elseif idx >= length(R_sorted)
        return R_sorted[end] - R_sorted[end-1]
    else
        return 0.5 * (R_sorted[idx+1] - R_sorted[idx-1])
    end
end

# ============================================================
# §4  HALO PHYSICS
# ============================================================
function rho_interp(rv, halo)
    r    = abs(rv[1])
    rhos = halo[:rho_s]
    rs   = halo[:r_s]
    x    = r / max(rs, 1e-30)
    halo[:type] === :nfw &&
        return rhos / (x * (1 + x)^2 + 1e-30)
    halo[:type] === :cored &&
        return rhos / ((1 + x) * (1 + x^2) + 1e-30)
    halo[:type] === :einasto && begin
        α = halo[:alpha]          # curvature parameter
        return rhos * exp(-2/α * (x^α - 1))
    end
    error("Unknown halo type: $(halo[:type])")
end

@inline function normalize_stellar_model(stellar_model)
    stellar_model === nothing && return nothing
    out = Dict{Symbol,Any}()
    for (k, v) in stellar_model
        ks = k isa Symbol ? k : Symbol(String(k))
        out[ks] = v
    end
    return out
end


@inline function stellar_model_sig(stellar_model)
    stellar_model === nothing && return UInt(0)
    return hash((get(stellar_model, :type, nothing), get(stellar_model, :Ltot, nothing), get(stellar_model, :a_pc, nothing)))
end

@inline function stellar_Menc_plummer(r::Float64, ML::Float64, Ltot::Float64, a::Float64)
    rr = max(r, 1e-30)
    Mtot = ML * Ltot * Msun
    return Mtot * rr^3 / (rr^2 + a^2)^(3/2)
end

@inline function stellar_Phi_plummer(r::Float64, ML::Float64, Ltot::Float64, a::Float64)
    rr = max(r, 1e-30)
    Mtot = ML * Ltot * Msun
    return -G * Mtot / sqrt(rr^2 + a^2)
end

function halo_from_theta(rho_s, r_s, MBH, ML; halo_type="nfw", alpha=nothing, stellar_model=nothing)
    ht = Symbol(lowercase(String(halo_type)))
    h = Dict(
        :rho_s => f64(rho_s) * Msun / pc^3,
        :r_s   => f64(r_s) * pc,
        :rs    => f64(r_s) * pc,
        :MBH   => f64(MBH) * Msun,
        :ML    => f64(ML),
        :type  => ht,
        :rmin  => 1e-6 * f64(r_s) * pc,
    )
    stellar_model !== nothing && (h[:stellar_model] = normalize_stellar_model(stellar_model))
    if ht === :einasto
        h[:alpha] = isnothing(alpha) ? 0.18 : f64(alpha)
    end
    return h
end

function tables_spherical(R, nlegup, halo, rhofn)
    halo=normalize_halo(halo); n=length(R)
    rho=similar(R); tabv=zeros(n); tabfr=zeros(n); Menc=zeros(n)

    @inbounds for i in eachindex(R)
        v=rhofn((R[i],0.0),halo)
        rho[i]=isfinite(v) ? v : 0.0
    end

    @inbounds for i in 2:n
        dr=R[i]-R[i-1]
        Menc[i]=Menc[i-1]+0.5*dr*(R[i]^2*rho[i]+R[i-1]^2*rho[i-1])
    end
    Menc .*= 4*pi

    J=zeros(n)
    @inbounds for i in (n-1):-1:1
        dr=R[i+1]-R[i]
        J[i]=J[i+1]+0.5*dr*(R[i+1]*rho[i+1] + R[i]*rho[i])
    end
    J .*= 4*pi

    @inbounds for i in eachindex(R)
        r=max(R[i],1e-30)
        tabv[i]  = -G*(Menc[i]/r + J[i])
        tabfr[i] = -(G*Menc[i])/(r*r)
    end

    tabv, tabfr, Menc
end

function make_potential_force_funcs(halo, R, nlegup, tabv, tabfr, Menc)
    halo = normalize_halo(halo)
    MBH  = f64(halo[:MBH])
    ML   = haskey(halo, :ML) ? f64(halo[:ML]) : 0.0
    rmin = f64(halo[:rmin])

    stellar_model = get(halo, :stellar_model, nothing)
    has_stars = stellar_model !== nothing && ML > 0.0

    rlgmin = log10(f64(R[1]))
    rlgmax = log10(f64(R[end]))
    np = length(R)
    rlgmax > rlgmin || error("Degenerate R grid")

    @inline function interp(arr, rr)
        r = max(f64(rr), rmin)
        lr = log10(r)
        x = (lr - rlgmin) * (np - 1) / (rlgmax - rlgmin)
        x = clamp(x, 0.0, np - 1.0)
        i0 = Int(floor(x)) + 1
        i1 = min(i0 + 1, np)
        t = x - (i0 - 1)
        (1 - t) * arr[i0] + t * arr[i1]
    end

    @inline function Mstar_enc(rr)
        if !has_stars
            return 0.0
        end
        stype = Symbol(lowercase(String(stellar_model[:type])))
        if stype === :plummer
            Ltot = f64(stellar_model[:Ltot])
            a    = f64(stellar_model[:a_pc]) * pc
            return stellar_Menc_plummer(rr, ML, Ltot, a)
        else
            error("Unknown stellar model type: $(stellar_model[:type])")
        end
    end

    @inline function Phistar(rr)
        if !has_stars
            return 0.0
        end
        stype = Symbol(lowercase(String(stellar_model[:type])))
        if stype === :plummer
            Ltot = f64(stellar_model[:Ltot])
            a    = f64(stellar_model[:a_pc]) * pc
            return stellar_Phi_plummer(rr, ML, Ltot, a)
        else
            error("Unknown stellar model type: $(stellar_model[:type])")
        end
    end

    pot(r, mu=0.0) = begin
        rr = max(abs(f64(r)), rmin)
        Ph   = interp(tabv, rr)
        Pbh  = MBH > 0 ? (-G * MBH / rr) : 0.0
        Pstr = has_stars ? Phistar(rr) : 0.0
        Ph + Pbh + Pstr
    end

    frc(r, mu=0.0) = begin
        rr = max(abs(f64(r)), rmin)
        frh  = interp(tabfr, rr)
        frbh = MBH > 0 ? (-G * MBH / (rr * rr)) : 0.0
        frst = has_stars ? (-G * Mstar_enc(rr) / (rr * rr)) : 0.0
        frh + frbh + frst, 0.0
    end

    pot, frc, R
end

function build_halo_context(rho_s, r_s, MBH, ML, halo_type; stellar_model=nothing, nR=DEFAULT_NR, rmax_factor=DEFAULT_RMAX_FACTOR)
    halo = halo_from_theta(rho_s, r_s, MBH, ML; halo_type=halo_type, stellar_model=stellar_model)
    R = build_R_halo_physical(nR; rmin=halo[:rmin], rmax=rmax_factor * halo[:rs])
    tabv, tabfr, Menc = tables_spherical(R, 1, halo, rho_interp)
    pot, frc, _ = make_potential_force_funcs(halo, R, 1, tabv, tabfr, Menc)
    HaloContext(halo, f64.(R), tabv, tabfr, Menc, pot, frc)
end

function get_halo_context(rho_s, r_s, MBH, ML, halo_type; stellar_model=nothing, nR=DEFAULT_NR, rmax_factor=DEFAULT_RMAX_FACTOR)
    ht = Symbol(lowercase(String(halo_type)))
    sig = stellar_model_sig(stellar_model)
    key = (_quant(f64(rho_s)), _quant(f64(r_s)), _quant(f64(MBH)), _quant(f64(ML)), sig, ht, nR, _quant(f64(rmax_factor)))
    lock(_HALO_LOCK)
    ctx = get(_HALO_CTX_CACHE, key, nothing)
    unlock(_HALO_LOCK)
    ctx !== nothing && return ctx
    newctx = build_halo_context(rho_s, r_s, MBH, ML, ht; stellar_model=stellar_model, nR=nR, rmax_factor=rmax_factor)
    lock(_HALO_LOCK)
    ctx = get(_HALO_CTX_CACHE, key, nothing)
    if ctx === nothing
        _HALO_CTX_CACHE[key] = newctx
        ctx = newctx
    end
    unlock(_HALO_LOCK)
    return ctx
end

# ============================================================
# §5  ORBIT INTEGRATION
# ============================================================
@inline function derivs(s::SVector{4,Float64}, Lz::Float64, frc, R)
    r,theta,vr,vtheta=s
    !(isfinite(r)&&isfinite(theta)&&isfinite(vr)&&isfinite(vtheta)) && return SVector(0.0,0.0,0.0,0.0)
    r_safe=max(abs(r),1e-12)
    st,ct=_sincos_safe(theta)
    r_tab=clamp(r_safe,R[1],R[end])
    fr,_=frc(r_tab,st)
    !isfinite(fr) && return SVector(0.0,0.0,0.0,0.0)
    dr=vr
    dtheta=vtheta/r_safe
    dvr=(vtheta*vtheta)/r_safe + (Lz*Lz)/(r_safe^3*st*st) + fr
    dvtheta=(Lz*Lz)*ct/(r_safe^3*st^3) - (vr*vtheta)/r_safe
    SVector(dr,dtheta,dvr,dvtheta)
end
function launch_orbit_apocenter(; rapo::Float64, theta0::Float64, Lz_frac::Float64, pot, frc, r0_frac::Float64=DEFAULT_R0_FRAC, dt_frac::Float64=DEFAULT_DT_FRAC, dt_floor::Float64=DEFAULT_DT_FLOOR, debug::Bool=true)
    ss = _ssin(theta0)
    if !(isfinite(ss) && abs(ss) > EPS_SIN)
        return (nothing, 0.0, 0.0, 0.0, :reject_sin)
    end
    frs, _ = frc(rapo, ss)
    if !(isfinite(frs) && isfinite(rapo) && rapo > 0.0)
        return (nothing, 0.0, 0.0, 0.0, :reject_force)
    end
    r_in  = rapo * (1 - BRACKET_FRAC)
    r_out = rapo * (1 + BRACKET_FRAC)
    fr_in,  _ = frc(r_in,  ss)
    fr_out, _ = frc(r_out, ss)
    fr_scale = max(abs(frs), abs(fr_in), abs(fr_out), EPS_FORCE)
    fr_tol   = max(EPS_FORCE, REL_FORCE * fr_scale)
    if frs > fr_tol
        return debug ?
            ((rapo, theta0, ss, frs, fr_tol, fr_scale), 0.0, 0.0, 0.0, :reject_force) :
            (nothing, 0.0, 0.0, 0.0, :reject_force)
    end
    vc2 = (-frs) * rapo
    if vc2 <= 0.0
        vc2 = fr_tol * rapo
    end
    vc = sqrt(vc2)

    if !(isfinite(vc) && vc > EPS_VEL)
        return debug ?
            ((rapo, theta0, ss, frs, vc, EPS_VEL), 0.0, 0.0, 0.0, :reject_vc) :
            (nothing, 0.0, 0.0, 0.0, :reject_vc)
    end
    Lz = Lz_frac * rapo * vc
    Papo = pot(rapo, ss)
    if !isfinite(Papo)
        return (nothing, 0.0, 0.0, vc, :reject_pot)
    end
    E = Papo + (Lz^2) / (2 * rapo^2 * ss^2)
    r0 = r0_frac * rapo
    P0 = pot(r0, ss)
    if !isfinite(P0)
        return (nothing, 0.0, E, vc, :reject_pot0)
    end
    arg = 2 * (E - P0) - (Lz^2) / (r0^2 * ss^2)
    if !(isfinite(arg) && arg > -EPS_ARG)
        return debug ?
            ((rapo, theta0, Lz, arg), Lz, E, vc, :reject_turning) :
            (nothing, Lz, E, vc, :reject_turning)
    end
    vr0 = -sqrt(max(arg, 0.0))
    Om  = abs(vc / r0)
    dt  = dt_frac / max(Om, dt_floor)
    return ((r0, theta0, dt, vr0, 0.0), Lz, E, vc, :ok)
end
function integrate_orbit_rk4(; ic, xLz, orbit_ctx, nsteps=DEFAULT_NSTEPS, stop_rmin_factor=DEFAULT_STOP_RMIN_FACTOR)
    halo = orbit_ctx.halo
    rmin_stop = stop_rmin_factor * f64(halo[:rmin])
    r0      = f64(ic[1])
    theta0  = f64(ic[2])
    dt      = f64(ic[3])
    vr0     = length(ic)>=4 ? f64(ic[4]) : 0.0
    vtheta0 = length(ic)>=5 ? f64(ic[5]) : 0.0
    state = SVector(r0,theta0,vr0,vtheta0)
    ns = Int(nsteps)
    r      = Vector{Float64}(undef, ns)
    vr     = Vector{Float64}(undef, ns)
    theta  = Vector{Float64}(undef, ns)
    vtheta = Vector{Float64}(undef, ns)
    rmax_stop = 10.0 * f64(orbit_ctx.R_pos[end])
    actual = 0
    @inbounds for step in 1:ns
        !all(isfinite,state) && break
        rr = state[1]
        tr = state[2]
        (rr<=rmin_stop || rr>=rmax_stop || abs(tr)>1e6) && break
        actual += 1
        r[actual]      = rr
        vr[actual]     = state[3]
        theta[actual]  = tr
        vtheta[actual] = state[4]
        k1 = derivs(state, xLz, orbit_ctx.frc, orbit_ctx.R_pos)
        k2 = derivs(state + 0.5*dt*k1, xLz, orbit_ctx.frc, orbit_ctx.R_pos)
        k3 = derivs(state + 0.5*dt*k2, xLz, orbit_ctx.frc, orbit_ctx.R_pos)
        k4 = derivs(state + dt*k3, xLz, orbit_ctx.frc, orbit_ctx.R_pos)
        state += (dt/6.0)*(k1 + 2k2 + 2k3 + k4)
        state = SVector(state[1], clamp(state[2],1e-6,pi-1e-6), state[3], state[4])
    end
    resize!(r, actual); resize!(vr, actual); resize!(theta, actual); resize!(vtheta, actual)
    return r, vr, theta, vtheta
end

# ============================================================
# §6  LIKELIHOOD & WEIGHT SOLVER
# ============================================================
function stellar_log_likelihood_jl(A::Matrix{Float64}, w::Vector{Float64}, verr::Vector{Float64}; rv_mask::Union{Vector{Bool},Nothing}=nothing, lambda_occ::Float64=0.0, Nocc::Int=0, eps::Float64=DEFAULT_LOG_EPS, sigma_floor_mps::Float64=DEFAULT_SIGMA_FLOOR_MPS)
    Nstar = size(A, 1) - Nocc
    p = A * w
    ll = 0.0
    if Nstar > 0
        mask = rv_mask !== nothing ? rv_mask : convert(Vector{Bool}, trues(Nstar))
        @inbounds for i in 1:Nstar
            mask[i] || continue
            sig = max(verr[i], sigma_floor_mps)
            pv  = max(p[i] * sig, eps)
            ll += log(pv)
        end
    end
    if Nocc > 0
        @inbounds for i in (Nstar+1):(Nstar+Nocc)
            ll += lambda_occ * log(max(p[i], eps))
        end
    end
    return ll
end

function stellar_log_likelihood_parts_jl(
    A::Matrix{Float64},
    w::Vector{Float64},
    verr::Vector{Float64};
    rv_mask::Union{Vector{Bool},Nothing}=nothing,
    lambda_occ::Float64=0.0,
    Nocc::Int=0,
    eps::Float64=DEFAULT_LOG_EPS,
    sigma_floor_mps::Float64=DEFAULT_SIGMA_FLOOR_MPS,
)
    Nstar = size(A, 1) - Nocc
    p = A * w

    ll_rows = zeros(Float64, max(Nstar, 0))
    ll_occ = 0.0

    if Nstar > 0
        mask = rv_mask !== nothing ? rv_mask : convert(Vector{Bool}, trues(Nstar))

        @inbounds for i in 1:Nstar
            if mask[i]
                sig = max(verr[i], sigma_floor_mps)
                pv  = max(p[i] * sig, eps)
                ll_rows[i] = log(pv)
            end
        end
    end

    if Nocc > 0
        @inbounds for i in (Nstar + 1):(Nstar + Nocc)
            ll_occ += lambda_occ * log(max(p[i], eps))
        end
    end

    return ll_rows, ll_occ
end

function solve_weights_stellar_jl(A::Matrix{Float64}, verr::Vector{Float64}; rv_mask::Union{Vector{Bool},Nothing}=nothing, Nocc::Int=0, lambda_occ::Float64=0.0, alpha::Float64=DEFAULT_ALPHA, maxiter::Int=DEFAULT_MAXITER, eps::Float64=DEFAULT_LOG_EPS, sigma_floor_mps::Float64=DEFAULT_SIGMA_FLOOR_MPS, seed::UInt=UInt(0))
    Nstar  = size(A, 1) - Nocc
    Norbit = size(A, 2)
    Norbit <= 0 && return (zeros(Float64, 0), false)
    alpha_eff = alpha * (Norbit / 200.0) * (max(Nstar, 1) / 90.0)
    rng = MersenneTwister(seed)
    w0  = rand(rng, Norbit) .+ 1e-3
    w0 ./= sum(w0)
    mask = rv_mask !== nothing ? rv_mask : trues(Nstar)
    function obj(w::Vector{Float64})
        ll = stellar_log_likelihood_jl(A, w, verr; rv_mask=mask, lambda_occ=lambda_occ, Nocc=Nocc, eps=eps, sigma_floor_mps=sigma_floor_mps)
        return -ll + alpha_eff * dot(w, w)
    end
    lo = zeros(Float64, Norbit)
    hi = fill(Inf,      Norbit)
    local result
    try
        result = Optim.optimize(obj, lo, hi, w0, Optim.Fminbox(Optim.LBFGS()),
            Optim.Options(iterations=maxiter, f_calls_limit=20000, g_tol=1e-6, show_trace=false, allow_f_increases=false))
    catch
        return (zeros(Float64, Norbit), false)
    end
    w = max.(Optim.minimizer(result), 0.0)
    s = sum(w)
    (!isfinite(s) || s <= 0.0) && return (zeros(Float64, Norbit), false)
    return (w ./ s, true)
end

mass_enclosed_two_radii(rin, rout, rho_s, r_s, MBH, ML, halo_type; stellar_model=nothing) = begin
    ctx = get_halo_context(rho_s, r_s, MBH, ML, halo_type; stellar_model=stellar_model)
    r1 = max(rin, ctx.halo[:rmin])
    r2 = max(rout, 1.001 * r1)
    fr1, _ = ctx.frc(r1, 0.0)
    fr2, _ = ctx.frc(r2, 0.0)
    (-r1 * r1 * fr1 / G, -r2 * r2 * fr2 / G)
end
# ============================================================
# Resolution / refinement support
# ============================================================

mutable struct ShellRefinementState
    chi2::Float64
    support::Float64
    reff::Float64
    score_history::Vector{Float64}
    status::Symbol
end

function effective_rank(A::AbstractMatrix; atol::Float64=1e-12)
    s = svdvals(A)
    isempty(s) && return 0.0
    smax = maximum(s)
    smax <= 0.0 && return 0.0
    sk = s[s .> (atol * smax)]
    isempty(sk) && return 0.0
    den = sum(abs2, sk)
    den <= 0.0 && return 0.0
    return (sum(sk)^2) / den
end

function shell_support(w::AbstractVector{<:Real}, cols_k::AbstractVector{<:Integer})
    isempty(cols_k) && return 0.0
    s = 0.0
    @inbounds for c in cols_k
        s += float(w[c])
    end
    return s
end

function shell_score(chi2_old::Real, chi2_new::Real, support_old::Real, support_new::Real, reff_old::Real, reff_new::Real; wchi::Float64=1.0, wsupport::Float64=1.0, wreff::Float64=1.0, eps::Float64=1e-12)
    dchi = max((float(chi2_old) - float(chi2_new)) / max(abs(float(chi2_old)), eps), 0.0)
    dsupport = abs(float(support_new) - float(support_old)) / max(float(support_old) + float(support_new), eps)
    dreff = max((float(reff_new) - float(reff_old)) / max(abs(float(reff_old)), eps), 0.0)
    support_penalty = max((float(support_old) - float(support_new)) / max(abs(float(support_old)), eps), 0.0)
    return wchi*dchi + wsupport*dsupport + wreff*dreff + support_penalty
end

function is_saturated(history::AbstractVector{<:Real}; tau::Float64=1e-3, patience::Int=3)
    n = length(history)
    n < patience && return false
    return all(@view(history[n-patience+1:n]) .< tau)
end

function update_shell_state!(st::ShellRefinementState, Ak_new::AbstractMatrix, chi2_new::Real, support_new::Real; wchi::Float64=1.0, wsupport::Float64=1.0, wreff::Float64=1.0, tau::Float64=1e-3, patience::Int=3, atol::Float64=1e-12)
    reff_new = effective_rank(Ak_new; atol=atol)
    score = shell_score(st.chi2, chi2_new, st.support, support_new, st.reff, reff_new; wchi=wchi, wsupport=wsupport, wreff=wreff)
    push!(st.score_history, score)
    st.chi2 = float(chi2_new)
    st.support = float(support_new)
    st.reff = reff_new
    st.status = is_saturated(st.score_history; tau=tau, patience=patience) ? :frozen : :active
    return score, reff_new, st.status
end


























# ################################################################################################################################################################################################################################################
# ################################################################################################################################################################################################################################################
# §DEAD  NOT ON THE HOT PATH — kept for back-compat / future
#        use.  Nothing below here is called by
#        evaluate_batch_theta or build_A_matrix_hybrid.
# ################################################################################################################################################################################################################################################
# ################################################################################################################################################################################################################################################

# --- Caches used only by the old sigma2 A-matrix path ---
const _DATA_LOCK = ReentrantLock()
const _RADIAL_EDGE_CACHE = Dict{UInt64,Vector{Float64}}()
const _ORBIT_TEMPLATE_CACHE = Dict{UInt64,Tuple{Vector{Float64},Vector{Float64}}}()
const _dbg_orbit_count=Ref(0)

reset_orbit_cache() = (lock(_HALO_LOCK); empty!(_HALO_CTX_CACHE); unlock(_HALO_LOCK); nothing)



function _orbit_template_key(shells::Vector{Float64}, Lfrac::NTuple)
    return hash((shells, Lfrac))
end
@inline function _dataset_key(r_centers_m::Vector{Float64})
    return hash(r_centers_m)
end
function get_radial_edges(r_centers_m::Vector{Float64})
    key = _dataset_key(r_centers_m)
    lock(_DATA_LOCK)
    edges = get(_RADIAL_EDGE_CACHE, key, nothing)
    unlock(_DATA_LOCK)
    edges !== nothing && return edges
    edges = Vector{Float64}(undef, length(r_centers_m)+1)
    edges[2:end-1] .= 0.5 .* (r_centers_m[1:end-1] .+ r_centers_m[2:end])
    edges[1] = 0.0
    edges[end] = Inf
    lock(_DATA_LOCK)
    _RADIAL_EDGE_CACHE[key] = edges
    unlock(_DATA_LOCK)
    return edges
end
function get_orbit_template(shells::Vector{Float64}, Lfrac::NTuple{N,Float64}) where N
    key = _orbit_template_key(shells, Lfrac)
    lock(_DATA_LOCK)
    tpl = get(_ORBIT_TEMPLATE_CACHE, key, nothing)
    unlock(_DATA_LOCK)
    tpl !== nothing && return tpl
    rapos = Float64[]
    lvals = Float64[]
    for r in shells
        rr = f64(r)
        for lf in Lfrac
            push!(rapos, rr)
            push!(lvals, f64(lf))
        end
    end
    tpl = (rapos, lvals)
    lock(_DATA_LOCK)
    _ORBIT_TEMPLATE_CACHE[key] = tpl
    unlock(_DATA_LOCK)
    return tpl
end

function orbit_to_sigma2_profile(; r_arr, th_arr, vr_arr, xLz, dt_orb, sini, r_centers_m, edges)
    nb=length(r_centers_m)
    w=zeros(nb); v=zeros(nb); v2=zeros(nb)
    _dbg_orbit_count[] += 1
    do_debug = _dbg_orbit_count[] <= 5
    hits = do_debug ? zeros(Int,nb) : nothing
    rmin_seen = do_debug ? Inf : 0.0
    rmax_seen = do_debug ? -Inf : 0.0

    sini=clamp01(f64(sini))
    phi=0.0
    @inbounds for i in eachindex(r_arr)
        rr=max(f64(r_arr[i]),1e-12)
        th=clamp(f64(th_arr[i]),1e-6,pi-1e-6)
        vr0=f64(vr_arr[i])
        ss,_=_sincos_safe(th)
        vphi=f64(xLz)/(rr*ss)
        cp,sp=cos(phi),sin(phi)
        vlos=sini*(vr0*cp - vphi*sp)
        phi+=f64(xLz)/(rr*rr*ss*ss)*f64(dt_orb)
        j=searchsortedfirst(edges,rr)-1
        (j<1 || j>nb) && continue
        w[j]+=1; v[j]+=vlos; v2[j]+=vlos^2
        if do_debug
            hits[j]+=1
            rmin_seen=min(rmin_seen,rr)
            rmax_seen=max(rmax_seen,rr)
        end
    end

    do_debug && println("BIN_HITS=",hits," rmin/max=",(rmin_seen,rmax_seen)," centers[1/end]=",(r_centers_m[1],r_centers_m[end]))
    sig2=zeros(nb)
    @inbounds for j in 1:nb
        w[j]>0 && (mv=v[j]/w[j]; sig2[j]=v2[j]/w[j]-mv^2)
    end
    sig2
end

function build_A_matrix_from_ctx(ctx::HaloContext, th0, dt, r_centers_m, valid, sini, nsteps)
    shells        = r_centers_m[valid]
    Ndat          = length(shells)
    theta0        = f64(th0[1])
    theta_hash    = hash((ctx.halo[:rho_s], ctx.halo[:r_s], ctx.halo[:MBH]))
    rng           = MersenneTwister(UInt(theta_hash))
    dt_frac_orbit = DEFAULT_DT_FRAC
    Lfrac         = DEFAULT_LFRAC
    rapos, lvals  = get_orbit_template(shells, Lfrac)
    Norbit        = length(rapos)
    A             = zeros(Float64, Ndat, Norbit)
    edges         = get_radial_edges(r_centers_m)
    orbit_ctx     = (frc = ctx.frc, R_pos = ctx.R, halo = ctx.halo)
    col           = 1
    @inbounds for idx in eachindex(rapos)
        rapo = rapos[idx]
        lf   = lvals[idx]

        if !(isfinite(rapo) && rapo > 0.0)
            col += 1
            continue
        end
        r0_frac = 0.95 + 0.04 * rand(rng)
        ic, Lz0, E0, vc, st = launch_orbit_apocenter(rapo=rapo, theta0=theta0, Lz_frac=lf, pot=ctx.pot, frc=ctx.frc, r0_frac=r0_frac, dt_frac=dt_frac_orbit)
        if st != :ok
            col += 1
            continue
        end
        r, vr, theta, vtheta = integrate_orbit_rk4(ic=ic, xLz=Lz0, orbit_ctx=orbit_ctx, nsteps=nsteps)
        isempty(r) && (col += 1; continue)
        dt_orb = f64(ic[3])
        sig2 = orbit_to_sigma2_profile(r_arr=r, th_arr=theta, vr_arr=vr, xLz=Lz0, dt_orb=dt_orb, sini=sini, r_centers_m=r_centers_m, edges=edges)
        A[:, col] .= sig2[valid]
        col += 1
    end

    return A
end

function build_A_matrix_julia(r0_unused, th0, dt, Etot_unused, xLz_unused, r_centers_m, valid, sini, nsteps, rho_s, r_s, MBH, ML, halo_type)
    build_A_matrix_julia(th0, dt, r_centers_m, valid, sini, nsteps, rho_s, r_s, MBH, ML, halo_type)
end

build_A_matrix_julia(th0, dt, r_centers_m, valid, sini, nsteps, rho_s, r_s, MBH, ML, halo_type) =
    build_A_matrix_from_ctx(get_halo_context(rho_s, r_s, MBH, ML, halo_type), th0, dt, r_centers_m, valid, sini, nsteps)

ospm_runcheck(theta, args...) =
    build_A_matrix_julia(args..., theta[1], theta[2], length(theta)>2 ? theta[3] : 0.0, length(theta)>3 ? theta[4] : 1.0)

# Back-compat wrapper: treats all stars as having vlos.
function build_A_matrix_stellar(Norbit::Int, R_star_m::Vector{Float64}, v_star_mps::Vector{Float64}, verr_star_mps::Vector{Float64}, sini::Float64, rho_s::Float64, r_s::Float64, MBH::Float64, ML::Float64, halo_type::String; stellar_model=nothing, nsteps::Int=DEFAULT_NSTEPS, Lfrac::NTuple{5,Float64}=DEFAULT_LFRAC, dt_frac_orbit::Float64=DEFAULT_DT_FRAC, dR_frac::Float64=DEFAULT_DR_FRAC, dR_floor_frac::Float64=DEFAULT_DR_FLOOR_FRAC, dR_floor_pc::Float64=DEFAULT_DR_FLOOR_PC, Nbins_occ::Int=DEFAULT_NBINS_OCC, return_occ::Bool=true, max_attempts_factor::Int=DEFAULT_MAX_ATTEMPTS, diag::Bool=false)
    has_vlos = trues(length(R_star_m))
    build_A_matrix_hybrid(Norbit, R_star_m, has_vlos, v_star_mps, verr_star_mps, sini, rho_s, r_s, MBH, ML, halo_type; stellar_model=stellar_model, nsteps=nsteps, Lfrac=Lfrac, dt_frac_orbit=dt_frac_orbit, dR_frac=dR_frac, Nbins_occ=Nbins_occ, return_occ=return_occ, max_attempts_factor=max_attempts_factor, diag=diag)
end

# --- WLS / NNLS (not used by current pipeline) ---
const _TARGET_LOCK = ReentrantLock()
const _TARGET_D  = Ref{Vector{Float64}}(Float64[])
const _TARGET_W2 = Ref{Vector{Float64}}(Float64[])   # w2 = 1/sigma^2
function set_target_wls!(d::AbstractVector{<:Real}, sigma::AbstractVector{<:Real})
    dd = Float64.(d)
    ss = Float64.(sigma)
    length(dd) == length(ss) || error("set_target_wls!: d and sigma must match")
    any(!isfinite, dd) && error("set_target_wls!: non-finite d")
    any(x->(!(isfinite(x) && x>0.0)), ss) && error("set_target_wls!: sigma must be finite and > 0")
    w2 = similar(ss)
    @inbounds for i in eachindex(ss)
        w2[i] = 1.0 / (ss[i]*ss[i])
    end
    lock(_TARGET_LOCK)
    _TARGET_D[]  = dd
    _TARGET_W2[] = w2
    unlock(_TARGET_LOCK)
    return nothing
end
function _get_target_wls()
    lock(_TARGET_LOCK)
    d  = _TARGET_D[]
    w2 = _TARGET_W2[]
    unlock(_TARGET_LOCK)
    isempty(d) && error("WLS target not initialized. Call set_target_wls! first.")
    return d, w2
end
@inline function _proj_nn!(x::Vector{Float64})
    @inbounds for i in eachindex(x)
        x[i] = x[i] < 0.0 ? 0.0 : x[i]
    end
    return x
end
function _lipschitz_wls(A::Matrix{Float64}, w2::Vector{Float64}; iters::Int=12)
    m, n = size(A)
    v  = fill(1.0/sqrt(n), n)
    Av = zeros(Float64, m)
    t  = zeros(Float64, m)
    g  = zeros(Float64, n)
    for _ in 1:iters
        mul!(Av, A, v)
        @inbounds for i in 1:m
            t[i] = w2[i] * Av[i]
        end
        mul!(g, transpose(A), t)
        ng = norm(g)
        ng > 0.0 || break
        @inbounds for j in 1:n
            v[j] = g[j] / ng
        end
    end
    mul!(Av, A, v)
    @inbounds for i in 1:m
        t[i] = w2[i] * Av[i]
    end
    mul!(g, transpose(A), t)
    L = dot(v, g)
    return max(L, 1e-30)
end
function _nnls_wls_pg!(A::Matrix{Float64}, d::Vector{Float64}, w2::Vector{Float64};
    max_iter::Int=400, tol::Float64=1e-8)
    m, n = size(A)
    x    = zeros(Float64, n)
    xnew = similar(x)
    g    = similar(x)
    Ax   = zeros(Float64, m)
    r    = zeros(Float64, m)
    t    = zeros(Float64, m)
    L = _lipschitz_wls(A, w2)
    α = 1.0 / L
    for _ in 1:max_iter
        mul!(Ax, A, x)
        @inbounds for i in 1:m
            r[i] = Ax[i] - d[i]
            t[i] = w2[i] * r[i]
        end
        mul!(g, transpose(A), t)
        @inbounds for j in 1:n
            xnew[j] = x[j] - α * g[j]
        end
        _proj_nn!(xnew)
        if norm(xnew .- x) <= tol * max(1.0, norm(x))
            x .= xnew
            break
        end
        x .= xnew
    end
    return x
end
