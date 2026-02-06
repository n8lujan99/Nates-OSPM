module OSPMPhysicsSpherical
@info "OSPMPhysicsSpherical loaded from" @__FILE__

using LinearAlgebra, StaticArrays, Statistics, Random, Base.Threads

export build_R_halo_physical, rho_interp, halo_from_theta, tables_spherical,
       make_potential_force_funcs, integrate_orbit_rk4, orbit_to_sigma2_profile,
       build_A_matrix_julia, ospm_runcheck, build_A_matrix_stellar, build_A_matrix_hybrid,
       mass_enclosed_two_radii, reset_orbit_cache, evaluate_batch_theta, NTHREADS

GC.enable(true)

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
const _dbg_orbit_count=Ref(0)

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

const _HALO_CTX_CACHE = Dict{Tuple{Float64,Float64,Float64,Symbol,Int,Float64},HaloContext}()
const _HALO_LOCK = ReentrantLock()

reset_orbit_cache() = (lock(_HALO_LOCK); empty!(_HALO_CTX_CACHE); unlock(_HALO_LOCK); nothing)

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

###########################################
# Dark matter halo part that needs to be extended later to test multiple halo types
#########################
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


function halo_from_theta(rho_s, r_s, MBH; halo_type="nfw", alpha=nothing)
    ht = Symbol(lowercase(String(halo_type)))
    h = Dict(
        :rho_s => f64(rho_s) * Msun / pc^3,
        :r_s   => f64(r_s)   * pc,
        :rs    => f64(r_s)   * pc,
        :MBH   => f64(MBH)   * Msun,
        :type  => ht,
        :rmin  => 1e-6 * f64(r_s) * pc
    )
    # Optional shape parameter, only used for alt halo types like Einasto. Default value is 0.18 if not provided.
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

    # Outer-shell term for potential:
    # Phi(r) = -G[ M(<r)/r + ∫_r^∞ 4pi r' rho(r') dr' ]
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
    halo=normalize_halo(halo); MBH=f64(halo[:MBH]); rmin=f64(halo[:rmin])
    rlgmin=log10(f64(R[1])); rlgmax=log10(f64(R[end])); np=length(R)
    rlgmax>rlgmin || error("Degenerate R grid")

    @inline function interp(arr, rr)
        r=max(f64(rr),rmin)
        lr=log10(r)
        x=(lr-rlgmin)*(np-1)/(rlgmax-rlgmin)
        x=clamp(x,0.0,np-1.0)
        i0=Int(floor(x))+1
        i1=min(i0+1,np)
        t=x-(i0-1)
        (1-t)*arr[i0] + t*arr[i1]
    end

    pot(r,mu=0.0)=begin
        rr=max(abs(f64(r)),rmin)
        Ph=interp(tabv,rr)
        Pbh=MBH>0 ? (-G*MBH/rr) : 0.0
        Ph+Pbh
    end

    frc(r,mu=0.0)=begin
        rr=max(abs(f64(r)),rmin)
        frh=interp(tabfr,rr)
        frbh=MBH>0 ? (-G*MBH/(rr*rr)) : 0.0
        frh+frbh, 0.0
    end

    pot, frc, R
end

function build_halo_context(rho_s, r_s, MBH, halo_type; nR=256, rmax_factor=300.0)
    halo=halo_from_theta(rho_s,r_s,MBH; halo_type=halo_type)
    R=build_R_halo_physical(nR; rmin=halo[:rmin], rmax=rmax_factor*halo[:rs])
    tabv,tabfr,Menc=tables_spherical(R,1,halo,rho_interp)
    pot,frc,_=make_potential_force_funcs(halo,R,1,tabv,tabfr,Menc)
    HaloContext(halo,f64.(R),tabv,tabfr,Menc,pot,frc)
end

function get_halo_context(rho_s, r_s, MBH, halo_type; nR=256, rmax_factor=300.0)
    ht=Symbol(lowercase(String(halo_type)))
    key=(f64(rho_s),f64(r_s),f64(MBH),ht,nR,f64(rmax_factor))

    lock(_HALO_LOCK); ctx=get(_HALO_CTX_CACHE,key,nothing); unlock(_HALO_LOCK)
    ctx!==nothing && return ctx

    newctx=build_halo_context(rho_s,r_s,MBH,ht;nR=nR,rmax_factor=rmax_factor)

    lock(_HALO_LOCK)
    ctx=get(_HALO_CTX_CACHE,key,nothing)
    if ctx===nothing
        _HALO_CTX_CACHE[key]=newctx
        ctx=newctx
    end
    unlock(_HALO_LOCK)

    ctx
end

mass_enclosed_two_radii(rin,rout,rho_s,r_s,MBH,halo_type)=begin
    ctx=get_halo_context(rho_s,r_s,MBH,halo_type)
    r1=max(rin,ctx.halo[:rmin]); r2=max(rout,1.001*r1)
    fr1,_=ctx.frc(r1,0.0); fr2,_=ctx.frc(r2,0.0)
    (-r1*r1*fr1/G, -r2*r2*fr2/G)
end

function build_A_matrix_julia(r0_unused, th0, dt, Etot_unused, xLz_unused, r_centers_m, valid, sini, nsteps, rho_s, r_s, MBH, halo_type)
    build_A_matrix_julia(th0, dt, r_centers_m, valid, sini, nsteps, rho_s, r_s, MBH, halo_type)
end

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

function integrate_orbit_rk4(; ic, xLz, ctx, nsteps=4000, stop_rmin_factor=1.001)
    halo=ctx.halo
    rmin_stop=stop_rmin_factor*f64(halo[:rmin])
    r0=f64(ic[1]); theta0=f64(ic[2]); dt=f64(ic[3])
    vr0=length(ic)>=4 ? f64(ic[4]) : 0.0
    vtheta0=length(ic)>=5 ? f64(ic[5]) : 0.0
    state=SVector(r0,theta0,vr0,vtheta0)
    r=Float64[]; vr=Float64[]; theta=Float64[]
    rmax_stop=10.0*f64(ctx.R_pos[end])

    @inbounds for step in 1:Int(nsteps)
        !all(isfinite,state) && break
        rr=state[1]; tr=state[2]
        (rr<=rmin_stop || rr>=rmax_stop || abs(tr)>1e6) && break
        push!(r,rr); push!(vr,state[3]); push!(theta,tr)
        k1=derivs(state,xLz,ctx.frc,ctx.R_pos)
        k2=derivs(state+0.5*dt*k1,xLz,ctx.frc,ctx.R_pos)
        k3=derivs(state+0.5*dt*k2,xLz,ctx.frc,ctx.R_pos)
        k4=derivs(state+dt*k3,xLz,ctx.frc,ctx.R_pos)
        state += (dt/6.0)*(k1+2k2+2k3+k4)
        state = SVector(state[1], clamp(state[2],1e-6,pi-1e-6), state[3], state[4])
    end
    r,vr,theta
end


function orbit_to_sigma2_profile(; r_arr, th_arr, vr_arr, xLz, r_centers_m, edges, sini)
    nb=length(r_centers_m)
    w=zeros(nb); v=zeros(nb); v2=zeros(nb)
    _dbg_orbit_count[] += 1
    do_debug = _dbg_orbit_count[] <= 5
    hits = do_debug ? zeros(Int,nb) : nothing
    rmin_seen = do_debug ? Inf : 0.0
    rmax_seen = do_debug ? -Inf : 0.0

    sini=f64(sini); sini=clamp01(sini)
    cosi=sqrt(max(1.0 - sini*sini, 0.0))

    @inbounds for i in eachindex(r_arr)
        rr=max(f64(r_arr[i]),1e-12)
        th=clamp(f64(th_arr[i]),1e-6,pi-1e-6)
        vr0=f64(vr_arr[i])
        ss,_=_sincos_safe(th)
        vphi=f64(xLz)/(rr*ss)
        vlos=sini*vr0 + cosi*vphi
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

function launch_orbit_apocenter(; rapo::Float64, theta0::Float64, Lz_frac::Float64,
    pot, frc, r0_frac::Float64=0.98, dt_frac::Float64=0.01, dt_floor::Float64=1e-30,
    debug::Bool=true)
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
    # inward is negative; allow a small ambiguous band
    if frs > fr_tol
        return debug ?
            ((rapo, theta0, ss, frs, fr_tol, fr_scale), 0.0, 0.0, 0.0, :reject_force) :
            (nothing, 0.0, 0.0, 0.0, :reject_force)
    end
    # Option 2: regularize vc if force is mushy
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


# ======================================================================
# HYBRID A-MATRIX:
#   - rows 1:Nvlos   : velocity likelihood for stars with has_vlos=true
#   - rows Nvlos+1:end : occupancy constraint (all stars), optional
# ======================================================================
function build_A_matrix_hybrid(
        Norbit::Int,
        R_star_m::Vector{Float64},
        has_vlos::AbstractVector{Bool},
        v_star_mps::Vector{Float64},
        verr_star_mps::Vector{Float64},
        sini::Float64, rho_s::Float64, r_s::Float64, MBH::Float64, halo_type::String;
        nsteps::Int=4000, Lfrac::NTuple{5,Float64}=(0.05,0.2,0.4,0.7,1.0), dt_frac_orbit::Float64=0.01,
        dR_frac::Float64=0.05, dR_floor_frac::Float64=0.01, dR_floor_pc::Float64=0.0,
        Nbins_occ::Int=6, return_occ::Bool=true, max_attempts_factor::Int=60,
        diag::Bool=false)

    GC.enable(true)

    Nstar=length(R_star_m)
    @assert length(has_vlos)==Nstar
    @assert length(v_star_mps)==Nstar
    @assert length(verr_star_mps)==Nstar
    Nstar==0 && return zeros(Float64,0,Norbit)

    # index map of vlos stars
    vlos_idx = Int[]
    sizehint!(vlos_idx, Nstar)
    @inbounds for i in 1:Nstar
        has_vlos[i] && push!(vlos_idx, i)
    end
    Nvlos = length(vlos_idx)

    ctx=get_halo_context(rho_s,r_s,MBH,halo_type)

    sini=clamp01(f64(sini))
    cosi=sqrt(max(1.0 - sini*sini, 0.0))

    Rmin=minimum(R_star_m); Rmax=maximum(R_star_m)
    !(isfinite(Rmin)&&isfinite(Rmax)&&Rmax>Rmin) && return zeros(Float64, Nvlos + (return_occ ? Nbins_occ : 0), Norbit)

    shells=sort(copy(R_star_m))
    dnn=Float64[]
    if length(shells)>=3
        resize!(dnn,length(shells)-1)
        @inbounds for i in 1:length(shells)-1
            dnn[i]=shells[i+1]-shells[i]
        end
    end
    med_dnn=(length(dnn)>0) ? median(abs.(dnn)) : 0.0

    auto_floor=max(dR_floor_frac*max(med_dnn,0.0), 0.02*Rmin)
    abs_floor=(dR_floor_pc>0) ? (dR_floor_pc*pc) : 0.0
    dR_floor=max(auto_floor,abs_floor,1e-12)

    occ_edges=collect(range(Rmin,Rmax; length=Nbins_occ+1))

    A_vlos = zeros(Float64, Nvlos, Norbit)
    A_occ  = zeros(Float64, Nbins_occ, Norbit)

    orbit_ctx=(frc=ctx.frc, R_pos=ctx.R, halo=ctx.halo)
    rng=MersenneTwister(0x5eed1234)
    theta0=f64(pi/2)

    col=1; idx=1; max_attempts=max_attempts_factor*Norbit; attempts=0
    while col<=Norbit && attempts<max_attempts
        attempts+=1
        rapo=f64(shells[idx])
        idx=(idx==length(shells)) ? 1 : (idx+1)
        !(isfinite(rapo)&&rapo>0) && continue

        lf=Lfrac[1 + (attempts % length(Lfrac))]
        r0_frac=0.95 + 0.04*rand(rng)

        ic,Lz0,E0,vc,st=launch_orbit_apocenter(rapo=rapo,theta0=theta0,Lz_frac=f64(lf),
                                              pot=ctx.pot,frc=ctx.frc,r0_frac=r0_frac,dt_frac=dt_frac_orbit)
        st!=:ok && continue

        r,vr,theta=integrate_orbit_rk4(ic=ic,xLz=Lz0,ctx=orbit_ctx,nsteps=nsteps)
        isempty(r) && continue

        s_arr=similar(r); vphi_arr=similar(r)
        @inbounds for i in eachindex(r)
            si=_ssin(f64(theta[i])); ri=f64(r[i])
            s_arr[i]=ri*si
            vphi_arr[i]=f64(Lz0)/(ri*si)
        end

        Nhits=length(s_arr)

        if return_occ && Nhits>0
            @inbounds for j in 1:Nbins_occ
                lo=occ_edges[j]; hi=occ_edges[j+1]; cnt=0
                @inbounds for x in s_arr
                    cnt += (x>=lo && x<hi) ? 1 : 0
                end
                A_occ[j,col]=cnt/Nhits
            end
        end

        vlos_arr = @. sini*f64(vr) + cosi*vphi_arr

        # sort by s for window queries
        pidx=sortperm(s_arr)
        s_sorted=s_arr[pidx]
        vlos_sorted=vlos_arr[pidx]
        Nhits=length(s_sorted)

        inv_sqrt2pi = inv(sqrt(2*pi))

        if Nvlos>0
            @inbounds for row in 1:Nvlos
                istar = vlos_idx[row]
                Ri=f64(R_star_m[istar])
                vi=f64(v_star_mps[istar])
                si=f64(verr_star_mps[istar])
                s2=si*si
                if !(isfinite(s2) && s2>0)
                    A_vlos[row,col]=0.0
                    continue
                end
                dR=max(dR_frac*Ri,dR_floor)
                lo=Ri-dR; hi=Ri+dR
                j0=searchsortedfirst(s_sorted,lo)
                j1=searchsortedlast(s_sorted,hi)
                if j1 < j0
                    A_vlos[row,col]=0.0
                    continue
                end
                p=0.0; nhit=j1-j0+1
                inv_norm = inv_sqrt2pi/si
                @inbounds for j in j0:j1
                    dv=vi - vlos_sorted[j]
                    p += exp(-0.5*(dv*dv)/s2)
                end
                A_vlos[row,col] = (p*inv_norm)/nhit
            end
        end

        col += 1
    end

    col<=Norbit && println("WARNING: build_A_matrix_hybrid filled ",col-1," / ",Norbit," after attempts=",attempts)

    A = return_occ ? vcat(A_vlos, A_occ) : A_vlos
    GC.gc()
    diag ? (A, Dict{String,Any}("filled"=>col-1,"attempts"=>attempts,"Nvlos"=>Nvlos,"Nbins_occ"=>(return_occ ? Nbins_occ : 0))) : A
end

# Back-compat: old name. Treats all stars as having vlos.
function build_A_matrix_stellar(Norbit::Int, R_star_m::Vector{Float64}, v_star_mps::Vector{Float64}, verr_star_mps::Vector{Float64},
        sini::Float64, rho_s::Float64, r_s::Float64, MBH::Float64, halo_type::String;
        nsteps::Int=4000, Lfrac::NTuple{5,Float64}=(0.05,0.2,0.4,0.7,1.0), dt_frac_orbit::Float64=0.01,
        dR_frac::Float64=0.05, dR_floor_frac::Float64=0.01, dR_floor_pc::Float64=0.0,
        Nbins_occ::Int=6, return_occ::Bool=true, max_attempts_factor::Int=60,
        diag::Bool=false)

    has_vlos = trues(length(R_star_m))
    build_A_matrix_hybrid(Norbit, R_star_m, has_vlos, v_star_mps, verr_star_mps,
                          sini, rho_s, r_s, MBH, halo_type;
                          nsteps=nsteps, Lfrac=Lfrac, dt_frac_orbit=dt_frac_orbit,
                          dR_frac=dR_frac, dR_floor_frac=dR_floor_frac, dR_floor_pc=dR_floor_pc,
                          Nbins_occ=Nbins_occ, return_occ=return_occ, max_attempts_factor=max_attempts_factor,
                          diag=diag)
end

function build_A_matrix_from_ctx(ctx::HaloContext, th0, dt, r_centers_m, valid, sini, nsteps)
    shells=r_centers_m[valid]; Ndat=length(shells); theta0=f64(th0[1]); rng=MersenneTwister(0xdeadbeef)
    dt_frac_orbit=0.01; Lfrac=(0.05,0.2,0.4,0.7,1.0); Norbit=length(shells)*length(Lfrac)
    A=zeros(Float64,Ndat,Norbit)

    edges=similar(r_centers_m,length(r_centers_m)+1)
    edges[2:end-1].=0.5.*(r_centers_m[1:end-1].+r_centers_m[2:end])
    edges[1]=0.0; edges[end]=Inf

    kept=0; rej_turn=0; rej_energy=0; rej_force=0; dErel_max=0.0; dErel_bad=0; col=1
    orbit_ctx=(frc=ctx.frc, R_pos=ctx.R, halo=ctx.halo)

    @inbounds for rapo0 in shells
        rapo=f64(rapo0)
        if !(isfinite(rapo)&&rapo>0); col += length(Lfrac); continue; end

        for lf in Lfrac
            r0_frac=0.95 + 0.04*rand(rng)
            ic,Lz0,E0,vc,st=launch_orbit_apocenter(rapo=rapo,theta0=theta0,Lz_frac=f64(lf),
                                                  pot=ctx.pot,frc=ctx.frc,r0_frac=r0_frac,dt_frac=dt_frac_orbit)

            if st!=:ok
                st==:reject_turning ? (rej_turn+=1) : (st==:reject_force || st==:reject_vc ? (rej_force+=1) : (rej_energy+=1))
                col += 1; continue
            end

            r,vr,theta=integrate_orbit_rk4(ic=ic,xLz=Lz0,ctx=orbit_ctx,nsteps=nsteps)
            isempty(r) && (rej_energy+=1; col+=1; continue)

            sig2=orbit_to_sigma2_profile(r_arr=r, th_arr=theta, vr_arr=vr, xLz=Lz0, r_centers_m=r_centers_m, edges=edges, sini=sini)
            A[:,col] .= sig2[valid]
            kept += 1

            if kept <= 8
                Emin=Inf; Emax=-Inf; Esum=0.0; nE=0
                @inbounds for i in eachindex(r)
                    si=_ssin(theta[i])
                    Ei=0.5*(vr[i]^2 + (Lz0/(r[i]*si))^2) + ctx.pot(r[i],si)
                    if isfinite(Ei); Emin=min(Emin,Ei); Emax=max(Emax,Ei); Esum+=Ei; nE+=1; end
                end
                if nE>0
                    Emean=Esum/nE
                    dErel=(Emax-Emin)/max(abs(Emean),1e-300)
                    dErel_max=max(dErel_max,dErel)
                    dErel_bad += (dErel>1e-3) ? 1 : 0
                    println("dE/<E>=",dErel," rapo=",rapo," Lfrac=",lf," n=",length(r))
                end
            end

            col += 1
        end
    end

    println("orbit columns kept = ",kept," / ",Norbit," | reject_turning=",rej_turn," reject_energy=",rej_energy," reject_force=",rej_force,
            " | dErel_max(first few)=",dErel_max," dErel_bad(first few)=",dErel_bad)
    A
end

build_A_matrix_julia(th0, dt, r_centers_m, valid, sini, nsteps, rho_s, r_s, MBH, halo_type) =
    build_A_matrix_from_ctx(get_halo_context(rho_s,r_s,MBH,halo_type), th0, dt, r_centers_m, valid, sini, nsteps)

ospm_runcheck(theta, args...) =
    build_A_matrix_julia(args..., theta[1], theta[2], length(theta)>2 ? theta[3] : 0.0)

function evaluate_batch_theta(
    thetas::AbstractMatrix{<:Real},
    r_centers_m::Vector{Float64},
    valid::Vector{Bool},
    sini::Float64,
    nsteps::Int,
    halo_type::String
)
    nrow, ncol = size(thetas)
    dim, nbatch = nrow <= ncol ? (nrow, ncol) : (ncol, nrow)

    getth(i) = nrow <= ncol ?
        (@inbounds ntuple(j->Float64(thetas[j,i]), dim)) :
        (@inbounds ntuple(j->Float64(thetas[i,j]), dim))

    status = Vector{Int}(undef, nbatch)
    chi2   = Vector{Float64}(undef, nbatch)

    Threads.@threads for i in 1:nbatch
        th = getth(i)
        try
            rho_s = th[1]; r_s = th[2]; MBH = dim>=3 ? th[3] : 0.0
            A = build_A_matrix_julia((pi/2,), 0.0, r_centers_m, valid, sini, nsteps, rho_s, r_s, MBH, halo_type)
            v = sum(abs2, A)
            if isfinite(v)
                status[i]=0
                chi2[i]=v
            else
                status[i]=2
                chi2[i]=Inf
            end
        catch
            status[i]=3
            chi2[i]=Inf
        end
    end

    status, chi2
end

end
