module OSPMPhysicsSpherical

using LinearAlgebra, StaticArrays, Statistics, Random, Base.Threads

export build_R_halo_physical, rho_interp, halo_from_theta, tables_spherical,
       make_potential_force_funcs, integrate_orbit_rk4, orbit_to_sigma2_profile,
       build_A_matrix_julia, ospm_runcheck, build_A_matrix_stellar,
       mass_enclosed_two_radii, reset_orbit_cache, evaluate_batch_theta, NTHREADS

GC.enable(true)

const NTHREADS = Threads.nthreads()
const G    = 6.67430e-11
const c    = 2.99792458e8
const pc   = 3.0856775814913673e16
const Msun = 1.98847e30

@inline f64(x)=Float64(x)
@inline safe_sign(x)=x>0 ? 1.0 : (x<0 ? -1.0 : 0.0)
@inline _ssin(θ::Float64)=begin s=sin(θ); abs(s)>1e-12 ? s : safe_sign(s)*1e-12 end
@inline function _sincos_safe(θ::Float64); s,c=sincos(θ); abs(s)>1e-12 ? (s,c) : (safe_sign(s)*1e-12,c) end

# ---- small HPC guard helpers ----

@inline clamp01(x::Float64) = x < 0 ? 0.0 : (x > 1 ? 1.0 : x)

# Sum of squares that does not hide NaNs and does not overflow silently if you can avoid it.
# Returns (isfinite::Bool, value::Float64)
function safe_sum_abs2(A::AbstractArray{<:Real})
    s = 0.0
    @inbounds for x in A
        fx = f64(x)
        if !isfinite(fx)
            return (false, Inf)
        end
        # scale-safe accumulation would be nicer, but this is already better than silent Inf.
        s += fx*fx
        if !isfinite(s)
            return (false, Inf)
        end
    end
    return (true, s)
end

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

function rho_interp(rv, halo)
    r=abs(rv[1]); ρs=halo[:rho_s]; rs=halo[:r_s]; x=r/max(rs,1e-30)
    halo[:type]===:nfw   && return ρs/(x*(1+x)^2 + 1e-30)
    halo[:type]===:cored && return ρs/((1+x)*(1+x^2) + 1e-30)
    error("Unknown halo type")
end

function halo_from_theta(rho_s, r_s, MBH; halo_type="nfw")
    ht=Symbol(lowercase(String(halo_type)))
    Dict(
        :rho_s => f64(rho_s)*Msun/pc^3,
        :r_s   => f64(r_s)*pc,
        :rs    => f64(r_s)*pc,
        :MBH   => f64(MBH)*Msun,
        :type  => ht,
        :rmin  => 1e-6*f64(r_s)*pc
    )
end

function tables_spherical(R, nlegup, halo, ρfn)
    halo=normalize_halo(halo); n=length(R)
    ρ=similar(R); tabv=zeros(n); tabfr=zeros(n); Menc=zeros(n)

    @inbounds for i in eachindex(R)
        v=ρfn((R[i],0.0),halo)
        ρ[i]=isfinite(v) ? v : 0.0
    end

    @inbounds for i in 2:n
        dr=R[i]-R[i-1]
        Menc[i]=Menc[i-1]+0.5*dr*(R[i]^2*ρ[i]+R[i-1]^2*ρ[i-1])
    end
    Menc .*= 4π

    J=zeros(n)
    @inbounds for i in (n-1):-1:1
        dr=R[i+1]-R[i]
        J[i]=J[i+1]+0.5*dr*(R[i+1]*ρ[i+1] + R[i]*ρ[i])
    end
    J .*= 4π

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

    pot(r,μ=0.0)=begin
        rr=max(abs(f64(r)),rmin)
        Φh=interp(tabv,rr)
        Φbh=MBH>0 ? (-G*MBH/rr) : 0.0
        Φh+Φbh
    end

    frc(r,μ=0.0)=begin
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
    r,θ,vr,vθ=s
    !(isfinite(r)&&isfinite(θ)&&isfinite(vr)&&isfinite(vθ)) && return SVector(0.0,0.0,0.0,0.0)
    r_safe=max(abs(r),1e-12)
    sθ,cθ=_sincos_safe(θ)
    r_tab=clamp(r_safe,R[1],R[end])
    fr,_=frc(r_tab,sθ)
    !isfinite(fr) && return SVector(0.0,0.0,0.0,0.0)
    dr=vr
    dθ=vθ/r_safe
    dvr=(vθ*vθ)/r_safe + (Lz*Lz)/(r_safe^3*sθ*sθ) + fr
    dvθ=(Lz*Lz)*cθ/(r_safe^3*sθ^3) - (vr*vθ)/r_safe
    SVector(dr,dθ,dvr,dvθ)
end

function integrate_orbit_rk4(; ic, xLz, ctx, nsteps=4000, stop_rmin_factor=1.001)
    # ctx must be a named tuple with fields: frc, R_pos, halo
    @assert hasproperty(ctx, :halo)
    @assert hasproperty(ctx, :frc)
    @assert hasproperty(ctx, :R_pos)

    halo=ctx.halo
    rmin_stop=stop_rmin_factor*f64(halo[:rmin])

    r0=f64(ic[1]); θ0=f64(ic[2]); dt=f64(ic[3])
    dt > 0.0 || return (Float64[], Float64[], Float64[])  # never integrate with dt<=0

    vr0=length(ic)>=4 ? f64(ic[4]) : 0.0
    vθ0=length(ic)>=5 ? f64(ic[5]) : 0.0
    state=SVector(r0,θ0,vr0,vθ0)

    r=Float64[]; vr=Float64[]; θ=Float64[]
    rmax_stop=10.0*f64(ctx.R_pos[end])

    @inbounds for step in 1:Int(nsteps)
        !all(isfinite,state) && break
        rr=state[1]; θr=state[2]
        (rr<=rmin_stop || rr>=rmax_stop || abs(θr)>1e6) && break

        push!(r,rr); push!(vr,state[3]); push!(θ,θr)

        k1=derivs(state,xLz,ctx.frc,ctx.R_pos)
        k2=derivs(state+0.5*dt*k1,xLz,ctx.frc,ctx.R_pos)
        k3=derivs(state+0.5*dt*k2,xLz,ctx.frc,ctx.R_pos)
        k4=derivs(state+dt*k3,xLz,ctx.frc,ctx.R_pos)

        state += (dt/6.0)*(k1+2k2+2k3+k4)
        state = SVector(state[1], clamp(state[2],1e-6,π-1e-6), state[3], state[4])
    end

    r,vr,θ
end

# thread-safe debug counter
const _dbg_orbit_count = Atomic{Int}(0)

function orbit_to_sigma2_profile(; r_arr, th_arr, vr_arr, xLz, r_centers_m, edges, sini)
    nb=length(r_centers_m)
    w=zeros(nb); v=zeros(nb); v2=zeros(nb)
    k = atomic_add!(_dbg_orbit_count, 1) + 1
    do_debug = k <= 5

    sini = clamp01(f64(sini))
    cosi = sqrt(max(1.0 - sini*sini, 0.0))

    hits = do_debug ? zeros(Int,nb) : nothing
    rmin_seen = do_debug ? Inf : 0.0
    rmax_seen = do_debug ? -Inf : 0.0

    @inbounds for i in eachindex(r_arr)
        rr=max(f64(r_arr[i]),1e-12)
        th=clamp(f64(th_arr[i]),1e-6,π-1e-6)
        vr0=f64(vr_arr[i])
        ss,_=_sincos_safe(th)
        vφ=f64(xLz)/(rr*ss)
        vlos=sini*vr0 + cosi*vφ

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

    σ2=zeros(nb)
    @inbounds for j in 1:nb
        w[j]>0 && (μv=v[j]/w[j]; σ2[j]=v2[j]/w[j]-μv^2)
    end
    σ2
end

# ... keep your launch_orbit_apocenter and build_A_matrix_stellar as-is ...

function build_A_matrix_from_ctx(ctx::HaloContext, th0, dt, r_centers_m, valid, sini, nsteps)
    # dt is not used here; keep it for API compatibility but never allow nonsense to leak
    shells=r_centers_m[valid]; Ndat=length(shells)
    θ0=f64(th0[1])
    rng=MersenneTwister(0xdeadbeef)

    sini = clamp01(f64(sini))

    dt_frac_orbit=0.01
    Lfrac=(0.05,0.2,0.4,0.7,1.0)
    Norbit=length(shells)*length(Lfrac)
    A=zeros(Float64,Ndat,Norbit)

    edges=similar(r_centers_m,length(r_centers_m)+1)
    edges[2:end-1].=0.5.*(r_centers_m[1:end-1].+r_centers_m[2:end])
    edges[1]=0.0; edges[end]=Inf

    kept=0; rej_turn=0; rej_energy=0; rej_force=0; col=1
    orbit_ctx=(frc=ctx.frc, R_pos=ctx.R, halo=ctx.halo)

    @inbounds for rapo0 in shells
        rapo=f64(rapo0)
        if !(isfinite(rapo)&&rapo>0)
            col += length(Lfrac)
            continue
        end

        for lf in Lfrac
            r0_frac=0.95 + 0.04*rand(rng)
            ic,Lz0,E0,vc,st=launch_orbit_apocenter(rapo=rapo,θ0=θ0,Lz_frac=f64(lf),pot=ctx.pot,frc=ctx.frc,r0_frac=r0_frac,dt_frac=dt_frac_orbit)

            if st!=:ok
                st==:reject_turning ? (rej_turn+=1) : (st==:reject_force || st==:reject_vc ? (rej_force+=1) : (rej_energy+=1))
                col += 1
                continue
            end

            r,vr,θ=integrate_orbit_rk4(ic=ic,xLz=Lz0,ctx=orbit_ctx,nsteps=nsteps)
            if isempty(r)
                rej_energy += 1
                col += 1
                continue
            end

            σ2=orbit_to_sigma2_profile(r_arr=r, th_arr=θ, vr_arr=vr, xLz=Lz0, r_centers_m=r_centers_m, edges=edges, sini=sini)
            A[:,col] .= σ2[valid]
            kept += 1
            col += 1
        end
    end

    println("orbit columns kept = ",kept," / ",Norbit," | reject_turning=",rej_turn," reject_energy=",rej_energy," reject_force=",rej_force)
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

    getθ(i) = nrow <= ncol ?
        (@inbounds ntuple(j->Float64(thetas[j,i]), dim)) :
        (@inbounds ntuple(j->Float64(thetas[i,j]), dim))

    status = fill(3, nbatch)          # default = exception
    chi2   = fill(Inf, nbatch)

    sini = clamp01(f64(sini))

    Threads.@threads for i in 1:nbatch
        θ = getθ(i)
        try
            rho_s = θ[1]; r_s = θ[2]; MBH = dim ≥ 3 ? θ[3] : 0.0

            # Never pass dt=0.0 to anything orbit-like.
            dt_dummy = 1.0

            A = build_A_matrix_julia((π/2,), dt_dummy, r_centers_m, valid, sini, nsteps, rho_s, r_s, MBH, halo_type)

            ok, v = safe_sum_abs2(A)
            if ok
                status[i] = 0
                chi2[i]   = v
            else
                status[i] = 2
                chi2[i]   = Inf
            end
        catch e
            status[i] = 3
            chi2[i]   = Inf
            @warn "evaluate_batch_theta failed" i=i θ=θ err=e
        end
    end

    status, chi2
end

end
