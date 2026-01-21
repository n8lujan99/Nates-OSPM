module OSPMPhysicsSpherical
using LinearAlgebra, StaticArrays, Statistics, Random, Base.Threads

export build_R_halo_physical, rho_interp, halo_from_theta, tables_spherical, make_potential_force_funcs,
       integrate_orbit_rk4, orbit_to_sigma2_profile, build_A_matrix_julia, ospm_runcheck, build_A_matrix_stellar,
       mass_enclosed_two_radii, reset_orbit_cache, evaluate_batch_theta, NTHREADS

const NTHREADS = Threads.nthreads()
const G=6.67430e-11; const c=2.99792458e8; const pc=3.0856775814913673e16; const Msun=1.98847e30
@inline f64(x)=Float64(x)
@inline safe_sign(x)=x>0 ? 1.0 : (x<0 ? -1.0 : 0.0)
@inline _ssin(θ::Float64)=begin s=sin(θ); abs(s)>1e-12 ? s : safe_sign(s)*1e-12 end
@inline function _sincos_safe(θ::Float64); s,c=sincos(θ); abs(s)>1e-12 ? (s,c) : (safe_sign(s)*1e-12,c) end

struct HaloContext
    halo::Dict{Symbol,Any}
    R::Vector{Float64}
    tabv::Vector{Float64}
    tabfr::Vector{Float64}
    pot::Function
    frc::Function
end

const _HALO_CTX_CACHE = Dict{Tuple{Float64,Float64,Float64,Symbol,Int,Float64},HaloContext}()
reset_orbit_cache() = (empty!(_HALO_CTX_CACHE); @info "OSPM: halo context cache cleared"; nothing)

@inline function normalize_halo(halo)
    h=Dict{Symbol,Any}(); for (k,v) in halo; h[k isa Symbol ? k : Symbol(String(k))]=v; end
    if haskey(h,:type) && !(h[:type] isa Symbol); h[:type]=Symbol(lowercase(String(h[:type]))); end
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
    rs_SI=f64(r_s)*pc
    rho_SI=f64(rho_s)*Msun/pc^3
    MBH_SI=f64(MBH)*Msun
    Dict(:rho_s=>rho_SI,:r_s=>rs_SI,:rs=>rs_SI,:MBH=>MBH_SI,:type=>ht,:rmin=>1e-6*rs_SI,:cslope=>(ht===:nfw ? 1.0 : 0.0))
end

function tables_spherical(R, nlegup, halo, ρfn; extend_factor=10.0)
    halo=normalize_halo(halo); n=length(R); tabv=zeros(n); tabfr=zeros(n); ρ=similar(R)
    @inbounds for i in eachindex(R); v=ρfn((R[i],0.0),halo); ρ[i]=isfinite(v) ? v : 0.0; end
    Menc=zeros(n); @inbounds for i in 2:n; dr=R[i]-R[i-1]; Menc[i]=Menc[i-1]+0.5*dr*(R[i]^2*ρ[i]+R[i-1]^2*ρ[i-1]); end; Menc .*= 4π
    J=zeros(n);    @inbounds for i in (n-1):-1:1; dr=R[i+1]-R[i]; J[i]=J[i+1]+0.5*dr*(R[i+1]*ρ[i+1]+R[i]*ρ[i]); end; J .*= 4π
    @inbounds for i in eachindex(R); r=max(R[i],1e-30); tabv[i]=-G*(Menc[i]/r + J[i]); end
    tabfr[2:end-1] .= -(tabv[3:end].-tabv[1:end-2])./(R[3:end].-R[1:end-2]); tabfr[1]=tabfr[2]; tabfr[end]=tabfr[end-1]
    tabv, tabfr, nothing
end

function make_potential_force_funcs(halo, R, nlegup, tabv, tabfr, tabfv)
    halo=normalize_halo(halo); MBH=f64(halo[:MBH]); rmin=f64(halo[:rmin])
    rlgmin=log10(f64(R[1])); rlgmax=log10(f64(R[end])); np=length(R); rlgmax>rlgmin || error("Degenerate R grid: logR range collapsed")
    function interp(arr, r)
        rr=max(f64(r),rmin); lr=log10(rr); x=(lr-rlgmin)*(np-1)/(rlgmax-rlgmin); x=clamp(x,0.0,np-1.0)
        i0=Int(floor(x))+1; i1=min(i0+1,np); t=x-(i0-1); (1-t)*arr[i0] + t*arr[i1]
    end
    pot(r,μ=0.0)=begin rr=max(abs(f64(r)),rmin); Φh=interp(tabv,rr); Φbh=MBH>0 ? (-G*MBH/rr) : 0.0; Φh+Φbh end
    frc(r,μ=0.0)=begin rr=max(abs(f64(r)),rmin); frh=interp(tabfr,rr); frbh=MBH>0 ? (-G*MBH/(rr*rr)) : 0.0; frh+frbh,0.0 end
    pot, frc, R
end

function build_halo_context(rho_s, r_s, MBH, halo_type; nR=256, rmax_factor=300.0)
    halo=halo_from_theta(rho_s,r_s,MBH; halo_type=halo_type)
    R=build_R_halo_physical(nR; rmin=halo[:rmin], rmax=rmax_factor*halo[:rs])
    tabv,tabfr,_=tables_spherical(R,1,halo,rho_interp); pot,frc,_=make_potential_force_funcs(halo,R,1,tabv,tabfr,nothing)
    HaloContext(halo,f64.(R),tabv,tabfr,pot,frc)
end

function get_halo_context(rho_s, r_s, MBH, halo_type; nR::Int=256, rmax_factor::Float64=300.0)
    ht=Symbol(lowercase(String(halo_type))); key=(f64(rho_s),f64(r_s),f64(MBH),ht,nR,f64(rmax_factor))
    ctx=get(_HALO_CTX_CACHE,key,nothing); ctx!==nothing && return ctx
    ctx=build_halo_context(rho_s,r_s,MBH,ht; nR=nR, rmax_factor=rmax_factor); _HALO_CTX_CACHE[key]=ctx; ctx
end

function mass_enclosed_two_radii(r_in::Float64, r_out::Float64, rho_s::Float64, r_s::Float64, MBH::Float64, halo_type::String)
    ctx=get_halo_context(rho_s,r_s,MBH,halo_type); rin=max(r_in,f64(ctx.halo[:rmin])); rout=max(r_out,rin*1.001)
    fr_in,_=ctx.frc(rin,0.0); fr_out,_=ctx.frc(rout,0.0); (-rin*rin*fr_in/G, -rout*rout*fr_out/G)
end

function build_A_matrix_julia(r0_unused, th0, dt, Etot_unused, xLz_unused, r_centers_m, valid, sini, nsteps, rho_s, r_s, MBH, halo_type)
    build_A_matrix_julia(th0, dt, r_centers_m, valid, sini, nsteps, rho_s, r_s, MBH, halo_type)
end

@inline function derivs(s::SVector{4,Float64}, Lz::Float64, frc, R)
    r,θ,vr,vθ=s
    !(isfinite(r)&&isfinite(θ)&&isfinite(vr)&&isfinite(vθ)) && return SVector(0.0,0.0,0.0,0.0)
    r_safe=max(abs(r),1e-12); sθ,cθ=_sincos_safe(θ); r_tab=clamp(r_safe,R[1],R[end]); fr,_=frc(r_tab,sθ); !isfinite(fr) && return SVector(0.0,0.0,0.0,0.0)
    dr=vr; dθ=vθ/r_safe
    dvr=(vθ*vθ)/r_safe + (Lz*Lz)/(r_safe^3*sθ*sθ) + fr
    dvθ=(Lz*Lz)*cθ/(r_safe^3*sθ^3) - (vr*vθ)/r_safe
    SVector(dr,dθ,dvr,dvθ)
end

function integrate_orbit_rk4(; ic, xLz, ctx, nsteps=4000, stop_rmin_factor=1.001)
    halo=ctx.halo; rmin_stop=stop_rmin_factor*f64(halo[:rmin])
    r0=f64(ic[1]); θ0=f64(ic[2]); dt=f64(ic[3]); vr0=length(ic)>=4 ? f64(ic[4]) : 0.0; vθ0=length(ic)>=5 ? f64(ic[5]) : 0.0
    state=SVector(r0,θ0,vr0,vθ0); r=Float64[]; vr=Float64[]; θ=Float64[]; rmax_stop=10.0*f64(ctx.R_pos[end])
    @inbounds for step in 1:Int(nsteps)
        !all(isfinite,state) && break
        rr=state[1]; θr=state[2]
        (rr<=rmin_stop || rr>=rmax_stop || abs(θr)>1e6) && break
        push!(r,rr); push!(vr,state[3]); push!(θ,θr)
        k1=derivs(state,xLz,ctx.frc,ctx.R_pos); k2=derivs(state+0.5*dt*k1,xLz,ctx.frc,ctx.R_pos)
        k3=derivs(state+0.5*dt*k2,xLz,ctx.frc,ctx.R_pos); k4=derivs(state+dt*k3,xLz,ctx.frc,ctx.R_pos)
        state += (dt/6.0)*(k1+2k2+2k3+k4)
        state = SVector(state[1], clamp(state[2],1e-6,π-1e-6), state[3], state[4])
    end
    r,vr,θ
end

const _dbg_orbit_count=Ref(0)
function orbit_to_sigma2_profile(; r_arr, th_arr, vr_arr, xLz, r_centers_m, edges, sini)
    nb=length(r_centers_m); w=zeros(nb); v=zeros(nb); v2=zeros(nb); _dbg_orbit_count[]+=1
    do_debug=_dbg_orbit_count[]<=5; hits=do_debug ? zeros(Int,nb) : nothing; rmin_seen=do_debug ? Inf : 0.0; rmax_seen=do_debug ? -Inf : 0.0
    @inbounds for i in eachindex(r_arr)
        rr=max(f64(r_arr[i]),1e-12); th=clamp(f64(th_arr[i]),1e-6,π-1e-6); vr0=f64(vr_arr[i])
        ss,_=_sincos_safe(th); vφ=f64(xLz)/(rr*ss); vlos=f64(sini)*vr0 + sqrt(1-f64(sini)^2)*vφ
        j=searchsortedfirst(edges,rr)-1; (j<1 || j>nb) && continue
        w[j]+=1; v[j]+=vlos; v2[j]+=vlos^2
        if do_debug; hits[j]+=1; rmin_seen=min(rmin_seen,rr); rmax_seen=max(rmax_seen,rr); end
    end
    do_debug && println("BIN_HITS=",hits," rmin/max=",(rmin_seen,rmax_seen)," centers[1/end]=",(r_centers_m[1],r_centers_m[end]))
    σ2=zeros(nb); @inbounds for j in 1:nb; w[j]>0 && (μv=v[j]/w[j]; σ2[j]=v2[j]/w[j]-μv^2); end
    σ2
end

function launch_orbit_apocenter(; rapo::Float64, θ0::Float64, Lz_frac::Float64, pot, frc, r0_frac::Float64=0.98, dt_frac::Float64=0.01, dt_floor::Float64=1e-30, debug::Bool=false)
    ss=_ssin(θ0); frs,_=frc(rapo,ss); !(isfinite(frs)&&isfinite(rapo)&&rapo>0) && return (nothing,0.0,0.0,0.0,:reject_force)
    vc=sqrt(abs(frs)*rapo); !(isfinite(vc)&&vc>0) && return (nothing,0.0,0.0,0.0,:reject_vc)
    Lz=Lz_frac*rapo*vc; Φapo=pot(rapo,ss); !isfinite(Φapo) && return (nothing,0.0,0.0,vc,:reject_pot)
    E=Φapo + (Lz^2)/(2*rapo^2*ss^2); r0=r0_frac*rapo; Φ0=pot(r0,ss); !isfinite(Φ0) && return (nothing,0.0,E,vc,:reject_pot0)
    arg=2*(E-Φ0) - (Lz^2)/(r0^2*ss^2)
    if !(isfinite(arg) && arg>0)
        debug ? ((rapo,θ0,Lz,arg),Lz,E,vc,:reject_turning) : (nothing,Lz,E,vc,:reject_turning)
    else
        vr0=-sqrt(arg); Ω=abs(vc/r0); dt=dt_frac/max(Ω,1e-30); ((r0,θ0,dt,vr0,0.0),Lz,E,vc,:ok)
    end
end

function build_A_matrix_stellar(Norbit::Int, R_star_m::Vector{Float64}, v_star_mps::Vector{Float64}, verr_star_mps::Vector{Float64},
        sini::Float64, rho_s::Float64, r_s::Float64, MBH::Float64, halo_type::String;
        nsteps::Int=4000, Lfrac::NTuple{5,Float64}=(0.05,0.2,0.4,0.7,1.0), dt_frac_orbit::Float64=0.01,
        ΔR_frac::Float64=0.05, ΔR_floor_frac::Float64=0.01, ΔR_floor_pc::Float64=0.0,
        Nbins_occ::Int=6, return_occ::Bool=true, max_attempts_factor::Int=60)
    Nstar=length(R_star_m); @assert length(v_star_mps)==Nstar; @assert length(verr_star_mps)==Nstar
    Nstar==0 && return zeros(Float64,0,Norbit)
    ctx=get_halo_context(rho_s,r_s,MBH,halo_type); cosi=sqrt(max(1.0-sini*sini,0.0))
    Rmin=minimum(R_star_m); Rmax=maximum(R_star_m); !(isfinite(Rmin)&&isfinite(Rmax)&&Rmax>Rmin) && return zeros(Float64,Nstar,Norbit)
    shells=sort(copy(R_star_m))
    dnn=Float64[]; if length(shells)>=3; resize!(dnn,length(shells)-1); @inbounds for i in 1:length(shells)-1; dnn[i]=shells[i+1]-shells[i]; end; end
    med_dnn=(length(dnn)>0) ? median(abs.(dnn)) : 0.0
    auto_floor=max(ΔR_floor_frac*max(med_dnn,0.0), 0.02*Rmin); abs_floor=(ΔR_floor_pc>0) ? (ΔR_floor_pc*pc) : 0.0
    ΔR_floor=max(auto_floor,abs_floor,1e-12); occ_edges=collect(range(Rmin,Rmax; length=Nbins_occ+1))
    A_star=zeros(Float64,Nstar,Norbit); A_occ=zeros(Float64,Nbins_occ,Norbit)
    orbit_ctx=(frc=ctx.frc, R_pos=ctx.R, halo=ctx.halo); rng=MersenneTwister(0x5eed1234); θ0=f64(π/2)
    col=1; idx=1; max_attempts=max_attempts_factor*Norbit; attempts=0
    while col<=Norbit && attempts<max_attempts
        attempts+=1; rapo=f64(shells[idx]); idx=(idx==length(shells)) ? 1 : (idx+1); !(isfinite(rapo)&&rapo>0) && continue
        lf=Lfrac[1 + (attempts % length(Lfrac))]; r0_frac=0.95 + 0.04*rand(rng)
        ic,Lz0,E0,vc,st=launch_orbit_apocenter(rapo=rapo,θ0=θ0,Lz_frac=f64(lf),pot=ctx.pot,frc=ctx.frc,r0_frac=r0_frac,dt_frac=dt_frac_orbit)
        st!=:ok && continue
        r,vr,θ=integrate_orbit_rk4(ic=ic,xLz=Lz0,ctx=orbit_ctx,nsteps=nsteps); isempty(r) && continue
        s_arr=similar(r); vphi_arr=similar(r)
        @inbounds for i in eachindex(r); si=_ssin(f64(θ[i])); ri=f64(r[i]); s_arr[i]=ri*si; vphi_arr[i]=f64(Lz0)/(ri*si); end
        Nhits=length(s_arr)
        if return_occ && Nhits>0
            @inbounds for j in 1:Nbins_occ
                lo=occ_edges[j]; hi=occ_edges[j+1]; cnt=0
                @inbounds for x in s_arr; cnt += (x≥lo && x<hi) ? 1 : 0; end
                A_occ[j,col]=cnt/Nhits
            end
        end
        vlos_arr=@. sini*f64(vr) + cosi*vphi_arr
        @inbounds for istar in 1:Nstar
            Ri=f64(R_star_m[istar]); vi=f64(v_star_mps[istar]); σi=f64(verr_star_mps[istar]); σ2=σi*σi
            if !(isfinite(σ2)&&σ2>0); A_star[istar,col]=0.0; continue; end
            ΔR=max(ΔR_frac*Ri,ΔR_floor); p=0.0; nhit=0
            @inbounds for j in eachindex(s_arr)
                if abs(s_arr[j]-Ri) < ΔR
                    dv=vi-vlos_arr[j]; p += exp(-0.5*(dv*dv)/σ2); nhit += 1
                end
            end
            A_star[istar,col] = (nhit==0) ? 0.0 : (p/(sqrt(2π*σ2)*nhit))
        end
        col += 1
    end
    col<=Norbit && println("WARNING: build_A_matrix_stellar filled ",col-1," / ",Norbit," after attempts=",attempts)
    return_occ ? vcat(A_star,A_occ) : A_star
end

function build_A_matrix_from_ctx(ctx::HaloContext, th0, dt, r_centers_m, valid, sini, nsteps)
    shells=r_centers_m[valid]; Ndat=length(shells); θ0=f64(th0[1]); rng=MersenneTwister(0xdeadbeef)
    dt_frac_orbit=0.01; Lfrac=(0.05,0.2,0.4,0.7,1.0); Norbit=length(shells)*length(Lfrac)
    A=zeros(Float64,Ndat,Norbit)
    edges=similar(r_centers_m,length(r_centers_m)+1); edges[2:end-1].=0.5.*(r_centers_m[1:end-1].+r_centers_m[2:end]); edges[1]=0.0; edges[end]=Inf
    kept=0; rej_turn=0; rej_energy=0; rej_force=0; dErel_max=0.0; dErel_bad=0; col=1
    orbit_ctx=(frc=ctx.frc, R_pos=ctx.R, halo=ctx.halo)
    @inbounds for rapo0 in shells
        rapo=f64(rapo0); if !(isfinite(rapo)&&rapo>0); col += length(Lfrac); continue; end
        for lf in Lfrac
            r0_frac=0.95 + 0.04*rand(rng)
            ic,Lz0,E0,vc,st=launch_orbit_apocenter(rapo=rapo,θ0=θ0,Lz_frac=f64(lf),pot=ctx.pot,frc=ctx.frc,r0_frac=r0_frac,dt_frac=dt_frac_orbit)
            if st!=:ok
                st==:reject_turning ? (rej_turn+=1) : (st==:reject_force || st==:reject_vc ? (rej_force+=1) : (rej_energy+=1))
                col += 1; continue
            end
            r,vr,θ=integrate_orbit_rk4(ic=ic,xLz=Lz0,ctx=orbit_ctx,nsteps=nsteps)
            isempty(r) && (rej_energy+=1; col+=1; continue)
            σ2=orbit_to_sigma2_profile(r_arr=r, th_arr=θ, vr_arr=vr, xLz=Lz0, r_centers_m=r_centers_m, edges=edges, sini=sini)
            A[:,col] .= σ2[valid]; kept += 1
            if kept <= 8
                Emin=Inf; Emax=-Inf; Esum=0.0; nE=0
                @inbounds for i in eachindex(r)
                    si=_ssin(θ[i]); Ei=0.5*(vr[i]^2 + (Lz0/(r[i]*si))^2) + ctx.pot(r[i],si)
                    if isfinite(Ei); Emin=min(Emin,Ei); Emax=max(Emax,Ei); Esum+=Ei; nE+=1; end
                end
                if nE>0
                    Emean=Esum/nE; dErel=(Emax-Emin)/max(abs(Emean),1e-300); dErel_max=max(dErel_max,dErel); dErel_bad += (dErel>1e-3) ? 1 : 0
                    println("ΔE/<E>=",dErel," rapo=",rapo," Lfrac=",lf," n=",length(r))
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

ospm_runcheck(theta, args...) = build_A_matrix_julia(args..., theta[1], theta[2], length(theta)>2 ? theta[3] : 0.0)

# ============================================================
# Public batch entry for Python
# ============================================================

# status: 0 pass | 1 orbit_fail | 2 numeric_fail | 3 unknown_fail
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

    status = Vector{Int}(undef, nbatch)
    chi2   = Vector{Float64}(undef, nbatch)

    Threads.@threads for i in 1:nbatch
        θ = getθ(i)
        try
            rho_s = θ[1]; r_s = θ[2]; MBH = dim ≥ 3 ? θ[3] : 0.0
            A = build_A_matrix_julia(nothing, nothing, r_centers_m, valid, sini, nsteps, rho_s, r_s, MBH, halo_type)
            v = sum(abs2, A)
            if isfinite(v); status[i]=0; chi2[i]=v
            else; status[i]=2; chi2[i]=Inf
            end
        catch
            status[i]=3; chi2[i]=Inf
        end
    end

    status, chi2
end

end
