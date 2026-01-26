# OSPM/Physics/OSPM_Physics.py
# contract + cache authority + diagnostics hooks (compressed)

import os, numpy as np

USE_JULIA=os.environ.get("OSPM_USE_JULIA","0").strip().lower() in ("1","true","yes")
_JL_READY=False; _Main=None

pc=3.085677581e16; kms=1.0e3; Msun=1.98847e30; G=6.67430e-11; c=2.99792458e8
_LAST_SIG=None

def _jl_init():
    global _JL_READY,_Main
    if _JL_READY: return
    if not USE_JULIA: raise RuntimeError("OSPM_USE_JULIA is not enabled")
    from julia.api import Julia
    Julia(compiled_modules=False)
    from julia import Main as _JuliaMain
    _Main=_JuliaMain
    here=os.path.dirname(os.path.abspath(__file__))
    jl_path=os.path.join(here,"OSPM_Physics_Spherical.jl")
    if not os.path.exists(jl_path): raise FileNotFoundError(f"Julia backend file not found: {jl_path}")
    _Main.include(jl_path)
    if not hasattr(_Main,"OSPMPhysicsSpherical"): raise RuntimeError("Julia loaded but module OSPMPhysicsSpherical not found in Main")
    _JL_READY=True

def _theta_sig(theta,halo_type):
    t=np.asarray(theta,float).ravel()
    if t.size<2: raise ValueError("theta must have at least [rho_s, r_s]")
    rho_s=float(t[0]); r_s=float(t[1]); MBH=float(t[2]) if t.size>=3 else 0.0
    ht=str(halo_type).strip().lower()
    return (rho_s,r_s,MBH,ht)

def assert_theta_contract(theta,*,halo_type,bounds=None,require_mbh=True):
    t=np.asarray(theta,float).ravel()
    if t.size<2: raise ValueError("theta too short")
    if require_mbh and t.size<3: raise ValueError("theta missing MBH (expects [rho_s,r_s,MBH])")
    if not np.all(np.isfinite(t[:3])): raise ValueError("theta has non-finite values")
    rho_s=float(t[0]); r_s=float(t[1]); MBH=float(t[2]) if t.size>=3 else 0.0
    if rho_s<=0 or r_s<=0: raise ValueError("theta expects rho_s>0 and r_s>0")
    if MBH<0: raise ValueError("theta expects MBH>=0")
    ht=str(halo_type).strip().lower()
    if ht not in ("nfw","cored"): raise ValueError(f"Unknown halo_type: {halo_type}")
    if bounds is not None:
        b=np.asarray(bounds,float)
        if b.shape[0]<3: raise ValueError("bounds must cover at least rho_s,r_s,MBH")
        for i,x in enumerate((rho_s,r_s,MBH)):
            lo,hi=float(b[i,0]),float(b[i,1])
            if not (lo<=x<=hi): raise ValueError(f"theta out of bounds at i={i}: {x} not in [{lo},{hi}]")
    return (rho_s,r_s,MBH,ht)

def reset_orbit_cache_julia():
    if not USE_JULIA: raise RuntimeError("Julia required")
    _jl_init(); _Main.OSPMPhysicsSpherical.reset_orbit_cache(); return None

def maybe_reset_orbit_cache(theta,halo_type):
    global _LAST_SIG
    sig=_theta_sig(theta,halo_type)
    if _LAST_SIG==sig: return
    _LAST_SIG=sig
    reset_orbit_cache_julia()

def mass_enclosed_two_radii_julia(*,r_in_m,r_out_m,theta,halo_type):
    if not USE_JULIA: raise RuntimeError("Julia required")
    _jl_init()
    rho_s,r_s,MBH,ht=assert_theta_contract(theta,halo_type=halo_type,require_mbh=True)
    Min,Mout=_Main.OSPMPhysicsSpherical.mass_enclosed_two_radii(float(r_in_m),float(r_out_m),float(rho_s),float(r_s),float(MBH),str(ht))
    return float(Min),float(Mout)

def make_inclination(inclination_deg: float):
    inc=np.radians(float(inclination_deg))
    sini=float(np.sin(inc)); cosi=float(np.cos(inc))
    edge_on=float(inclination_deg)>=85.0
    return sini,cosi,edge_on

def halo_from_theta_astro(theta,halo_type="nfw"):
    rho_s,r_s,MBH,ht=assert_theta_contract(theta,halo_type=halo_type,require_mbh=True)
    return {"rho_s":rho_s,"r_s":r_s,"MBH":MBH,"type":ht}

def build_dynamics_context(*,theta,halo_type,**_ignored):
    if not USE_JULIA: raise RuntimeError("Python backend is disabled. Set OSPM_USE_JULIA=1.")
    return {"halo":halo_from_theta_astro(theta,halo_type=halo_type)}

def halo_kwargs_from_ctx(ctx):
    halo=ctx["halo"] if isinstance(ctx,dict) else ctx.halo
    for k in ("rho_s","r_s","MBH","type"):
        if k not in halo: raise KeyError(f"halo missing required key '{k}'")
    return dict(rho_s=float(halo["rho_s"]),r_s=float(halo["r_s"]),MBH=float(halo["MBH"]),halo_type=str(halo["type"]))

def build_A_matrix_stellar_julia(*,R_star_m,v_star_mps,verr_star_mps,sini,Norbit,theta,halo_type,return_occ=False,Nbins_occ=6,diag=False):
    if not USE_JULIA: raise RuntimeError("Stellar mode requires Julia")
    _jl_init()
    rho_s,r_s,MBH,ht=assert_theta_contract(theta,halo_type=halo_type,require_mbh=True)
    maybe_reset_orbit_cache((rho_s,r_s,MBH),ht)
    out=_Main.OSPMPhysicsSpherical.build_A_matrix_stellar(
        int(Norbit),
        np.asarray(R_star_m,float),
        np.asarray(v_star_mps,float),
        np.asarray(verr_star_mps,float),
        float(sini),
        float(rho_s),float(r_s),float(MBH),
        str(ht),
        return_occ=bool(return_occ),
        Nbins_occ=int(Nbins_occ),
        diag=bool(diag),
    )
    if diag:
        A,meta=out
        return np.asarray(A,float), dict(meta)
    return np.asarray(out,float)

def build_A_matrix(obs,ctx,*,return_occ=False,Nbins_occ=6,diag=False):
    mode=str(getattr(obs,"mode","")).strip().lower()
    if mode!="stellar": raise RuntimeError("build_A_matrix only supports obs.mode=='stellar'")
    hk=halo_kwargs_from_ctx(ctx)
    theta=[hk["rho_s"],hk["r_s"],hk["MBH"]]
    halo_type=hk["halo_type"]

    vv=np.asarray(getattr(obs,"valid_vlos",None),bool) if hasattr(obs,"valid_vlos") else None
    if vv is None:
        has=np.asarray(getattr(obs,"has_vlos",None),bool) if hasattr(obs,"has_vlos") else None
        if has is None: raise RuntimeError("obs missing valid_vlos (or has_vlos)")
        vv=has & np.isfinite(obs.v_star_mps) & np.isfinite(obs.verr_star_mps) & (np.asarray(obs.verr_star_mps,float)>0)

    R=np.asarray(obs.R_star_m,float)[vv]
    v=np.asarray(obs.v_star_mps,float)[vv]
    ve=np.asarray(obs.verr_star_mps,float)[vv]
    if len(R)==0: raise RuntimeError("No RV-valid stars to build A-matrix")

    return build_A_matrix_stellar_julia(
        R_star_m=R,v_star_mps=v,verr_star_mps=ve,
        sini=float(obs.sini),Norbit=int(obs.Norbit),
        theta=theta,halo_type=halo_type,
        return_occ=bool(return_occ),Nbins_occ=int(Nbins_occ),
        diag=bool(diag),
    )

def build_A_matrix_from_theta(obs,theta,*,halo_type="nfw",return_occ=True,Nbins_occ=6,diag=False):
    ctx=build_dynamics_context(theta=theta,halo_type=halo_type)
    return build_A_matrix(obs,ctx,return_occ=bool(return_occ),Nbins_occ=int(Nbins_occ),diag=bool(diag))

def rho_interp(*args,**kwargs):
    raise RuntimeError("rho_interp is legacy-only. It should never be called in Julia mode.")