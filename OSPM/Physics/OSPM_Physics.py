import os
import numpy as np

USE_JULIA=os.environ.get("OSPM_USE_JULIA","0").strip().lower() in ("1","true","yes")

_JL_READY=False
_Main=None

pc=3.085677581e16
kms=1.0e3
Msun=1.98847e30
G=6.67430e-11
c=2.99792458e8

_ORBIT_LIB_SIG=None

def _potential_signature(theta,halo_type):
    rho_s,r_s,MBH=map(float,theta[:3])
    ht=str(halo_type).strip().lower()
    return (ht,float(np.round(rho_s,12)),float(np.round(r_s,12)),float(np.round(MBH,12)))

def _maybe_reset_orbit_library(theta,halo_type):
    global _ORBIT_LIB_SIG
    sig=_potential_signature(theta,halo_type)
    if _ORBIT_LIB_SIG==sig: return
    _ORBIT_LIB_SIG=sig
    if (_Main is None) or (not hasattr(_Main,"OSPMPhysicsSpherical")):
        raise RuntimeError("Julia backend not initialized")
    if hasattr(_Main.OSPMPhysicsSpherical,"reset_orbit_cache"):
        _Main.OSPMPhysicsSpherical.reset_orbit_cache()
        return
    raise RuntimeError("Julia backend missing reset_orbit_cache(). Add it to OSPM_Physics_Spherical.jl.")

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
    if not hasattr(_Main,"OSPMPhysicsSpherical"):
        raise RuntimeError("Julia loaded but module OSPMPhysicsSpherical not found in Main")
    _JL_READY=True

def mass_enclosed_two_radii_julia(*,r_in_m,r_out_m,theta,halo_type):
    if not USE_JULIA: raise RuntimeError("Julia required")
    _jl_init()
    rho_s,r_s,MBH=map(float,theta[:3])
    Min,Mout=_Main.OSPMPhysicsSpherical.mass_enclosed_two_radii(
        float(r_in_m),float(r_out_m),float(rho_s),float(r_s),float(MBH),str(halo_type)
    )
    return float(Min),float(Mout)

def make_inclination(inclination_deg: float):
    inc=np.radians(float(inclination_deg))
    sini=float(np.sin(inc)); cosi=float(np.cos(inc))
    edge_on=float(inclination_deg)>=85.0
    return sini,cosi,edge_on

def make_rng(seed: int):
    return np.random.default_rng(int(seed))

def halo_from_theta_astro(theta,halo_type="nfw"):
    theta=np.asarray(theta,float)
    rho_s=float(theta[0]); r_s=float(theta[1])
    MBH=float(theta[2]) if len(theta)>2 else 0.0
    ht=str(halo_type).strip().lower()
    return {"rho_s":rho_s,"r_s":r_s,"MBH":MBH,"type":ht}

def build_dynamics_context(*,theta,halo_type,**_ignored):
    if not USE_JULIA: raise RuntimeError("Python backend is disabled. Set OSPM_USE_JULIA=1.")
    return {"halo":halo_from_theta_astro(theta,halo_type=halo_type)}

def halo_kwargs_from_ctx(ctx):
    halo=ctx["halo"] if isinstance(ctx,dict) else ctx.halo
    for k in ("rho_s","r_s","MBH","type"):
        if k not in halo: raise KeyError(f"halo missing required key '{k}'")
    return dict(rho_s=float(halo["rho_s"]),r_s=float(halo["r_s"]),MBH=float(halo["MBH"]),halo_type=str(halo["type"]))

def build_A_matrix_stellar_julia(*,R_star_m,v_star_mps,verr_star_mps,sini,Norbit,theta,halo_type,return_occ=False,Nbins_occ=6):
    if not USE_JULIA: raise RuntimeError("Stellar mode requires Julia")
    _jl_init()
    _maybe_reset_orbit_library(theta,halo_type)
    rho_s,r_s,MBH=map(float,theta[:3])
    A=_Main.OSPMPhysicsSpherical.build_A_matrix_stellar(
        int(Norbit),
        np.asarray(R_star_m,float),
        np.asarray(v_star_mps,float),
        np.asarray(verr_star_mps,float),
        float(sini),
        float(rho_s),float(r_s),float(MBH),
        str(halo_type),
        return_occ=bool(return_occ),
        Nbins_occ=int(Nbins_occ),
    )
    return np.asarray(A,float)

def build_A_matrix(obs,ctx,*,return_occ=False,Nbins_occ=6):
    mode=str(getattr(obs,"mode","")).strip().lower()
    if mode!="stellar":
        raise RuntimeError("build_A_matrix currently only supports obs.mode=='stellar'")

    hk=halo_kwargs_from_ctx(ctx)
    theta=[hk["rho_s"],hk["r_s"],hk["MBH"]]
    halo_type=hk["halo_type"]

    vv=np.asarray(getattr(obs,"valid_vlos",None),bool) if hasattr(obs,"valid_vlos") else None
    if vv is None:
        has=np.asarray(getattr(obs,"has_vlos",None),bool) if hasattr(obs,"has_vlos") else None
        if has is None:
            raise RuntimeError("obs missing valid_vlos (or has_vlos) for stellar mode")
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
    )

def rho_interp(*args,**kwargs):
    raise RuntimeError("rho_interp is legacy-only. It should never be called in Julia mode.")
