from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from ..physics.OSPM_Physics import build_A_matrix_stellar_julia

def _as_1d_float(x, name: str) -> np.ndarray:
    a=np.asarray(x,dtype=float).reshape(-1)
    if a.size==0: raise ValueError(f"{name} is empty")
    if not np.all(np.isfinite(a)): raise ValueError(f"{name} contains non-finite values")
    return a

def stellar_log_likelihood(A:np.ndarray,w:np.ndarray,*,verr_star_mps:np.ndarray,Nstar:int,Nocc:int,
    lambda_occ:float=1.0,eps:float=1e-300,sigma_floor_mps:float=2e3)->float:
    A=np.asarray(A,float); w=_as_1d_float(w,"w")
    if A.ndim!=2: raise ValueError("A must be 2D")
    if A.shape[1]!=w.size: raise ValueError(f"A.shape[1]={A.shape[1]} != len(w)={w.size}")
    Nrow=A.shape[0]
    if Nstar<0 or Nocc<0: raise ValueError("Nstar and Nocc must be >= 0")
    if Nstar+Nocc!=Nrow: raise ValueError(f"Nstar+Nocc={Nstar+Nocc} does not match A rows={Nrow}")
    p=A@w; p=np.maximum(p,float(eps))
    if Nstar>0:
        p_star=p[:Nstar]
        sigma=np.asarray(verr_star_mps,float); sigma=np.maximum(sigma,float(sigma_floor_mps))
        if len(sigma)!=Nstar: raise ValueError("verr_star_mps length must match Nstar")
        p_star=np.maximum(p_star*sigma,float(eps))
        logL_star=float(np.sum(np.log(p_star)))
    else:
        logL_star=0.0
    if Nocc==0: return logL_star
    p_occ=np.maximum(p[Nstar:],float(eps))
    logL_occ=float(np.sum(np.log(p_occ)))
    return logL_star+float(lambda_occ)*logL_occ

def solve_weights_stellar(A:np.ndarray,*,verr_star_mps:np.ndarray,Nstar:int,Nocc:int,lambda_occ:float=1.0,
    alpha:float=1e-3,eps:float=1e-300,maxiter:int=500,maxfun:int=20000,p0_floor:float=1e-15)->np.ndarray:
    A=np.asarray(A,float); Norbit=int(A.shape[1])
    alpha_eff=float(alpha)*(Norbit/200.0)*(max(Nstar,1)/90.0)
    w0=np.random.rand(Norbit); w0+=1e-3; w0/=w0.sum()
    p0=A@w0; scale=1.0/np.maximum(p0,float(p0_floor))
    A_eff=A*scale[:,None]
    bounds=[(0.0,None)]*Norbit
    def objective(w:np.ndarray)->float:
        ll=stellar_log_likelihood(A_eff,w,verr_star_mps=verr_star_mps,Nstar=Nstar,Nocc=Nocc,lambda_occ=lambda_occ,eps=eps)
        return -ll+alpha_eff*float(np.dot(w,w))
    res=minimize(objective,w0,method="L-BFGS-B",bounds=bounds,options={"maxiter":int(maxiter),"maxfun":int(maxfun)})
    if (not res.success) or (res.x is None):
        res=minimize(objective,w0,method="SLSQP",bounds=bounds,options={"maxiter":int(maxiter*3),"ftol":1e-9,"disp":False})
        if (not res.success) or (res.x is None):
            p0_raw=A@w0; p0_eff=A_eff@w0; msg=getattr(res,"message","unknown")
            raise RuntimeError(
                f"Weight solve failed: {msg} | "
                f"p0raw_min={p0_raw.min():.3e} p0raw_max={p0_raw.max():.3e} "
                f"bad0raw={(p0_raw<=1e-200).sum()}/{len(p0_raw)} | "
                f"p0eff_min={p0_eff.min():.3e} p0eff_max={p0_eff.max():.3e} "
                f"bad0eff={(p0_eff<=1e-200).sum()}/{len(p0_eff)}"
            )
    w=np.asarray(res.x,float); w[w<0]=0.0
    s=float(w.sum())
    if (not np.isfinite(s)) or (s<=0): w=np.full(Norbit,1.0/Norbit)
    else: w/=s
    return w

def solve_ospm_theta_stellar(theta, obs, *, halo_type: str="nfw"):
    rho_s,r_s,MBH=map(float,theta)

    R_all=_as_1d_float(obs.R_star_m,"obs.R_star_m")
    v_all=_as_1d_float(obs.v_star_mps,"obs.v_star_mps")
    ve_all=_as_1d_float(obs.verr_star_mps,"obs.verr_star_mps")
    hv=np.asarray(getattr(obs,"has_vlos",None),bool) if hasattr(obs,"has_vlos") else None
    vv=np.asarray(getattr(obs,"valid_vlos",None),bool) if hasattr(obs,"valid_vlos") else None

    if not (len(R_all)==len(v_all)==len(ve_all)): raise ValueError("obs star arrays length mismatch")
    if vv is None:
        if hv is None: vv=np.isfinite(v_all)&np.isfinite(ve_all)&(ve_all>0)
        else: vv=hv&np.isfinite(v_all)&np.isfinite(ve_all)&(ve_all>0)
    if len(vv)!=len(R_all): raise ValueError("obs.valid_vlos length mismatch")

    R=R_all[vv]; v=v_all[vv]; ve=ve_all[vv]
    Nstar=int(len(R))

    Norbit=int(obs.Norbit)
    if Norbit<1: raise ValueError("obs.Norbit must be >= 1")
    sini=float(obs.sini)
    if not np.isfinite(sini): raise ValueError("obs.sini is non-finite")

    Nocc=int(getattr(obs,"Nocc",0))
    lambda_occ=float(getattr(obs,"lambda_occ",1.0))

    if Nstar==0:
        raise RuntimeError("No valid RV stars. Velocity likelihood has zero constraints. Enable occupancy-only mode if desired.")

    A=build_A_matrix_stellar_julia(
        R_star_m=R,v_star_mps=v,verr_star_mps=ve,
        sini=sini,Norbit=Norbit,theta=[rho_s,r_s,MBH],
        halo_type=halo_type,return_occ=False
    )
    A=np.asarray(A,float)
    if A.ndim!=2 or A.shape[1]!=Norbit: raise ValueError(f"Unexpected A shape: {A.shape}, expected (*,{Norbit})")

    if A.shape[0]==Nstar: Nocc_eff=0
    elif A.shape[0]==Nstar+Nocc: Nocc_eff=Nocc
    else: raise ValueError(f"A rows={A.shape[0]} do not match Nstar={Nstar} or Nstar+Nocc={Nstar+Nocc}")

    w=solve_weights_stellar(
        A,verr_star_mps=ve,Nstar=Nstar,Nocc=Nocc_eff,lambda_occ=lambda_occ,
        alpha=float(getattr(obs,"alpha_w",1e-2)),maxiter=int(getattr(obs,"w_maxiter",500))
    )
    ll=stellar_log_likelihood(A,w,verr_star_mps=ve,Nstar=Nstar,Nocc=Nocc_eff,lambda_occ=lambda_occ)
    chi2_like=-2.0*ll
    nu=max(Nstar-3,1)
    chi2_red=float(chi2_like)/float(nu)
    model=A@w
    return float(chi2_red),w,model
