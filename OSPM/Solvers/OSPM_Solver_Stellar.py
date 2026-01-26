from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from ..Physics.OSPM_Physics import build_A_matrix_stellar_julia
from ..Physics.OSPM_PhysicsEngine import chi2_resolution_penalty, mass_slope_penalty

def _as_1d_float(x,name):
    a=np.asarray(x,float).reshape(-1)
    if a.size==0 or not np.all(np.isfinite(a)): raise ValueError(f"{name} invalid")
    return a

def stellar_log_likelihood(A,w,*,verr_star_mps,rv_mask,Nstar,Nocc,lambda_occ=1.0,eps=1e-300,sigma_floor_mps=2e3):
    A=np.asarray(A,float); w=_as_1d_float(w,"w")
    if A.ndim!=2 or A.shape[1]!=w.size or A.shape[0]!=Nstar+Nocc: raise ValueError("A shape mismatch")
    p=np.maximum(A@w,eps)
    ll=0.0
    if Nstar>0:
        rv_mask=np.asarray(rv_mask,bool)
        if rv_mask.size!=Nstar: raise ValueError("rv_mask mismatch")
        if rv_mask.any():
            sig=np.asarray(verr_star_mps,float)[rv_mask]
            sig=np.maximum(np.nan_to_num(sig,np.inf),sigma_floor_mps)
            ll+=float(np.sum(np.log(np.maximum(p[:Nstar][rv_mask]*sig,eps))))
    if Nocc>0: ll+=float(lambda_occ)*float(np.sum(np.log(np.maximum(p[Nstar:],eps))))
    return ll

def solve_weights_stellar(A,*,verr_star_mps,rv_mask,Nstar,Nocc,lambda_occ=1.0,alpha=1e-3,eps=1e-300,maxiter=500,maxfun=20000,p0_floor=1e-15):
    A=np.asarray(A,float); Norbit=A.shape[1]
    alpha*=Norbit/200.0*max(Nstar,1)/90.0
    w0=np.random.rand(Norbit)+1e-3; w0/=w0.sum()
    rv_mask=np.asarray(rv_mask,bool)
    active=np.concatenate([rv_mask,np.ones(Nocc,bool)]) if Nocc>0 else rv_mask
    p0=A@w0; scale=1.0/np.maximum(p0[active],p0_floor)
    Aeff=A.copy(); Aeff[active]*=scale[:,None]
    bounds=[(0.0,None)]*Norbit

    def obj(w):
        ll=stellar_log_likelihood(Aeff,w,verr_star_mps=verr_star_mps,rv_mask=rv_mask,Nstar=Nstar,Nocc=Nocc,lambda_occ=lambda_occ,eps=eps)
        return -ll+alpha*np.dot(w,w)

    res=minimize(obj,w0,method="L-BFGS-B",bounds=bounds,options={"maxiter":maxiter,"maxfun":maxfun})
    if not res.success:
        res=minimize(obj,w0,method="SLSQP",bounds=bounds,options={"maxiter":3*maxiter,"ftol":1e-9})
        if not res.success: raise RuntimeError(f"Weight solve failed: {res.message}")
    w=np.maximum(res.x,0.0); s=w.sum()
    return w/s if s>0 else np.full(Norbit,1.0/Norbit)

def solve_ospm_theta_stellar(theta,obs,*,halo_type="nfw"):
    rho_s,r_s,MBH=map(float,theta)
    R=np.asarray(obs.R_star_m,float); v=np.asarray(obs.v_star_mps,float); ve=np.asarray(obs.verr_star_mps,float)
    if R.size==0 or not np.all(np.isfinite(R)): raise ValueError("R invalid")
    rv=np.asarray(obs.valid_vlos,bool)
    if rv.size!=R.size: raise ValueError("rv mask mismatch")
    Norbit=int(obs.Norbit); sini=float(obs.sini)
    Nocc=int(getattr(obs,"Nocc",0)); lambda_occ=float(getattr(obs,"lambda_occ",1.0))

    A=build_A_matrix_stellar_julia(R_star_m=R,v_star_mps=v,verr_star_mps=ve,sini=sini,Norbit=Norbit,theta=[rho_s,r_s,MBH],halo_type=halo_type,return_occ=True)
    A=np.asarray(A,float)
    if A.ndim!=2 or A.shape[1]!=Norbit or A.shape[0]!=R.size+Nocc: raise ValueError("A malformed")

    w=solve_weights_stellar(A,verr_star_mps=ve,rv_mask=rv,Nstar=R.size,Nocc=Nocc,lambda_occ=lambda_occ,alpha=float(getattr(obs,"alpha_w",1e-2)),maxiter=int(getattr(obs,"w_maxiter",500)))
    ll=stellar_log_likelihood(A,w,verr_star_mps=ve,rv_mask=rv,Nstar=R.size,Nocc=Nocc,lambda_occ=lambda_occ)
    chi2_like=-2.0*ll

    pen1,_,_,_=chi2_resolution_penalty(MBH_msun=MBH,R_star_m=R,v_star_mps=v,strength=float(getattr(obs,"PEN_SPHERE_STRENGTH",2.0)),power=float(getattr(obs,"PEN_SPHERE_POWER",2.0)))
    pen2,_,_,_,_,_=mass_slope_penalty(theta=[rho_s,r_s,MBH],halo_type=halo_type,R_star_m=R,strength=float(getattr(obs,"PEN_SLOPE_STRENGTH",0.5)))
    chi2=chi2_like+pen1+pen2

    nu=max(int(rv.sum())-3,1)
    return float(chi2/nu), w, A@w