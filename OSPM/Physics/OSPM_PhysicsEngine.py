# OSPM/Physics/OSPM_PhysicsEngine.py
# reward shaping (compressed) + optional physics diagnostics pass-through

import numpy as np
from . import OSPM_Physics as phys
from scipy.optimize import nnls

G=6.67430e-11
Msun=1.98847e30
_PRINT_EVERY=200
_print_counter=0

def _robust_sigma(v):
    v=np.asarray(v,float); v=v[np.isfinite(v)]
    if v.size<10: return float(np.std(v)) if v.size else np.nan
    med=np.median(v); mad=np.median(np.abs(v-med))
    return float(1.4826*mad)

def _inner_resolution_scale(R):
    R=np.sort(np.asarray(R,float)); R=R[np.isfinite(R)]
    if R.size<10: return np.nan
    dR=np.diff(R); dR=dR[np.isfinite(dR)&(dR>0)]
    if dR.size==0: return np.nan
    return float(max(np.median(dR),0.05*R[0]))

def chi2_resolution_penalty(MBH_msun,R_star_m,v_star_mps,strength=2.0,power=2.0):
    sig=_robust_sigma(v_star_mps); rmin=_inner_resolution_scale(R_star_m)
    if not np.isfinite(sig) or sig<=0: return 0.0,np.nan,np.nan,np.nan
    if not np.isfinite(rmin) or rmin<=0: return 0.0,np.nan,np.nan,np.nan
    if not np.isfinite(MBH_msun) or MBH_msun<=0: return 0.0,np.nan,rmin,np.nan
    r_soi=G*(float(MBH_msun)*Msun)/(sig*sig)
    x=r_soi/rmin
    if x<=1.0: return 0.0,r_soi,rmin,x
    return float(strength*(x-1.0)**power), r_soi,rmin,x

def mass_slope_penalty(theta,halo_type,R_star_m,strength=0.5):
    R=np.asarray(R_star_m,float); R=R[np.isfinite(R)]
    if R.size<20: return 0.0,np.nan,np.nan,np.nan,np.nan,np.nan
    r_in=float(np.quantile(R,0.10)); r_out=float(np.quantile(R,0.60))
    if not (np.isfinite(r_in) and np.isfinite(r_out) and r_out>r_in>0):
        return 0.0,np.nan,r_in,r_out,np.nan,np.nan
    Min,Mout=phys.mass_enclosed_two_radii_julia(r_in_m=r_in,r_out_m=r_out,theta=theta,halo_type=halo_type)
    Min=float(Min); Mout=float(Mout)
    if not (np.isfinite(Min) and np.isfinite(Mout)):
        return 0.0,np.nan,r_in,r_out,Min,Mout
    if Min<=0 or Mout<=0 or Mout<=Min:
        return float(strength*10.0), np.nan,r_in,r_out,Min,Mout
    alpha=float(np.log(Mout/Min)/np.log(r_out/r_in))
    lo=0.3; hi=2.8
    if alpha<lo: return float(strength*(lo-alpha)**2), alpha,r_in,r_out,Min,Mout
    if alpha>hi: return float(strength*(alpha-hi)**2), alpha,r_in,r_out,Min,Mout
    return 0.0,alpha,r_in,r_out,Min,Mout

def _chi2_from_A(A, *, b=None, w=None):
    A=np.asarray(A,float)
    if A.ndim!=2 or A.size==0: return None
    if not np.all(np.isfinite(A)): return None
    nrow,ncol=A.shape
    if b is None:
        b=np.zeros(nrow,float)
    else:
        b=np.asarray(b,float).ravel()
        if b.size!=nrow: return None
        if not np.all(np.isfinite(b)): return None
    if w is None:
        w=np.ones(nrow,float)
    else:
        w=np.asarray(w,float).ravel()
        if w.size!=nrow: return None
        if not np.all(np.isfinite(w)): return None
        w=np.where(w>0,w,0.0)
    sw=np.sqrt(w)
    Aw=A*sw[:,None]
    bw=b*sw
    col_norm=np.linalg.norm(Aw,axis=0)
    good=col_norm>0
    if not np.any(good): return None
    A2=Aw[:,good]/col_norm[good]
    b2=bw
    w2,_=nnls(A2,b2)
    wfull=np.zeros(ncol,float)
    wfull[good]=w2/col_norm[good]
    r=A@wfull-b
    if not np.all(np.isfinite(r)): return None
    return float(np.sum((sw*r)*(sw*r)))

def wrap_physics_engine(base_engine, *, obs, halo_type, config=None):
    cfg=config or {}
    s1=float(cfg.get("PEN_SPHERE_STRENGTH",2.0))
    p1=float(cfg.get("PEN_SPHERE_POWER",2.0))
    s2=float(cfg.get("PEN_SLOPE_STRENGTH",0.5))
    print_every=int(cfg.get("PRINT_EVERY",_PRINT_EVERY))

    R_star_m=getattr(obs,"R_star_m",None)
    v_star_mps=getattr(obs,"v_star_mps",None)
    verr_star_mps=getattr(obs,"verr_star_mps",None)
    if R_star_m is None or v_star_mps is None:
        R_star_m=getattr(obs,"R_m",None)
        v_star_mps=getattr(obs,"v_mps",None)
        verr_star_mps=getattr(obs,"verr_mps",None)
    if R_star_m is None or v_star_mps is None:
        raise AttributeError("obs must expose R_star_m and v_star_mps for penalties")

    v_star_mps=np.asarray(v_star_mps,float)
    if verr_star_mps is not None: verr_star_mps=np.asarray(verr_star_mps,float)

    def engine(theta):
        global _print_counter
        _print_counter+=1
        try: A=base_engine(theta,return_A=True)
        except TypeError: A=base_engine(theta)

        if (_print_counter%50)==0: print("[PHYS] base_engine type:",type(A))

        if isinstance(A,(float,int,np.floating,np.integer)):
            chi2=float(A)
        else:
            if A is None: return float("inf")
            A=np.asarray(A,float)
            if A.ndim!=2 or A.size==0: return float("inf")

            nrow=A.shape[0]
            b=v_star_mps[:nrow]
            if b.size!=nrow: return float("inf")

            w=None
            if verr_star_mps is not None and verr_star_mps.size>=nrow:
                ve=verr_star_mps[:nrow]
                ok=np.isfinite(ve)&(ve>0)
                if np.any(ok):
                    w=np.zeros(nrow,float)
                    w[ok]=1.0/(ve[ok]*ve[ok])

            chi2=_chi2_from_A(A,b=b,w=w)
            if chi2 is None: return float("inf")

        MBH=float(theta[2]) if len(theta)>2 else 0.0
        pen1,r_soi,rmin,ratio=chi2_resolution_penalty(MBH_msun=MBH,R_star_m=R_star_m,v_star_mps=v_star_mps,strength=s1,power=p1)
        pen2,alpha,r_in,r_out,Min,Mout=mass_slope_penalty(theta=theta,halo_type=halo_type,R_star_m=R_star_m,strength=s2)
        chi2_tot=float(chi2+pen1+pen2)

        if print_every>0 and (_print_counter%print_every==0):
            print(f"[PHYS] chi2={chi2:10.4f} +pen_res={pen1:9.4f} +pen_slope={pen2:9.4f} =>chi2_tot={chi2_tot:10.4f} MBH={MBH:9.3e}")
            if np.isfinite(ratio): print(f"[SOI]  r_soi={r_soi:9.3e} r_min={rmin:9.3e} ratio={ratio:6.2f}")
            if np.isfinite(alpha): print(f"[SLOPE] alpha={alpha:6.3f} r_in={r_in:9.3e} r_out={r_out:9.3e} Min={Min:9.3e} Mout={Mout:9.3e}")
            elif pen2>0 or (np.isfinite(Min) and np.isfinite(Mout)):
                print(f"[SLOPE] alpha=nan r_in={r_in:9.3e} r_out={r_out:9.3e} Min={Min:9.3e} Mout={Mout:9.3e}")

        return chi2_tot

    return engine