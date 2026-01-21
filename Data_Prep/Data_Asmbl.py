#Data_Asmbl.py
import numpy as np, pandas as pd
import astropy.units as u
import astropy.coordinates as coord

FIELD={ "ra":["ra","RA","RA_ICRS","RAJ2000","RAdeg","RA_deg","RAJ2000_deg"],
        "dec":["dec","DEC","DE_ICRS","DEJ2000","DEdeg","Dec_deg","DEJ2000_deg"],
        "vlos":["vlos","RV","HRV","VHELIO_AVG","Vlos","Vlos_kms","radial_velocity","radial_velocity_kms"],
        "vlos_err":["vlos_err","e_RV","e_HRV","VERR","e_Vlos","radial_velocity_error","radial_velocity_err"],
        }

def _pick(df,names):
    for n in names:
        if n in df.columns: return n
    return None

def _to_float(x):
    return pd.to_numeric(x,errors="coerce")

def _try_sexagesimal(ra,dec):
    sc=coord.SkyCoord(ra.astype(str),dec.astype(str),unit=(u.hourangle,u.deg),frame="icrs")
    return sc.ra.deg,sc.dec.deg

def assemble_source(df,*,src,v_sign=+1.0,parse_sexagesimal=True):
    if df is None or len(df)==0: return pd.DataFrame()
    df=df.copy()
    c_ra=_pick(df,FIELD["ra"]); c_dec=_pick(df,FIELD["dec"])
    c_v=_pick(df,FIELD["vlos"]); c_ev=_pick(df,FIELD["vlos_err"])
    out=pd.DataFrame()
    if c_ra is not None: out["ra"]=df[c_ra]
    if c_dec is not None: out["dec"]=df[c_dec]
    if c_v is not None: out["vlos"]=df[c_v]
    if c_ev is not None: out["vlos_err"]=df[c_ev]
    if "ra" not in out or "dec" not in out: return pd.DataFrame()
    ra=_to_float(out["ra"]); dec=_to_float(out["dec"])
    if parse_sexagesimal:
        bad=(~np.isfinite(ra.to_numpy()))|(~np.isfinite(dec.to_numpy()))
        if bad.any():
            try:
                ra2,dec2=_try_sexagesimal(out["ra"],out["dec"])
                ra=pd.Series(ra2,index=out.index); dec=pd.Series(dec2,index=out.index)
            except Exception:
                pass
    out["ra"]=ra; out["dec"]=dec
    if "vlos" in out: out["vlos"]=v_sign*_to_float(out["vlos"])
    if "vlos_err" in out: out["vlos_err"]=np.abs(_to_float(out["vlos_err"]))
    out["src"]=src
    keep=["ra","dec","vlos","vlos_err","src"]
    return out[[c for c in keep if c in out.columns]]
