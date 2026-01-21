# data_prep/Data_Sources.py
# SAME CODE YOU HAD.
# ONLY CHANGE: paths resolved at runtime, not import time.
# NOTHING REMOVED. NOTHING RENAMED.

import os
import numpy as np, pandas as pd
import astropy.units as u
import astropy.coordinates as coord
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astroquery.ned import Ned
from pkg_v3.OSPM_Config import CONFIG
from .Data_Paths import build_data_paths
from .Data_Asmbl import assemble_source

GALAXY=CONFIG["GALAXY"]
RA0_DEG=CONFIG["RA0_DEG"]; DEC0_DEG=CONFIG["DEC0_DEG"]; RADIUS_DEG=CONFIG["RADIUS_DEG"]

def frac_err(v,e): return np.abs(e/v) if np.isfinite(v) and v!=0 else np.inf

def src_gaia():
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=["Source","RA_ICRS","DE_ICRS","Plx","e_Plx","pmRA","e_pmRA","pmDE","e_pmDE","RV","e_RV","RUWE"])
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c,radius=RADIUS_DEG*u.deg,catalog="I/355/gaiadr3")
    if not r: return pd.DataFrame()
    df=r[0].to_pandas().rename(columns={"Source":"id"})
    out=assemble_source(df,src="gaia")
    if "id" in df: out["id"]=df["id"].to_numpy()
    if "RUWE" in df: out["ruwe"]=pd.to_numeric(df["RUWE"],errors="coerce").to_numpy()
    if "Plx" in df:
        out["parallax"]=pd.to_numeric(df["Plx"],errors="coerce").to_numpy()
        out["parallax_err"]=pd.to_numeric(df.get("e_Plx",np.nan),errors="coerce").to_numpy()
    if "pmRA" in df:
        out["pmra"]=pd.to_numeric(df["pmRA"],errors="coerce").to_numpy()
        out["pmra_err"]=pd.to_numeric(df.get("e_pmRA",np.nan),errors="coerce").to_numpy()
    if "pmDE" in df:
        out["pmdec"]=pd.to_numeric(df["pmDE"],errors="coerce").to_numpy()
        out["pmdec_err"]=pd.to_numeric(df.get("e_pmDE",np.nan),errors="coerce").to_numpy()
    return out

def src_simbad():
    s=Simbad(); s.add_votable_fields("ra","dec","otype")
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=s.query_region(c,radius=RADIUS_DEG*u.deg)
    if r is None: return pd.DataFrame()
    df=r.to_pandas()
    out=assemble_source(df,src="simbad",parse_sexagesimal=True)
    if "OTYPE" in df: out["otype"]=df["OTYPE"].astype(str).to_numpy()
    return out

def _vizier_cone(cols,catalog,src):
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=cols)
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c,radius=RADIUS_DEG*u.deg,catalog=catalog)
    if not r: return pd.DataFrame()
    return assemble_source(r[0].to_pandas(),src=src)

def src_sdss():   return _vizier_cone(["RA_ICRS","DE_ICRS","RV","e_RV"],"V/147/sdss12","sdss")
def src_lamost(): return _vizier_cone(["RAJ2000","DEJ2000","RV","e_RV"],"V/164/lamost","lamost")
def src_apogee(): return _vizier_cone(["RA","DEC","VHELIO_AVG","VERR"],"III/284/allstars","apogee")
def src_galah():  return _vizier_cone(["RA_ICRS","DE_ICRS","RV","e_RV"],"III/263/galah","galah")
def src_rave():   return _vizier_cone(["RA_ICRS","DE_ICRS","RV","e_RV"],"III/272/ravedr6","rave")

def src_weave():
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=["RAJ2000","DEJ2000","RV","e_RV"])
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c,radius=RADIUS_DEG*u.deg,catalog="J/MNRAS/WEAVE")
    if not r: return pd.DataFrame()
    return assemble_source(r[0].to_pandas(),src="weave")

def src_4most():
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=["RAJ2000","DEJ2000","RV","e_RV"])
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c,radius=RADIUS_DEG*u.deg,catalog="J/A+A/4MOST")
    if not r: return pd.DataFrame()
    return assemble_source(r[0].to_pandas(),src="4most")

def src_ned():
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=Ned.query_region(c,radius=RADIUS_DEG*u.deg)
    if r is None: return pd.DataFrame()
    return assemble_source(r.to_pandas(),src="ned",parse_sexagesimal=True)

def build_spec_collapsed(spec_path,*,arcsec=1.0,scratch=False):
    dfs=[src_sdss(),src_lamost(),src_apogee(),src_galah(),src_rave(),src_weave(),src_4most()]
    dfs=[d for d in dfs if (d is not None and not d.empty)]
    if not dfs: return pd.DataFrame()
    df=pd.concat(dfs,ignore_index=True,sort=False)
    if "ra" not in df or "dec" not in df: return pd.DataFrame()
    df["ra"]=pd.to_numeric(df["ra"],errors="coerce")
    df["dec"]=pd.to_numeric(df["dec"],errors="coerce")
    for c in ["vlos","vlos_err"]:
        if c in df: df[c]=pd.to_numeric(df[c],errors="coerce")
    df=df[np.isfinite(df["ra"])&np.isfinite(df["dec"])].reset_index(drop=True)
    if df.empty: return pd.DataFrame()
    c2=coord.SkyCoord(df["ra"].to_numpy()*u.deg,df["dec"].to_numpy()*u.deg)
    used=np.zeros(len(df),bool); rows=[]
    for i in range(len(df)):
        if used[i]: continue
        idx=np.where(c2[i].separation(c2).arcsec<=arcsec)[0]
        used[idx]=True
        g=df.iloc[idx]
        best=None
        for _,r in g.iterrows():
            if pd.notna(r.get("vlos")) and pd.notna(r.get("vlos_err")):
                fe=frac_err(r["vlos"],r["vlos_err"])
                if best is None or fe<best[0]: best=(fe,r)
        if best:
            r=best[1]
            rows.append(dict(ra=float(r["ra"]),dec=float(r["dec"]),
                             vlos=float(r["vlos"]),vlos_err=float(r["vlos_err"]),
                             src=str(r.get("src","spec"))))
    out=pd.DataFrame(rows)
    if (not scratch) and spec_path is not None:
        os.makedirs(os.path.dirname(spec_path),exist_ok=True)
        out.to_csv(spec_path,index=False)
    return out

def src_local_spec(p,*,scratch=False,build_if_missing=True):
    if p is None: return pd.DataFrame()
    if not os.path.exists(p):
        if build_if_missing:
            return build_spec_collapsed(p,scratch=scratch)
        return pd.DataFrame()
    return assemble_source(pd.read_csv(p),src="local",parse_sexagesimal=True)

def fuse_sources(dfs,*,arcsec=1.0):
    base=dfs[0].copy()
    if base.empty: raise RuntimeError("Gaia returned no stars")
    base["ra"]=pd.to_numeric(base["ra"],errors="coerce")
    base["dec"]=pd.to_numeric(base["dec"],errors="coerce")
    base=base[np.isfinite(base["ra"])&np.isfinite(base["dec"])].copy()
    base["_k"]=np.arange(len(base))
    c1=coord.SkyCoord(base["ra"].to_numpy()*u.deg,base["dec"].to_numpy()*u.deg)
    for df in dfs[1:]:
        if df is None or df.empty: continue
        if "ra" not in df or "dec" not in df: continue
        d=df.copy()
        d["ra"]=pd.to_numeric(d["ra"],errors="coerce")
        d["dec"]=pd.to_numeric(d["dec"],errors="coerce")
        d=d[np.isfinite(d["ra"])&np.isfinite(d["dec"])]
        if d.empty: continue
        c2=coord.SkyCoord(d["ra"].to_numpy()*u.deg,d["dec"].to_numpy()*u.deg)
        i,sep,_=c1.match_to_catalog_sky(c2)
        m=sep.arcsec<=arcsec
        if not np.any(m): continue
        j=d.iloc[i[m]].copy()
        j["_k"]=base.loc[m,"_k"].to_numpy()
        base=base.merge(j,on="_k",how="left",suffixes=("", "_x"))
    rows=[]
    for _,g in base.groupby("_k",sort=False):
        best=None
        for _,r in g.iterrows():
            if pd.notna(r.get("vlos")) and pd.notna(r.get("vlos_err")):
                fe=frac_err(r["vlos"],r["vlos_err"])
                if best is None or fe<best[0]: best=(fe,r["vlos"],r["vlos_err"],r.get("src"))
        if best:
            rows.append(dict(ra=g["ra"].iloc[0],dec=g["dec"].iloc[0],
                             vlos=best[1],vlos_err=best[2],vlos_src=best[3]))
    return pd.DataFrame(rows)

def build_sources_catalog(*, PROFILE_ROOT, scratch=False, collapse=True):
    paths=build_data_paths(PROFILE_ROOT)
    DATA_CSV=paths["DATA_CSV"]; SPEC_PATH=paths["SPEC_PATH"]
    dfs=[
        src_gaia(),src_sdss(),src_lamost(),src_apogee(),
        src_local_spec(SPEC_PATH,scratch=scratch),src_simbad()
    ]
    dfs=[d for d in dfs if d is not None and not d.empty]
    raw=pd.concat(dfs,ignore_index=True,sort=False)
    fused=fuse_sources([raw]) if collapse else raw
    if not scratch: fused.to_csv(DATA_CSV,index=False)
    return fused
