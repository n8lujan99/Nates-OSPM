# Holds all the data source functions as they are structurally similar

import pandas as pd
import astropy.units as u
from astropy import coordinates as coord
from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astroquery.ned import Ned
from pkg_v3.OSPM_Config import CONFIG

RA0_DEG=CONFIG["RA0_DEG"]; DEC0_DEG=CONFIG["DEC0_DEG"]; RADIUS_DEG=CONFIG["RADIUS_DEG"]

# ---------------- Gaia (authoritative spine)
def src_gaia():
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=["Source","RA_ICRS","DE_ICRS","Plx","e_Plx", "pmRA","e_pmRA","pmDE","e_pmDE","RV","e_RV","RUWE"])
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c,radius=RADIUS_DEG*u.deg,catalog="I/355/gaiadr3")
    if not r: return pd.DataFrame()
    df=r[0].to_pandas().rename(columns={"Source":"id","RA_ICRS":"ra","DE_ICRS":"dec", "Plx":"parallax", "e_Plx":"parallax_err",
        "pmRA":"pmra","e_pmRA":"pmra_err", "pmDE":"pmdec","e_pmDE":"pmdec_err", "RV":"vlos","e_RV":"vlos_err"})
    df["src"]="gaia"
    return df

# ---------------- SIMBAD (verification only)
def src_simbad():
    s=Simbad(); s.add_votable_fields("ra","dec","otype")
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=s.query_region(c,radius=RADIUS_DEG*u.deg)
    if r is None: return pd.DataFrame()
    df=r.to_pandas().rename(columns={"RA":"ra","DEC":"dec"})
    df["src"]="simbad"
    return df

# ---------------- Spectroscopic surveys
def src_sdss():
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=["RA_ICRS","DE_ICRS","RV","e_RV"])
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c,radius = RADIUS_DEG*u.deg,catalog="V/147/sdss12")
    if not r: return pd.DataFrame()
    df=r[0].to_pandas().rename(columns={"RA_ICRS":"ra","DE_ICRS":"dec","RV":"vlos","e_RV":"vlos_err"})
    df["src"]="sdss"
    return df

def src_lamost():
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=["RAJ2000","DEJ2000","RV","e_RV"])
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c,radius = RADIUS_DEG*u.deg,catalog="V/164/lamost")
    if not r: return pd.DataFrame()
    df=r[0].to_pandas().rename(columns={ "RAJ2000":"ra","DEJ2000":"dec", "RV":"vlos","e_RV":"vlos_err"})
    df["src"]="lamost"
    return df

def src_apogee():
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=["RA","DEC","VHELIO_AVG","VERR"])
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c, radius = RADIUS_DEG*u.deg,catalog="III/284/allstars")
    if not r: return pd.DataFrame()
    df=r[0].to_pandas().rename(columns={"RA":"ra","DEC":"dec", "VHELIO_AVG":"vlos","VERR":"vlos_err"})
    df["src"]="apogee"
    return df

#############################################################################################
# ---------------- Additional auxiliary sources (photometry, background veto, etc.) -----------------
#############################################################################################

# -------- GALAH (optical, high-res, southern)
def src_galah():
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=["RA_ICRS","DE_ICRS","RV","e_RV"])
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c,radius=RADIUS_DEG*u.deg,catalog="III/263/galah")
    if not r: return pd.DataFrame()
    df=r[0].to_pandas().rename(columns={ "RA_ICRS":"ra","DE_ICRS":"dec", "RV":"vlos","e_RV":"vlos_err"})
    df["src"]="galah"; return df

# -------- RAVE (legacy southern RVs)
def src_rave():
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=["RA_ICRS","DE_ICRS","RV","e_RV"])
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c,radius=RADIUS_DEG*u.deg,catalog="III/272/ravedr6")
    if not r: return pd.DataFrame()
    df=r[0].to_pandas().rename(columns={ "RA_ICRS":"ra","DE_ICRS":"dec", "RV":"vlos","e_RV":"vlos_err"})
    df["src"]="rave"; return df

# -------- WEAVE (early public DRs; catalog id may evolve)
def src_weave():
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=["RAJ2000","DEJ2000","RV","e_RV"])
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c,radius=RADIUS_DEG*u.deg,catalog="J/MNRAS/WEAVE")  # placeholder DR
    if not r: return pd.DataFrame()
    df=r[0].to_pandas().rename(columns={ "RAJ2000":"ra","DEJ2000":"dec", "RV":"vlos","e_RV":"vlos_err"})
    df["src"]="weave"; return df

# -------- 4MOST (early public DRs; catalog id may evolve)
def src_4most():
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=["RAJ2000","DEJ2000","RV","e_RV"])
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c,radius=RADIUS_DEG*u.deg,catalog="J/A+A/4MOST")  # placeholder DR
    if not r: return pd.DataFrame()
    df=r[0].to_pandas().rename(columns={ "RAJ2000":"ra","DEJ2000":"dec", "RV":"vlos","e_RV":"vlos_err"})
    df["src"]="4most"; return df

# -------- DECaLS / Legacy Survey (photometry only; no RV)
def src_decals():
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=["RAJ2000","DEJ2000"])
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c,radius=RADIUS_DEG*u.deg,catalog="II/371/ls_dr10")
    if not r: return pd.DataFrame()
    df=r[0].to_pandas().rename(columns={ "RAJ2000":"ra","DEJ2000":"dec"})
    df["src"]="decals"; return df

# -------- Pan-STARRS1 (photometry only)
def src_ps1():
    Vizier.ROW_LIMIT=-1
    v=Vizier(columns=["RAJ2000","DEJ2000"])
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=v.query_region(c,radius=RADIUS_DEG*u.deg,catalog="II/349/ps1")
    if not r: return pd.DataFrame()
    df=r[0].to_pandas().rename(columns={"RAJ2000":"ra","DEJ2000":"dec"})
    df["src"]="ps1"; return df

# -------- NED (background galaxy veto)
def src_ned():
    c=coord.SkyCoord(RA0_DEG*u.deg,DEC0_DEG*u.deg)
    r=Ned.query_region(c,radius=RADIUS_DEG*u.deg)
    if r is None: return pd.DataFrame()
    df=r.to_pandas()
    if "RAJ2000" in df and "DEJ2000" in df: df=df.rename(columns={"RAJ2000":"ra","DEJ2000":"dec"})
    if "RA" in df and "DEC" in df: df=df.rename(columns={"RA":"ra","DEC":"dec"})
    if "HRV" in df:   df=df.rename(columns={"HRV":"vlos"})
    if "e_HRV" in df: df=df.rename(columns={"e_HRV":"vlos_err"})
    df["src"]="ned"
    keep=[c for c in ["ra","dec","vlos","vlos_err","src"] if c in df.columns]
    return df[keep]
