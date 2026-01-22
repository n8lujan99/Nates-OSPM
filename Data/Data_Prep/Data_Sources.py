# Data_Prep/Data_Sources.py
# Source catalog construction (Gaia + spectroscopy + databases)
# Geometry is injected at runtime. No solver-layer imports.

import os
import numpy as np
import pandas as pd

import astropy.units as u
import astropy.coordinates as coord

from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
from astroquery.ned import Ned

from .Data_Paths import build_data_paths
from .Data_Asmbl import assemble_source


# ------------------------------------------------------------
def frac_err(v, e):
    return np.abs(e / v) if np.isfinite(v) and v != 0 else np.inf


# ------------------------------------------------------------
# Cone-search helpers (geometry injected)
# ------------------------------------------------------------
def src_gaia(*, ra0_deg, dec0_deg, radius_deg):
    Vizier.ROW_LIMIT = -1
    v = Vizier(columns=[
        "Source", "RA_ICRS", "DE_ICRS",
        "Plx", "e_Plx",
        "pmRA", "e_pmRA",
        "pmDE", "e_pmDE",
        "RV", "e_RV",
        "RUWE"
    ])
    c = coord.SkyCoord(ra0_deg * u.deg, dec0_deg * u.deg)
    r = v.query_region(c, radius=radius_deg * u.deg, catalog="I/355/gaiadr3")
    if not r:
        return pd.DataFrame()

    df = r[0].to_pandas().rename(columns={"Source": "id"})
    out = assemble_source(df, src="gaia")

    if "id" in df:
        out["id"] = df["id"].to_numpy()
    if "RUWE" in df:
        out["ruwe"] = pd.to_numeric(df["RUWE"], errors="coerce").to_numpy()
    if "Plx" in df:
        out["parallax"] = pd.to_numeric(df["Plx"], errors="coerce").to_numpy()
        out["parallax_err"] = pd.to_numeric(df.get("e_Plx", np.nan), errors="coerce").to_numpy()
    if "pmRA" in df:
        out["pmra"] = pd.to_numeric(df["pmRA"], errors="coerce").to_numpy()
        out["pmra_err"] = pd.to_numeric(df.get("e_pmRA", np.nan), errors="coerce").to_numpy()
    if "pmDE" in df:
        out["pmdec"] = pd.to_numeric(df["pmDE"], errors="coerce").to_numpy()
        out["pmdec_err"] = pd.to_numeric(df.get("e_pmDE", np.nan), errors="coerce").to_numpy()

    return out


def src_simbad(*, ra0_deg, dec0_deg, radius_deg):
    s = Simbad()
    s.add_votable_fields("ra", "dec", "otype")
    c = coord.SkyCoord(ra0_deg * u.deg, dec0_deg * u.deg)
    r = s.query_region(c, radius=radius_deg * u.deg)
    if r is None:
        return pd.DataFrame()

    df = r.to_pandas()
    out = assemble_source(df, src="simbad", parse_sexagesimal=True)
    if "OTYPE" in df:
        out["otype"] = df["OTYPE"].astype(str).to_numpy()
    return out


def _vizier_cone(cols, catalog, src, *, ra0_deg, dec0_deg, radius_deg):
    Vizier.ROW_LIMIT = -1
    v = Vizier(columns=cols)
    c = coord.SkyCoord(ra0_deg * u.deg, dec0_deg * u.deg)
    r = v.query_region(c, radius=radius_deg * u.deg, catalog=catalog)
    if not r:
        return pd.DataFrame()
    return assemble_source(r[0].to_pandas(), src=src)


def src_sdss(*, ra0_deg, dec0_deg, radius_deg):
    return _vizier_cone(
        ["RA_ICRS", "DE_ICRS", "RV", "e_RV"],
        "V/147/sdss12",
        "sdss",
        ra0_deg=ra0_deg, dec0_deg=dec0_deg, radius_deg=radius_deg
    )


def src_lamost(*, ra0_deg, dec0_deg, radius_deg):
    return _vizier_cone(
        ["RAJ2000", "DEJ2000", "RV", "e_RV"],
        "V/164/lamost",
        "lamost",
        ra0_deg=ra0_deg, dec0_deg=dec0_deg, radius_deg=radius_deg
    )


def src_apogee(*, ra0_deg, dec0_deg, radius_deg):
    return _vizier_cone(
        ["RA", "DEC", "VHELIO_AVG", "VERR"],
        "III/284/allstars",
        "apogee",
        ra0_deg=ra0_deg, dec0_deg=dec0_deg, radius_deg=radius_deg
    )


def src_galah(*, ra0_deg, dec0_deg, radius_deg):
    return _vizier_cone(
        ["RA_ICRS", "DE_ICRS", "RV", "e_RV"],
        "III/263/galah",
        "galah",
        ra0_deg=ra0_deg, dec0_deg=dec0_deg, radius_deg=radius_deg
    )


def src_rave(*, ra0_deg, dec0_deg, radius_deg):
    return _vizier_cone(
        ["RA_ICRS", "DE_ICRS", "RV", "e_RV"],
        "III/272/ravedr6",
        "rave",
        ra0_deg=ra0_deg, dec0_deg=dec0_deg, radius_deg=radius_deg
    )


def src_ned(*, ra0_deg, dec0_deg, radius_deg):
    c = coord.SkyCoord(ra0_deg * u.deg, dec0_deg * u.deg)
    r = Ned.query_region(c, radius=radius_deg * u.deg)
    if r is None:
        return pd.DataFrame()
    return assemble_source(r.to_pandas(), src="ned", parse_sexagesimal=True)


# ------------------------------------------------------------
# Spectroscopy handling (unchanged)
# ------------------------------------------------------------
def build_spec_collapsed(spec_path, *, arcsec=1.0, scratch=False):
    dfs = [
        src_sdss, src_lamost, src_apogee,
        src_galah, src_rave
    ]
    rows = []
    for fn in dfs:
        try:
            d = fn()
            if d is not None and not d.empty:
                rows.append(d)
        except Exception:
            pass

    if not rows:
        return pd.DataFrame()

    df = pd.concat(rows, ignore_index=True, sort=False)
    if "ra" not in df or "dec" not in df:
        return pd.DataFrame()

    df["ra"] = pd.to_numeric(df["ra"], errors="coerce")
    df["dec"] = pd.to_numeric(df["dec"], errors="coerce")
    for c in ["vlos", "vlos_err"]:
        if c in df:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    df = df[np.isfinite(df["ra"]) & np.isfinite(df["dec"])].reset_index(drop=True)
    if df.empty:
        return pd.DataFrame()

    c2 = coord.SkyCoord(df["ra"].to_numpy() * u.deg, df["dec"].to_numpy() * u.deg)
    used = np.zeros(len(df), bool)

    out_rows = []
    for i in range(len(df)):
        if used[i]:
            continue
        idx = np.where(c2[i].separation(c2).arcsec <= arcsec)[0]
        used[idx] = True
        g = df.iloc[idx]

        best = None
        for _, r in g.iterrows():
            if pd.notna(r.get("vlos")) and pd.notna(r.get("vlos_err")):
                fe = frac_err(r["vlos"], r["vlos_err"])
                if best is None or fe < best[0]:
                    best = (fe, r)

        if best:
            r = best[1]
            out_rows.append({
                "ra": float(r["ra"]),
                "dec": float(r["dec"]),
                "vlos": float(r["vlos"]),
                "vlos_err": float(r["vlos_err"]),
                "src": str(r.get("src", "spec")),
            })

    out = pd.DataFrame(out_rows)
    if (not scratch) and spec_path is not None:
        os.makedirs(os.path.dirname(spec_path), exist_ok=True)
        out.to_csv(spec_path, index=False)

    return out


def src_local_spec(p, *, scratch=False, build_if_missing=True):
    if p is None:
        return pd.DataFrame()
    if not os.path.exists(p):
        if build_if_missing:
            return build_spec_collapsed(p, scratch=scratch)
        return pd.DataFrame()
    return assemble_source(pd.read_csv(p), src="local", parse_sexagesimal=True)


# ------------------------------------------------------------
# Final catalog builder (geometry injected here)
# ------------------------------------------------------------
def build_sources_catalog(*, PROFILE_ROOT, CONFIG, scratch=False, collapse=True):

    paths = build_data_paths(PROFILE_ROOT)
    DATA_CSV = paths["DATA_CSV"]
    SPEC_PATH = paths["SPEC_PATH"]

    ra0 = CONFIG["RA0_DEG"]
    dec0 = CONFIG["DEC0_DEG"]
    radius = CONFIG["RADIUS_DEG"]

    dfs = [
        src_gaia(ra0_deg=ra0, dec0_deg=dec0, radius_deg=radius),
        src_sdss(ra0_deg=ra0, dec0_deg=dec0, radius_deg=radius),
        src_lamost(ra0_deg=ra0, dec0_deg=dec0, radius_deg=radius),
        src_apogee(ra0_deg=ra0, dec0_deg=dec0, radius_deg=radius),
        src_local_spec(SPEC_PATH, scratch=scratch),
        src_simbad(ra0_deg=ra0, dec0_deg=dec0, radius_deg=radius),
    ]

    dfs = [d for d in dfs if d is not None and not d.empty]
    raw = pd.concat(dfs, ignore_index=True, sort=False) if dfs else pd.DataFrame()

    if raw.empty:
        return raw

    if not scratch:
        raw.to_csv(DATA_CSV, index=False)

    return raw
