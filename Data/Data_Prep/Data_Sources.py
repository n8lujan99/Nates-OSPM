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
from .Spec_Registry import SPEC_REGISTRY
import re

# ------------------------------------------------------------
def frac_err(v, e):
    return np.abs(e / v) if np.isfinite(v) and v != 0 else np.inf
# ------------------------------------------------------------
# Cone-search helpers (geometry injected)
# ------------------------------------------------------------
def src_gaia(*, ra0_deg, dec0_deg, radius_deg):
    Vizier.ROW_LIMIT = -1
    v = Vizier(columns=[  "Source", "RA_ICRS", "DE_ICRS", "Plx", "e_Plx",
        "pmRA", "e_pmRA", "pmDE", "e_pmDE", "RV", "e_RV", "RUWE"])
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
    # SDSS on Vizier only has redshifts (zsp), not stellar RVs.
    # DR16 (V/154/sdss16) is the same — no RV column.
    # Kept for completeness; returns empty for stellar work.
    return _vizier_cone(
        ["RA_ICRS", "DE_ICRS", "RV", "e_RV"],
        "V/154/sdss16",
        "sdss",
        ra0_deg=ra0_deg, dec0_deg=dec0_deg, radius_deg=radius_deg)


def src_lamost(*, ra0_deg, dec0_deg, radius_deg):
    # V/164/dr5 is the latest LAMOST release on Vizier (DR5, ~9M spectra)
    return _vizier_cone(
        ["RAJ2000", "DEJ2000", "HRV", "e_HRV"],
        "V/164/dr5",
        "lamost",
        ra0_deg=ra0_deg, dec0_deg=dec0_deg, radius_deg=radius_deg
    )


def src_apogee(*, ra0_deg, dec0_deg, radius_deg):
    # III/284/allstars — APOGEE-2 DR16 (473K stars, still latest on Vizier)
    return _vizier_cone(
        ["RA", "DEC", "VHELIO_AVG", "VERR"],
        "III/284/allstars",
        "apogee",
        ra0_deg=ra0_deg, dec0_deg=dec0_deg, radius_deg=radius_deg
    )


def src_galah(*, ra0_deg, dec0_deg, radius_deg):
    # GALAH DR4 (J/A+A/703/A104, 917K stars, 2025)
    return _vizier_cone(
        ["RAJ2000", "DEJ2000", "RV1", "e_RV1"],
        "J/A+A/703/A104/catalog",
        "galah",
        ra0_deg=ra0_deg, dec0_deg=dec0_deg, radius_deg=radius_deg
    )


def src_rave(*, ra0_deg, dec0_deg, radius_deg):
    # RAVE DR6 — final release
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
# Registry-based spectroscopy (dedicated studies per galaxy)
# ------------------------------------------------------------
def src_registry(galaxy_name, *, ra0_deg, dec0_deg, radius_deg):
    """Query all Vizier catalogs registered for *galaxy_name*."""
    entries = SPEC_REGISTRY.get(galaxy_name, [])
    if not entries:
        return pd.DataFrame()

    frames = []
    for spec in entries:
        try:
            cat = spec["catalog"]
            v = Vizier(columns=["**"], row_limit=-1)
            c = coord.SkyCoord(ra0_deg * u.deg, dec0_deg * u.deg)
            r = v.query_region(c, radius=radius_deg * u.deg, catalog=cat)
            if not r:
                print(f"[registry] {spec['ref']}: no data")
                continue

            df = r[0].to_pandas()

            # Filter by target pattern if needed
            tcol = spec.get("target_col")
            tpat = spec.get("target_pat")
            if tcol and tpat and tcol in df.columns:
                mask = df[tcol].astype(str).str.contains(tpat, flags=re.IGNORECASE, na=False)
                df = df[mask].reset_index(drop=True)
                if df.empty:
                    print(f"[registry] {spec['ref']}: 0 rows after target filter")
                    continue

            # Build output with explicit column mapping
            out = pd.DataFrame()
            ra_c = spec["ra_col"]
            dec_c = spec["dec_col"]
            v_c = spec["v_col"]
            ve_c = spec.get("v_err_col")

            # Try specified columns, fall back to Vizier computed _RA/_DE
            if ra_c in df.columns:
                out["ra"] = pd.to_numeric(df[ra_c], errors="coerce")
            elif "_RA" in df.columns:
                out["ra"] = pd.to_numeric(df["_RA"], errors="coerce")
            if dec_c in df.columns:
                out["dec"] = pd.to_numeric(df[dec_c], errors="coerce")
            elif "_DE" in df.columns:
                out["dec"] = pd.to_numeric(df["_DE"], errors="coerce")
            if v_c in df.columns:
                out["vlos"] = pd.to_numeric(df[v_c], errors="coerce")
            if ve_c and ve_c in df.columns:
                out["vlos_err"] = pd.to_numeric(df[ve_c], errors="coerce").abs()

            out["src"] = spec["ref"]
            out = out.dropna(subset=["ra", "dec"]).reset_index(drop=True)
            if not out.empty:
                print(f"[registry] {spec['ref']}: {len(out)} rows, "
                      f"{out['vlos'].notna().sum()} with vlos")
                frames.append(out)

        except Exception as e:
            print(f"[registry] {spec['ref']} failed: {e}")

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


# ------------------------------------------------------------
# Spectroscopy handling (survey-based)
# ------------------------------------------------------------
def build_spec_collapsed(spec_path, *, ra0_deg, dec0_deg, radius_deg, arcsec=1.0, scratch=False):
    dfs = [
        src_sdss, src_lamost, src_apogee,
        src_galah, src_rave
    ]
    rows = []
    for fn in dfs:
        try:
            d = fn(ra0_deg=ra0_deg, dec0_deg=dec0_deg, radius_deg=radius_deg)
            if d is not None and not d.empty:
                rows.append(d)
                print(f"[spec] {fn.__name__}: {len(d)} rows")
        except Exception as e:
            print(f"[spec] {fn.__name__} failed: {e}")

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


def src_local_spec(p, *, ra0_deg, dec0_deg, radius_deg, scratch=False, build_if_missing=True):
    if p is None:
        return pd.DataFrame()
    if not os.path.exists(p):
        if build_if_missing:
            return build_spec_collapsed(p, ra0_deg=ra0_deg, dec0_deg=dec0_deg,
                                        radius_deg=radius_deg, scratch=scratch)
        return pd.DataFrame()
    return assemble_source(pd.read_csv(p), src="local", parse_sexagesimal=True)


# ------------------------------------------------------------
# Final catalog builder (geometry injected here)
# ------------------------------------------------------------
def build_sources_catalog(*, PROFILE_ROOT, CONFIG, scratch=False, collapse=True, run_label="default"):

    paths = build_data_paths(PROFILE_ROOT, run_label=run_label)
    DATA_CSV  = paths["DATA_CSV"]
    SPEC_PATH = paths["SPEC_PATH"]

    ra0    = CONFIG["RA0_DEG"]
    dec0   = CONFIG["DEC0_DEG"]
    radius = CONFIG["RADIUS_DEG"]

    galaxy_name = CONFIG.get("GALAXY", "")

    # ------------------------------------------------------------
    # 1a. Registry-based spectroscopy (dedicated studies)
    # ------------------------------------------------------------
    registry_spec = src_registry(galaxy_name, ra0_deg=ra0, dec0_deg=dec0,
                                  radius_deg=radius)

    # ------------------------------------------------------------
    # 1b. Local spectroscopy (hand-curated, authoritative)
    # ------------------------------------------------------------
    spec = src_local_spec(SPEC_PATH, ra0_deg=ra0, dec0_deg=dec0,
                          radius_deg=radius, scratch=scratch)

    dfs = []
    has_spec = False
    if spec is not None and not spec.empty:
        dfs.append(spec)
        has_spec = True
    if registry_spec is not None and not registry_spec.empty:
        dfs.append(registry_spec)
        has_spec = True

    # ------------------------------------------------------------
    # 2. Load cone-search catalogs (surveys)
    # ------------------------------------------------------------
    dfs += [
        src_gaia(ra0_deg=ra0, dec0_deg=dec0, radius_deg=radius),
        src_sdss(ra0_deg=ra0, dec0_deg=dec0, radius_deg=radius),
        src_lamost(ra0_deg=ra0, dec0_deg=dec0, radius_deg=radius),
        src_apogee(ra0_deg=ra0, dec0_deg=dec0, radius_deg=radius),
        src_simbad(ra0_deg=ra0, dec0_deg=dec0, radius_deg=radius),
    ]

    dfs = [d for d in dfs if d is not None and not d.empty]
    if not dfs:
        return pd.DataFrame()

    raw = pd.concat(dfs, ignore_index=True, sort=False)

    # ------------------------------------------------------------
    # 3. If spectroscopy exists, suppress Gaia-only RV-less rows
    # ------------------------------------------------------------
    if has_spec and "vlos" in raw:
        raw = raw[
            ~(raw["src"].eq("gaia") & raw["vlos"].isna())
        ].reset_index(drop=True)

    # ------------------------------------------------------------
    # 4. Persist
    # ------------------------------------------------------------
    if not scratch:
        raw.to_csv(DATA_CSV, index=False)

    return raw
