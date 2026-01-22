# data_prep/DATA_NEW_GAL.py
# Galaxy data preparation only
# Running: python -m data_prep.DATA_NEW_GAL

from importlib import import_module

from .Data_Paths import build_data_paths, ensure_data_dirs
from .Data_Sources import build_sources_catalog
from .Data_Preprocess import (
    preprocess_stars_for_ospm,
    quality_mask_from_config,
    PREPROCESS_MODE,
)

# -------------------------------------------------------------------
# CONFIG LOADING
# -------------------------------------------------------------------
# This script runs ONE galaxy per environment.
# The active galaxy is defined by which config module is imported.

# Example expected module:
# Data.Gal_Profiles.Draco.Draco_OSPM_Config

CFG = import_module("Data.Gal_Profiles.Draco.Draco_OSPM_Config")
CONFIG = CFG.CONFIG
PROFILE_ROOT = CFG.PROFILE_ROOT


# -------------------------------------------------------------------
def main(*, run_label="default", write=True):

    # Paths
    paths = build_data_paths(PROFILE_ROOT, run_label=run_label)
    ensure_data_dirs(paths)

    # Raw catalog
    df_raw = build_sources_catalog(PROFILE_ROOT=PROFILE_ROOT, scratch=not write)
    print(f"[DATA] raw sources: {len(df_raw)} rows")

    # Quality mask
    qmask = (
        quality_mask_from_config(CONFIG)
        if PREPROCESS_MODE == "regular"
        else None
    )

    # Geometry + cuts
    df_clean = preprocess_stars_for_ospm(
        df_raw,
        ra_col="ra",
        dec_col="dec",
        v_col="vlos",
        v_err_col="vlos_err",
        ra0_deg=CONFIG["RA0_DEG"],
        dec0_deg=CONFIG["DEC0_DEG"],
        distance_pc=CONFIG["DISTANCE_PC"],
        r_max_pc=CONFIG.get("R_MAX_STARS_PC"),
        quality_mask=qmask,
        bin_mode=CONFIG.get("BIN_MODE", "circular"),
        pa_deg=CONFIG.get("PA_DEG", 0.0),
        axis_ratio_q=CONFIG.get("AXIS_RATIO_Q", 1.0),
    )

    print(f"[DATA] preprocessed stars: {len(df_clean)} rows")

    if write:
        df_clean.to_csv(paths["STAR_CSV"], index=False)
        print(f"[OK] wrote STAR_CSV -> {paths['STAR_CSV']}")

    return df_clean


if __name__ == "__main__":
    main(run_label="default", write=True)
