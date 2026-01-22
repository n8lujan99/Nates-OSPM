# Data/Tools/generate_galaxy_profiles.py
# Authoritative generator for galaxy geometry + OSPM config

from pathlib import Path
import pandas as pd
import math

# ------------------------------------------------------------------
# Paths (repo-true)
# ------------------------------------------------------------------
DATA_ROOT = Path(__file__).resolve().parents[1]      # .../Data
GAL_FILE  = DATA_ROOT / "Galaxies.txt"
GAL_ROOT  = DATA_ROOT / "Gal_Profiles"

def is_nan(x):
    return x != x
def get_axis_ratio(row):
    return getattr(row, "axis_ratio_q", 1.0)

# ------------------------------------------------------------------
# center.txt
# ------------------------------------------------------------------
def write_center_txt(gdir, row):
    path = gdir / "center.txt"
    with open(path, "w") as f:
        f.write(f"RA0_DEG {row.ra_deg}\n")
        f.write(f"DEC0_DEG {row.dec_deg}\n")
        f.write(f"DISTANCE_PC {row.distance_kpc * 1000.0}\n")
        f.write(f"R_HALF_LIGHT_PC {row.r_half_light_pc}\n")
        q = get_axis_ratio(row)
        f.write(f"AXIS_RATIO_Q {q}\n")
        if not is_nan(row.pa_deg):
            f.write(f"PA_DEG {row.pa_deg}\n")

# ------------------------------------------------------------------
# Galaxy-specific OSPM config
# ------------------------------------------------------------------
def write_config_py(gdir, row):
    gname = row.name
    cfg_path = gdir / f"{gname}_OSPM_Config.py"
    pa_val = "None" if is_nan(row.pa_deg) else row.pa_deg
    rmax   = 5.0 * row.r_half_light_pc

    with open(cfg_path, "w") as f:
        f.write(f'''# OSPM_Config_Center — {gname}
# Only place that should have Galaxy-specific configuration variables

from pathlib import Path
from Data.Data_Prep.Data_Paths import build_data_paths

Galaxy = "{gname}"

OSPM_ROOT    = Path(__file__).resolve().parent
PROFILE_ROOT = OSPM_ROOT

CONFIG = {{
    # =========================================================
    # Galaxy geometry (declared, never fitted)
    # =========================================================
    "RA0_DEG":          {row.ra_deg},
    "DEC0_DEG":         {row.dec_deg},
    "DISTANCE_PC":      {row.distance_kpc * 1000.0},

    "PA_DEG":           {pa_val},
    "AXIS_RATIO_Q":     {get_axis_ratio(row)},
    "R_HALF_LIGHT_PC":  {row.r_half_light_pc},
    "R_MAX_STARS_PC":   {rmax},

    "INCLINATION_DEG":  90.0,

    # =========================================================
    # Parameter space
    # =========================================================
    "PARAMETER_NAMES": ["rho_s", "r_s", "MBH"],
    "INITIAL_THETA":   [1.0, 500.0, 0.0],
    "THETA_BOUNDS": [
        (1e-1, 1e6),
        (200, 1e5),
        (1.0, 2e6),
    ],

    ################################################################
    # THE FOLLOWING IS ONLY FOR DATA GENERATION — DO NOT CHANGE
    ################################################################
    "RADIUS_DEG":   0.6,
    "RUWE_MAX":     1.4,
    "PAR_SNR_MIN":  5.0,

    "STAR_R_COL":    "r_pc",
    "STAR_V_COL":    "radial_velocity",
    "STAR_VERR_COL": "radial_velocity_error",
    "RA_COL":        "RA_deg",
    "DEC_COL":       "Dec_deg",
    "VLOS_COL":      "radial_velocity",

    ################################################################
    # GLOBAL — IDENTICAL FOR ALL GALAXIES
    ################################################################
    "MODE":      "stellar",
    "GALAXY":    Galaxy,
    "HALO_TYPE": "nfw",

    "NORBIT": 1000,

    "BINNING": {{
        "MIN_BINS":         5,
        "N_TARGET_CIRC":    5,
        "MIN_PER_BIN_CIRC": 6,
    }},

    "PEN_SPHERE_STRENGTH": 200,
    "PEN_SPHERE_POWER":   2.0,
    "PEN_SLOPE_STRENGTH": 5000,

    "REQUIRE_COLUMNS": ["rho_s","r_s","MBH","chi2","reward","status","proposal_id"],
    "ALLOWED_STATUSES": ["todo","seed","pass","orbit_fail","numeric_fail","unknown_fail","forbidden"],
    "FILL_DEFAULT_STATUS": "todo",

    "BATCH_SIZE": 250,
    "MAX_RUNS":   100000,

    "N_WORKERS": 4,

    "G":    6.67430e-11,
    "Msun": 1.98847e30,

    **build_data_paths(PROFILE_ROOT),
}}
''')

# ------------------------------------------------------------------
def main():
    df = pd.read_csv(GAL_FILE)
    for row in df.itertuples(index=False):
        if row.name == "Draco":
            print("[skip] Draco (explicitly skipped)")
            continue
        gdir = GAL_ROOT / row.name
        if not gdir.exists():
            print(f"[skip] {row.name} (no directory)")
            continue
        write_center_txt(gdir, row)
        write_config_py(gdir, row)
        print(f"[ok] {row.name}")
if __name__ == "__main__":
    main()
