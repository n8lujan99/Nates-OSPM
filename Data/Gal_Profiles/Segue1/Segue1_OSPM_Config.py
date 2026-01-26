# OSPM_Config_Center — Segue1
# Only place that should have Galaxy-specific configuration variables

from pathlib import Path
from Data.Data_Prep.Data_Paths import build_data_paths

Galaxy = "Segue1"

OSPM_ROOT    = Path(__file__).resolve().parent
PROFILE_ROOT = OSPM_ROOT

CONFIG = {
    # Galaxy identity
    "GALAXY": Galaxy,

    # Geometry (declared, never fitted)
    "RA0_DEG":         151.7667,
    "DEC0_DEG":        16.0819,
    "DISTANCE_PC":     23000.0,
    "PA_DEG":          None,
    "AXIS_RATIO_Q":    1.0,
    "R_HALF_LIGHT_PC": 29.0,
    "R_MAX_STARS_PC":  145.0,
    "INCLINATION_DEG": 90.0,

    # Parameter space
    "PARAMETER_NAMES": ["rho_s", "r_s", "MBH"],
    "INITIAL_THETA":   [0.1, 30.0, 0.0],
    "THETA_BOUNDS": [
        (1e-3, 10.0),
        (1.0, 300),
        (0, 1e6),
    ],

    # Data-generation metadata
    "RADIUS_DEG":  0.6,
    "RUWE_MAX":    1.4,
    "PAR_SNR_MIN": 5.0,

    # Column authority
    "STAR_R_COL":    "r_pc",
    "STAR_V_COL":    "vlos",
    "STAR_VERR_COL": "vlos_err",
    "RA_COL":        "ra",
    "DEC_COL":       "dec",
    "VLOS_COL":      "vlos",

    # Paths
    **build_data_paths(PROFILE_ROOT),
}