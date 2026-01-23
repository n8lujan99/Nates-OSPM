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
    "INITIAL_THETA":   [1.0, 500.0, 0.0],
    "THETA_BOUNDS": [
        (1e-1, 1e6),
        (20.0, 1e5),
        (1.0, 2e6),
    ],

    # Data-generation metadata
    "RADIUS_DEG":  0.6,
    "RUWE_MAX":    1.4,
    "PAR_SNR_MIN": 5.0,

    # Column authority
    "STAR_R_COL":    "r_pc",
    "STAR_V_COL":    "radial_velocity",
    "STAR_VERR_COL": "radial_velocity_error",
    "RA_COL":        "RA_deg",
    "DEC_COL":       "Dec_deg",
    "VLOS_COL":      "radial_velocity",

    # Paths
    **build_data_paths(PROFILE_ROOT),
}