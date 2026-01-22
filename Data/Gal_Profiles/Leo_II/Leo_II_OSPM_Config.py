# OSPM_Config_Center — Leo_II
# Only place that should have Galaxy-specific configuration variables

from pathlib import Path
from Data.Data_Prep.Data_Paths import build_data_paths

Galaxy = "Leo_II"

OSPM_ROOT    = Path(__file__).resolve().parent
PROFILE_ROOT = OSPM_ROOT

CONFIG = {
    # =========================================================
    # Galaxy geometry (declared, never fitted)
    # =========================================================
    "RA0_DEG":          168.3708,
    "DEC0_DEG":         22.1517,
    "DISTANCE_PC":      233000.0,

    "PA_DEG":           13.0,
    "AXIS_RATIO_Q":     1.0,
    "R_HALF_LIGHT_PC":  176,
    "R_MAX_STARS_PC":   880.0,

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

    "BINNING": {
        "MIN_BINS":         5,
        "N_TARGET_CIRC":    5,
        "MIN_PER_BIN_CIRC": 6,
    },

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
}
