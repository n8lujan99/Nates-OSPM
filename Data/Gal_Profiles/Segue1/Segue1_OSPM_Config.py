# OSPM_Config_Center — Segue1
# Only place that should have Galaxy-specific configuration variables
# so all other modules remain galaxy-agnostic

from pathlib import Path
from Data.Data_Prep.Data_Paths import build_data_paths

Galaxy = "Segue1"

# ---------------------------------------------------------
# Profile root (this directory contains default/, data.csv, etc.)
# ---------------------------------------------------------
PROFILE_ROOT = Path(__file__).resolve().parent
if not PROFILE_ROOT.exists():
    raise FileNotFoundError(f"PROFILE_ROOT does not exist: {PROFILE_ROOT}")

CONFIG = {
    # =========================================================
    # Parallelization
    # =========================================================
    "N_WORKERS": 6,   # local default; override via launcher on HPC

    # =========================================================
    # Identity
    # =========================================================
    "MODE":        "stellar",
    "GALAXY":      Galaxy,
    "HALO_TYPE":   "nfw",

    # =========================================================
    # Galaxy geometry (declared, never fitted)
    # =========================================================
    "RA0_DEG":          151.7667,
    "DEC0_DEG":         16.0819,
    "DISTANCE_PC":      23000.0,

    # Morphology (Segue 1 treated as spherical)
    "PA_DEG":           None,
    "AXIS_RATIO_Q":     1.0,
    "R_HALF_LIGHT_PC":  29.0,
    "R_MAX_STARS_PC":   145.0,

    # Viewing geometry
    "INCLINATION_DEG":  90.0,

    # =========================================================
    # Data harvesting & quality
    # =========================================================
    "RADIUS_DEG":   0.6,
    "RUWE_MAX":     1.4,
    "PAR_SNR_MIN":  5.0,

    # Column authority (Segue 1 conventions)
    "STAR_R_COL":      "r_pc",
    "STAR_V_COL":      "vlos",
    "STAR_VERR_COL":   "vlos_err",
    "RA_COL":          "ra",
    "DEC_COL":         "dec",
    "VLOS_COL":        "vlos",

    # =========================================================
    # OSPM numerical setup
    # =========================================================

    "NORBIT":              1000,   # RESOLUTION
    "BINNING": {
        "MIN_BINS":            5,
        "N_TARGET_CIRC":       5,
        "MIN_PER_BIN_CIRC":    6,
    },

    # =========================================================
    # Parameter space
    # =========================================================
    "PARAMETER_NAMES": ["rho_s", "r_s", "MBH"],
    "INITIAL_THETA":   [0.1, 2.0, 50000],
    "THETA_BOUNDS": [
        (1e-4, 100.0),     # rho_s
        (5, 5000.0),    # r_s
        (0.0, 1e6),      # MBH
    ],

    # Penalties
    "PEN_SPHERE_STRENGTH": 200,
    "PEN_SPHERE_POWER":    2.0,
    "PEN_SLOPE_STRENGTH":  5000,

    # =========================================================
    # Physical domain (solver only)
    # =========================================================
    "MIN_DISTANCE":             5e-4,
    "MAX_DISTANCE":             2e3,
    "R_GRID_POINTS":            256,
    "POTENTIAL_EXTENT":         6.0,
    "BH_MIN_RADIUS_MULTIPLIER": 2.0,

    # =========================================================
    # Deck semantics
    # =========================================================
    "REQUIRE_COLUMNS": ["rho_s", "r_s", "MBH", "chi2", "reward", "status", "proposal_id"],
    "ALLOWED_STATUSES": [ "todo", "seed", "pass", "orbit_fail", "numeric_fail", "unknown_fail", "forbidden" ],
    "FILL_DEFAULT_STATUS": "todo",

    # =========================================================
    # Sampling & control
    # =========================================================
    "BATCH_SIZE":          64,
    "MIN_BATCH_SIZE":      32,
    "MAX_BATCH_SIZE":      256,
    "_PRINT_EVERY":        10,
    "_print_counter":      4,

    # =========================================================
    # AI / learning
    # =========================================================
    "AI_START_AFTER":        500,
    "MIN_TRAIN_POINTS":      500,
    "TRAIN_WINDOW":          2000,
    "AI_NOISE_INIT":         0.30,
    "AI_NOISE_MIN":          0.02,
    "AI_NOISE_TAU":          5000,
    "AI_MIN_DISTINCT_PASS":  800,
    "RESET_INTERVAL":        10000,
    "AI_DEBUG_EVERY":        200,
    "AI_SNAPSHOT_EVERY":     2000,
    "FLAT_WINDOW":           200,
    "FLAT_THRESHOLD":        1e-6,
    "FLAT_PATIENCE":         10,
    "AI_RESET_ON_FLAT":      True,

    # =========================================================
    # Termination
    # =========================================================
    "MAX_RUNS":              100000,
    "STOP_NO_IMPROVEMENT":   1000,
    "IMPROVEMENT_EPSILON":   1e-6,
    "LOG_INTERVAL":          10,

    # =========================================================
    # Constants
    # =========================================================
    "G":    6.67430e-11,
    "Msun": 1.98847e30,

    # =========================================================
    # Paths (authoritative)
    # =========================================================
    **build_data_paths(PROFILE_ROOT),
}
