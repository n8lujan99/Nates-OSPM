# OSPM_Config_Center — Draco
# Only place that should have Galaxy-specific configuration variables
# that way other modules can remain Galaxy-agnostic
from pathlib import Path
import pathlib
from data_prep.Data_Paths import build_data_paths


Galaxy = "Draco"


# Paths need to be specified in Data_Paths.py sticking to the following convention:
OSPM_ROOT    = Path(__file__).resolve().parent
PROFILE_ROOT = OSPM_ROOT / "data" / "profiles" / Galaxy


CONFIG = {
    # =========================================================
    # Parallelization 
    # =========================================================
    "N_WORKERS": 4, # Equal to the number of CPU cores available 
    # Local machine: 4 cores
    # HPC: ## cores

    # =========================================================
    # Identity
    # =========================================================
    "MODE":        "stellar",
    "GALAXY":      Galaxy,
    "HALO_TYPE":   "nfw",
    
    # =========================================================
    # Galaxy geometry (declared, never fitted)
    # =========================================================
    "RA0_DEG":          260.0517,
    "DEC0_DEG":         57.9153,
    "DISTANCE_PC":      76000.0,

    # Morphology (literature / policy)
    "PA_DEG":           90.0,
    "AXIS_RATIO_Q":     0.70,
    "R_HALF_LIGHT_PC":  221.0,
    "R_MAX_STARS_PC":   1500.0,

    # Viewing geometry
    "INCLINATION_DEG":  78.0,

    # =========================================================
    # Data harvesting & quality
    # =========================================================
    "RADIUS_DEG":   0.6,
    "RUWE_MAX":     1.4,
    "PAR_SNR_MIN":  5.0,

    # Column authority
    "STAR_R_COL":      "r_pc",
    "STAR_V_COL":      "radial_velocity",
    "STAR_VERR_COL":   "radial_velocity_error",
    "RA_COL":          "RA_deg",
    "DEC_COL":         "Dec_deg",
    "VLOS_COL":        "radial_velocity",

    # =========================================================
    # OSPM numerical setup
    # =========================================================
    "NORBIT":              1000, # of orbits to sample per evaluation (Karl's default: 10,000) 10,000 is a lot and bogging down the system. 
    "BINNING": {
        "MIN_BINS":            5,
        "N_TARGET_CIRC":       5,
        "MIN_PER_BIN_CIRC":    6,
    },

    # =========================================================
    # Parameter space
    # =========================================================
    "PARAMETER_NAMES": ["rho_s", "r_s", "MBH"],
    "INITIAL_THETA":   [1.0, 500.0, 0.0],   
    "THETA_BOUNDS": [
        (1e-1, 1e6),     # rho_s
        (200, 1e5),      # r_s
        (1.0, 2e6),      # MBH
    ],

    # Penalties
    "PEN_SPHERE_STRENGTH": 200,
    "PEN_SPHERE_POWER":    2.0,
    "PEN_SLOPE_STRENGTH":  5000,

    # =========================================================
    # Physical domain (solver only)
    # =========================================================
    "MIN_DISTANCE":             1e-6,
    "MAX_DISTANCE":             5e3,
    "R_GRID_POINTS":            256,
    "POTENTIAL_EXTENT":         10.0,
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
    "BATCH_SIZE":          250,
    "MIN_BATCH_SIZE":      8,
    "MAX_BATCH_SIZE":      256,
    "_PRINT_EVERY":        200,
    "_print_counter":      0,

    # =========================================================
    # AI / learning
    # =========================================================
    "AI_START_AFTER":        250,
    "MIN_TRAIN_POINTS":      800,
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
    "FLAT_PATIENCE":         3,
    "AI_RESET_ON_FLAT":      False,

    # =========================================================
    # Termination
    # =========================================================
    "MAX_RUNS":              100000,
    "STOP_NO_IMPROVEMENT":   1000,
    "IMPROVEMENT_EPSILON":   1e-6,
    "LOG_INTERVAL":          500,

    # =========================================================
    # Constants
    # =========================================================
    "G":    6.67430e-11,
    "Msun": 1.98847e30,
    **build_data_paths(PROFILE_ROOT),
}
