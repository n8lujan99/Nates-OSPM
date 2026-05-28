# OSPM_Config — Segue1_v3
# Probe run using v3 data (registry-sourced spectroscopy).
# Reduced NORBIT for fast structural exploration.

from pathlib import Path
from Data.Data_Prep.Data_Paths import build_data_paths
from OSPM.AI_defaults import CONFIG as AI_DEFAULTS

Galaxy = "Segue1_v3"

# Data lives in the main Segue1 profile under v3/
SEGUE1_ROOT  = Path(__file__).resolve().parents[1] / "Segue1"
PROFILE_ROOT = SEGUE1_ROOT     # build_data_paths reads from here

CONFIG = {
    # =========================================================
    # Parallelization
    # =========================================================
    "N_WORKERS": 80,

    # =========================================================
    # Identity
    # =========================================================
    "MODE":        "stellar",
    "GALAXY":      "Segue1",   # physics identity stays Segue1
    "HALO_TYPE":   "nfw",

    # =========================================================
    # Galaxy geometry (declared, never fitted)
    # =========================================================
    "RA0_DEG":          151.7667,
    "DEC0_DEG":         16.0819,
    "DISTANCE_PC":      23000.0,

    "PA_DEG":           None,
    "AXIS_RATIO_Q":     1.0,
    "R_HALF_LIGHT_PC":  29.0,
    "R_MAX_STARS_PC":   120.0,

    "INCLINATION_DEG":  90.0,

    # =========================================================
    # Data harvesting & quality
    # =========================================================
    "RADIUS_DEG":   0.6,
    "RUWE_MAX":     1.4,
    "PAR_SNR_MIN":  5.0,

    "STAR_R_COL":      "r_pc",
    "STAR_V_COL":      "vlos",
    "STAR_VERR_COL":   "vlos_err",
    "RA_COL":          "ra",
    "DEC_COL":         "dec",
    "VLOS_COL":        "vlos",

    # =========================================================
    # OSPM numerical setup — reduced for probe run
    # =========================================================
    "NORBIT": 2500,

    "BINNING": {
        "MIN_BINS":            4,
        "N_TARGET_CIRC":       5,
        "MIN_PER_BIN_CIRC":    3,
    },
    "OBSERVABLES": {
        "NBINS_OCC": 6,
        "LAMBDA_OCC": 0.3,
    },

    # =========================================================
    # Parameter space (same as production)
    # =========================================================
    "PARAMETER_NAMES": ["rho_s", "r_s", "MBH", "ML"],
    "INITIAL_THETA":   [37.37, 1759.0, 1.2e5, 2],
    "THETA_BOUNDS": [
        (0.0, 100.0),     # rho_s
        (300, 5000.0),    # r_s
        (0.0, 1.5e6),      # MBH
        (0.2, 2.0),       # ML
    ],

    "PEN_SPHERE_STRENGTH": 2500,
    "PEN_SPHERE_POWER":    2.0,
    "PEN_SLOPE_STRENGTH":  5000,

    # =========================================================
    # Physical domain
    # =========================================================
    "MIN_DISTANCE":             5e-4,
    "MAX_DISTANCE":             2e3,
    "R_GRID_POINTS":            256,
    "POTENTIAL_EXTENT":         6.0,
    "BH_MIN_RADIUS_MULTIPLIER": 2.0,

    # =========================================================
    # Deck semantics
    # =========================================================
    "REQUIRE_COLUMNS": [ "rho_s", "r_s", "MBH", "ML", "chi2", "chi2_inner", "chi2_outer", "chi2_occ", "N_inner", "N_outer", 
                        "N_nonzero_weights", "effective_N_orbits", "max_weight_fraction", "reward", "status", "proposal_id", "refine_passes"],
    "ALLOWED_STATUSES": ["todo", "seed", "pass", "orbit_fail", "numeric_fail", "unknown_fail", "forbidden"],
    "FILL_DEFAULT_STATUS": "todo",

    # =========================================================
    # Sampling & control
    # =========================================================
    "BATCH_SIZE":          80, # 80
    "MIN_BATCH_SIZE":      80, # 80
    "MAX_BATCH_SIZE":      256,
    "_PRINT_EVERY":        10,
    "_print_counter":      1,
    "R_INNER_DIAG_PC":     30.0,
    # =========================================================
    # AI / learning
    # =========================================================
    "AI_START_AFTER":        600,
    "EXPLORE_FRACTION":      0.25,
    "MIN_TRAIN_POINTS":      300,
    "TRAIN_WINDOW":          500,
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
    "MAX_RUNS":              10000, # 10000
    "STOP_NO_IMPROVEMENT":   1000,
    "IMPROVEMENT_EPSILON":   1e-6,
    "LOG_INTERVAL":          10,

    # =========================================================
    # Constants
    # =========================================================
    "G":    6.67430e-11,
    "Msun": 1.98847e30,

    # =========================================================
    # Paths — read v3 Walker-cleaned data, write to v3 deck
    # =========================================================
    **build_data_paths(SEGUE1_ROOT, run_label="v3"),
    "DATA_CSV": str(SEGUE1_ROOT / "v3" / "segue1_walker2023_clean.csv"),
    "CSV_PATH": str(SEGUE1_ROOT / "v3" / "daemon_deck_v3_diag.csv"),
}
print("[CONFIG] Segue1_v3 probe run | NORBIT=2500 | data:", CONFIG["DATA_CSV"])