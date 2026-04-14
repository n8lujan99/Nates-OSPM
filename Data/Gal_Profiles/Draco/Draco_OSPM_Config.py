# OSPM_Config_Center — Draco
# !LOCAL VERSION! NOT FOR ARC OR CLUSTERS
# Just so I have the option to run local tests that wont affect the cluster data 

from pathlib import Path
import pathlib
import os
import multiprocessing as mp
from Data.Data_Prep.Data_Paths import build_data_paths

Galaxy = "Draco"
# Run mode toggle
LOCAL_DEBUG = False   # True = local responsive run, False = full run / HPC-style

PROFILE_ROOT = Path(__file__).resolve().parent
if not PROFILE_ROOT.exists():
    raise FileNotFoundError(f"PROFILE_ROOT does not exist: {PROFILE_ROOT}")

def detect_workers():
    slurm = os.getenv("SLURM_CPUS_PER_TASK")
    if slurm and slurm.isdigit():
        return int(slurm)
    return mp.cpu_count()

WORKERS = detect_workers()


# Mode-dependent knobs
NORBIT = 1000 if LOCAL_DEBUG else 2500
BATCH_SIZE = 40 if LOCAL_DEBUG else 90
MIN_BATCH_SIZE = 40 if LOCAL_DEBUG else 90
MAX_BATCH_SIZE = 120 if LOCAL_DEBUG else 270
CHUNK_SIZE = 30 if LOCAL_DEBUG else 90
LOG_INTERVAL = 1 if LOCAL_DEBUG else 10
PROF_EVERY = 2 if LOCAL_DEBUG else 20
EVAL_TIMEOUT_S = 20.0 if LOCAL_DEBUG else 600.0


CONFIG = {
    "N_WORKERS": WORKERS,

    "MODE":        "stellar",
    "GALAXY":      Galaxy,
    "HALO_TYPE":   "nfw",
    
    "RA0_DEG":          260.0517,
    "DEC0_DEG":         57.9153,
    "DISTANCE_PC":      76000.0,

    "PA_DEG":           90.0,
    "AXIS_RATIO_Q":     0.70,
    "R_HALF_LIGHT_PC":  221.0,
    "R_MAX_STARS_PC":   1500.0,

    "INCLINATION_DEG":  78.0,

    "RADIUS_DEG":   0.6,
    "RUWE_MAX":     1.4,
    "PAR_SNR_MIN":  5.0,

    "STAR_R_COL":      "r_pc",
    "STAR_V_COL":      "vlos",
    "STAR_VERR_COL":   "vlos_err",
    "RA_COL":          "ra",
    "DEC_COL":         "dec",
    "VLOS_COL":        "vlos",

    "NORBIT":              NORBIT, 
    "BINNING": {
        "MIN_BINS":            5,
        "N_TARGET_CIRC":       5,
        "MIN_PER_BIN_CIRC":    6,
    },
    
    "OBSERVABLES": {
    "NBINS_OCC": 6,
    "LAMBDA_OCC": 0.3,
    },

    # =========================================================
    # Parameter space — UPDATED (Draco tuned box)
    "PARAMETER_NAMES": ["rho_s", "r_s", "MBH", "ML"],
    "INITIAL_THETA": [1, 1800.0, 9e5, 2.0],
    "THETA_BOUNDS": [
       (0.01, 15.0),      # rho_s  → high-density basin only
       (0, 15000.0),    # r_s    → keep basin width
       (0, 5e6),    # MBH    → upper cluster region
       (0.2, 2.0),  # ML     → mass-to-light ratio
    ],
    "PEN_SPHERE_STRENGTH": 200,
    "PEN_SPHERE_POWER":    2.0,
    "PEN_SLOPE_STRENGTH":  5000,

    "MIN_DISTANCE":             1e-6,
    "MAX_DISTANCE":             5e3,
    "R_GRID_POINTS":            256,
    "POTENTIAL_EXTENT":         10.0,
    "BH_MIN_RADIUS_MULTIPLIER": 2.0,

    "REQUIRE_COLUMNS": ["rho_s", "r_s", "MBH", "ML", "chi2", "reward", "status", "proposal_id" ],
    "ALLOWED_STATUSES": [ "todo", "seed", "pass", "orbit_fail", "numeric_fail", "unknown_fail", "forbidden" ],
    "FILL_DEFAULT_STATUS": "todo",

    "BATCH_SIZE":          BATCH_SIZE,
    "MIN_BATCH_SIZE":      MIN_BATCH_SIZE,
    "MAX_BATCH_SIZE":      MAX_BATCH_SIZE,
    "CHUNK_SIZE":          CHUNK_SIZE,
    "_PRINT_EVERY":        10,
    "_print_counter":      0,

    "AI_START_AFTER":        400,
    "MIN_TRAIN_POINTS":      300,
    "TRAIN_WINDOW":          3000,
    "AI_NOISE_INIT":         0.30,
    "AI_NOISE_MIN":          0.02,
    "AI_NOISE_TAU":          8000,
    "AI_MIN_DISTINCT_PASS":  500,
    "RESET_INTERVAL":        10000,
    "AI_DEBUG_EVERY":        200,
    "AI_SNAPSHOT_EVERY":     2000,
    "FLAT_WINDOW":           300,
    "FLAT_THRESHOLD":        1e-6,
    "FLAT_PATIENCE":         5,
    "AI_RESET_ON_FLAT":      False,

    "MAX_RUNS":              100000,
    "STOP_NO_IMPROVEMENT":   5000,
    "IMPROVEMENT_EPSILON":   1e-6,
    "LOG_INTERVAL":          1,

    "G":    6.67430e-11,
    "Msun": 1.98847e30,
    **build_data_paths(PROFILE_ROOT),
    "CSV_PATH": str(PROFILE_ROOT / "default" / "daemon_deck_clusterM2500L.csv"), # This creates the new local deck this is the only difference
}
print("[CONFIG] CSV_PATH =", CONFIG["CSV_PATH"])
