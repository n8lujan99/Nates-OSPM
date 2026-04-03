# OSPM_Config — Draco_v3
# v3 data (Walker+2023 registry spectroscopy, 2896 vlos stars)

from pathlib import Path
import os
import multiprocessing as mp
from Data.Data_Prep.Data_Paths import build_data_paths

Galaxy = "Draco_v3"

DRACO_ROOT   = Path(__file__).resolve().parents[1] / "Draco"
PROFILE_ROOT = DRACO_ROOT

def detect_workers():
    slurm = os.getenv("SLURM_CPUS_PER_TASK")
    if slurm and slurm.isdigit():
        return int(slurm)
    return mp.cpu_count()

WORKERS = detect_workers()

CONFIG = {
    "N_WORKERS": WORKERS,

    "MODE":        "stellar",
    "GALAXY":      "Draco",
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

    "NORBIT": 2500,

    "BINNING": {
        "MIN_BINS":            5,
        "N_TARGET_CIRC":       5,
        "MIN_PER_BIN_CIRC":    6,
    },

    "OBSERVABLES": {
        "NBINS_OCC": 6,
        "LAMBDA_OCC": 0.3,
    },

    "PARAMETER_NAMES": ["rho_s", "r_s", "MBH"],
    "INITIAL_THETA": [3, 1800.0, 9e5],
    "THETA_BOUNDS": [
        (2.0, 15.0),
        (1800.0, 8000.0),
        (5.5e5, 5e6),
    ],

    "PEN_SPHERE_STRENGTH": 200,
    "PEN_SPHERE_POWER":    2.0,
    "PEN_SLOPE_STRENGTH":  5000,

    "MIN_DISTANCE":             1e-6,
    "MAX_DISTANCE":             5e3,
    "R_GRID_POINTS":            256,
    "POTENTIAL_EXTENT":         10.0,
    "BH_MIN_RADIUS_MULTIPLIER": 2.0,

    "REQUIRE_COLUMNS": ["rho_s", "r_s", "MBH", "chi2", "reward", "status", "proposal_id"],
    "ALLOWED_STATUSES": ["todo", "seed", "pass", "orbit_fail", "numeric_fail", "unknown_fail", "forbidden"],
    "FILL_DEFAULT_STATUS": "todo",

    "BATCH_SIZE":          WORKERS,
    "MIN_BATCH_SIZE":      16,
    "MAX_BATCH_SIZE":      400,
    "_PRINT_EVERY":        1,
    "_print_counter":      0,

    "AI_START_AFTER":        200,
    "MIN_TRAIN_POINTS":      200,
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

    **build_data_paths(DRACO_ROOT, run_label="v3"),
    "CSV_PATH": str(DRACO_ROOT / "v3" / "daemon_deck_v3.csv"),
}
print("[CONFIG] Draco_v3 | NORBIT=2500 | data:", CONFIG["DATA_CSV"])
