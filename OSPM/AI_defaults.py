# OSPM_Config_Center — Segue1
# Canonical, structurally complete galaxy configuration
# Matches Draco surface exactly. Values may be adjusted later.

from pathlib import Path


AI_DEFAULTS = {
    # Mode
    "MODE":      "stellar",
    "HALO_TYPE": "nfw",

    # Parallelism
    "N_WORKERS": 4,

    # Orbit / binning
    "NORBIT": 800, # default for local machine cluster will be 10000+
    "BINNING": {
        "MIN_BINS":         5,
        "N_TARGET_CIRC":    5,
        "MIN_PER_BIN_CIRC": 6,
    },

    # Physical domain
    "MIN_DISTANCE":             1e-6,
    "MAX_DISTANCE":             5e3,
    "R_GRID_POINTS":            256,
    "POTENTIAL_EXTENT":         10.0,
    "BH_MIN_RADIUS_MULTIPLIER": 2.0,

    # Penalties
    "PEN_SPHERE_STRENGTH": 200,
    "PEN_SPHERE_POWER":    2.0,
    "PEN_SLOPE_STRENGTH":  5000,

    # Deck semantics
    "REQUIRE_COLUMNS": [
        "rho_s","r_s","MBH","chi2","reward","status","proposal_id"
    ],
    "ALLOWED_STATUSES": [
        "todo","seed","pass","orbit_fail",
        "numeric_fail","unknown_fail","forbidden"
    ],
    "FILL_DEFAULT_STATUS": "todo",

    # Sampling / termination
    "BATCH_SIZE":            250,
    "MIN_BATCH_SIZE":        8,
    "MAX_BATCH_SIZE":        256,
    "MAX_RUNS":              100000,
    "STOP_NO_IMPROVEMENT":   1000,
    "IMPROVEMENT_EPSILON":   1e-6,
    "LOG_INTERVAL":          500,

    # AI control
    "AI_START_AFTER":       250,
    "MIN_TRAIN_POINTS":     800,
    "TRAIN_WINDOW":         2000,
    "AI_NOISE_INIT":        0.30,
    "AI_NOISE_MIN":         0.02,
    "AI_NOISE_TAU":         5000,
    "AI_MIN_DISTINCT_PASS": 800,
    "RESET_INTERVAL":       10000,
    "FLAT_WINDOW":          200,
    "FLAT_THRESHOLD":       1e-6,
    "FLAT_PATIENCE":        3,
    "AI_RESET_ON_FLAT":     False,

    # Constants
    "G":    6.67430e-11,
    "Msun": 1.98847e30,
}
