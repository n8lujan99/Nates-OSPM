from .AI_defaults import AI_DEFAULTS
from Data.Gal_Profiles.Segue1 import Segue1_OSPM_Config
from Data.Gal_Profiles.Draco  import Draco_OSPM_Config
from Data.Gal_Profiles.Carina import Carina_OSPM_Config

_GALAXY_MAP = {
    "Segue1": Segue1_OSPM_Config,
    "Draco":  Draco_OSPM_Config,
    "Carina": Carina_OSPM_Config,
}

def load_config(galaxy: str):
    if galaxy not in _GALAXY_MAP:
        raise KeyError(f"Unknown galaxy: {galaxy}")

    cfg = {**AI_DEFAULTS, **_GALAXY_MAP[galaxy].CONFIG}

    required = [
        "GALAXY","MODE","HALO_TYPE",
        "MIN_DISTANCE","MAX_DISTANCE",
        "NORBIT","BATCH_SIZE","MAX_RUNS",
    ]
    missing = [k for k in required if k not in cfg]
    if missing:
        raise KeyError(f"CONFIG missing required keys: {missing}")

    return cfg
