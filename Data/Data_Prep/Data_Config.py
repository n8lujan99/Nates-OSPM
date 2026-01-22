# Get_All_Gal_Data.py
# Utils/run_all_galaxies.py

import os
import importlib

Galaxy = os.environ.get("OSPM_GALAXY")
if Galaxy is None:
    raise RuntimeError(
        "OSPM_GALAXY not set. Example:\n"
        "  export OSPM_GALAXY=Carina"
    )

modname = f"Data.Gal_Profiles.{Galaxy}.{Galaxy}_OSPM_Config"

try:
    mod = importlib.import_module(modname)
except ModuleNotFoundError as e:
    raise RuntimeError(
        f"Could not load galaxy config module:\n"
        f"  {modname}"
    ) from e

CONFIG = mod.CONFIG
PROFILE_ROOT = mod.PROFILE_ROOT
