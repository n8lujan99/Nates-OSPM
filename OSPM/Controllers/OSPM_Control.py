# OSPM_Control.py
# Control center for OSPM configuration
# Pure orchestration. No physics. No data logic.

import os
import datetime
import pandas as pd
from ..Gal_Registry import load_config

def get_active_galaxy():
    gal = os.environ.get("OSPM_GALAXY")
    if not gal:
        raise RuntimeError(
            "OSPM_GALAXY not set. Example:\n"
            "export OSPM_GALAXY=Segue1"
        )
    return gal

def load_active_config():
    galaxy = get_active_galaxy()
    return load_config(galaxy)

def ensure_deck(ctrl):
    path = ctrl["CSV_PATH"]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path):
        return
    row = dict(zip(ctrl["PARAMETER_NAMES"], ctrl["INITIAL_THETA"]))
    row.update({
        "chi2": float("inf"),
        "reward": None,
        "status": ctrl["FILL_DEFAULT_STATUS"],
    })
    pd.DataFrame([row]).to_csv(path, index=False)

def build_runtime(ctrl):
    ts  = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    pid = os.getpid()
    wid = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    rt = dict(ctrl)
    rt["RUN_ID"] = f"{ts}_pid{pid}" + (f"_w{wid}" if wid else "")
    rt["WORKER_ID"] = wid
    rt["RANDOM_SEED_EFFECTIVE"] = rt.get("RANDOM_SEED", 123456789) + wid
    return rt
