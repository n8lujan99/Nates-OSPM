# Control center for OSPM configuration
# STAYS in PYTHON
import importlib.util
import pathlib
import os
import pandas as pd
import datetime
from .OSPM_Config import CONFIG 
# NO Variables defined here (except SEED but that should stay the same)

def ensure_deck(ctrl):
    path = ctrl["CSV_PATH"]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    if os.path.exists(path):
        return
    row = dict(zip(ctrl["PARAMETER_NAMES"], ctrl["INITIAL_THETA"]))
    row.update({ "chi2": float("inf"), "reward": None, "status": ctrl["FILL_DEFAULT_STATUS"]}) 
    pd.DataFrame([row]).to_csv(path, index=False)
def load_config(galaxy: str):
    if galaxy != CONFIG.get("GALAXY"):
        raise RuntimeError(f"No config for galaxy '{galaxy}'")
    return CONFIG
def build_runtime(ctrl):
    rt = dict(ctrl)
    ts  = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    pid = os.getpid()
    wid = int(os.environ.get("SLURM_ARRAY_TASK_ID", "0"))
    rt["RUN_ID"] = f"{ts}_pid{pid}" + (f"_w{wid}" if wid else "")
    rt["WORKER_ID"] = wid
    rt["RANDOM_SEED_EFFECTIVE"] = rt.get("RANDOM_SEED", 123456789) + wid
    return rt