# Data/Get_All_Gal_Data.py
# Run data-prep for every galaxy profile

import subprocess
import os
import sys
from pathlib import Path

GAL_ROOT = Path("Data/Gal_Profiles")
SKIP = {"Draco", "Segue1", "segue1", "Segue1_v3"}   # already have v3 data / not a real galaxy

def main(run_label="default"):
    galaxies = sorted(d.name for d in GAL_ROOT.iterdir() if d.is_dir())

    for gal in galaxies:
        if gal in SKIP:
            print(f"[skip] {gal}")
            continue

        print(f"[run] {gal}")

        env = dict(os.environ)
        env["OSPM_GALAXY"] = gal

        result = subprocess.run(
            ["python", "-m", "Data.Data_Prep.DATA_NEW_GAL", run_label],
                env=env,
        )
        if result.returncode != 0:
            print(f"[FAIL] {gal}")

if __name__ == "__main__":
    label = sys.argv[1] if len(sys.argv) > 1 else "default"
    main(run_label=label)

