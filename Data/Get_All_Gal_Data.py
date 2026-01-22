# Data/Get_All_Gal_Data.py
# Run data-prep for every galaxy profile

import subprocess
import os
from pathlib import Path

GAL_ROOT = Path("Data/Gal_Profiles")
SKIP = {"Draco"}   # keep your explicit skip

def main():
    galaxies = sorted(d.name for d in GAL_ROOT.iterdir() if d.is_dir())

    for gal in galaxies:
        if gal in SKIP:
            print(f"[skip] {gal}")
            continue

        print(f"[run] {gal}")

        env = dict(os.environ)
        env["OSPM_GALAXY"] = gal

        result = subprocess.run(
            ["python", "-m", "Data.Data_Prep.DATA_NEW_GAL"],
                env=env,
        )
        if result.returncode != 0:
            print(f"[FAIL] {gal}")

if __name__ == "__main__":
    main()

