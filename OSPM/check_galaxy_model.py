# check_galaxy_model.py
#
# Galaxy-level sanity check for OSPM.
# Verifies:
# - star.csv integrity
# - observables construction
# - Julia A-matrix build
# - single-theta solver evaluation
#
# This is a REQUIRED preflight for every galaxy.
import inspect
import OSPM.Solvers.OSPM_Solver_Stellar as S

print("[CHECK] solver file:", S.__file__)
print("[CHECK] solve_ospm_theta_stellar starts at line:",
      inspect.getsourcelines(S.solve_ospm_theta_stellar)[1])

import sys
import numpy as np
import pandas as pd

from OSPM.Gal_Registry import load_galaxy
from OSPM.Observables.OSPM_Observables_Stellar import OSPMObservablesStellar
from OSPM.Solvers.OSPM_Solver_Stellar import solve_ospm_theta_stellar



def main(galaxy):

    print(f"[CHECK] Galaxy = {galaxy}")

    gal = load_galaxy(galaxy)

    if gal.get("error"):
        raise RuntimeError(f"Failed to load galaxy: {gal['error']}")

    star_df = gal["stars"]
    cfg     = gal["config"]

    if star_df is None:
        raise RuntimeError("star.csv not found")

    print("[CHECK] star.csv loaded")
    print("  Nstar =", len(star_df))
    print("  columns =", list(star_df.columns))

    if "r_pc" not in star_df.columns:
        raise RuntimeError("Missing r_pc column")

    if "has_vlos" not in star_df.columns:
        raise RuntimeError("Missing has_vlos column")

    n_has = int(star_df["has_vlos"].sum())
    n_fin = int(star_df["vlos"].notna().sum()) if "vlos" in star_df.columns else 0

    print(f"  has_vlos = {n_has}")
    print(f"  finite vlos = {n_fin}")

    # --------------------------------------------------
    # Build observables
    # --------------------------------------------------

    Norbit = 8

    obs = OSPMObservablesStellar.from_star_table(
        gal["profile_root"] / "default" / "star.csv",
        inclination_deg=cfg["INCLINATION_DEG"],
        Norbit=Norbit,
    )

    print("[CHECK] Observables built")
    print("  mode =", obs.mode)
    print("  Nstar =", obs.Nstar)
    print("  Nstar_vlos =", obs.Nstar_vlos)

    if obs.Nstar == 0:
        raise RuntimeError("No valid stars after geometry filter")

    # --------------------------------------------------
    # Single-theta solver test
    # --------------------------------------------------

    theta = cfg["INITIAL_THETA"]
    halo_type = cfg.get("HALO_TYPE", "nfw")

    print("[CHECK] Solver test")
    print("  theta =", theta)
    print("  halo_type =", halo_type)

    chi2, w, meta = solve_ospm_theta_stellar(theta, obs, halo_type=halo_type)

    if not np.isfinite(chi2):
        raise RuntimeError("chi2 is not finite")

    if not np.isfinite(w).all():
        raise RuntimeError("weights contain NaN or Inf")

    ws = w.sum()
    if abs(ws - 1.0) > 1e-6:
        raise RuntimeError(f"weights not normalized: sum={ws}")

    print("[CHECK] Solver OK")
    print("  chi2 =", chi2)
    print("  weights sum =", ws)
    print("  meta =", meta)

    print(f"[CHECK] {galaxy} PASSED")


if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("Usage: python check_galaxy_model.py <GalaxyName>")
        sys.exit(1)

    main(sys.argv[1])
