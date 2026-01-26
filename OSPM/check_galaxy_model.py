# check_galaxy_model.py
#
# Galaxy-level sanity check for OSPM.
# Uses the ACTIVE galaxy only.
#
# Verifies:
# - star.csv integrity
# - observables construction
# - single-theta solver evaluation

import inspect
import numpy as np
import pandas as pd

from OSPM.load_config import load_config, get_profile_root
from OSPM.Observables.OSPM_Observables_Stellar import OSPMObservablesStellar
from OSPM.Solvers.OSPM_Solver_Stellar import solve_ospm_theta_stellar
import OSPM.Solvers.OSPM_Solver_Stellar as S


def main():

    cfg = load_config()
    galaxy = cfg["GALAXY"]
    profile_root = get_profile_root()

    print(f"[CHECK] Galaxy = {galaxy}")
    print("[CHECK] solver file:", S.__file__)
    print("[CHECK] solve_ospm_theta_stellar starts at line:",
          inspect.getsourcelines(S.solve_ospm_theta_stellar)[1])

    star_csv = profile_root / "default" / "star.csv"

    if not star_csv.exists():
        raise RuntimeError("star.csv not found")

    star_df = pd.read_csv(star_csv)

    print("[CHECK] star.csv loaded")
    print("  Nstar =", len(star_df))
    print("  columns =", list(star_df.columns))

    for col in ("r_pc", "has_vlos"):
        if col not in star_df.columns:
            raise RuntimeError(f"Missing {col} column")

    n_has = int(star_df["has_vlos"].sum())
    n_fin = int(star_df["vlos"].notna().sum()) if "vlos" in star_df.columns else 0

    print(f"  has_vlos = {n_has}")
    print(f"  finite vlos = {n_fin}")

    # --------------------------------------------------
    # Build observables
    # --------------------------------------------------

    Norbit = min(8, cfg["NORBIT"])

    obs = OSPMObservablesStellar.from_star_table(
        star_csv,
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
    main()
