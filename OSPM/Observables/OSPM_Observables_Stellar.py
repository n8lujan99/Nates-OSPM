"""
OSPM_Observables_Stellar
Star-level observable container for OSPM.
Holds per-star data only. No binning. No physics.

Invariants:
- one row per star
- units fixed at construction
- mode == "stellar"
"""

import numpy as np
from ..physics.OSPM_Physics import pc, kms, make_inclination


class OSPMObservablesStellar:
    def __init__(
        self,
        *,
        R_star_pc,
        v_star_kms,
        verr_star_kms,
        inclination_deg,
        Norbit,
    ):
        self.mode = "stellar"

        R_star_pc     = np.asarray(R_star_pc, float)
        v_star_kms    = np.asarray(v_star_kms, float)
        verr_star_kms = np.asarray(verr_star_kms, float)

        if not (len(R_star_pc) == len(v_star_kms) == len(verr_star_kms)):
            raise ValueError("Star arrays must have equal length")

        valid = (
            np.isfinite(R_star_pc)
            & np.isfinite(v_star_kms)
            & np.isfinite(verr_star_kms)
            & (R_star_pc > 0)
            & (verr_star_kms > 0)
        )

        if not np.any(valid):
            raise RuntimeError("No valid stars after filtering")

        self.R_star_pc     = R_star_pc[valid]
        self.v_star_kms    = v_star_kms[valid]
        self.verr_star_kms = verr_star_kms[valid]

        self.R_star_m     = self.R_star_pc * pc
        self.v_star_mps   = self.v_star_kms * kms
        self.verr_star_mps = self.verr_star_kms * kms

        self.sini, self.cosi, self.edge_on = make_inclination(inclination_deg)

        self.Norbit = int(Norbit)
        self.Nstar  = len(self.R_star_m)

        # Occupancy likelihood control (used only if Julia returns stacked A)
        self.Nocc        = 6
        self.lambda_occ  = 0.3

    @classmethod
    def from_star_table(
        cls,
        csv_path,
        *,
        r_col="r_pc",
        v_col="vlos",
        verr_col="vlos_err",
        inclination_deg,
        Norbit,
    ):
        import pandas as pd

        df = pd.read_csv(csv_path)

        missing = [c for c in (r_col, v_col, verr_col) if c not in df.columns]
        if missing:
            raise KeyError(f"Missing required columns in star table: {missing}")

        return cls(
            R_star_pc     = df[r_col].values,
            v_star_kms    = df[v_col].values,
            verr_star_kms = df[verr_col].values,
            inclination_deg = inclination_deg,
            Norbit          = Norbit,
        )
