"""
OSPM_Observables_Stellar
Star-level observable container for OSPM.
One row per star. No binning.
This file is the observational interpretation layer.
The raw catalog is not modified. The data file remains the measured input from
the source catalog. If the config declares a pressure-supported system, the raw
line-of-sight velocity column is passed directly into OSPM. If the config declares
a non-pressure-supported system, this file applies the configured galaxy-motion
correction in memory before building the velocity array used by the model.
This is not where gravitational physics is added. Things the stars physically
feel, such as the dark matter halo, stellar mass profile, or central black hole,
belong in the physics/potential layer. This file only decides how the observed
velocity column should be interpreted before the orbit-superposition fit sees it.
"""

import numpy as np
from ..Physics.OSPM_Physics import pc, kms, make_inclination
def _apply_motion_model(df, *, v_col, config=None):
    v_raw = np.asarray(df[v_col].values, float)
    v_model = np.zeros_like(v_raw, dtype=float)
    if config is None:
        return v_raw, v_raw, v_model
    dyn_mode = config.get("DYNAMICAL_MODE", "pressure_supported")
    if dyn_mode in ("pressure_supported", "spherical_pressure_supported"):
        return v_raw, v_raw, v_model
    if dyn_mode != "non_pressure_supported":
        raise ValueError(f"Unknown DYNAMICAL_MODE: {dyn_mode}")
    motion = config.get("MOTION_MODEL", None)
    if motion is None:
        raise KeyError("DYNAMICAL_MODE='non_pressure_supported' requires MOTION_MODEL")
    mode = motion.get("mode", None)
    if mode == "systemic":
        raw_v_col = motion.get("raw_v_col", v_col)
        if raw_v_col not in df.columns:
            raise KeyError(f"MOTION_MODEL raw_v_col not found in star table: {raw_v_col}")
        v_raw = np.asarray(df[raw_v_col].values, float)
        v_model = np.full_like(v_raw, float(motion["v_sys_kms"]), dtype=float)
        v_used = v_raw - v_model
        return v_used, v_raw, v_model
    raise ValueError(f"Unknown MOTION_MODEL mode: {mode}")

class OSPMObservablesStellar:
    def __init__(self, *, R_star_pc, v_star_kms, verr_star_kms, has_vlos=None, inclination_deg,
        Norbit, stellar_model=None, dynamical_mode=None, motion_model=None, v_star_raw_kms=None, v_motion_model_kms=None):
        self.mode = "stellar"
        self.dynamical_mode = dynamical_mode
        self.motion_model = motion_model
        if stellar_model is not None and not isinstance(stellar_model, dict):
            raise TypeError("stellar_model must be a dict or None")
        self.stellar_model = stellar_model
        R = np.asarray(R_star_pc, float)
        v = np.asarray(v_star_kms, float)
        ve = np.asarray(verr_star_kms, float)
        vraw = None if v_star_raw_kms is None else np.asarray(v_star_raw_kms, float)
        vmod = None if v_motion_model_kms is None else np.asarray(v_motion_model_kms, float)
        hv = (np.isfinite(v) & np.isfinite(ve)) if has_vlos is None else np.asarray(has_vlos, bool)
        if not (len(R) == len(v) == len(ve) == len(hv)):
            raise ValueError("Star arrays must have equal length")
        g = np.isfinite(R) & (R > 0)
        if not np.any(g):
            raise RuntimeError("No valid stars after geometric filtering")
        self.R_star_pc = R[g]
        self.v_star_kms = v[g]
        self.verr_star_kms = ve[g]
        self.has_vlos = hv[g]
        self.v_star_raw_kms = None if vraw is None else vraw[g]
        self.v_motion_model_kms = None if vmod is None else vmod[g]
        self.motion_applied = (
            dynamical_mode == "non_pressure_supported"
            and motion_model is not None
        )
        self.valid_vlos = (
            self.has_vlos
            & np.isfinite(self.v_star_kms)
            & np.isfinite(self.verr_star_kms)
            & (self.verr_star_kms > 0)
        )
        self.R_star_m = self.R_star_pc * pc
        self.v_star_mps = self.v_star_kms * kms
        self.verr_star_mps = self.verr_star_kms * kms
        self.sini, self.cosi, self.edge_on = make_inclination(inclination_deg)
        self.Norbit = int(Norbit)
        self.Nstar = len(self.R_star_m)
        self.Nstar_vlos = int(self.valid_vlos.sum())
        self.Nocc = 6
        self.lambda_occ = 0.3

    @classmethod
    def from_star_table( cls, csv_path, *, r_col="r_pc", v_col="vlos", verr_col="vlos_err",
        inclination_deg, Norbit, stellar_model=None, config=None,):
        import pandas as pd
        df = pd.read_csv(csv_path)
        needed = [r_col, v_col, verr_col]
        if config is not None and config.get("DYNAMICAL_MODE") == "non_pressure_supported":
            motion = config.get("MOTION_MODEL", {})
            raw_v_col = motion.get("raw_v_col", v_col)
            if raw_v_col not in needed:
                needed.append(raw_v_col)
        missing = [c for c in needed if c not in df.columns]
        if missing:
            raise KeyError(f"Missing required columns in star table: {missing}")
        v_used, v_raw, v_model = _apply_motion_model(df, v_col=v_col, config=config)
        return cls( R_star_pc=df[r_col].values, v_star_kms=v_used, verr_star_kms=df[verr_col].values, inclination_deg=inclination_deg,
            Norbit=Norbit, stellar_model=stellar_model, dynamical_mode=None if config is None else config.get("DYNAMICAL_MODE"),
            motion_model=None if config is None else config.get("MOTION_MODEL"), v_star_raw_kms=v_raw, v_motion_model_kms=v_model,
        )