# halobounds.py
#
# Generic helper for galaxy-specific NFW halo proposal bounds.
# This is NOT a grader / veto after the fact.
# It defines the admissible halo region once, then provides:
#   - sample_pair(): direct sampling from that region
#   - project_pair(): deterministic projection of a proposed pair back onto that region
#
# First-pass rule:
#   1) The NFW turnover scale may not hide below the inner resolved stellar scale
#      -> r_s >= r_cut
#   2) The halo cusp strength may not exceed the densest halo at that smallest resolved scale
#      -> mu = rho_s * r_s <= mu_max = rho_raw_hi * r_s_min
#      -> rho_s <= min(rho_raw_hi, mu_max / r_s)
#
# This gives a simple curved admissible region:
#   rho_s_max(r_s) = min(raw rho upper bound, mu_max / r_s)
#
# The geometric mean of the first two stellar radii is used so the rule is tied
# to the inner resolved interval rather than to a single tracer.

from __future__ import annotations
import math
import os
from dataclasses import dataclass
from typing import Optional, Tuple
import numpy as np
import pandas as pd


def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))

@dataclass
class HaloBounds:
    
    galaxy: str
    halo_type: str
    rho_lo_raw: float
    rho_hi_raw: float
    rs_lo_raw: float
    rs_hi_raw: float
    r1: float
    r2: float
    r_cut: float
    rs_min: float
    rs_max: float
    mu_max: float
    radius_path: Optional[str] = None
    radius_col: str = "r_pc"

    @classmethod
    def from_config(cls, config: dict) -> "HaloBounds":
        halo_type = str(config.get("HALO_TYPE", "nfw")).lower()
        if halo_type != "nfw":
            raise ValueError(f"HaloBounds only supports HALO_TYPE='nfw', got {halo_type!r}")

        galaxy = str(config.get("GALAXY", "unknown"))
        params = list(config["PARAMETER_NAMES"])
        bounds = list(config["THETA_BOUNDS"])

        try:
            i_rho = params.index("rho_s")
            i_rs = params.index("r_s")
        except ValueError as e:
            raise KeyError("PARAMETER_NAMES must include 'rho_s' and 'r_s'") from e

        rho_lo_raw, rho_hi_raw = map(float, bounds[i_rho])
        rs_lo_raw, rs_hi_raw = map(float, bounds[i_rs])

        radius_col = str(config.get("STAR_R_COL", "r_pc"))
        radius_path = _resolve_radius_path(config)

        r1, r2 = _load_inner_tracers(
            config=config,
            radius_path=radius_path,
            radius_col=radius_col,
        )

        # Inner resolved scale from first two tracers
        # Geometric mean is less arbitrary than a hand-tuned weighted average.
        r_cut = math.sqrt(r1 * r2)

        # Halo turnover must stay outside the inner resolved interval
        rs_min = max(rs_lo_raw, r_cut)
        rs_max = rs_hi_raw
        
        # First-pass curved cap:
        # the strongest allowed inner cusp is the densest halo at the smallest resolved scale
        # allowed by the raw config.
        mu_max = rho_hi_raw * rs_min

        hb = cls(
            galaxy=galaxy,
            halo_type=halo_type,
            rho_lo_raw=rho_lo_raw,
            rho_hi_raw=rho_hi_raw,
            rs_lo_raw=rs_lo_raw,
            rs_hi_raw=rs_hi_raw,
            r1=r1,
            r2=r2,
            r_cut=r_cut,
            rs_min=rs_min,
            rs_max=rs_max,
            mu_max=mu_max,
            radius_path=radius_path,
            radius_col=radius_col,
        )
        hb._sanity_check()
        return hb

    def _sanity_check(self) -> None:
        if not np.isfinite(self.r1) or not np.isfinite(self.r2):
            raise ValueError("Inner tracer radii are not finite")
        if self.r1 <= 0 or self.r2 <= 0:
            raise ValueError("Inner tracer radii must be positive")
        if self.rs_min <= 0:
            raise ValueError("rs_min must be positive")
        if self.rs_max <= self.rs_min:
            raise ValueError(
                f"Invalid r_s range after halo bounds construction: "
                f"rs_min={self.rs_min}, rs_max={self.rs_max}"
            )
        if self.rho_hi_raw <= self.rho_lo_raw:
            raise ValueError(
                f"Invalid rho_s raw range: rho_lo_raw={self.rho_lo_raw}, rho_hi_raw={self.rho_hi_raw}"
            )

    def rho_max(self, r_s: float) -> float:
        """Curved upper boundary for rho_s at a given r_s."""
        r_s = max(float(r_s), self.rs_min)
        return min(self.rho_hi_raw, self.mu_max / r_s)

    def valid_pair(self, rho_s: float, r_s: float) -> bool:
        """Membership check for the admissible region. Mostly for debugging."""
        if r_s < self.rs_min or r_s > self.rs_max:
            return False
        if rho_s < self.rho_lo_raw:
            return False
        if rho_s > self.rho_max(r_s):
            return False
        return True

    def sample_pair(self, rng=np.random) -> Tuple[float, float]:
        """
        Sample directly from the admissible NFW region.
        We draw r_s first, then rho_s under the curve.
        """
        r_s = float(rng.uniform(self.rs_min, self.rs_max))
        rho_hi = self.rho_max(r_s)
        rho_s = float(rng.uniform(self.rho_lo_raw, rho_hi))
        return rho_s, r_s

    def project_pair(self, rho_s: float, r_s: float) -> Tuple[float, float]:
        """
        Deterministically map a raw proposal back into the admissible region.
        This is for AI proposals, so the policy is not ignored.
        """
        r_s = clamp(float(r_s), self.rs_min, self.rs_max)
        rho_hi = self.rho_max(r_s)
        rho_s = clamp(float(rho_s), self.rho_lo_raw, rho_hi)
        return rho_s, r_s

    def describe(self) -> dict:
        return {
            "galaxy": self.galaxy,
            "halo_type": self.halo_type,
            "radius_path": self.radius_path,
            "radius_col": self.radius_col,
            "inner_tracer_1_pc": self.r1,
            "inner_tracer_2_pc": self.r2,
            "r_cut_pc": self.r_cut,
            "rs_min_pc": self.rs_min,
            "rs_max_pc": self.rs_max,
            "rho_lo_raw": self.rho_lo_raw,
            "rho_hi_raw": self.rho_hi_raw,
            "mu_max": self.mu_max,
        }


def _resolve_radius_path(config: dict) -> Optional[str]:
    """
    Try a few likely config keys first.
    Then fall back to the profile directory inferred from CSV_PATH.
    """
    candidates = [
        "STAR_PATH",
        "STARS_PATH",
        "STAR_FILE",
        "DATA_PATH",
        "OBS_PATH",
        "MEMBER_PATH",
    ]

    for key in candidates:
        p = config.get(key)
        if isinstance(p, str) and os.path.isfile(p):
            return p

    csv_path = config.get("CSV_PATH")
    if isinstance(csv_path, str) and csv_path:
        profile_dir = os.path.dirname(csv_path)
        guessed = [
            os.path.join(profile_dir, "data.csv"),
            os.path.join(profile_dir, "star.csv"),
            os.path.join(profile_dir, "spec.csv"),
        ]
        for p in guessed:
            if os.path.isfile(p):
                return p

        galaxy_dir = os.path.dirname(profile_dir)
        guessed2 = [
            os.path.join(galaxy_dir, "default", "data.csv"),
            os.path.join(galaxy_dir, "default", "star.csv"),
            os.path.join(galaxy_dir, "v3", "data.csv"),
            os.path.join(galaxy_dir, "v3", "spec.csv"),
        ]
        for p in guessed2:
            if os.path.isfile(p):
                return p

    return None


def _load_inner_tracers(config: dict, radius_path: Optional[str], radius_col: str) -> Tuple[float, float]:
    """
    Load the first two positive finite stellar radii.
    Fallback uses half-light radius if a file or column cannot be read.
    """
    if radius_path is not None:
        try:
            df = pd.read_csv(radius_path)
            if radius_col not in df.columns:
                raise KeyError(f"Radius column {radius_col!r} not found in {radius_path}")

            r = pd.to_numeric(df[radius_col], errors="coerce").values.astype(float)
            r = r[np.isfinite(r) & (r > 0)]
            r = np.sort(r)

            if len(r) >= 2:
                return float(r[0]), float(r[1])
            if len(r) == 1:
                return float(r[0]), float(r[0])
        except Exception as e:
            print(f"[HaloBounds] warning: failed to read inner tracers from {radius_path}: {e}")

    # Fallback if no usable stellar file was found/read.
    r_half = float(config.get("R_HALF_LIGHT_PC", 1.0))
    if not np.isfinite(r_half) or r_half <= 0:
        r_half = 1.0

    # Conservative fallback: use a small inner interval around half-light scale.
    r1 = 0.75 * r_half
    r2 = 1.25 * r_half
    print(
        "[HaloBounds] warning: using fallback inner tracer radii "
        f"from R_HALF_LIGHT_PC={r_half:.6g} pc -> r1={r1:.6g}, r2={r2:.6g}"
    )
    return float(r1), float(r2)