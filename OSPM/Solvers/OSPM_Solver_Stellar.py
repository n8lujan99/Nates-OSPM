from __future__ import annotations
import time

import numpy as np
from scipy.optimize import minimize

from ..Physics import OSPM_Physics as phys
from ..Physics.OSPM_Physics import build_A_matrix_stellar_julia
from ..Physics.OSPM_PhysicsEngine import chi2_resolution_penalty, mass_slope_penalty


# Utilities
def _as_1d_float(x, name: str) -> np.ndarray:
    a = np.asarray(x, float).reshape(-1)
    if a.size == 0 or not np.all(np.isfinite(a)):
        raise ValueError(f"{name} invalid")
    return a


def _safe_float(x, default=np.inf) -> float:
    try:
        y = float(x)
        return y if np.isfinite(y) else float(default)
    except Exception:
        return float(default)


# Likelihood
def stellar_log_likelihood(A, w, *, verr_star_mps, rv_mask, Nstar: int, Nocc: int, lambda_occ: float = 1.0, eps: float = 1e-300, sigma_floor_mps: float = 2e3) -> float:
    A = np.asarray(A, float)
    w = _as_1d_float(w, "w")

    if A.ndim != 2 or A.shape[1] != w.size or A.shape[0] != Nstar + Nocc:
        raise ValueError("A shape mismatch")

    p = np.maximum(A @ w, eps)
    ll = 0.0

    # Velocity likelihood
    if Nstar > 0:
        rv_mask = np.asarray(rv_mask, bool)
        if rv_mask.size != Nstar:
            raise ValueError("rv_mask mismatch")

        if rv_mask.any():
            sig = np.asarray(verr_star_mps, float)[rv_mask]
            sig = np.maximum(np.nan_to_num(sig, np.inf), sigma_floor_mps)
            ll += float(np.sum(np.log(np.maximum(p[:Nstar][rv_mask] * sig, eps))))

    # Occupancy likelihood
    if Nocc > 0:
        ll += float(lambda_occ) * float(np.sum(np.log(np.maximum(p[Nstar:], eps))))

    return float(ll)


# Weight Solver
def solve_weights_stellar(A, *, verr_star_mps, rv_mask, Nstar: int, Nocc: int, lambda_occ: float = 1.0, alpha: float = 1e-3, eps: float = 1e-300, maxiter: int = 500, maxfun: int = 20000, p0_floor: float = 1e-15, seed: int | None = None):
    """
    Returns (w_normed, meta) on success, (None, meta) on failure.
    meta includes solver status and timing.
    """
    t0 = time.perf_counter()

    A = np.asarray(A, float)
    if A.ndim != 2 or A.shape[0] != Nstar + Nocc:
        return None, {"ok": False, "reason": "A_shape", "t_solve": 0.0}

    Norbit = A.shape[1]
    if Norbit <= 0:
        return None, {"ok": False, "reason": "Norb0", "t_solve": 0.0}

    # Scale regularization with problem size
    alpha_eff = float(alpha) * (Norbit / 200.0) * (max(Nstar, 1) / 90.0)

    rng = np.random.default_rng(seed)
    w0 = rng.random(Norbit) + 1e-3
    w0 = w0 / w0.sum()

    rv_mask = np.asarray(rv_mask, bool)
    if rv_mask.size != Nstar:
        return None, {"ok": False, "reason": "rv_mask_mismatch", "t_solve": 0.0}

    active = np.concatenate([rv_mask, np.ones(Nocc, bool)]) if Nocc > 0 else rv_mask

    # Rescale active rows to avoid catastrophic scaling
    p0 = A @ w0
    Aeff = A.copy()
    if active.any():
        scale = 1.0 / np.maximum(p0[active], p0_floor)
        Aeff[active] *= scale[:, None]

    bounds = [(0.0, None)] * Norbit

    def obj(w):
        ll = stellar_log_likelihood(Aeff, w, verr_star_mps=verr_star_mps, rv_mask=rv_mask, Nstar=Nstar, Nocc=Nocc, lambda_occ=lambda_occ, eps=eps)
        return -ll + alpha_eff * float(np.dot(w, w))

    # Try L-BFGS-B first
    try:
        res = minimize(obj, w0, method="L-BFGS-B", bounds=bounds, options={"maxiter": int(maxiter), "maxfun": int(maxfun)})
        method = "L-BFGS-B"
    except Exception as e:
        t1 = time.perf_counter()
        return None, {"ok": False, "reason": f"lbfgsb_exc:{type(e).__name__}", "t_solve": t1 - t0}

    # Fallback: SLSQP
    if not getattr(res, "success", False):
        try:
            res2 = minimize(obj, w0, method="SLSQP", bounds=bounds, options={"maxiter": int(3 * maxiter), "ftol": 1e-9})
            if getattr(res2, "success", False):
                res = res2
                method = "SLSQP"
        except Exception as e:
            t1 = time.perf_counter()
            return None, {
                "ok": False,
                "reason": f"slsqp_exc:{type(e).__name__}",
                "t_solve": t1 - t0,
                "method": "SLSQP",
            }

    t1 = time.perf_counter()

    if not getattr(res, "success", False):
        msg = getattr(res, "message", "fail")
        return None, {
            "ok": False,
            "reason": f"no_converge:{msg}",
            "t_solve": t1 - t0,
            "method": method,
            "nit": int(getattr(res, "nit", -1)),
            "nfev": int(getattr(res, "nfev", -1)),
        }

    w = np.maximum(np.asarray(res.x, float), 0.0)
    s = float(w.sum())
    if not np.isfinite(s) or s <= 0.0:
        return None, {
            "ok": False,
            "reason": "bad_sum",
            "t_solve": t1 - t0,
            "method": method,
        }

    return (w / s), {
        "ok": True,
        "t_solve": t1 - t0,
        "method": method,
        "nit": int(getattr(res, "nit", -1)),
        "nfev": int(getattr(res, "nfev", -1)),
    }


# Theta Solve
def solve_ospm_theta_stellar(theta, obs, *, halo_type: str = "nfw", diag: bool = False):
    """
    Returns:
      - chi2_red (float)
      - w (np.ndarray | None)
      - model (np.ndarray | None)  # A@w
    If diag=True, returns a 4th item: meta dict with timings and failure stage.
    """

    meta = {"ok": False, "stage": None}
    t_all0 = time.perf_counter()
    out = (float(np.inf), None, None)

    def fail(stage: str, **extra):
        if diag:
            meta.update(stage=stage, t_total=time.perf_counter() - t_all0, **extra)
            return (*out, meta)
        return out

    try:
        rho_s, r_s, MBH = map(float, theta)
    except Exception:
        return fail("theta_parse")

    R = np.asarray(obs.R_star_m, float)
    v = np.asarray(obs.v_star_mps, float)
    ve = np.asarray(obs.verr_star_mps, float)

    if R.size == 0 or not np.all(np.isfinite(R)):
        return fail("obs_R")

    rv = np.asarray(obs.valid_vlos, bool)
    if rv.size != R.size:
        return fail("obs_rv_mask")

    Norbit = int(obs.Norbit)
    sini = float(obs.sini)
    Nocc = int(getattr(obs, "Nocc", 0))
    lambda_occ = float(getattr(obs, "lambda_occ", 1.0))

    # Build orbit matrix (Julia)
    tA0 = time.perf_counter()
    try:
        A = build_A_matrix_stellar_julia(R_star_m=R, v_star_mps=v, verr_star_mps=ve, sini=sini, Norbit=Norbit, theta=[rho_s, r_s, MBH], halo_type=halo_type, return_occ=True)
    except Exception as e:
        return fail("build_A_exc", exc=type(e).__name__, t_buildA=time.perf_counter() - tA0)
    tA1 = time.perf_counter()

    A = np.asarray(A, float)
    if A.ndim != 2 or A.shape[1] != Norbit or A.shape[0] != R.size + Nocc:
        return fail("A_malformed", A_shape=tuple(A.shape), t_buildA=tA1 - tA0)

    # Solve weights (SciPy)
    tW0 = time.perf_counter()
    w, wmeta = solve_weights_stellar(A, verr_star_mps=ve, rv_mask=rv, Nstar=R.size, Nocc=Nocc, lambda_occ=lambda_occ, alpha=float(getattr(obs, "alpha_w", 1e-2)), maxiter=int(getattr(obs, "w_maxiter", 500)), seed=None)
    tW1 = time.perf_counter()

    if w is None:
        return fail("weights_fail", t_buildA=tA1 - tA0, t_weights=tW1 - tW0, wmeta=wmeta)

    # Likelihood
    try:
        ll = stellar_log_likelihood(A, w, verr_star_mps=ve, rv_mask=rv, Nstar=R.size, Nocc=Nocc, lambda_occ=lambda_occ)
    except Exception as e:
        return fail("ll_exc", exc=type(e).__name__, t_buildA=tA1 - tA0, t_weights=tW1 - tW0)

    chi2_like = -2.0 * float(ll)

    # Penalties
    pen1 = _safe_float(chi2_resolution_penalty(MBH_msun=MBH, R_star_m=R, v_star_mps=v, strength=float(getattr(obs, "PEN_SPHERE_STRENGTH", 2.0)), power=float(getattr(obs, "PEN_SPHERE_POWER", 2.0)))[0], default=0.0)
    pen2 = _safe_float(mass_slope_penalty(theta=[rho_s, r_s, MBH], halo_type=halo_type, R_star_m=R, strength=float(getattr(obs, "PEN_SLOPE_STRENGTH", 0.5)))[0], default=0.0)

    chi2 = chi2_like + pen1 + pen2
    nu = max(int(rv.sum()) - 3, 1)
    chi2_red = float(chi2 / nu)
    model = A @ w
    r_grid = np.logspace(np.log10(np.min(R)), np.log10(np.max(R)), 50)
    M_r = np.array([phys.mass_enclosed_two_radii_julia(r_in_m=1e-6, r_out_m=r, theta=[rho_s, r_s, MBH], halo_type=halo_type)[0] for r in r_grid], dtype=float)
    Vcirc = np.sqrt(np.maximum(phys.G * M_r / r_grid, 0.0))
    L_r = np.array([np.sum(R <= r) for r in r_grid], dtype=float)
    diag_data = {"r": r_grid, "M": M_r, "Vcirc": Vcirc, "ML": M_r / np.maximum(L_r, 1e-10)}

    if diag:
        meta |= {"ok": True, "stage": "ok", "t_buildA": tA1 - tA0, "t_weights": tW1 - tW0, "t_total": time.perf_counter() - t_all0, "wmeta": wmeta, "chi2_like": float(chi2_like), "pen_res": float(pen1), "pen_slope": float(pen2), **diag_data}
        return chi2_red, w, model, meta

    return chi2_red, w, model, diag_data
