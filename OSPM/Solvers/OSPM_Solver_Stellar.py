from __future__ import annotations
import numpy as np
from scipy.optimize import minimize
from ..Physics.OSPM_Physics import build_A_matrix_stellar_julia

def _as_1d_float(x, name: str) -> np.ndarray:
    a = np.asarray(x, dtype=float).reshape(-1)
    if a.size == 0:
        raise ValueError(f"{name} is empty")
    if not np.all(np.isfinite(a)):
        raise ValueError(f"{name} contains non-finite values")
    return a

def stellar_log_likelihood(
    A: np.ndarray,
    w: np.ndarray,
    *,
    verr_star_mps: np.ndarray,
    rv_mask: np.ndarray,
    Nstar: int,
    Nocc: int,
    lambda_occ: float = 1.0,
    eps: float = 1e-300,
    sigma_floor_mps: float = 2e3
) -> float:
    A = np.asarray(A, float)
    w = _as_1d_float(w, "w")
    if A.ndim != 2:
        raise ValueError("A must be 2D")
    if A.shape[1] != w.size:
        raise ValueError(f"A.shape[1]={A.shape[1]} != len(w)={w.size}")
    if Nstar + Nocc != A.shape[0]:
        raise ValueError(f"Nstar+Nocc={Nstar+Nocc} does not match A rows={A.shape[0]}")
    p = A @ w
    p = np.maximum(p, float(eps))
    # --- RV likelihood on subset only ---
    logL_star = 0.0
    if Nstar > 0:
        rv_mask = np.asarray(rv_mask, bool).reshape(-1)
        if rv_mask.size != Nstar:
            raise ValueError("rv_mask length must match Nstar")
        if rv_mask.any():
            p_star = p[:Nstar][rv_mask]
            sigma  = np.asarray(verr_star_mps, float).reshape(-1)[rv_mask]
            sigma  = np.nan_to_num(sigma, nan=np.inf, posinf=np.inf, neginf=np.inf)
            sigma  = np.maximum(sigma, float(sigma_floor_mps))
            p_star = np.maximum(p_star * sigma, float(eps))
            logL_star = float(np.sum(np.log(p_star)))
    # --- occupancy likelihood always on occupancy rows ---
    if Nocc == 0:
        return logL_star
    p_occ = np.maximum(p[Nstar:], float(eps))
    logL_occ = float(np.sum(np.log(p_occ)))
    return logL_star + float(lambda_occ) * logL_occ

def solve_weights_stellar(
    A: np.ndarray,
    *,
    verr_star_mps: np.ndarray,
    rv_mask: np.ndarray,
    Nstar: int,
    Nocc: int,
    lambda_occ: float = 1.0,
    alpha: float = 1e-3,
    eps: float = 1e-300,
    maxiter: int = 500,
    maxfun: int = 20000,
    p0_floor: float = 1e-15
) -> np.ndarray:
    A = np.asarray(A, float)
    Norbit = int(A.shape[1])
    alpha_eff = float(alpha) * (Norbit / 200.0) * (max(Nstar, 1) / 90.0)
    w0 = np.random.rand(Norbit)
    w0 += 1e-3
    w0 /= w0.sum()
    # Scale rows using only "active" rows:
    # RV rows where rv_mask True plus all occupancy rows.
    rv_mask = np.asarray(rv_mask, bool).reshape(-1)
    if rv_mask.size != Nstar:
        raise ValueError("rv_mask length must match Nstar")
    
    active_rows = np.concatenate([rv_mask, np.ones(int(Nocc), dtype=bool)]) if (Nocc > 0) else rv_mask.copy()
    
    if active_rows.size != A.shape[0]:
        raise ValueError("active_rows length mismatch")
    
    p0 = A @ w0
    p0_active = p0[active_rows]
    scale_active = 1.0 / np.maximum(p0_active, float(p0_floor))
    A_eff = A.copy()
    A_eff[active_rows, :] *= scale_active[:, None]
    bounds = [(0.0, None)] * Norbit

    def objective(w: np.ndarray) -> float:
        ll = stellar_log_likelihood(
            A_eff, w,
            verr_star_mps=verr_star_mps,
            rv_mask=rv_mask,
            Nstar=Nstar, Nocc=Nocc,
            lambda_occ=lambda_occ,
            eps=eps
        )
        return -ll + alpha_eff * float(np.dot(w, w))

    res = minimize(objective, w0, method="L-BFGS-B", bounds=bounds,
                   options={"maxiter": int(maxiter), "maxfun": int(maxfun)})

    if (not res.success) or (res.x is None):
        res = minimize(objective, w0, method="SLSQP", bounds=bounds,
                       options={"maxiter": int(maxiter * 3), "ftol": 1e-9, "disp": False})
        if (not res.success) or (res.x is None):
            msg = getattr(res, "message", "unknown")
            raise RuntimeError(f"Weight solve failed: {msg}")

    w = np.asarray(res.x, float)
    w[w < 0] = 0.0
    s = float(w.sum())
    if (not np.isfinite(s)) or (s <= 0):
        w[:] = 1.0 / Norbit
    else:
        w /= s
    return w

def solve_ospm_theta_stellar(theta, obs, *, halo_type: str = "nfw"):
    rho_s, r_s, MBH = map(float, theta)

    # --- all stars always (geometry + RV) ---
    R_all  = np.asarray(obs.R_star_m, float)
    v_all  = np.asarray(obs.v_star_mps, float)
    ve_all = np.asarray(obs.verr_star_mps, float)

    if not np.all(np.isfinite(R_all)):
        raise ValueError("obs.R_star_m contains non-finite values")
    if R_all.size == 0:
        raise ValueError("obs.R_star_m is empty")

    rv = np.asarray(getattr(obs, "valid_vlos", None), bool)
    if rv is None or rv.size != R_all.size:
        raise ValueError("obs.valid_vlos missing or length mismatch")

    # Enforce finiteness only where RV is valid
    if rv.any():
        if not np.all(np.isfinite(v_all[rv])):
            raise ValueError("v_star_mps contains non-finite values for valid_vlos stars")
        if not np.all(np.isfinite(ve_all[rv])):
            raise ValueError("verr_star_mps contains non-finite values for valid_vlos stars")

    Norbit = int(obs.Norbit)
    if Norbit < 1:
        raise ValueError("obs.Norbit must be >= 1")

    sini = float(obs.sini)
    if not np.isfinite(sini):
        raise ValueError("obs.sini is non-finite")

    # Occupancy knobs
    Nocc = int(getattr(obs, "Nocc", 0))
    lambda_occ = float(getattr(obs, "lambda_occ", 1.0))

    # Build A for ALL stars. Julia already zeros out NaN/invalid sigma rows.
    # return_occ=True gives star rows + occupancy rows.
    A = build_A_matrix_stellar_julia(
        R_star_m=R_all, v_star_mps=v_all, verr_star_mps=ve_all,
        sini=sini, Norbit=Norbit, theta=[rho_s, r_s, MBH],
        halo_type=halo_type, return_occ=True
    )
    A = np.asarray(A, float)
    if A.ndim != 2 or A.shape[1] != Norbit:
        raise ValueError(f"Unexpected A shape: {A.shape}, expected (*,{Norbit})")

    Nstar = int(R_all.size)
    if A.shape[0] != Nstar + Nocc:
        raise ValueError(f"A rows={A.shape[0]} do not match Nstar+Nocc={Nstar+Nocc}")

    # Solve weights with RV subset + occupancy. Geometry-only stars still influence occupancy.
    w = solve_weights_stellar(
        A,
        verr_star_mps=ve_all,
        rv_mask=rv,
        Nstar=Nstar,
        Nocc=Nocc,
        lambda_occ=lambda_occ,
        alpha=float(getattr(obs, "alpha_w", 1e-2)),
        maxiter=int(getattr(obs, "w_maxiter", 500))
    )

    ll = stellar_log_likelihood(
        A, w,
        verr_star_mps=ve_all,
        rv_mask=rv,
        Nstar=Nstar,
        Nocc=Nocc,
        lambda_occ=lambda_occ
    )

    chi2_like = -2.0 * ll
    # For mixed capability, dof should use Nrv, not Nstar.
    Nrv = int(rv.sum())
    nu = max(Nrv - 3, 1)
    chi2_red = float(chi2_like) / float(nu)

    model = A @ w
    return float(chi2_red), w, model
