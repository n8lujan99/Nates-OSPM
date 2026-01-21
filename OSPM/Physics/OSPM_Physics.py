import os
import numpy as np
from scipy.optimize import nnls

USE_JULIA = os.environ.get("OSPM_USE_JULIA", "0").strip().lower() in ("1", "true", "yes")

_JL_READY = False
_Main = None

# ---------------- Constants ----------------
pc   = 3.085677581e16
kms  = 1.0e3
Msun = 1.98847e30
G    = 6.67430e-11
c    = 2.99792458e8

# ============================================================
# Orbit-library cache invalidation (Python gate)
# ============================================================

_ORBIT_LIB_SIG = None

def _potential_signature(theta, halo_type):
    # theta is astro units: rho_s [Msun/pc^3], r_s [pc], MBH [Msun]
    rho_s, r_s, MBH = map(float, theta[:3])
    ht = str(halo_type).strip().lower()
    # round to avoid thrash from tiny float jitter
    # keep this aggressive enough to treat "same proposal" as same sig
    rho_s_r = float(np.round(rho_s, 12))
    r_s_r   = float(np.round(r_s,   12))
    MBH_r   = float(np.round(MBH,   12))
    return (ht, rho_s_r, r_s_r, MBH_r)

def _maybe_reset_orbit_library(theta, halo_type):
    global _ORBIT_LIB_SIG
    sig = _potential_signature(theta, halo_type)
    if _ORBIT_LIB_SIG == sig:
        return
    _ORBIT_LIB_SIG = sig
    if hasattr(_Main.OSPMPhysicsSpherical, "reset_orbit_cache"):
        _Main.OSPMPhysicsSpherical.reset_orbit_cache()
    else:
        # hard fail.
        raise RuntimeError( "Julia backend missing reset_orbit_cache(). "
                            "Add it to OSPM_Physics_Spherical.jl.")

# ---------------- Mass enclosed (Julia) ----------------
def mass_enclosed_two_radii_julia(*, r_in_m, r_out_m, theta, halo_type):
    if not USE_JULIA:
        raise RuntimeError("Julia required")
    _jl_init()
    rho_s, r_s, MBH = map(float, theta)
    Min, Mout = _Main.OSPMPhysicsSpherical.mass_enclosed_two_radii( float(r_in_m), float(r_out_m),
        float(rho_s), float(r_s), float(MBH), str(halo_type) )
    return float(Min), float(Mout)

print("IMPORTING OSPM_Physics FROM:", __file__)

def _jl_init():
    global _JL_READY, _Main
    if _JL_READY:
        return
    if not USE_JULIA:
        raise RuntimeError("OSPM_USE_JULIA is not enabled")
    from julia.api import Julia
    jl = Julia(compiled_modules=False)
    from julia import Main as _JuliaMain
    _Main = _JuliaMain
    here = os.path.dirname(os.path.abspath(__file__))
    jl_path = os.path.join(here, "OSPM_Physics_Spherical.jl")
    if not os.path.exists(jl_path):
        raise FileNotFoundError(f"Julia backend file not found: {jl_path}")
    _Main.include(jl_path)
    # Hard sanity check
    if not hasattr(_Main, "OSPMPhysicsSpherical"):
        raise RuntimeError( "Julia loaded but module OSPMPhysicsSpherical not found in Main" )
    _JL_READY = True

# ---------------- Cheap helpers (Python) ----------------
def make_inclination(inclination_deg: float):
    inc = np.radians(float(inclination_deg))
    sini = float(np.sin(inc))
    cosi = float(np.cos(inc))
    edge_on = float(inclination_deg) >= 85.0
    return sini, cosi, edge_on

def make_rng(seed: int):
    return np.random.default_rng(int(seed))

# ---------------- Halo / context ----------------
def halo_from_theta_astro(theta, halo_type="nfw"):

    theta = np.asarray(theta, float)
    rho_s = float(theta[0])
    r_s   = float(theta[1])
    MBH   = float(theta[2]) if len(theta) > 2 else 0.0
    ht    = str(halo_type).lower()
    return {"rho_s": rho_s, "r_s": r_s, "MBH": MBH, "type": ht}

def build_dynamics_context(*, theta, halo_type, **_ignored):
    if not USE_JULIA:
        raise RuntimeError("Python backend is disabled. Set OSPM_USE_JULIA=1.")
    halo_py = halo_from_theta_astro(theta, halo_type=halo_type)
    return {"halo": halo_py}

def halo_kwargs_from_ctx(ctx):
    halo = ctx["halo"] if isinstance(ctx, dict) else ctx.halo
    for k in ("rho_s", "r_s", "MBH", "type"):
        if k not in halo:
            raise KeyError(f"halo missing required key '{k}'")
    return dict( rho_s=float(halo["rho_s"]), r_s=float(halo["r_s"]), MBH=float(halo["MBH"]), halo_type=str(halo["type"]))

# ---------------- A-matrix (Julia) ----------------
def extract_ic_arrays(obs):
    Norbit = int(obs.Norbit)
    r0   = np.empty(Norbit, dtype=float)
    th0  = np.empty(Norbit, dtype=float)
    dt   = np.empty(Norbit, dtype=float)
    Etot = np.empty(Norbit, dtype=float)
    xLz  = np.empty(Norbit, dtype=float)
    for k in range(Norbit):
        r0[k], th0[k], dt[k], Etot[k], xLz[k] = obs.initial_conditions(k)
    return r0, th0, dt, Etot, xLz

def build_A_matrix_stellar_julia( *, R_star_m, v_star_mps, verr_star_mps, sini, Norbit,
                            theta, halo_type, return_occ=False, Nbins_occ=6):
    if not USE_JULIA:
        raise RuntimeError("Stellar mode requires Julia")
    _jl_init()
    _maybe_reset_orbit_library(theta, halo_type)
    rho_s, r_s, MBH = map(float, theta[:3])
    A = _Main.OSPMPhysicsSpherical.build_A_matrix_stellar( int(Norbit), np.asarray(R_star_m, float),
        np.asarray(v_star_mps, float), np.asarray(verr_star_mps, float), float(sini), float(rho_s),
        float(r_s), float(MBH), str(halo_type), return_occ=bool(return_occ), Nbins_occ=int(Nbins_occ) )
    A = np.asarray(A, float)
    return A

def build_A_matrix(obs, ctx, Norbit=None):
    r0, th0, dt, Etot, xLz = extract_ic_arrays(obs)  # r0/Etot/xLz currently unused by Julia
    hk = halo_kwargs_from_ctx(ctx)
    return build_A_matrix_stellar_julia( th0=th0, dt=dt, r_centers_m=obs.r_centers_m,
        valid=obs.valid, sini=float(obs.sini), nsteps=int(obs.nsteps_orbit), **hk )

def rho_interp(*args, **kwargs):
    raise RuntimeError("rho_interp is legacy-only. It should never be called in Julia mode.")
