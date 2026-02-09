# OSPM/Physics/OSPM_Physics.py
# Julia-backed stellar dynamics (LS6-safe, hardened)

import os, sys, numpy as np, subprocess

# --- Pin Python for PyCall ---
os.environ["PYTHON"] = sys.executable

USE_JULIA = os.environ.get("OSPM_USE_JULIA","0").strip().lower() in ("1","true","yes")
_JL_READY = False
_Main = None

print("[PY] OSPM_Physics imported from:", __file__)

# --- physical constants ---
pc   = 3.085677581e16
kms  = 1.0e3
Msun = 1.98847e30
G    = 6.67430e-11
c    = 2.99792458e8

_LAST_SIG = None


# ------------------------------------------------------------
# helper: locate gcc libstdc++ (for GLIBCXX_3.4.26)
# ------------------------------------------------------------
def _gcc_libstdcpp_dir():
    try:
        p = subprocess.check_output(
            ["g++", "-print-file-name=libstdc++.so.6"],
            text=True
        ).strip()
        if os.path.isabs(p) and os.path.exists(p):
            return os.path.dirname(p)
    except Exception:
        pass
    return None


# ------------------------------------------------------------
# Julia initialization (this is the critical section)
# ------------------------------------------------------------
def _jl_init():
    global _JL_READY, _Main
    if _JL_READY:
        return
    if not USE_JULIA:
        raise RuntimeError("OSPM_USE_JULIA is not enabled")

    # --- Julia executable (fixed, no guessing) ---
    julia_exe = os.environ.get(
        "OSPM_JULIA_EXE",
        "/scratch/09909/n8lujan99/julia-1.10.10/bin/julia"
    )
    if not os.path.exists(julia_exe):
        raise FileNotFoundError(f"Julia executable not found: {julia_exe}")

    # --- Always use scratch depot ---
    os.environ["JULIA_DEPOT_PATH"] = "/scratch/09909/n8lujan99/.julia"

    from julia.api import Julia

    # --- snapshot environment ---
    _ld = os.environ.get("LD_LIBRARY_PATH")
    _lp = os.environ.get("LD_PRELOAD")

    # --- thin LD_LIBRARY_PATH: only libstdc++ ---
    gcc_lib = _gcc_libstdcpp_dir()

    try:
        if gcc_lib:
            os.environ["LD_LIBRARY_PATH"] = gcc_lib
        else:
            os.environ.pop("LD_LIBRARY_PATH", None)

        os.environ.pop("LD_PRELOAD", None)

        # --- start Julia ---
        Julia(runtime=julia_exe, compiled_modules=False)

    finally:
        # restore original environment for Python
        if _ld is not None:
            os.environ["LD_LIBRARY_PATH"] = _ld
        else:
            os.environ.pop("LD_LIBRARY_PATH", None)

        if _lp is not None:
            os.environ["LD_PRELOAD"] = _lp
        else:
            os.environ.pop("LD_PRELOAD", None)

    from julia import Main as _JuliaMain
    _Main = _JuliaMain

    # --- load backend ---
    here = os.path.dirname(os.path.abspath(__file__))
    jl_path = os.path.join(here, "OSPM_Physics_Spherical.jl")
    if not os.path.exists(jl_path):
        raise FileNotFoundError(f"Julia backend file not found: {jl_path}")

    if not hasattr(_Main, "OSPMPhysicsSpherical"):
        _Main.include(jl_path)

    if not hasattr(_Main, "OSPMPhysicsSpherical"):
        raise RuntimeError("OSPMPhysicsSpherical failed to load")

    _JL_READY = True


# ------------------------------------------------------------
# parameter contract / cache
# ------------------------------------------------------------
def _theta_sig(theta, halo_type):
    t = np.asarray(theta, float).ravel()
    rho_s = float(t[0])
    r_s   = float(t[1])
    MBH   = float(t[2]) if t.size >= 3 else 0.0
    return (rho_s, r_s, MBH, str(halo_type).lower())


def assert_theta_contract(theta, *, halo_type, bounds=None, require_mbh=True):
    t = np.asarray(theta, float).ravel()
    if t.size < 2:
        raise ValueError("theta too short")
    if require_mbh and t.size < 3:
        raise ValueError("theta missing MBH")

    rho_s = float(t[0])
    r_s   = float(t[1])
    MBH   = float(t[2]) if t.size >= 3 else 0.0

    if rho_s <= 0 or r_s <= 0:
        raise ValueError("rho_s and r_s must be > 0")
    if MBH < 0:
        raise ValueError("MBH must be >= 0")

    ht = str(halo_type).strip().lower()
    if ht not in ("nfw", "cored"):
        raise ValueError(f"Unknown halo_type: {halo_type}")

    return rho_s, r_s, MBH, ht


def reset_orbit_cache_julia():
    _jl_init()
    _Main.OSPMPhysicsSpherical.reset_orbit_cache()


def maybe_reset_orbit_cache(theta, halo_type):
    global _LAST_SIG
    sig = _theta_sig(theta, halo_type)
    if sig != _LAST_SIG:
        _LAST_SIG = sig
        reset_orbit_cache_julia()


# ------------------------------------------------------------
# public Julia-backed API
# ------------------------------------------------------------
def build_A_matrix_stellar_julia(
    *,
    R_star_m,
    v_star_mps,
    verr_star_mps,
    sini,
    Norbit,
    theta,
    halo_type,
    return_occ=False,
    Nbins_occ=6,
    diag=False,
):
    _jl_init()

    rho_s, r_s, MBH, ht = assert_theta_contract(
        theta, halo_type=halo_type, require_mbh=True
    )

    maybe_reset_orbit_cache((rho_s, r_s, MBH), ht)

    out = _Main.OSPMPhysicsSpherical.build_A_matrix_stellar(
        int(Norbit),
        np.asarray(R_star_m, float),
        np.asarray(v_star_mps, float),
        np.asarray(verr_star_mps, float),
        float(sini),
        float(rho_s),
        float(r_s),
        float(MBH),
        str(ht),
        return_occ=bool(return_occ),
        Nbins_occ=int(Nbins_occ),
        diag=bool(diag),
    )

    if diag:
        A, meta = out
        return np.asarray(A, float), dict(meta)

    return np.asarray(out, float)
