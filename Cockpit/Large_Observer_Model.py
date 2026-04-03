# Large Observer Model (LOM) for OSPM


import os
import time
import json
import shutil
from datetime import datetime

import numpy as np
import pandas as pd

REPO_ROOT = "/home/nate/research/Nates-OSPM"
COCKPIT_ROOT = f"{REPO_ROOT}/Cockpit"
SHARED_ROOT = f"{COCKPIT_ROOT}/shared"
PARAMS = ["rho_s", "r_s", "MBH"]

LOM_MODE = os.environ.get("LOM_MODE", "local").strip().lower()
IS_HPC = LOM_MODE in {"hpc", "cluster", "headless"}
IS_LOOP = LOM_MODE in {"hpc", "cluster", "headless", "dev"}

LOCAL_SHOW_PLOTS = (
    os.environ.get("LOM_SHOW_PLOTS", "0").strip().lower() in {"1", "true", "yes", "on"}
) and (not IS_HPC)

PLOT_ENABLED = (
    os.environ.get("LOM_MAKE_PLOTS", "1").strip().lower() not in {"0", "false", "no", "off"}
) and (not IS_HPC)

SUPPORT_FRAC = 0.20
SUPPORT_MIN = 50
MIN_POINTS = 100
MAX_GP_POINTS = 800

if PLOT_ENABLED:
    import matplotlib
    if os.environ.get("DISPLAY", "") == "":
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
else:
    plt = None

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C


# =========================
# PATHS
# =========================
def get_paths():
    with open(f"{REPO_ROOT}/which_galaxy") as f:
        galaxy = f.read().strip()

    profile = f"{REPO_ROOT}/Data/Gal_Profiles/{galaxy}"
    deck = f"{profile}/default/daemon_deck_local.csv"
    config = f"{profile}/{galaxy}_OSPM_Config.py"
    gal_dir = f"{COCKPIT_ROOT}/{galaxy}"

    return galaxy, deck, config, gal_dir


# =========================
# LOAD
# =========================
def load_deck(path):
    df = pd.read_csv(path)
    df = df[df["status"].str.startswith("pass")].copy()
    df["label"] = df["status"].str.replace("pass_", "", regex=False)
    df = df[np.isfinite(df["chi2"])].copy()
    return df


# =========================
# SUPPORT CLOUD
# =========================
def compute_support(df):
    n = max(SUPPORT_MIN, int(SUPPORT_FRAC * len(df)))
    n = min(n, len(df))
    return df.nsmallest(n, "chi2").copy()


# =========================
# GP (secondary hint only)
# =========================
def fit_gp(df):
    gp_df = df
    if len(gp_df) > MAX_GP_POINTS:
        gp_df = gp_df.sample(MAX_GP_POINTS, random_state=42)

    X = gp_df[PARAMS].values.astype(float).copy()

    for j, p in enumerate(PARAMS):
        if np.all(X[:, j] > 0):
            X[:, j] = np.log10(X[:, j])

    y = gp_df["chi2"].values.astype(float)

    kernel = C(1.0) * RBF(length_scale=np.ones(len(PARAMS)))
    gp = GaussianProcessRegressor(
        kernel=kernel,
        alpha=1e-2,
        optimizer=None,
        normalize_y=True,
    )
    gp.fit(X, y)
    return gp


def curvature_proxy(gp, center):
    x0 = center.astype(float).copy()
    for j in range(len(x0)):
        if x0[j] > 0:
            x0[j] = np.log10(x0[j])

    eps = 1e-2
    base = gp.predict(x0.reshape(1, -1))[0]
    curv = []

    for i in range(len(x0)):
        d = np.zeros_like(x0)
        d[i] = eps
        up = gp.predict((x0 + d).reshape(1, -1))[0]
        dn = gp.predict((x0 - d).reshape(1, -1))[0]
        curv.append(float(up + dn - 2.0 * base))

    return curv


def classify_shape(curv, tol=1e-6):
    s = []
    for x in curv:
        if abs(x) <= tol:
            s.append(0)
        elif x > 0:
            s.append(1)
        else:
            s.append(-1)

    if all(v > 0 for v in s):
        return "basin_like"
    if all(v < 0 for v in s):
        return "peak_like"
    if any(v == 0 for v in s):
        return "flat_or_unresolved"
    return "mixed_or_ridge_like"


# =========================
# METRICS
# =========================
def boundary_pressure(best, df):
    out = {}
    for p in PARAMS:
        lo = float(df[p].min())
        hi = float(df[p].max())
        span = hi - lo if hi > lo else 1.0

        out[p] = {
            "low": float(np.mean(best[p] < lo + 0.10 * span)),
            "high": float(np.mean(best[p] > hi - 0.10 * span)),
        }
    return out


def safe_corr(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3:
        return np.nan
    if np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def build_evidence(df, best):
    lo = best[PARAMS].quantile(0.16)
    hi = best[PARAMS].quantile(0.84)

    evidence = {
        "chi2": {
            "global_best": float(df["chi2"].min()),
            "global_median": float(df["chi2"].median()),
            "subset_best": float(best["chi2"].min()),
            "subset_median": float(best["chi2"].median()),
            "subset_std": float(best["chi2"].std()) if len(best) > 1 else 0.0,
        },
        "spread": {},
        "correlation_vs_chi2": {},
    }

    for p in PARAMS:
        evidence["spread"][p] = {
            "lo": float(lo[p]),
            "hi": float(hi[p]),
            "width": float(hi[p] - lo[p]),
        }
        evidence["correlation_vs_chi2"][p] = safe_corr(best[p].values, best["chi2"].values)

    return evidence, lo, hi


def build_notes(evidence, pressure, shape_hint, curvature, df, best, center, patch_lo, patch_hi):
    LOG_PARAMS = {"rho_s", "MBH"}
    notes = []

    # ── 1. Dataset quality ──────────────────────────────────────────────────
    n_total   = len(df)
    n_support = len(best)
    frac      = n_support / n_total
    g_best    = evidence["chi2"]["global_best"]
    s_best    = evidence["chi2"]["subset_best"]
    s_med     = evidence["chi2"]["subset_median"]

    notes.append(
        f"dataset: {n_total} pass rows — top {frac*100:.1f}% ({n_support} models) form the support cloud"
    )
    notes.append(
        f"chi2: global_best={g_best:.4g}  subset_best={s_best:.4g}  subset_median={s_med:.4g}"
    )
    if evidence["chi2"]["subset_std"] < 1e-3:
        notes.append("chi2 is nearly flat inside the support cloud — model is insensitive in this region")

    # ── 2. Landscape shape from GP ───────────────────────────────────────────
    shape_interp = {
        "basin_like":          "GP sees a bowl — genuine minimum, bounds look reliable",
        "peak_like":           "GP sees a ridge/peak — low-chi2 region may be a narrow crest, not a basin",
        "flat_or_unresolved":  "GP surface is flat — chi2 is not constraining here yet, sample more broadly",
        "mixed_or_ridge_like": "GP shows mixed curvature — some parameters constrained, others may be degenerate",
        "not evaluated":       "GP was not run this cycle",
        "gp_failed":           "GP fitting failed — check data range or increase sample size",
    }
    notes.append(f"GP landscape: {shape_hint} — {shape_interp.get(shape_hint, shape_hint)}")

    if curvature is not None:
        for i, p in enumerate(PARAMS):
            c = curvature[i]
            if abs(c) < 1e-6:
                desc = "flat — no curvature detected"
            elif c > 0:
                desc = f"concave up ({c:.3e}) — well-constrained direction"
            else:
                desc = f"concave down ({c:.3e}) — ridge or boundary artifact, treat with caution"
            notes.append(f"  GP curvature {p}: {desc}")

    # ── 3. Per-parameter status ──────────────────────────────────────────────
    for p in PARAMS:
        sp  = evidence["spread"][p]
        lo_v, hi_v = sp["lo"], sp["hi"]
        c   = center[p]
        corr = evidence["correlation_vs_chi2"].get(p, np.nan)

        if p in LOG_PARAMS and c > 0:
            rel = sp["width"] / c * 100
        else:
            rel = (sp["width"] / (abs(c) + 1e-30)) * 100

        notes.append(
            f"{p}: center={c:.4g}  16-84 range=[{lo_v:.4g}, {hi_v:.4g}]  spread~{rel:.1f}% of center"
        )

        if np.isfinite(corr):
            if abs(corr) > 0.60:
                direction = "lower values preferred" if corr < 0 else "higher values preferred"
                notes.append(f"  strong chi2 trend (r={corr:.2f}): {direction} — consider shifting bounds")
            elif abs(corr) > 0.30:
                direction = "lower" if corr < 0 else "higher"
                notes.append(f"  mild chi2 trend (r={corr:.2f}): {direction} values slightly preferred")

    # ── 4. Boundary pressure ─────────────────────────────────────────────────
    any_pressure = False
    for p in PARAMS:
        plo = pressure[p]["low"]
        phi = pressure[p]["high"]
        if plo > 0.30:
            notes.append(
                f"WARNING {p}: {plo*100:.0f}% of support cloud near LOWER edge"
                f" — true minimum may be outside bounds, expand lower limit"
            )
            any_pressure = True
        if phi > 0.30:
            notes.append(
                f"WARNING {p}: {phi*100:.0f}% of support cloud near UPPER edge"
                f" — true minimum may be outside bounds, expand upper limit"
            )
            any_pressure = True

    # ── 5. Suggested next bounds ─────────────────────────────────────────────
    notes.append("suggested bounds for next run (written to patch):")
    for p in PARAMS:
        notes.append(f"  {p}: [{patch_lo[p]:.4g}, {patch_hi[p]:.4g}]  center={center[p]:.4g}")

    # ── 6. Overall recommendation ────────────────────────────────────────────
    if any_pressure and shape_hint in {"basin_like", "mixed_or_ridge_like"}:
        notes.append("ACTION: expand bounds toward the pressured edges before tightening further")
    elif shape_hint == "basin_like" and not any_pressure:
        notes.append("ACTION: basin confirmed, no edge pressure — safe to tighten bounds or increase resolution")
    elif shape_hint in {"flat_or_unresolved", "not evaluated", "gp_failed"}:
        notes.append("ACTION: landscape unresolved — run more samples before narrowing bounds")
    elif shape_hint == "peak_like":
        notes.append("ACTION: ridge detected — verify this is physical and not a numerical artifact")
    else:
        notes.append("ACTION: continue sampling at current bounds")

    return notes


def compute_patch_bounds(df, best, center):
    """
    Suggest wider exploration bounds by expanding the support cloud's full
    min/max range outward by 60% of that span (log-space for log params),
    clipped to the full deck extent so we never suggest untested territory.
    """
    LOG_PARAMS = {"rho_s", "MBH"}
    lo_out = {}
    hi_out = {}

    for p in PARAMS:
        b_min = float(best[p].min())
        b_max = float(best[p].max())
        d_min = float(df[p].min())
        d_max = float(df[p].max())

        if p in LOG_PARAMS and b_min > 0 and b_max > b_min:
            log_lo   = np.log10(b_min)
            log_hi   = np.log10(b_max)
            log_span = log_hi - log_lo
            lo_v = max(10 ** (log_lo - 0.6 * log_span), d_min)
            hi_v = min(10 ** (log_hi + 0.6 * log_span), d_max)
        else:
            span = b_max - b_min if b_max > b_min else abs(center[p]) * 0.3 + 1.0
            lo_v = max(b_min - 0.6 * span, d_min)
            hi_v = min(b_max + 0.6 * span, d_max)

        lo_out[p] = float(lo_v)
        hi_out[p] = float(hi_v)

    return lo_out, hi_out


# =========================
# PATCH OUTPUT
# =========================
def _fmt(v):
    """Format a parameter value to a clean, readable number."""
    if v == 0:
        return "0.0"
    mag = abs(v)
    if 0.1 <= mag < 1e5:
        # fixed notation, drop trailing zeros
        s = f"{v:.4g}"
        return s if "." in s else s + ".0"
    return f"{v:.3e}"


def build_config_patch(galaxy, center, lo, hi):
    rc  = _fmt(center["rho_s"])
    rsc = _fmt(center["r_s"])
    mc  = _fmt(center["MBH"])
    rlo = _fmt(lo["rho_s"]);  rhi = _fmt(hi["rho_s"])
    rslo = _fmt(lo["r_s"]);   rshi = _fmt(hi["r_s"])
    mlo = _fmt(lo["MBH"]);    mhi = _fmt(hi["MBH"])

    block = (
        f"# --- LOM suggested region update for {galaxy} ---\n"
        f"# generated: {datetime.now().isoformat(timespec='seconds')}\n"
        f'"PARAMETER_NAMES": ["rho_s", "r_s", "MBH"],\n'
        f'"INITIAL_THETA": [{rc}, {rsc}, {mc}],\n'
        f'"THETA_BOUNDS": [\n'
        f"   ({rlo}, {rhi}),      # rho_s\n"
        f"   ({rslo}, {rshi}),    # r_s\n"
        f"   ({mlo}, {mhi}),      # MBH\n"
        f"],"
    )
    script = (
        "#!/bin/bash\n\n"
        f"cd {REPO_ROOT}\n"
        "cat <<'EOF' >> "
        f"Data/Gal_Profiles/{galaxy}/{galaxy}_OSPM_Config.py\n"
        f"{block}\nEOF\n"
    )
    return block, script


# =========================
# PLOTS (VERTICAL STRIP MODE)
# =========================
PLOT_DPI = 140

def _save_and_optionally_show(fig, outpath):
    fig.subplots_adjust(
        top=0.90,
        bottom=0.10,
        left=0.08,
        right=0.80,
        hspace=0.35
    )
    fig.savefig(outpath, dpi=PLOT_DPI, bbox_inches="tight")

    if LOCAL_SHOW_PLOTS:
        plt.show(block=False)
        plt.pause(0.001)

    plt.close(fig)


def _style_axes(ax):
    ax.grid(alpha=0.12)
    ax.tick_params(labelsize=8)


def _scatter_pair(ax, df, x, y, title=None, alpha_all=0.25):
    sc = ax.scatter(df[x], df[y], c=df["chi2"], s=6, alpha=alpha_all)

    ax.set_xlabel(x, fontsize=9)
    ax.set_ylabel(y, fontsize=9)

    if title:
        ax.set_title(title, fontsize=10)

    if x in {"rho_s", "MBH"}:
        ax.set_xscale("log")
    if y in {"rho_s", "MBH"}:
        ax.set_yscale("log")

    _style_axes(ax)
    return sc


def _draw_interpreted_overlay(ax, best, x, y, center, lo, hi):
    ax.scatter(best[x], best[y], c=best["chi2"], s=8, alpha=0.95, edgecolors="none")
    ax.scatter([center[x]], [center[y]], marker="x", s=80, linewidths=2)

    ax.axvline(lo[x], linestyle="--", linewidth=1.0)
    ax.axvline(hi[x], linestyle="--", linewidth=1.0)
    ax.axhline(lo[y], linestyle="--", linewidth=1.0)
    ax.axhline(hi[y], linestyle="--", linewidth=1.0)


# -------------------------
# GEOMETRY
# -------------------------
def plot_geometry_raw(df, outpath, zoom=False):
    pairs = [("rho_s", "r_s"), ("rho_s", "MBH"), ("r_s", "MBH")]

    fig, axes = plt.subplots(3, 1, figsize=(12, 24))

    last_sc = None
    for ax, (x, y) in zip(axes, pairs):
        last_sc = _scatter_pair(ax, df, x, y, title=f"{x} vs {y}")

    if zoom:
        for ax, (x, y) in zip(axes, pairs):
            for param, setter in ((x, ax.set_xlim), (y, ax.set_ylim)):
                pvals = df[param].values
                mask = np.isfinite(pvals)
                if param in {"rho_s", "MBH"}:
                    mask &= pvals > 0
                pvals = pvals[mask]
                if len(pvals) >= 2:
                    setter(np.quantile(pvals, 0.05), np.quantile(pvals, 0.95))

    cbar = fig.colorbar(last_sc, ax=axes, pad=0.01)
    cbar.set_label("chi2")

    fig.suptitle("Raw parameter space" + (" [zoom]" if zoom else ""), fontsize=11)
    _save_and_optionally_show(fig, outpath)


def plot_geometry_interpreted(df, best, center, lo, hi, outpath, zoom=False):
    pairs = [("rho_s", "r_s"), ("rho_s", "MBH"), ("r_s", "MBH")]

    fig, axes = plt.subplots(3, 1, figsize=(12, 24))

    last_sc = None
    for ax, (x, y) in zip(axes, pairs):
        last_sc = ax.scatter(df[x], df[y], c=df["chi2"], s=5, alpha=0.10)
        _draw_interpreted_overlay(ax, best, x, y, center, lo, hi)

        ax.set_title(f"{x} vs {y}", fontsize=10)
        ax.set_xlabel(x)
        ax.set_ylabel(y)

        if x in {"rho_s", "MBH"}:
            ax.set_xscale("log")
        if y in {"rho_s", "MBH"}:
            ax.set_yscale("log")

        _style_axes(ax)

    if zoom:
        for ax, (x, y) in zip(axes, pairs):
            for param, setter in ((x, ax.set_xlim), (y, ax.set_ylim)):
                lo_v, hi_v = lo[param], hi[param]
                if param in {"rho_s", "MBH"} and lo_v > 0 and hi_v > lo_v:
                    setter(lo_v * 0.3, hi_v * 3.0)
                else:
                    span = (hi_v - lo_v) if hi_v > lo_v else abs(hi_v) * 0.2 + 1.0
                    setter(lo_v - 0.3 * span, hi_v + 0.3 * span)

    cbar = fig.colorbar(last_sc, ax=axes, location="right", fraction=0.02, pad=0.04, label="chi2", aspect=20, shrink=0.8, ticks=[df["chi2"].min(), df["chi2"].max()], extend="both", extendfrac=0.1)
    cbar.set_ticklabels(["low", "high"])
    cbar.set_label("chi2")

    fig.suptitle("Interpreted region" + (" [y-zoom]" if zoom else ""), fontsize=11)
    _save_and_optionally_show(fig, outpath)


# -------------------------
# HISTOGRAMS
# -------------------------
def plot_hist_raw(df, outpath, zoom=False):
    fig, axes = plt.subplots(3, 1, figsize=(12, 24))

    for ax, p in zip(axes, PARAMS):
        ax.hist(df[p].values, bins=40, alpha=0.85)

        ax.set_title(p, fontsize=10)
        ax.set_xlabel(p)
        ax.set_ylabel("count")

        if p in {"rho_s", "MBH"}:
            ax.set_xscale("log")

        _style_axes(ax)

    if zoom:
        for ax, p in zip(axes, PARAMS):
            pvals = df[p].values
            mask = np.isfinite(pvals)
            if p in {"rho_s", "MBH"}:
                mask &= pvals > 0
            pvals = pvals[mask]
            if len(pvals) >= 2:
                ax.set_xlim(np.quantile(pvals, 0.05), np.quantile(pvals, 0.95))

    fig.suptitle("Raw distributions" + (" [zoom]" if zoom else ""), fontsize=11)
    _save_and_optionally_show(fig, outpath)


def plot_hist_interpreted(df, best, center, lo, hi, outpath, zoom=False):
    fig, axes = plt.subplots(3, 1, figsize=(12, 24))

    for ax, p in zip(axes, PARAMS):
        ax.hist(df[p].values, bins=40, alpha=0.25, label="all")
        ax.hist(best[p].values, bins=40, alpha=0.85, label="low-chi2")

        ax.axvline(center[p], linewidth=2)
        ax.axvline(lo[p], linestyle="--")
        ax.axvline(hi[p], linestyle="--")

        ax.set_title(p, fontsize=10)
        ax.set_xlabel(p)
        ax.set_ylabel("count")

        if p in {"rho_s", "MBH"}:
            ax.set_xscale("log")

        ax.legend(fontsize=8)
        _style_axes(ax)

    if zoom:
        for ax, p in zip(axes, PARAMS):
            lo_v, hi_v = lo[p], hi[p]
            if p in {"rho_s", "MBH"} and lo_v > 0 and hi_v > lo_v:
                ax.set_xlim(lo_v * 0.3, hi_v * 3.0)
            else:
                span = (hi_v - lo_v) if hi_v > lo_v else abs(hi_v) * 0.2 + 1.0
                ax.set_xlim(lo_v - 0.3 * span, hi_v + 0.3 * span)

    fig.suptitle("Interpreted distributions" + (" [zoom]" if zoom else ""), fontsize=11)
    _save_and_optionally_show(fig, outpath)


# -------------------------
# CHI2 TRENDS
# -------------------------
def plot_chi2_trends_raw(df, outpath, zoom=False):
    _LABEL_COLORS = {
        "full":        "black",
        "no_bh":       "red",
        "no_halo":     "blue",
        "bh_small":    "orange",
        "halo_small":  "green",
    }

    best_full_chi2 = df[df["label"] == "full"]["chi2"].min()
    if not np.isfinite(best_full_chi2):
        best_full_chi2 = df["chi2"].min()
    df = df.copy()
    df["dchi2"] = df["chi2"] - best_full_chi2

    fig, axes = plt.subplots(3, 1, figsize=(12, 24))

    for ax, p in zip(axes, PARAMS):
        for label, sub in df.groupby("label"):
            c = _LABEL_COLORS.get(label, "gray")
            ax.scatter(sub[p], sub["dchi2"], s=6, alpha=0.5, color=c, label=label)

        ax.axhline(0, linewidth=0.8, linestyle=":", color="black")
        ax.legend(fontsize=7)
        ax.set_title(p, fontsize=10)
        ax.set_xlabel(p)
        ax.set_ylabel("Δχ2 (relative to best full)")

        if p in {"rho_s", "MBH"}:
            ax.set_xscale("log")

        _style_axes(ax)

    if zoom:
        dchi2_vals = df["dchi2"].values
        dchi2_vals = dchi2_vals[np.isfinite(dchi2_vals)]
        q95_y = float(np.quantile(dchi2_vals, 0.95))
        pad_y = max(q95_y * 0.05, 1e-10)
        for ax, p in zip(axes, PARAMS):
            ax.set_ylim(-pad_y, q95_y + pad_y)
            pvals = df[p].values
            mask = np.isfinite(pvals)
            if p in {"rho_s", "MBH"}:
                mask &= pvals > 0
            pvals = pvals[mask]
            if len(pvals) >= 2:
                ax.set_xlim(np.quantile(pvals, 0.05), np.quantile(pvals, 0.95))

    fig.suptitle("Raw Δχ2 trends by label" + (" [zoom]" if zoom else ""), fontsize=11)
    _save_and_optionally_show(fig, outpath)


def plot_chi2_trends_interpreted(df, best, center, lo, hi, outpath, zoom=False):
    fig, axes = plt.subplots(3, 1, figsize=(12, 24))

    for ax, p in zip(axes, PARAMS):
        ax.scatter(df[p], df["chi2"], s=5, alpha=0.10, label="all")
        ax.scatter(best[p], best["chi2"], s=8, alpha=0.9, label="low-chi2")

        ax.axvline(center[p], linewidth=2)
        ax.axvline(lo[p], linestyle="--")
        ax.axvline(hi[p], linestyle="--")

        ax.set_title(p, fontsize=10)
        ax.set_xlabel(p)
        ax.set_ylabel("chi2")

        if p in {"rho_s", "MBH"}:
            ax.set_xscale("log")

        ax.legend(fontsize=8)
        _style_axes(ax)

    if zoom:
        chi2_min = float(best["chi2"].min())
        chi2_hi = float(best["chi2"].quantile(0.95))
        pad_y = (chi2_hi - chi2_min) * 0.1 + 1e-10
        for ax, p in zip(axes, PARAMS):
            ax.set_ylim(chi2_min - pad_y, chi2_hi + pad_y)
            lo_v, hi_v = lo[p], hi[p]
            if p in {"rho_s", "MBH"} and lo_v > 0 and hi_v > lo_v:
                ax.set_xlim(lo_v * 0.3, hi_v * 3.0)
            else:
                span = (hi_v - lo_v) if hi_v > lo_v else abs(hi_v) * 0.2 + 1.0
                ax.set_xlim(lo_v - 0.3 * span, hi_v + 0.3 * span)

    fig.suptitle("Interpreted chi2 trends" + (" [zoom]" if zoom else ""), fontsize=11)
    _save_and_optionally_show(fig, outpath)


# -------------------------
# GP SLICES
# -------------------------
def plot_gp_slices(df, center, outpath, zoom=False):
    try:
        gp = fit_gp(df)
    except Exception as e:
        return {"status": "skipped", "reason": f"GP fit failed: {e}"}

    fig, axes = plt.subplots(3, 1, figsize=(12, 24))
    center_vec = np.array([center[p] for p in PARAMS], dtype=float)
    yhats = []

    for i, p in enumerate(PARAMS):
        ax = axes[i]
        vals = df[p].values
        lo = np.quantile(vals, 0.05)
        hi = np.quantile(vals, 0.95)

        if p in {"rho_s", "MBH"} and lo > 0 and hi > lo:
            grid = np.logspace(np.log10(lo), np.log10(hi), 200)
            ax.set_xscale("log")
        else:
            grid = np.linspace(lo, hi, 200)

        X = np.tile(center_vec, (len(grid), 1))
        X[:, i] = grid

        Xp = X.copy()
        for j, q in enumerate(PARAMS):
            if np.all(df[q].values > 0):
                Xp[:, j] = np.log10(Xp[:, j])

        yhat = gp.predict(Xp)
        yhats.append(yhat)
        ax.plot(grid, yhat, linewidth=1.5)
        ax.axvline(center[p], linestyle="--", linewidth=1.0)

        ax.set_title(f"GP slice: chi2 vs {p}", fontsize=10)
        ax.set_xlabel(p)
        ax.set_ylabel("predicted chi2")
        _style_axes(ax)

    if zoom:
        for ax, yhat in zip(axes, yhats):
            pad = (yhat.max() - yhat.min()) * 0.1 + 1e-10
            ax.set_ylim(float(yhat.min()) - pad, float(yhat.max()) + pad)

    fig.suptitle("GP slices" + (" [y-zoom]" if zoom else ""), fontsize=11)
    _save_and_optionally_show(fig, outpath)
    return {"status": "ok"}


# -------------------------
# DIRS AND PLOT WRITER
# -------------------------
def ensure_galaxy_dirs(gal_dir):
    os.makedirs(gal_dir, exist_ok=True)
    os.makedirs(f"{gal_dir}/raw", exist_ok=True)
    os.makedirs(f"{gal_dir}/interpreted", exist_ok=True)
    os.makedirs(f"{gal_dir}/state", exist_ok=True)
    os.makedirs(f"{gal_dir}/patch", exist_ok=True)


def refresh_plot_dirs(gal_dir):
    for sub in ("raw", "interpreted"):
        p = f"{gal_dir}/{sub}"
        if os.path.exists(p):
            shutil.rmtree(p)
        os.makedirs(p, exist_ok=True)


def make_plots(gal_dir, df, best, center, lo, hi, meta):
    refresh_plot_dirs(gal_dir)

    plot_geometry_raw(df, f"{gal_dir}/raw/geometry.png")
    plot_geometry_raw(df, f"{gal_dir}/raw/geometry_zoom.png", zoom=True)
    plot_hist_raw(df, f"{gal_dir}/raw/hist.png")
    plot_hist_raw(df, f"{gal_dir}/raw/hist_zoom.png", zoom=True)
    plot_chi2_trends_raw(df, f"{gal_dir}/raw/chi2.png")
    plot_chi2_trends_raw(df, f"{gal_dir}/raw/chi2_zoom.png", zoom=True)

    plot_geometry_interpreted(df, best, center, lo, hi, f"{gal_dir}/interpreted/geometry.png")
    plot_geometry_interpreted(df, best, center, lo, hi, f"{gal_dir}/interpreted/geometry_zoom.png", zoom=True)
    plot_hist_interpreted(df, best, center, lo, hi, f"{gal_dir}/interpreted/hist.png")
    plot_hist_interpreted(df, best, center, lo, hi, f"{gal_dir}/interpreted/hist_zoom.png", zoom=True)
    plot_chi2_trends_interpreted(df, best, center, lo, hi, f"{gal_dir}/interpreted/chi2.png")
    plot_chi2_trends_interpreted(df, best, center, lo, hi, f"{gal_dir}/interpreted/chi2_zoom.png", zoom=True)

    gp_meta = plot_gp_slices(df, center, f"{gal_dir}/interpreted/gp_slices.png")
    plot_gp_slices(df, center, f"{gal_dir}/interpreted/gp_slices_zoom.png", zoom=True)
    meta["gp_plot"] = gp_meta

    with open(f"{gal_dir}/state/meta.json", "w") as f:
        json.dump(meta, f, indent=2)
# =========================
# MAIN ANALYSIS
# =========================
def analyze(deck_path, galaxy, config_path=None):
    df = load_deck(deck_path)

    if len(df) < MIN_POINTS:
        return {
            "status": "insufficient_data",
            "galaxy": galaxy,
            "n_points": len(df),
            "notes": ["not enough pass points yet for a large-scale cloud summary"],
        }, None, None, None, None, None

    best = compute_support(df)

    center_s = best[PARAMS].median()
    center = {p: float(center_s[p]) for p in PARAMS}

    evidence, lo_s, hi_s = build_evidence(df, best)
    lo = {p: float(lo_s[p]) for p in PARAMS}
    hi = {p: float(hi_s[p]) for p in PARAMS}

    pressure = boundary_pressure(best, df)

    shape_hint = "not evaluated"
    curvature = None
    try:
        gp = fit_gp(df)
        curvature = curvature_proxy(gp, np.array([center[p] for p in PARAMS], dtype=float))
        shape_hint = classify_shape(curvature)
    except Exception:
        shape_hint = "gp_failed"
        curvature = None

    patch_lo, patch_hi = compute_patch_bounds(df, best, center)

    notes = build_notes(evidence, pressure, shape_hint, curvature, df, best, center, patch_lo, patch_hi)
    patch_text, patch_script = build_config_patch(galaxy, center, patch_lo, patch_hi)

    result = {
        "status": "ok",
        "galaxy": galaxy,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "n_points": int(len(df)),
        "support_points": int(len(best)),
        "support_fraction": float(len(best) / len(df)),
        "center": center,
        "range_16_84": {p: [lo[p], hi[p]] for p in PARAMS},
        "boundary_pressure": pressure,
        "evidence": evidence,
        "gp_shape_hint": shape_hint,
        "gp_curvature_proxy": curvature,
        "notes": notes,
        "config_patch_text": patch_text,
        "config_patch_cmd": patch_script,
        "deck_path": deck_path,
        "config_path": config_path,
        "lom_mode": LOM_MODE,
    }

    return result, df, best, center, lo, hi


# =========================
# WRITE OUTPUTS
# =========================
def write_outputs(galaxy, gal_dir, result):
    with open(f"{gal_dir}/state/well_state.json", "w") as f:
        json.dump(result, f, indent=2)
    if "config_patch_text" in result:
        with open(f"{gal_dir}/patch/config_patch.txt", "w") as f:
            f.write(result["config_patch_text"] + "\n")
    if "config_patch_cmd" in result:
        patch_path = f"{gal_dir}/patch/config_patch.sh"
        with open(patch_path, "w") as f:
            f.write(result["config_patch_cmd"])
        os.chmod(patch_path, 0o755)


def update_shared(galaxy):
    os.makedirs(SHARED_ROOT, exist_ok=True)

    with open(f"{SHARED_ROOT}/latest_galaxy", "w") as f:
        f.write(galaxy)

    with open(f"{SHARED_ROOT}/latest_update.json", "w") as f:
        json.dump(
            {
                "galaxy": galaxy,
                "timestamp": datetime.now().isoformat(timespec="seconds"),
                "epoch": time.time(),
                "lom_mode": LOM_MODE,
            },
            f,
            indent=2,
        )


def run_once():
    os.makedirs(COCKPIT_ROOT, exist_ok=True)
    os.makedirs(SHARED_ROOT, exist_ok=True)

    galaxy, deck, config, gal_dir = get_paths()
    ensure_galaxy_dirs(gal_dir)

    print(f"[LOM] snapshot | galaxy={galaxy}")

    try:
        result, df, best, center, lo, hi = analyze(deck, galaxy, config_path=config)

        if result["status"] == "ok":
            if PLOT_ENABLED:
                meta = {
                    "galaxy": galaxy,
                    "timestamp": result["timestamp"],
                    "n_points": result["n_points"],
                    "support_points": result["support_points"],
                    "support_fraction": result["support_fraction"],
                    "lom_mode": LOM_MODE,
                }
                make_plots(gal_dir, df, best, center, lo, hi, meta)

            write_outputs(galaxy, gal_dir, result)
            update_shared(galaxy)

            print(
                f"[LOM] snapshot complete | "
                f"rho_s~{center['rho_s']:.3e} "
                f"r_s~{center['r_s']:.3e} "
                f"MBH~{center['MBH']:.3e}"
            )
        else:
            write_outputs(galaxy, gal_dir, result)
            update_shared(galaxy)
            print(f"[LOM] insufficient data | n={result.get('n_points')}")

    except Exception as e:
        print(f"[LOM ERROR] {e}")


# =========================
# LOOP
# =========================
def run_loop(interval=60):
    os.makedirs(COCKPIT_ROOT, exist_ok=True)
    os.makedirs(SHARED_ROOT, exist_ok=True)

    last_counts = {}
    last_mtimes = {}

    while True:
        print(f"[LOM] tick {time.strftime('%H:%M:%S')} mode={LOM_MODE}")

        galaxy, deck, config, gal_dir = get_paths()
        ensure_galaxy_dirs(gal_dir)

        try:
            deck_mtime = os.path.getmtime(deck) if os.path.exists(deck) else None
            prev_mtime = last_mtimes.get(galaxy, None)

            result, df, best, center, lo, hi = analyze(deck, galaxy, config_path=config)

            if result["status"] == "ok":
                current_n = result["n_points"]
                prev_n = last_counts.get(galaxy, None)

                if prev_n is None:
                    print(f"[LOM] {galaxy} initial load with n={current_n}")
                elif current_n == prev_n and deck_mtime == prev_mtime:
                    print(f"[LOM] {galaxy} | no new data | n={current_n}")
                else:
                    delta = current_n - prev_n if prev_n is not None else current_n
                    print(f"[LOM] {galaxy} | +{delta} pass rows | n={current_n}")

                if PLOT_ENABLED:
                    meta = {
                        "galaxy": galaxy,
                        "timestamp": result["timestamp"],
                        "n_points": result["n_points"],
                        "support_points": result["support_points"],
                        "support_fraction": result["support_fraction"],
                        "lom_mode": LOM_MODE,
                    }
                    make_plots(gal_dir, df, best, center, lo, hi, meta)
                    print(f"[LOM] plots refreshed: {gal_dir}")

                write_outputs(galaxy, gal_dir, result)
                update_shared(galaxy)

                print(
                    f"[LOM] {galaxy} | "
                    f"rho_s~{center['rho_s']:.3e} "
                    f"r_s~{center['r_s']:.3e} "
                    f"MBH~{center['MBH']:.3e}"
                )
                for line in result["notes"]:
                    print(f"  - {line}")

                print("\n[LOM] config patch script:")
                print(f"{gal_dir}/patch/config_patch.sh")

                last_counts[galaxy] = current_n
                last_mtimes[galaxy] = deck_mtime

            else:
                print(f"[LOM] {galaxy} | {result['status']} | n={result.get('n_points')}")
                write_outputs(galaxy, gal_dir, result)
                update_shared(galaxy)

        except Exception as e:
            print(f"[LOM ERROR] {e}")

        time.sleep(interval)


# =========================
# ENTRY
# =========================
if __name__ == "__main__":
    if IS_LOOP:
        run_loop()
    else:
        run_once()