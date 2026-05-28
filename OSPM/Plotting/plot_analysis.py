import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse, glob, os, re
from pathlib import Path

plt.style.use("dark_background")

# --------------------------------------------------
# Resolve repo root
# --------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[2]

# --------------------------------------------------
# Resolve galaxy from repo root if present
# --------------------------------------------------
WHICH = REPO_ROOT / "which_galaxy"

if WHICH.exists():
    GALAXY = WHICH.read_text().strip()
    if not GALAXY:
        GALAXY = "unknown"
else:
    GALAXY = "unknown"

# --------------------------------------------------
# Defaults
# --------------------------------------------------
DEFAULT_PROFILE_DIR = REPO_ROOT / "Data" / "Gal_Profiles" / GALAXY / "default"

# --------------------------------------------------
# Helpers
# --------------------------------------------------
def _timestamp_from_filename(path):
    m = re.search(r"(\d{8})_(\d{6})", os.path.basename(path))
    return int(m.group(1) + m.group(2)) if m else None


def find_latest_daemon_deck_csv(profile_dir, pattern="*daemon_deck*.csv"):
    profile_dir = Path(profile_dir)

    if not profile_dir.exists():
        raise FileNotFoundError(f"Profile directory not found:\n  {profile_dir}")

    candidates = glob.glob(str(profile_dir / pattern))

    if not candidates:
        raise FileNotFoundError(f"No daemon_deck CSV found in:\n  {profile_dir}")

    stamped = [(_timestamp_from_filename(p), p) for p in candidates]
    stamped = [x for x in stamped if x[0] is not None]

    if stamped:
        return max(stamped, key=lambda x: x[0])[1]

    return max(candidates, key=os.path.getmtime)


def resolve_path(path):
    path = Path(path)

    if path.is_absolute():
        return path

    return REPO_ROOT / path


def infer_galaxy_from_csv(csv_file):
    csv_file = Path(csv_file).resolve()

    try:
        rel = csv_file.relative_to(REPO_ROOT / "Data" / "Gal_Profiles")
        return rel.parts[0]
    except Exception:
        return GALAXY


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument(
        "csv_positional",
        nargs="?",
        default=None,
        help="Optional daemon deck CSV path."
    )

    p.add_argument(
        "--csv",
        default=None,
        help="Daemon deck CSV path. Can be absolute or relative to repo root."
    )

    p.add_argument(
        "--profile-dir",
        default=None,
        help="Profile directory containing a daemon deck. Can be absolute or relative to repo root."
    )

    p.add_argument(
        "--pattern",
        default="*daemon_deck*.csv",
        help="Glob pattern used when searching profile-dir."
    )

    p.add_argument(
        "--save",
        action="store_true",
        help="Save plots instead of showing them."
    )

    p.add_argument(
        "--out",
        default=None,
        help="Output image path for the main landscape. Diagnostic plot gets _diagnostics suffix."
    )

    p.add_argument(
        "--zoom-dchi",
        type=float,
        default=5.0,
        help="Zoom cut: keep rows with chi2 < best + zoom_dchi."
    )

    return p.parse_args()


def require_columns(df, cols):
    missing = [c for c in cols if c not in df.columns]

    if missing:
        raise KeyError(
            "Missing required column(s): "
            + ", ".join(missing)
            + "\nAvailable columns:\n  "
            + "\n  ".join(df.columns)
        )


def finite_positive_floor(values):
    values = pd.to_numeric(values, errors="coerce")
    good = values[np.isfinite(values) & (values > 0)]

    if len(good) == 0:
        return 1.0

    return good.min() / 3.0


def numeric_clean(df, cols):
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.replace([np.inf, -np.inf], np.nan)


def split_bh(df, mbh_floor):
    bh0 = df["MBH"] == 0
    blue = df[bh0].copy()
    teal = df[~bh0].copy()

    blue["MBH_plot"] = mbh_floor
    teal["MBH_plot"] = teal["MBH"]

    return blue, teal


def add_mbh_zero_marker(ax, mbh_floor):
    ax.axvline(mbh_floor, color="#4FA3FF", alpha=0.45, linestyle="--")
    ymin, ymax = ax.get_ylim()
    ax.text(
        mbh_floor,
        ymin,
        " MBH=0",
        color="#4FA3FF",
        fontsize=9,
        rotation=90,
        va="bottom",
    )


def scatter_split(ax, teal, blue, x, y, s_teal=10, s_blue=16, label=True):
    ax.scatter(
        teal[x],
        teal[y],
        s=s_teal,
        c="#39EB33",
        label="BH models" if label else None,
    )
    ax.scatter(
        blue[x],
        blue[y],
        s=s_blue,
        c="#4FA3FF",
        label="MBH = 0" if label else None,
    )


def set_log_if_positive(ax, values):
    values = pd.to_numeric(values, errors="coerce")
    good = values[np.isfinite(values) & (values > 0)]
    if len(good) > 0:
        ax.set_xscale("log")


def save_or_show(fig, args, csv_file, suffix, default_name):
    if args.save:
        if args.out and suffix == "landscape":
            out = resolve_path(args.out)
        else:
            out = Path(csv_file).parent / default_name

        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=250, bbox_inches="tight")
        print(f"[SAVE] {out}")
    else:
        plt.show()


# --------------------------------------------------
# Load data
# --------------------------------------------------
args = parse_args()

if args.csv:
    csv_file = resolve_path(args.csv)

elif args.csv_positional:
    csv_file = resolve_path(args.csv_positional)

else:
    profile_dir = resolve_path(args.profile_dir) if args.profile_dir else DEFAULT_PROFILE_DIR
    csv_file = Path(find_latest_daemon_deck_csv(profile_dir, args.pattern))

if not csv_file.exists():
    raise FileNotFoundError(f"CSV file not found:\n  {csv_file}")

galaxy_label = infer_galaxy_from_csv(csv_file)

print(f"[PLOT] galaxy = {galaxy_label}")
print(f"[PLOT] using  = {csv_file}")

df = pd.read_csv(csv_file)

required = ["chi2", "MBH", "rho_s", "r_s", "ML"]
require_columns(df, required)

diagnostic_cols = [
    "chi2_inner",
    "chi2_outer",
    "chi2_occ",
    "N_inner",
    "N_outer",
    "N_nonzero_weights",
    "effective_N_orbits",
    "max_weight_fraction",
]

all_numeric = required + [c for c in diagnostic_cols if c in df.columns]
df = numeric_clean(df, all_numeric)

df = df.dropna(subset=required)

if len(df) == 0:
    raise ValueError("No usable rows left after cleaning chi2, MBH, rho_s, r_s, and ML.")

best = df["chi2"].min()
zoom = df[df["chi2"] < best + args.zoom_dchi].copy()

if len(zoom) == 0:
    zoom = df.nsmallest(min(50, len(df)), "chi2").copy()

mbh_floor = finite_positive_floor(df["MBH"])

df["MBH_plot"] = df["MBH"]
df.loc[df["MBH"] <= 0, "MBH_plot"] = mbh_floor

zoom["MBH_plot"] = zoom["MBH"]
zoom.loc[zoom["MBH"] <= 0, "MBH_plot"] = mbh_floor

blue, teal = split_bh(df, mbh_floor)
zoom_blue, zoom_teal = split_bh(zoom, mbh_floor)

print(f"[INFO] rows      = {len(df)}")
print(f"[INFO] best chi2 = {best}")
print(f"[INFO] zoom rows = {len(zoom)}")
print(f"[INFO] MBH floor = {mbh_floor}")

# --------------------------------------------------
# Figure 1: total chi2 landscape, now including ML
# --------------------------------------------------
fig, ax = plt.subplots(2, 4, figsize=(24, 10))

params = [
    ("MBH_plot", r"$M_{\rm BH}$", "MBH"),
    ("rho_s", r"$\rho_s$", "rho_s"),
    ("r_s", r"$r_s$", "r_s"),
    ("ML", r"$M/L$", "ML"),
]

for j, (xplot, label, real_col) in enumerate(params):

    scatter_split(ax[0, j], teal, blue, xplot, "chi2", label=True)
    set_log_if_positive(ax[0, j], df[xplot])

    ax[0, j].set_xlabel(label)
    ax[0, j].set_ylabel(r"$\chi^2$")
    ax[0, j].set_title(f"{galaxy_label} : $\\chi^2$ vs {label}")

    scatter_split(ax[1, j], zoom_teal, zoom_blue, xplot, "chi2", s_teal=20, s_blue=24, label=False)
    set_log_if_positive(ax[1, j], zoom[xplot])

    ax[1, j].set_xlabel(label)
    ax[1, j].set_ylabel(r"$\chi^2$")
    ax[1, j].set_title(f"{galaxy_label} : zoomed {label}")

    if real_col == "MBH":
        add_mbh_zero_marker(ax[0, j], mbh_floor)
        add_mbh_zero_marker(ax[1, j], mbh_floor)

for a in ax[0]:
    y = df["chi2"].values
    ypad = 0.05 * (np.nanmax(y) - np.nanmin(y))
    if ypad == 0:
        ypad = 1.0
    a.set_ylim(np.nanmin(y) - ypad, np.nanmax(y) + ypad)
    a.legend(loc="best", fontsize=8)

z = zoom["chi2"].values
zpad = 0.05 * (np.nanmax(z) - np.nanmin(z))
if zpad == 0:
    zpad = 1.0

for a in ax[1]:
    a.set_ylim(np.nanmin(z) - zpad, np.nanmax(z) + zpad)

fig.suptitle(
    f"{galaxy_label} OSPM χ² landscape\n"
    f"source: {Path(csv_file).name}",
    fontsize=14,
)

plt.tight_layout()

save_or_show(
    fig,
    args,
    csv_file,
    "landscape",
    f"{Path(csv_file).stem}_chi2_landscape.png",
)

# --------------------------------------------------
# Figure 2: diagnostic panels
# --------------------------------------------------
has_diag = all(c in df.columns for c in [
    "chi2_inner",
    "chi2_outer",
    "effective_N_orbits",
    "max_weight_fraction",
])

if has_diag:
    diag_required = ["chi2_inner", "chi2_outer", "effective_N_orbits", "max_weight_fraction"]
    diag_df = df.dropna(subset=diag_required).copy()

    if len(diag_df) > 0:
        diag_df["MBH_plot"] = diag_df["MBH"]
        diag_df.loc[diag_df["MBH"] <= 0, "MBH_plot"] = mbh_floor
        diag_blue, diag_teal = split_bh(diag_df, mbh_floor)

        fig2, ax2 = plt.subplots(2, 4, figsize=(24, 10))

        diag_y = [
            ("chi2", r"total $\chi^2$", r"total $\chi^2$"),
            ("chi2_inner", r"inner $\chi^2$", r"inner stars"),
            ("chi2_outer", r"outer $\chi^2$", r"outer stars"),
            ("chi2_occ", r"occupancy diagnostic", r"occupancy"),
            ("effective_N_orbits", r"$N_{\rm eff}$ orbits", r"orbit-weight spread"),
            ("max_weight_fraction", r"max weight fraction", r"largest orbit weight"),
            ("N_nonzero_weights", r"nonzero orbit weights", r"active orbit count"),
            ("ML", r"$M/L$", r"stellar mass-to-light"),
        ]

        for k, (ycol, ylabel, title) in enumerate(diag_y):
            r = k // 4
            c = k % 4

            if ycol not in diag_df.columns:
                ax2[r, c].set_visible(False)
                continue

            scatter_split(ax2[r, c], diag_teal, diag_blue, "MBH_plot", ycol, s_teal=14, s_blue=22, label=(k == 0))
            set_log_if_positive(ax2[r, c], diag_df["MBH_plot"])
            ax2[r, c].set_xlabel(r"$M_{\rm BH}$")
            ax2[r, c].set_ylabel(ylabel)
            ax2[r, c].set_title(title)
            add_mbh_zero_marker(ax2[r, c], mbh_floor)

            if ycol == "max_weight_fraction":
                ax2[r, c].set_ylim(0.0, min(1.0, max(0.05, 1.1 * np.nanmax(diag_df[ycol]))))

            if k == 0:
                ax2[r, c].legend(loc="best", fontsize=8)

        fig2.suptitle(
            f"{galaxy_label} OSPM diagnostic landscape\n"
            f"source: {Path(csv_file).name}",
            fontsize=14,
        )

        plt.tight_layout()

        save_or_show(
            fig2,
            args,
            csv_file,
            "diagnostics",
            f"{Path(csv_file).stem}_diagnostics.png",
        )

    else:
        print("[SKIP] Diagnostic columns exist, but no finite diagnostic rows were found.")
else:
    print("[SKIP] Diagnostic plot skipped. Missing one or more diagnostic columns.")