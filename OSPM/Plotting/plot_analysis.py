import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse, glob, os, re
from pathlib import Path

from OSPM.load_config import load_config, get_profile_root

plt.style.use("dark_background")

# --------------------------------------------------
# Context
# --------------------------------------------------

cfg          = load_config()
GALAXY       = cfg["GALAXY"]
PROFILE_ROOT = get_profile_root()

DEFAULT_PROFILE_DIR = PROFILE_ROOT

# --------------------------------------------------
# Helpers
# --------------------------------------------------

def _timestamp_from_filename(path):
    m = re.search(r"(\d{8})_(\d{6})", os.path.basename(path))
    return int(m.group(1) + m.group(2)) if m else None

def find_latest_daemon_deck_csv(profile_dir, pattern="*daemon_deck*"):
    c = glob.glob(os.path.join(profile_dir, pattern))
    if not c:
        raise FileNotFoundError(f"No daemon_deck in {profile_dir}")
    ts = [(_timestamp_from_filename(p), p) for p in c]
    ts = [t for t in ts if t[0] is not None]
    return max(ts, key=lambda x: x[0])[1] if ts else max(c, key=os.path.getmtime)

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--profile-dir", default=DEFAULT_PROFILE_DIR)
    p.add_argument("--csv", default=None)
    p.add_argument("--pattern", default="*daemon_deck*")
    return p.parse_args()

# --------------------------------------------------
# Load data
# --------------------------------------------------

args = parse_args()
csv_file = args.csv if args.csv else find_latest_daemon_deck_csv(
    args.profile_dir, args.pattern
)

print(f"[PLOT] galaxy={GALAXY}")
print(f"[PLOT] using {csv_file}")

df = pd.read_csv(csv_file)

df["chi2"] = pd.to_numeric(df["chi2"], errors="coerce")
df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["chi2"])

for c in ("MBH", "rho_s", "r_s"):
    if c not in df:
        raise KeyError(f"Missing {c}")

# --------------------------------------------------
# Split + zoom
# --------------------------------------------------

bh0  = df["MBH"] == 0
blue = df[bh0]
teal = df[~bh0]

best = df["chi2"].min()
zoom = df[df["chi2"] < best * 1.2]

# --------------------------------------------------
# Plot
# --------------------------------------------------

fig, ax = plt.subplots(2, 3, figsize=(18, 10))

params = [
    ("MBH",   r"$M_{\rm BH}$"),
    ("rho_s", r"$\rho_s$"),
    ("r_s",   r"$r_s$"),
]

for j, (x, label) in enumerate(params):
    ax[0, j].scatter(teal[x], teal["chi2"], s=10, c="#39EB33")
    ax[0, j].scatter(blue[x], blue["chi2"], s=16, c="#4FA3FF")
    ax[0, j].set_xscale("log")
    ax[0, j].set_title(f"{GALAXY} : $\\chi^2$ vs {label}")

    ax[1, j].scatter(zoom[x], zoom["chi2"], s=20, c="#39EB33")
    ax[1, j].set_xscale("log")
    ax[1, j].set_title(f"{GALAXY} : zoomed {label}")

y = df["chi2"].values
ypad = 0.05 * (y.max() - y.min())
for a in ax[0]:
    a.set_ylim(y.min() - ypad, y.max() + ypad)

z = zoom["chi2"].values
zpad = 0.05 * (z.max() - z.min())
for a in ax[1]:
    a.set_ylim(z.min() - zpad, z.max() + zpad)

fig.suptitle(
    f"{GALAXY} OSPM χ² landscape\n"
    f"source: {Path(csv_file).name}",
    fontsize=14
)

plt.tight_layout()
plt.show()
