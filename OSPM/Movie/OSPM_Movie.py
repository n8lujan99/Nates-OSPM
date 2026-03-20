import os, sys, subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

G = 4.30091e-6
KMS_TO_KPC_MYR = 1.0227121650537077e-3
KMS2_TO_KPC2_MYR2 = KMS_TO_KPC_MYR**2

def load_best_theta():
    repo_root = Path(__file__).resolve().parents[2]
    galaxy = (repo_root / "which_galaxy").read_text().strip()
    gal_dir = repo_root / "Data" / "Gal_Profiles" / galaxy / "default"
    df = pd.read_csv(gal_dir / "daemon_deck.csv")
    row = df.loc[df["chi2"].idxmin()]
    return galaxy, float(row["rho_s"]), float(row["r_s"]), float(row["MBH"])

def nfw_Menc(r, rho_s, r_s):
    x = r / r_s
    return 4 * np.pi * rho_s * r_s**3 * (np.log(1 + x) - x / (1 + x))

def accel_xy_vec(x, y, rho_s, r_s, MBH):
    r = np.sqrt(x*x + y*y) + 1e-12
    x_nfw = r / r_s
    Menc = 4 * np.pi * rho_s * r_s**3 * (np.log(1 + x_nfw) - x_nfw / (1 + x_nfw))
    a = (-G * MBH / (r*r) - G * Menc / (r*r)) * KMS2_TO_KPC2_MYR2
    return a * (x / r), a * (y / r)

def integrate_orbits_vec(x, y, vx, vy, rho_s, r_s, MBH, dt, steps):
    N = len(x)
    traj_x = np.zeros((steps, N), dtype=np.float32)
    traj_y = np.zeros((steps, N), dtype=np.float32)
    ax, ay = accel_xy_vec(x, y, rho_s, r_s, MBH)
    vx += 0.5 * dt * ax;  vy += 0.5 * dt * ay
    for i in range(steps):
        x += dt * vx;  y += dt * vy
        ax, ay = accel_xy_vec(x, y, rho_s, r_s, MBH)
        vx += dt * ax;  vy += dt * ay
        traj_x[i] = x;  traj_y[i] = y
    return traj_x, traj_y

def main():
    galaxy, rho_s, r_s, MBH = load_best_theta()
    print("Galaxy:", galaxy)
    print("theta:", rho_s, r_s, MBH)
    np.random.seed(0)
    N_orbits = 650;  dt = 5e-4;  steps = 12000;  frame_skip = 5;  trail_len = 700
    r_min = 0.015;   r_max = 0.060
    r = np.linspace(r_min, r_max, N_orbits)
    v_circ = np.sqrt(G * (nfw_Menc(r, rho_s, r_s) + MBH) / r)
    v = v_circ * (0.75 + 0.45 * np.random.rand(N_orbits))
    theta = 2 * np.pi * np.random.rand(N_orbits)
    x  = (r * np.cos(theta)).astype(np.float64)
    y  = (r * np.sin(theta)).astype(np.float64)
    vx = (-v * np.sin(theta) * KMS_TO_KPC_MYR).astype(np.float64)
    vy = ( v * np.cos(theta) * KMS_TO_KPC_MYR).astype(np.float64)
    traj_x, traj_y = integrate_orbits_vec(x, y, vx, vy, rho_s, r_s, MBH, dt, steps)
    all_x = traj_x[::200].ravel();  all_y = traj_y[::200].ravel()
    lim = 0.8 * np.max(np.sqrt(all_x**2 + all_y**2))
    r_infl = (G * MBH) / (6.0**2)
    bins = np.linspace(r_min, r_max, 8)
    star_colors = np.clip(np.random.rand(N_orbits, 3) * 1.2, 0, 1)
    fig = plt.figure(figsize=(8, 8), facecolor="black", dpi=100)
    ax = plt.gca()
    ax.set_facecolor("black");  ax.set_xlim(-lim, lim);  ax.set_ylim(-lim, lim)
    ax.set_aspect("equal");     ax.axis("off")
    ax.plot(0, 0, "o", color="white", markersize=3)
    ax.add_patch(plt.Circle((0, 0), r_infl, color="red", fill=False, linestyle="--", linewidth=1.2))
    for b in bins:
        ax.add_patch(plt.Circle((0, 0), b, color=(0.7, 0.7, 1, 0.12), fill=False, linewidth=0.6))
    lc = LineCollection([], linewidths=0.8, antialiaseds=False)
    ax.add_collection(lc)
    out = Path(__file__).resolve().parents[2] / "outputs" / f"{galaxy}_orbits_v3.mp4"
    cmd = [
        "ffmpeg", "-y",
        "-f", "rawvideo", "-pix_fmt", "rgba", "-s", "800x800", "-r", "60", "-i", "-",
        "-vcodec", "h264_nvenc", "-preset", "p5", "-cq", "18", "-pix_fmt", "yuv420p",
        str(out)
    ]
    pipe = subprocess.Popen(cmd, stdin=subprocess.PIPE)
    canvas = fig.canvas
    for frame in range(0, steps, frame_skip):
        j0 = max(0, frame - trail_len);  j1 = frame
        xs = traj_x[j0:j1];  ys = traj_y[j0:j1]
        segs = [np.column_stack((xs[:, i], ys[:, i])) for i in range(N_orbits)]
        cols = [(*star_colors[i], 0.35) for i in range(N_orbits)]
        lc.set_segments(segs);  lc.set_color(cols)
        if j1 > 0:
            ax.scatter(traj_x[j1-1], traj_y[j1-1], s=3, c=star_colors, alpha=0.9, edgecolors="none")
        canvas.draw()
        pipe.stdin.write(np.frombuffer(canvas.buffer_rgba(), dtype=np.uint8).tobytes())
    pipe.stdin.close();  pipe.wait()
    print("Wrote:", out)
if __name__ == "__main__":
    main()