import os,sys,subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# ---------------------------------------------------------
# Units
# ---------------------------------------------------------
G = 4.30091e-6                         # kpc (km/s)^2 / Msun
KMS_TO_KPC_MYR = 1.0227121650537077e-3
KMS2_TO_KPC2_MYR2 = KMS_TO_KPC_MYR**2

# ---------------------------------------------------------
# Load best-fit halo parameters
# ---------------------------------------------------------
def load_best_theta():
    repo_root = Path(__file__).resolve().parents[2]
    which_file = repo_root/"which_galaxy"
    galaxy = which_file.read_text().strip()
    gal_dir = repo_root/"Data"/"Gal_Profiles"/galaxy/"default"
    deck_file = gal_dir/"daemon_deck.csv"
    print("Using deck:",deck_file)
    df = pd.read_csv(deck_file)
    row = df.loc[df["chi2"].idxmin()]
    return galaxy,float(row["rho_s"]),float(row["r_s"]),float(row["MBH"])

# ---------------------------------------------------------
# NFW enclosed mass
# ---------------------------------------------------------
def nfw_Menc(r,rho_s,r_s):
    x = r/r_s
    return 4*np.pi*rho_s*(r_s**3)*(np.log(1+x)-x/(1+x))

# ---------------------------------------------------------
# Vectorized acceleration
# ---------------------------------------------------------
def accel_xy_vec(x,y,rho_s,r_s,MBH):
    r = np.sqrt(x*x+y*y)+1e-12
    a_bh = -G*MBH/(r*r)
    x_nfw = r/r_s
    Menc = 4*np.pi*rho_s*(r_s**3)*(np.log(1+x_nfw)-x_nfw/(1+x_nfw))
    a_h = -G*Menc/(r*r)
    a = (a_bh + a_h) * KMS2_TO_KPC2_MYR2
    ax = a*(x/r)
    ay = a*(y/r)
    return ax,ay

# ---------------------------------------------------------
# Vectorized leapfrog integrator
# ---------------------------------------------------------
def integrate_orbits_vec(x,y,vx,vy,rho_s,r_s,MBH,dt,steps):
    N = len(x)
    traj_x = np.zeros((steps,N),dtype=np.float32)
    traj_y = np.zeros((steps,N),dtype=np.float32)
    ax,ay = accel_xy_vec(x,y,rho_s,r_s,MBH)
    vx = vx + 0.5*dt*ax
    vy = vy + 0.5*dt*ay
    for i in range(steps):
        x = x + dt*vx
        y = y + dt*vy
        ax,ay = accel_xy_vec(x,y,rho_s,r_s,MBH)
        vx = vx + dt*ax
        vy = vy + dt*ay
        traj_x[i] = x
        traj_y[i] = y
    vx = vx - 0.5*dt*ax
    vy = vy - 0.5*dt*ay
    return traj_x,traj_y

# ---------------------------------------------------------
# Colors (fast, stable per orbit)
# ---------------------------------------------------------
def hsv_to_rgb(h,s,v):
    # h,s,v in [0,1]
    i = np.floor(h*6).astype(int)
    f = h*6 - i
    p = v*(1-s)
    q = v*(1-f*s)
    t = v*(1-(1-f)*s)
    i = i % 6
    r = np.select([i==0,i==1,i==2,i==3,i==4,i==5],[v,q,p,p,t,v])
    g = np.select([i==0,i==1,i==2,i==3,i==4,i==5],[t,v,v,q,p,p])
    b = np.select([i==0,i==1,i==2,i==3,i==4,i==5],[p,p,t,v,v,q])
    return np.stack([r,g,b],axis=1)

# ---------------------------------------------------------
# Main
# ---------------------------------------------------------
def main():
    galaxy,rho_s,r_s,MBH = load_best_theta()
    print("Galaxy:",galaxy)
    print("Using theta:",rho_s,r_s,MBH)

    np.random.seed(0)

    # -----------------------------
    # Speed knobs
    # -----------------------------
    N_orbits   = 800          # drop this if you want a bigger jump
    dt         = 5e-4
    steps      = 20000
    frame_skip = 6            # bigger = faster encode + render
    trail_len  = 900          # only draw last trail_len samples per orbit

    # Optional: if you want *way* faster, do this instead:
    # steps = 12000
    # frame_skip = 8
    # trail_len = 600

    r_min = 0.015
    r_max = 0.060

    r = np.linspace(r_min,r_max,N_orbits)
    Menc = nfw_Menc(r,rho_s,r_s)
    v_circ = np.sqrt(G*(Menc+MBH)/r)
    v = v_circ*(0.8+0.4*np.random.rand(N_orbits))

    theta = 2*np.pi*np.random.rand(N_orbits)
    x = (r*np.cos(theta)).astype(np.float64)
    y = (r*np.sin(theta)).astype(np.float64)
    vx = (-v*np.sin(theta)*KMS_TO_KPC_MYR).astype(np.float64)
    vy = ( v*np.cos(theta)*KMS_TO_KPC_MYR).astype(np.float64)

    traj_x,traj_y = integrate_orbits_vec(x,y,vx,vy,rho_s,r_s,MBH,dt,steps)

    # Plot limits
    all_x = traj_x[::max(1,steps//2000)].ravel()
    all_y = traj_y[::max(1,steps//2000)].ravel()
    rmax = np.max(np.sqrt(all_x*all_x+all_y*all_y))
    lim = 1.1*rmax

    # Influence radius
    v_ref = 6.0
    r_infl = (G*MBH)/(v_ref*v_ref)

    # Figure
    fig = plt.figure(figsize=(8,8),facecolor="black",dpi=100)
    ax = plt.gca()
    ax.set_facecolor("black")
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)
    ax.set_aspect("equal")
    ax.axis("off")

    # Center + influence ring
    ax.plot(0,0,"o",color="white",markersize=3)
    infl = plt.Circle((0,0),r_infl,color="red",fill=False,linestyle="--",linewidth=1.0,alpha=0.8)
    ax.add_patch(infl)

    # Stable per-orbit colors
    h = (np.linspace(0,1,N_orbits,endpoint=False) + 0.17*np.random.rand(N_orbits)) % 1.0
    s = 0.70 + 0.25*np.random.rand(N_orbits)
    vcol = 0.85 + 0.15*np.random.rand(N_orbits)
    rgb = hsv_to_rgb(h,s,vcol)
    colors = np.concatenate([rgb, 0.28*np.ones((N_orbits,1))],axis=1)  # RGBA with global alpha

    # One artist instead of 800 Line2D objects
    lc = LineCollection([], linewidths=0.65, antialiaseds=False)
    lc.set_colors(colors)
    ax.add_collection(lc)

    out_mp4 = Path(__file__).resolve().parents[2]/"outputs"/f"{galaxy}_orbits_nvenc.mp4"
    cmd=[
        "ffmpeg","-y",
        "-f","rawvideo","-vcodec","rawvideo",
        "-pix_fmt","rgba","-s","800x800",
        "-r","60","-i","-",
        "-an",
        "-vcodec","h264_nvenc",
        "-preset","p5",
        "-cq","19",
        "-pix_fmt","yuv420p",
        str(out_mp4)
    ]
    print("Launching NVENC encode...")
    pipe = subprocess.Popen(cmd,stdin=subprocess.PIPE)

    canvas = fig.canvas

    # Render loop
    # We only build polylines for the last trail_len points. That caps work.
    for frame in range(0,steps,frame_skip):
        j0 = max(0, frame - trail_len)
        j1 = frame

        if j1 - j0 < 2:
            # nothing to draw yet
            canvas.draw()
            buf = np.frombuffer(canvas.buffer_rgba(),dtype=np.uint8)
            pipe.stdin.write(buf.tobytes())
            continue

        # Build per-orbit polyline arrays for LineCollection
        xs = traj_x[j0:j1]   # (L,N)
        ys = traj_y[j0:j1]

        # list of (L,2) arrays, one per orbit
        segs = [np.column_stack((xs[:,i],ys[:,i])) for i in range(N_orbits)]
        lc.set_segments(segs)

        canvas.draw()
        buf = np.frombuffer(canvas.buffer_rgba(),dtype=np.uint8)
        pipe.stdin.write(buf.tobytes())

    pipe.stdin.close()
    pipe.wait()
    print("Wrote:",out_mp4)

if __name__ == "__main__":
    main()