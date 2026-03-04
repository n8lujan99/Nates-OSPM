import os,sys,subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
G=4.30091e-6
def load_best_theta():
    repo_root=Path(__file__).resolve().parents[2]
    which_file=repo_root/"which_galaxy"
    if not which_file.exists():
        raise RuntimeError(f"which_galaxy not found at {which_file}")
    galaxy=which_file.read_text().strip()
    gal_dir=repo_root/"Data"/"Gal_Profiles"/galaxy/"default"
    if not gal_dir.exists():
        raise RuntimeError(f"Galaxy folder not found: {gal_dir}")
    deck_file=gal_dir/"daemon_deck.csv"
    if not deck_file.exists():
        raise RuntimeError(f"daemon_deck.csv not found: {deck_file}")
    print("Using deck:",deck_file)
    df=pd.read_csv(deck_file)
    for col in ("rho_s","r_s","MBH","chi2"):
        if col not in df.columns:
            raise RuntimeError(f"Missing column in daemon_deck.csv: {col}")
    row=df.loc[df["chi2"].idxmin()]
    return galaxy,float(row["rho_s"]),float(row["r_s"]),float(row["MBH"])
def nfw_Menc(r,rho_s,r_s):
    x=r/r_s
    return 4.0*np.pi*rho_s*(r_s**3)*(np.log(1.0+x)-x/(1.0+x))
def accel_xy(x,y,rho_s,r_s,MBH,eps=1e-6):
    r=np.sqrt(x*x+y*y)+eps
    a_bh=-G*MBH/(r*r)
    Menc=nfw_Menc(r,rho_s,r_s)
    a_h=-G*Menc/(r*r)
    a=a_bh+a_h
    return a*(x/r),a*(y/r)
def integrate_orbit(x0,y0,vx0,vy0,rho_s,r_s,MBH,dt=5e-4,steps=9000,r_stop=1e-4):
    s=np.array([x0,y0,vx0,vy0],dtype=float)
    traj=np.empty((steps,4),dtype=float)
    def rhs(state):
        x,y,vx,vy=state
        ax,ay=accel_xy(x,y,rho_s,r_s,MBH)
        return np.array([vx,vy,ax,ay],dtype=float)
    for i in range(steps):
        traj[i]=s
        r=np.hypot(s[0],s[1])
        if r<r_stop:
            traj=traj[:i+1]
            break
        k1=rhs(s)
        k2=rhs(s+0.5*dt*k1)
        k3=rhs(s+0.5*dt*k2)
        k4=rhs(s+dt*k3)
        s=s+(dt/6.0)*(k1+2*k2+2*k3+k4)
    return traj
def main():
    galaxy,rho_s,r_s,MBH=load_best_theta()
    print("Galaxy:",galaxy)
    print("Using theta:",rho_s,r_s,MBH)
    np.random.seed(0)
    N_orbits=800
    r_min=0.015
    r_max=0.060
    inits=[]
    for i in range(N_orbits):
        frac=i/(N_orbits-1)
        r=r_min+(r_max-r_min)*frac
        Menc=4.0*np.pi*rho_s*(r_s**3)*(np.log(1.0+r/r_s)-(r/r_s)/(1.0+r/r_s))
        v_circ=np.sqrt(G*(Menc+MBH)/r)
        v=v_circ*(0.8+0.4*np.random.rand())
        inits.append((r,0.0,0.0,v))
    dt=5e-4
    steps=50000
    trajs=[]
    maxlen=0
    for (x0,y0,vx0,vy0) in inits:
        tr=integrate_orbit(x0,y0,vx0,vy0,rho_s,r_s,MBH,dt=dt,steps=steps)
        trajs.append(tr)
        maxlen=max(maxlen,len(tr))
    for i in range(len(trajs)):
        last=trajs[i][-1]
        if len(trajs[i])<maxlen:
            pad=np.tile(last,(maxlen-len(trajs[i]),1))
            trajs[i]=np.vstack([trajs[i],pad])
    all_xy=np.vstack([tr[:,0:2] for tr in trajs])
    rmax=np.max(np.sqrt(all_xy[:,0]**2+all_xy[:,1]**2))
    lim=1.10*rmax
    v_ref=6.0
    r_infl=(G*MBH)/(v_ref*v_ref)
    fig=plt.figure(figsize=(8,8),facecolor="black")
    ax=plt.gca()
    ax.set_facecolor("black")
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.plot(0,0,"o",markersize=5,color="white")
    infl=plt.Circle((0,0),r_infl,color="red",fill=False,linestyle="--",linewidth=1.5)
    ax.add_patch(infl)
    lines=[]
    for i in range(N_orbits):
        ln,=ax.plot([],[],color="white",alpha=0.20,linewidth=0.5)
        lines.append(ln)
    out_mp4=Path(__file__).resolve().parents[2]/"outputs"/f"{galaxy}_orbits_nvenc.mp4"
    cmd=["ffmpeg","-y","-f","rawvideo","-vcodec","rawvideo","-pix_fmt","rgba","-s","800x800","-r","60","-i","-","-an","-vcodec","h264_nvenc","-preset","p5","-cq","19","-pix_fmt","yuv420p",str(out_mp4)]
    print("Launching NVENC encode...")
    pipe=subprocess.Popen(cmd,stdin=subprocess.PIPE)
    canvas=fig.canvas
    frame_skip=3
    for frame in range(0,maxlen,frame_skip):
        for i,ln in enumerate(lines):
            x=trajs[i][:frame,0]
            y=trajs[i][:frame,1]
            ln.set_data(x,y)
        canvas.draw()
        buf=np.frombuffer(canvas.buffer_rgba(),dtype=np.uint8)
        pipe.stdin.write(buf.tobytes())
    pipe.stdin.close()
    pipe.wait()
    print("Wrote:",out_mp4)
if __name__=="__main__":
    main()
