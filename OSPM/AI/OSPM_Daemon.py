# OSPM_Runner.py
# STAYS IN PYTHON FOREVER
# Parallelism lives in Julia, not here

import os, time, shutil, datetime
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from collections import deque

torch.set_num_threads(1)
torch.backends.cudnn.benchmark = False

try:
    from sklearn.preprocessing import StandardScaler
except Exception:
    StandardScaler = None

# ============================================================
# Small utilities
# ============================================================

def clamp(x, lo, hi): return max(lo, min(hi, x))
def clamp_vector(theta, bounds):
    return [clamp(theta[i], bounds[i][0], bounds[i][1]) for i in range(len(theta))]
def random_theta(bounds):
    return [np.random.uniform(lo, hi) for lo, hi in bounds]
def min_dist(theta, arr):
    if len(arr) == 0: return np.inf
    return np.linalg.norm(np.asarray(arr) - np.asarray(theta), axis=1).min()

class IdentityScaler:
    def fit(self, X): return self
    def transform(self, X): return X
    def fit_transform(self, X): return X

# ============================================================
# Models
# ============================================================

class Model(nn.Module):
    def __init__(self, dim, width=64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(dim, width), nn.ReLU(),
            nn.Linear(width, width), nn.ReLU(),
            nn.Linear(width, 1)
        )
    def forward(self, x): return self.net(x)

class Agent(nn.Module):
    def __init__(self, dim, hidden=128):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(dim, hidden), nn.ReLU(),
            nn.Linear(hidden, hidden), nn.ReLU(),
            nn.Linear(hidden, dim), nn.Tanh()
        )
    def forward(self, x): return self.net(x)
    @torch.no_grad()
    def act(self, x, noise):
        mu = self.forward(x)
        a = mu + noise * torch.randn_like(mu)
        return torch.clamp(a, -1.0, 1.0)

# ============================================================
# Deck
# ============================================================

class Deck:
    def __init__(self, config):
        self.config  = config
        self.path    = config["CSV_PATH"]
        self.cols    = config["REQUIRE_COLUMNS"]
        self.params  = config["PARAMETER_NAMES"]
        self.flush_n = int(config.get("CSV_FLUSH_INTERVAL", 50))
        self.allowed = config.get("ALLOWED_STATUSES", ["todo","pass","orbit_fail","numeric_fail","unknown_fail"])
        self._dirty  = 0
        self._load()

    def _load(self):
        os.makedirs(os.path.dirname(self.path), exist_ok=True)
        if os.path.exists(self.path):
            self.df = pd.read_csv(self.path)
        else:
            row = {k: np.nan for k in self.cols}
            for i,k in enumerate(self.params):
                row[k] = self.config["INITIAL_THETA"][i]
            row["status"] = "todo"
            self.df = pd.DataFrame([row])
            self.save()
        self.df = self.df[self.cols]

    def save(self): self.df.to_csv(self.path, index=False)

    def is_forbidden(self, theta, ndp=12):
        A = np.round(self.df[self.params].values, ndp)
        t = np.round(theta, ndp)
        mask = (A == t).all(axis=1)
        return mask.any() and (self.df.loc[mask,"status"] == "forbidden").any()

    def nearest_distance(self, theta, tol):
        A = self.df[self.params].values.astype(float)
        mask = np.all(np.abs(A - theta) < tol, axis=1)
        if not mask.any(): return np.inf
        return np.linalg.norm(A[mask] - theta, axis=1).min()

    def add(self, theta, chi2, reward, pid, status):
        row = {k: theta[i] for i,k in enumerate(self.params)}
        row |= dict(chi2=chi2, reward=reward, status=status, proposal_id=pid)
        self.df.loc[len(self.df)] = row
        self._dirty += 1
        if self._dirty >= self.flush_n:
            self.save(); self._dirty = 0

# ============================================================
# Physics wrapper
# ============================================================

class Corpo:
    def __init__(self, engine): self.engine = engine
    def eval(self, theta):
        try:
            chi2 = self.engine(theta)
            if not np.isfinite(chi2): return "numeric_fail", np.inf
            return "pass", float(chi2)
        except FloatingPointError: return "numeric_fail", np.inf
        except RuntimeError:       return "orbit_fail",  np.inf
        except Exception:          return "unknown_fail",np.inf

# ============================================================
# AI helpers
# ============================================================

class Fixer:
    def __init__(self, cfg):
        self.warmup = int(cfg.get("AI_START_AFTER",500))
        self.unlocked = False
    def unlock(self, deck, runner):
        if self.unlocked: return
        if (deck.df["status"]=="pass").sum() >= self.warmup:
            runner.enable_ai()
            self.unlocked = True
            print("[AI] unlocked")
    def reward(self, status, chi2):
        return -chi2 if status=="pass" else -1e9

class FlatDetector:
    def __init__(self, w, eps, p):
        self.w, self.eps, self.p = w, eps, p
        self.buf, self.cnt = [], 0
    def push(self, x):
        if not np.isfinite(x): return
        self.buf.append(x)
        if len(self.buf)>self.w: self.buf.pop(0)
        if len(self.buf)<self.w: self.cnt=0; return
        self.cnt = self.cnt+1 if np.std(self.buf)<self.eps else 0
    def flat(self): return self.cnt>=self.p

# ============================================================
# Runner
# ============================================================

class Runner:
    def __init__(self, cfg):
        self.cfg = cfg
        self.bounds = cfg["THETA_BOUNDS"]
        self.cols   = cfg["PARAMETER_NAMES"]
        self.dim    = len(self.cols)
        self.batch  = int(cfg["BATCH_SIZE"])
        self.min_d  = float(cfg["MIN_DISTANCE"])
        self.ai     = False
        self.model  = None
        self.agent  = None
        self.opt_m  = None
        self.opt_a  = None
        self.noise0 = float(cfg.get("AI_NOISE_INIT",0.3))
        self.noise1 = float(cfg.get("AI_NOISE_MIN",0.02))
        self.tau    = float(cfg.get("AI_NOISE_TAU",5000))
        self.step   = 0
        self.recent = deque(maxlen=5000)
        self.scaler = IdentityScaler() if StandardScaler is None else StandardScaler()
        self.scaled = False

    def enable_ai(self):
        self.model = Model(self.dim)
        self.agent = Agent(self.dim)
        self.opt_m = torch.optim.Adam(self.model.parameters(),1e-3)
        self.opt_a = torch.optim.Adam(self.agent.parameters(),1e-3)
        self.ai = True

    def _noise(self):
        if not self.ai: return self.noise0
        return max(self.noise1, self.noise0*np.exp(-self.step/self.tau))

    def propose(self, deck):
        out=[]
        while len(out)<self.batch:
            if self.ai:
                base = deck.df[self.cols].dropna().sample(1).values[0]
                x = torch.tensor(self.scaler.transform(base.reshape(1,-1)),dtype=torch.float32)
                a = self.agent.act(x,self._noise()).squeeze().numpy()
                theta=[clamp(base[i]+0.2*(hi-lo)*a[i],lo,hi)
                       for i,(lo,hi) in enumerate(self.bounds)]
            else:
                theta=random_theta(self.bounds)
            if deck.is_forbidden(theta): continue
            if min_dist(theta,self.recent)<self.min_d: continue
            if deck.nearest_distance(theta,self.min_d)<self.min_d: continue
            self.recent.append(theta); self.step+=1
            out.append((theta,self.step))
        return out

    def train(self, deck):
        if not self.ai: return
        df=deck.df[(deck.df.status=="pass") & np.isfinite(deck.df.reward)]
        if len(df)<200: return
        X=df[self.cols].values; y=df.reward.values.reshape(-1,1)
        if not self.scaled:
            self.scaler.fit(X); self.scaled=True
        Xt=torch.tensor(self.scaler.transform(X),dtype=torch.float32)
        yt=torch.tensor(y,dtype=torch.float32)
        loss=((self.model(Xt)-yt)**2).mean()
        self.opt_m.zero_grad(); loss.backward(); self.opt_m.step()

# ============================================================
# Daemon
# ============================================================

def run_daemon(config, physics_engine):

    deck   = Deck(config)
    runner = Runner(config)
    corpo  = Corpo(physics_engine)
    fixer  = Fixer(config)
    flat   = FlatDetector(
        config.get("FLAT_WINDOW",200),
        config.get("FLAT_THRESHOLD",1e-6),
        config.get("FLAT_PATIENCE",3)
    )

    start=time.time()
    runs,best=0,np.inf

    while runs<config["MAX_RUNS"]:
        props = runner.propose(deck)
        for theta,pid in props:
            status,chi2 = corpo.eval(theta)
            rwd = fixer.reward(status,chi2)
            deck.add(theta,chi2,rwd,pid,status)
            flat.push(chi2)
            fixer.unlock(deck,runner)
            if chi2<best: best=chi2
            runs+=1
            if flat.flat(): 
                deck.save()
                print(f"[Daemon] Terminating early due to flat region detected after {runs} runs.")
                return
        runner.train(deck)

    deck.save()
