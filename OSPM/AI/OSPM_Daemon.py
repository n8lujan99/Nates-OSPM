# OSPM_Daemon.py
# STAYS IN PYTHON FOREVER
# Parallelism lives in Julia, not here
import os, time, sys
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from collections import deque
torch.backends.cudnn.benchmark = False
try:
    from sklearn.preprocessing import StandardScaler
except Exception:
    StandardScaler = None
# ============================================================
# Small utilities
# ============================================================
def clamp(x, lo, hi): return max(lo, min(hi, x))
def random_theta(bounds): return [np.random.uniform(lo, hi) for lo, hi in bounds]
def min_dist(theta, arr):
    if len(arr)==0: return np.inf
    return np.linalg.norm(np.asarray(arr)-np.asarray(theta), axis=1).min()
class IdentityScaler:
    def fit(self, X): return self
    def transform(self, X): return X
# ============================================================
# Models
# ============================================================
class Model(nn.Module):
    def __init__(self, dim, width=64):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(dim, width), nn.ReLU(), nn.Linear(width, width), nn.ReLU(), nn.Linear(width, 1))
    def forward(self, x): return self.net(x)
class Agent(nn.Module):
    def __init__(self, dim, hidden=128):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(dim, hidden), nn.ReLU(), nn.Linear(hidden, hidden), nn.ReLU(), nn.Linear(hidden, dim), nn.Tanh())
    def forward(self, x): return self.net(x)
    @torch.no_grad()
    def act(self, x, noise):
        a = self.forward(x) + noise*torch.randn_like(x)
        return torch.clamp(a, -1.0, 1.0)
# ============================================================
# Deck
# ============================================================
class Deck:
    def __init__(self, config):
        self.config = config
        self.path   = config["CSV_PATH"]
        self.cols   = config["REQUIRE_COLUMNS"]
        self.params = config["PARAMETER_NAMES"]
        self.flush  = int(config.get("CSV_FLUSH_INTERVAL", 50))
        self._dirty = 0
        self._buf   = []
        self._pbuf  = []
        self._sbuf  = []
        self._load()
    def _load(self):
        dir_name = os.path.dirname(self.path)
        if dir_name:
            os.makedirs(dir_name, exist_ok=True)
        if os.path.exists(self.path):
            df = pd.read_csv(self.path)
        else:
            row = {k:np.nan for k in self.cols}
            for i,k in enumerate(self.params):
                row[k] = self.config["INITIAL_THETA"][i]
            row["status"]="todo"
            df = pd.DataFrame([row])
            df.to_csv(self.path,index=False)
        missing=[c for c in self.cols if c not in df.columns]
        if missing:
            raise KeyError(f"Deck missing required columns: {missing}")
        self.df = df[self.cols].copy()
        self._params_arr = self.df[self.params].values.astype(float)
        self._status_arr = self.df["status"].values.astype(str)
    def _flush_buf(self):
        if not self._buf:
            return
        new_df = pd.DataFrame(self._buf, columns=self.cols)
        self.df = pd.concat([self.df, new_df], ignore_index=True)
        self._params_arr = self.df[self.params].values.astype(float)
        self._status_arr = self.df["status"].values.astype(str)
        self._buf.clear()
        self._pbuf.clear()
        self._sbuf.clear()
    def save(self):
        self._flush_buf()
        self.df.to_csv(self.path,index=False)
        print(f"[Deck] saved {len(self.df)} rows → {self.path}", flush=True)
    def _all_params(self):
        if self._pbuf:
            return np.vstack([self._params_arr, np.array(self._pbuf)])
        return self._params_arr
    def _all_status(self):
        if self._sbuf:
            return np.concatenate([self._status_arr, np.array(self._sbuf)])
        return self._status_arr
    def is_forbidden(self, theta, ndp=12):
        A = np.round(self._all_params(), ndp)
        t = np.round(theta, ndp)
        m = (A==t).all(axis=1)
        if not m.any(): return False
        return (self._all_status()[m]=="forbidden").any()
    def nearest_distance(self, theta, tol):
        A = self._all_params()
        m=np.all(np.abs(A-theta)<tol, axis=1)
        if not m.any(): return np.inf
        return np.linalg.norm(A[m]-theta,axis=1).min()
    def add(self, theta, chi2, reward, pid, status, refine_passes=None):
        row_dict = {k: theta[i] for i, k in enumerate(self.params)}
        row_dict |= dict(chi2=chi2, reward=reward, status=status, proposal_id=pid, refine_passes=refine_passes)
        self._buf.append([row_dict.get(k) for k in self.cols])
        self._pbuf.append([theta[i] for i in range(len(self.params))])
        self._sbuf.append(status)
        self._dirty+=1
        if self._dirty>=self.flush:
            self._flush_buf()
            self.save(); self._dirty=0
# ============================================================
# Physics wrapper
# ============================================================
class Corpo:
    def __init__(self, engine):
        self.engine = engine
    def eval(self, theta):
        try:
            chi2=float(self.engine(theta))
            if not np.isfinite(chi2): return "numeric_fail",np.inf
            return "pass",chi2
        except FloatingPointError:
            return "numeric_fail",np.inf
        except RuntimeError as e:
            import traceback
            print("[Corpo] RuntimeError:", repr(e))
            traceback.print_exc()
            raise
        except Exception:
            raise
# ============================================================
# AI helpers
# ============================================================
class Fixer:
    def __init__(self,cfg):
        self.warmup=int(cfg.get("AI_START_AFTER",500))
        self.unlocked=False
    def unlock(self, deck, runner):
        if self.unlocked: return
        if deck.df.status.str.startswith("pass").sum()>=self.warmup:
            runner.enable_ai()
            self.unlocked=True
            print("[AI] unlocked")
    def reward(self,status,chi2):
        if status!="pass": return -1e6
        return -float(chi2)
class FlatDetector:
    def __init__(self,w,eps,p):
        self.w,self.eps,self.p=w,eps,p
        self.buf=deque(maxlen=w); self.cnt=0
    def push(self,x):
        if not np.isfinite(x): return
        self.buf.append(x)
        if len(self.buf)<self.w: self.cnt=0; return
        self.cnt=self.cnt+1 if np.std(self.buf)<self.eps and np.isfinite(x) else 0
    def flat(self): return self.cnt>=self.p
# ============================================================
# Runner
# ============================================================
class Runner:
    def __init__(self,cfg):
        self.cfg=cfg
        self.bounds=cfg["THETA_BOUNDS"]
        self.cols=cfg["PARAMETER_NAMES"]
        self.dim=len(self.cols)
        self.batch=int(cfg["BATCH_SIZE"])
        self.min_d=float(cfg["MIN_DISTANCE"])
        self.ai=False
        self.model=None; self.agent=None
        self.opt_m=None; self.opt_a=None
        self.noise0=float(cfg.get("AI_NOISE_INIT",0.3))
        self.noise1=float(cfg.get("AI_NOISE_MIN",0.02))
        self.tau=float(cfg.get("AI_NOISE_TAU",5000))
        self.step=0
        self.recent=deque(maxlen=5000)
        self.scaler=IdentityScaler() if StandardScaler is None else StandardScaler()
        self.scaled=False
        self.fill_mode = False
        self.fill_triggered = False
        self.explore_frac = float(cfg.get("EXPLORE_FRACTION", 0.0))
        
    def enable_ai(self):
        self.model=Model(self.dim)
        self.agent=Agent(self.dim)
        self.opt_m=torch.optim.Adam(self.model.parameters(),1e-3)
        self.opt_a=torch.optim.Adam(self.agent.parameters(),1e-3)
        self.ai=True
    def _noise(self):
        if not self.ai: return self.noise0
        return max(self.noise1, self.noise0*np.exp(-self.step/self.tau))
    def _base(self,deck):
        good=deck.df[deck.df.status.str.startswith("pass")]
        if len(good)>=10:
            if np.random.rand() < 0.15:
                return good[self.cols].sample(1).values[0]
            return good.nsmallest(min(len(good),500),"chi2")[self.cols].sample(1).values[0]
        return deck.df[self.cols].dropna().sample(1).values[0]
    def detect_basin(self, deck):
        good = deck.df[deck.df.status.str.startswith("pass")]
        if len(good) < 500:
            return False
        top = good.nsmallest(min(len(good), 200), "chi2")
        chi_std = top["chi2"].std()
        X = top[self.cols].values
        spread = np.std(X, axis=0)
        span = np.array([hi - lo for lo, hi in self.bounds])
        rel_spread = np.mean(spread / span)
        return (chi_std < 1.0) and (rel_spread < 0.15) and ("refine_passes" in good.columns and good["refine_passes"].fillna(0).median() >= self.cfg.get("MAX_REFINE", 0))
    
    def step_scale(self, deck):
        if not self.ai:
            return 0.2
        if not self.fill_mode:
            return 0.2
        good = deck.df[deck.df.status.str.startswith("pass")]
        top = good.nsmallest(min(len(good), 200), "chi2")
        X = top[self.cols].values
        spread = np.std(X, axis=0)
        span = np.array([hi - lo for lo, hi in self.bounds])
        rel = np.mean(spread / span)
        return clamp(0.01 + 0.2 * rel, 0.01, 0.05)


    def propose(self,deck):
        out=[]
        while len(out)<self.batch:
            if self.ai and not (self.explore_frac > 0 and np.random.rand() < self.explore_frac):
                if self.fill_mode:
                    good = deck.df[deck.df.status.str.startswith("pass")]
                    base = good.nsmallest(100,"chi2")[self.cols].sample(1).values[0]
                else:
                    base = self._base(deck)
                xb=self.scaler.transform(base.reshape(1,-1)) if self.scaled else base.reshape(1,-1)
                a=self.agent.act(torch.tensor(xb,dtype=torch.float32),self._noise()).numpy().squeeze()
                s = self.step_scale(deck)
                if self.fill_mode and self.step % 200 == 0:
                    print(f"[FillMode] step_scale={s:.4f}")
                theta=[clamp(base[i]+s*(hi-lo)*a[i],lo,hi)
                    for i,(lo,hi) in enumerate(self.bounds)]
            else:
                theta=random_theta(self.bounds)
            if deck.is_forbidden(theta): continue
            if not self.fill_mode:
                if min_dist(theta,self.recent)<self.min_d: continue
                if deck.nearest_distance(theta,self.min_d)<self.min_d: continue
            self.recent.append(theta)
            self.step+=1
            out.append((theta,self.step))
        return out
    def train(self,deck):
        if not self.ai: return
        df=deck.df[deck.df.status.str.startswith("pass") & np.isfinite(deck.df.reward)]
        if len(df)<200: return
        if len(df)>5000: df=df.tail(5000)
        X=df[self.cols].values
        y=df.reward.values.reshape(-1,1)
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
    from collections import defaultdict

    deck   = Deck(config)
    runner = Runner(config)
    corpo  = Corpo(physics_engine)
    fixer  = Fixer(config)
    flat   = FlatDetector(config.get("FLAT_WINDOW", 200), config.get("FLAT_THRESHOLD", 1e-6), config.get("FLAT_PATIENCE", 3))

    runs = 0
    best = np.inf

    t_acc  = defaultdict(float)
    t_cnt  = defaultdict(int)
    PROF_EVERY = int(config.get("PROF_EVERY", 25))

    obs       = getattr(physics_engine, "__wrapped_obs__", None)
    halo_type = getattr(physics_engine, "__halo_type__", config.get("HALO_TYPE", "nfw"))
    use_batch = obs is not None

    if use_batch:
        from juliacall import Main
        import juliacall
        jl_batch = Main.OSPMPhysicsSpherical.evaluate_batch_theta
        sini     = float(obs.sini)
        Norbit   = int(obs.Norbit)
        print(f"[Daemon] batch mode ON — Norbit={Norbit}, Nstar_vlos={obs.Nstar_vlos}")
    else:
        print("[Daemon] batch mode OFF — falling back to serial corpo.eval")

    while runs < config["MAX_RUNS"]:

        print(f"[Daemon] loop iter runs={runs}", flush=True)
        t0 = time.perf_counter()

        deck._flush_buf()
        base_props = runner.propose(deck)

        props = []
        for theta, pid in base_props:
            rho_s, r_s, MBH = theta
            variants = [
                ("full", theta),

                # clean isolation
                ("bh_only",   [0.0,   r_s, MBH]),
                ("halo_only", [rho_s, r_s, 0.0]),

                # true local perturbations (keep the other component!)
                ("bh_up",     [rho_s, r_s, MBH * 2.0]),
                ("bh_down",   [rho_s, r_s, MBH * 0.5]),

                ("halo_up",   [rho_s * 2.0, r_s, MBH]),
                ("halo_down", [rho_s * 0.5, r_s, MBH]),
            ]
            for label, tvar in variants:
                props.append((tvar, pid, label))

        t_acc["propose"] += time.perf_counter() - t0
        t_cnt["propose"] += 1

        print(f"[Daemon] proposing {len(props)} variants, starting eval...", flush=True)

        t0 = time.perf_counter()
        CHUNK = 80

        def _record(theta, pid, label, status, chi2, refine_passes):
            nonlocal best, runs
            base_reward = fixer.reward(status, chi2)

            if status == "pass":
                stability = 1.0 - (refine_passes / max(1, config.get("MAX_REFINE", 1)))
                reward = base_reward + 0.5 * stability
            else:
                reward = base_reward

            t_add = time.perf_counter()
            deck.add(theta, chi2, reward, pid, f"{status}_{label}", refine_passes=refine_passes)
            t_acc["add"] += time.perf_counter() - t_add
            t_cnt["add"] += 1

            flat.push(chi2 if refine_passes >= config.get("MAX_REFINE", 0) else np.inf)
            fixer.unlock(deck, runner)

            if chi2 < best:
                best = chi2

            runs += 1

            if not runner.fill_triggered:
                if runner.detect_basin(deck):
                    runner.fill_mode = True
                    runner.fill_triggered = True
                    print(f"[Daemon] Basin detected at run {runs} — switching to fill mode")

            if runs % PROF_EVERY == 0:
                def avg(k):
                    n = t_cnt[k]
                    return (t_acc[k] / n) if n else 0.0

                t_eval    = t_acc["eval"]
                per_batch = t_eval / max(t_cnt["propose"], 1)
                per_theta = t_eval / max(t_cnt["eval"], 1)

                print(f"[PROF] runs={runs} best={best:.4f} propose={avg('propose'):.4f}s eval/batch={per_batch:.4f}s eval/theta={per_theta:.4f}s add={avg('add'):.4f}s", flush=True)

                t_acc.clear()
                t_cnt.clear()

            if flat.flat():
                deck.save()
                print(f"[Daemon] Flat region detected after {runs} runs")
                return True

            return False

        if use_batch:
            thetas = [theta for theta, pid, label in props]
            stop = False

            for i in range(0, len(thetas), CHUNK):
                chunk_props  = props[i:i+CHUNK]
                chunk_thetas = thetas[i:i+CHUNK]
                theta_mat    = np.array(chunk_thetas, dtype=float).T

                chunk_t0 = time.perf_counter()

                try:
                    NBINS_OCC  = config["OBSERVABLES"]["NBINS_OCC"]
                    LAMBDA_OCC = config["OBSERVABLES"]["LAMBDA_OCC"]

                    status_code_vec, chi2_vec, refine_vec = jl_batch(
                        juliacall.convert(Main.Matrix[Main.Float64], theta_mat), juliacall.convert(Main.Vector[Main.Float64], obs.R_star_m), juliacall.convert(Main.Vector[Main.Bool], obs.valid_vlos),
                        juliacall.convert(Main.Vector[Main.Float64], obs.v_star_mps), juliacall.convert(Main.Vector[Main.Float64], obs.verr_star_mps),
                        sini, Norbit, halo_type, Nocc=NBINS_OCC, lambda_occ=LAMBDA_OCC, max_refine=config.get("MAX_REFINE", 0))

                    t_acc["eval"] += time.perf_counter() - chunk_t0
                    t_cnt["eval"] += len(chunk_thetas)

                    for j, (theta, pid, label) in enumerate(chunk_props):

                        chi2 = float(chi2_vec[j])
                        code = int(status_code_vec[j])
                        refine_passes = int(refine_vec[j])   # ← FIX

                        if code == 0 and np.isfinite(chi2):
                            status = "pass"
                        elif code == 1:
                            status = "orbit_fail"
                        elif code == 2:
                            status = "numeric_fail"
                        else:
                            status = "unknown_fail"

                        if not np.isfinite(chi2):
                            chi2 = np.inf

                        if _record(theta, pid, label, status, chi2, refine_passes):  # ← FIX
                            stop = True
                            break

                except Exception as e:
                    t_acc["eval"] += time.perf_counter() - chunk_t0
                    t_cnt["eval"] += len(chunk_thetas)
                    print(f"[Daemon] chunk failed (size={len(chunk_thetas)}): {e}", flush=True)

                    for theta, pid, label in chunk_props:
                        if _record(theta, pid, label, "numeric_fail", np.inf, 0):  # ← FIX
                            stop = True
                            break

                if stop:
                    break

        else:
            for theta, pid, label in props:
                chunk_t0 = time.perf_counter()
                status, chi2 = corpo.eval(theta)

                t_acc["eval"] += time.perf_counter() - chunk_t0
                t_cnt["eval"] += 1

                if _record(theta, pid, label, status, chi2, 0):  # ← FIX
                    stop = True
                    break

        if stop:
            return

        t0 = time.perf_counter()
        deck._flush_buf()
        runner.train(deck)
        t_acc["train"] += time.perf_counter() - t0
        t_cnt["train"] += 1

    deck.save()