import os
import sys
os.environ["JULIA_NUM_THREADS"] = "auto"  # Must be first, before ANY other imports

from OSPM.load_config import load_config
from .OSPM_Control import build_runtime
from .OSPM_MASTER import build_observables, solve_ospm_theta
from ..Physics.OSPM_PhysicsEngine import wrap_physics_engine

def build_physics_engine(config):
    obs = build_observables(config)
    print("torch imported?", "torch" in sys.modules)
    print("obs type:", type(obs))
    print("has R_star_m:", hasattr(obs, "R_star_m"))
    print("has v_star_mps:", hasattr(obs, "v_star_mps"))
    print("has verr_star_mps:", hasattr(obs, "verr_star_mps"))
    def base_engine(theta, *, return_A=False, **_ignored):
        chi2, w, model, diag_data = solve_ospm_theta(theta, obs, halo_type=config["HALO_TYPE"])
        if return_A:
            return model
        return float(chi2)
    return wrap_physics_engine(base_engine, obs=obs, halo_type=config["HALO_TYPE"], config=config)

def main():
    config = load_config()
    runtime = build_runtime(config)
    from ..Physics import OSPM_Physics as P
    P._jl_init()
    from juliacall import Main
    print("Julia threads seen by module:", Main.OSPMPhysicsSpherical.NTHREADS)
    physics_engine = build_physics_engine(config)
    print("torch imported after build?", "torch" in sys.modules)
    from .OSPM_API import OSPM_API
    api = OSPM_API(runtime)
    api.set_physics_engine(physics_engine)
    result = api.run()
    print(result)
if __name__ == "__main__":
    main()