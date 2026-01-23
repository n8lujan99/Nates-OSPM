from OSPM.load_config import load_config
from .OSPM_Control import build_runtime
from .OSPM_MASTER import build_observables, solve_ospm_theta
from .OSPM_API import OSPM_API
from ..Physics.OSPM_PhysicsEngine import wrap_physics_engine

import sys

def main():
    Galaxy = sys.argv[1]                 # ← single source of galaxy choice
    config = load_config(Galaxy)         # ← ONLY config entry point

    runtime = build_runtime(config)
    obs = build_observables(config)

    def base_engine(theta):
        chi2, _, _ = solve_ospm_theta(theta, obs, halo_type=config["HALO_TYPE"])
        return float(chi2)

    physics_engine = wrap_physics_engine(
        base_engine,
        obs=obs,
        halo_type=config["HALO_TYPE"],
        config=config
    )

    api = OSPM_API(runtime)
    api.set_physics_engine(physics_engine)

    result = api.run()
    print(result)

if __name__ == "__main__":
    main()


