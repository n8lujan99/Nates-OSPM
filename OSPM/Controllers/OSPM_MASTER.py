"""
OSPM_Master
Dispatch layer for OSPM.
Currently supports STELLAR mode only.
"""
from ..Observables.OSPM_Observables_Stellar import OSPMObservablesStellar
from ..Solvers.OSPM_Solver_Stellar import solve_ospm_theta_stellar
def build_observables(config):
    return OSPMObservablesStellar.from_star_table( config["DATA_CSV"], 
            inclination_deg=config["INCLINATION_DEG"], Norbit=config["NORBIT"])
def solve_ospm_theta(theta, obs, *, halo_type="nfw"):
    return solve_ospm_theta_stellar(theta, obs, halo_type=halo_type)
