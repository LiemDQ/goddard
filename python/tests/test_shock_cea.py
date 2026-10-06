"""Compare Goddard's shock solver against NASA CEA (the `cea` Python module).

The reference case is RP-1311 example 7: a shock tube filled with 0.05 H2 / 0.05 O2 / 0.9 Ar
(mole fractions) at 10 mmHg and 300 K, with incident shock speeds u1 from 1100 to 1400 m/s and
the shock reflected from the closed end. The unburned mixture is state 1: behind an equilibrium
shock it burns.

Both codes are given the same product species, every gas species of `nasa9_gas.yaml` whose
elements are among H, O and Ar, so the comparison is same-database.

Facts about `cea.ShockSolver`, verified directly against pycea before writing these tests:

* `ShockSolution` holds stations [1, 2, 5]. `P21`, `T21`, `rho12`, `P52`, `T52`, `rho52` are
  the jump ratios; `rho12` is rho1/rho2 while `rho52` is rho5/rho2.
* `v2` is the lab-frame velocity u_p of gas 2 and `velocity[2]` is the reflected shock speed
  W_R, so `u5_p_v2` = u_p + W_R is the speed of gas 2 relative to the reflected shock.
* `Mach[0]` is u1 over the frozen sound speed of the unburned mixture.
* With `reflected_frozen=True`, `sonic_velocity[2]` (and so `Mach[2]`) uses the incident
  frozen cp at the reflected temperature, so reflected Mach numbers are compared only for
  equilibrium reflected shocks.

Tests are skipped automatically if the `cea` package is not installed.
"""
import os
from functools import lru_cache

import numpy as np
import pytest

cea = pytest.importorskip("cea")

import goddard
from goddard import GasChemistry
from conftest import find_data_dir, assert_close_rel
from test_fac_cea import run_isolated

DATA_DIR = find_data_dir()
NASA9_GAS = os.path.join(DATA_DIR, "nasa9_gas.yaml")

MMHG = 133.322387415  # Pa
PRESSURE = 10.0 * MMHG
TEMPERATURE = 300.0
MOLE_FRACTIONS = {"H2": 0.05, "O2": 0.05, "Ar": 0.9}
SHOCK_SPEEDS = (1100.0, 1200.0, 1250.0, 1300.0, 1350.0, 1400.0)

# Same database, and both codes converge the jump conditions to ~1e-5 or better. Observed worst
# error 2.4e-4 (P2/P1, equilibrium, u1 = 1100 m/s, just above the ~1094 m/s CJ speed where the
# equilibrium jump is most sensitive), 8e-5 at every other speed.
RATIO_RTOL = 5e-4

# (incident, reflected) chemistry
CHEMISTRY_CASES = [
    (GasChemistry.EQUILIBRIUM, GasChemistry.EQUILIBRIUM),
    (GasChemistry.FROZEN, GasChemistry.FROZEN),
    (GasChemistry.FROZEN, GasChemistry.EQUILIBRIUM),
]


@lru_cache(maxsize=1)
def _product_species():
    """Every gas species of nasa9_gas.yaml whose elements are among H, O and Ar."""
    import cantera as ct
    species_list = ct.Species.list_from_file(NASA9_GAS)
    return sorted(sp.name for sp in species_list if set(sp.composition) <= {"H", "O", "Ar"})


def _cea_shock_tube(product_names, incident_frozen, reflected_frozen):
    """Solve RP-1311 example 7 with `cea.ShockSolver`, inside the forked child."""
    names = list(MOLE_FRACTIONS)
    reactants = cea.Mixture(names)
    products = cea.Mixture(product_names)
    solver = cea.ShockSolver(products, reactants=reactants)
    solution = cea.ShockSolution(solver, reflected=True)
    weights = reactants.moles_to_weights(np.array([MOLE_FRACTIONS[n] for n in names]))

    results = []
    for u1 in SHOCK_SPEEDS:
        solver.solve(solution, weights, TEMPERATURE, PRESSURE / 1e5, u1=u1, reflected=True,
                     incident_frozen=incident_frozen, reflected_frozen=reflected_frozen)
        results.append({
            "converged": bool(solution.converged),
            "Mach": [float(v) for v in solution.Mach],
            "P21": float(solution.P21),
            "T21": float(solution.T21),
            "rho12": float(solution.rho12),
            "P52": float(solution.P52),
            "T52": float(solution.T52),
            "rho52": float(solution.rho52),
            "u5_p_v2": float(solution.u5_p_v2),
        })
    return results


def _goddard_solver():
    gas = goddard.Gas(NASA9_GAS, "gas", GasChemistry.EQUILIBRIUM, set(_product_species()))
    gas.set_state_TPX(TEMPERATURE, PRESSURE, MOLE_FRACTIONS)
    return goddard.ShockSolver(gas)


@pytest.mark.parametrize("incident_chemistry,reflected_chemistry", CHEMISTRY_CASES,
                         ids=["equilibrium", "frozen", "frozen-equilibrium"])
def test_shock_tube_matches_cea(incident_chemistry, reflected_chemistry):
    reference = run_isolated(_cea_shock_tube, _product_species(),
                             incident_chemistry == GasChemistry.FROZEN,
                             reflected_chemistry == GasChemistry.FROZEN)
    solver = _goddard_solver()

    for u1, expected in zip(SHOCK_SPEEDS, reference):
        assert expected["converged"], f"CEA did not converge at u1={u1}"
        label = f"u1={u1:g}"
        result = solver.reflected_shock_from_velocity(u1, incident_chemistry, reflected_chemistry)
        assert result.valid, label
        incident, reflected = result.incident, result.reflected

        assert_close_rel(incident.mach_in, expected["Mach"][0], RATIO_RTOL, f"{label} Mach1")
        assert_close_rel(incident.mach_out, expected["Mach"][1], RATIO_RTOL, f"{label} Mach2")
        assert_close_rel(incident.static_pressure_ratio, expected["P21"], RATIO_RTOL, f"{label} P2/P1")
        assert_close_rel(incident.static_temperature_ratio, expected["T21"], RATIO_RTOL, f"{label} T2/T1")
        assert_close_rel(1.0 / incident.density_ratio, expected["rho12"], RATIO_RTOL, f"{label} rho1/rho2")

        assert_close_rel(reflected.static_pressure_ratio, expected["P52"], RATIO_RTOL, f"{label} P5/P2")
        assert_close_rel(reflected.static_temperature_ratio, expected["T52"], RATIO_RTOL, f"{label} T5/T2")
        assert_close_rel(reflected.density_ratio, expected["rho52"], RATIO_RTOL, f"{label} rho5/rho2")

        # Speed of gas 2 relative to the reflected shock, u_p + W_R = u_p rho52/(rho52 - 1)
        particle_velocity = u1 * (1.0 - 1.0 / incident.density_ratio)
        relative_speed = particle_velocity * reflected.density_ratio / (reflected.density_ratio - 1.0)
        assert_close_rel(relative_speed, expected["u5_p_v2"], RATIO_RTOL, f"{label} u5+v2")

        if reflected_chemistry == GasChemistry.EQUILIBRIUM:
            assert_close_rel(reflected.mach_out, expected["Mach"][2], RATIO_RTOL, f"{label} Mach5")
