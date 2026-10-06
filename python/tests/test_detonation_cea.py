"""Compare Goddard's detonation solver against NASA CEA (the `cea` Python module).

Chapman-Jouguet cases:

* RP-1311 example 6: stoichiometric H2/O2 at T1 = 298.15 and 500 K, P1 = 1 and 20 bar.
* Lean (phi = 0.5) and rich (phi = 2) H2/O2, and stoichiometric H2/air, at 298.15 K and 1 bar.
* The RP-1311 example 7 shock-tube mixture, 5% H2, 5% O2, 90% Ar at 300 K and 10 mmHg.

Overdriven cases: above the CJ speed (~1094 m/s) the equilibrium shocks of RP-1311 example 7
(1100-1400 m/s) are overdriven detonations, so `cea.ShockSolver` is also the reference for
`DetonationSolver.detonation_from_velocity(u1, OVERDRIVEN)`.

Both codes are given the same product species: every gas species of `nasa9_gas.yaml` whose
elements are among those of the reactants.

Facts about `cea.DetonationSolver` (cea 3.3.4), checked directly before writing these tests:

* `DetonationSolver.solve(solution, weights, T1, p1_bar)` solves the CJ detonation only;
  `frozen=True` does not converge.
* `Mach` is the CJ velocity over the frozen sound speed of the unburned gas (`sonic_velocity1`,
  with `gamma1` the frozen cp/cv), as `DetonationResult.mach_in`.
* `M_M1` is the ratio of the gas molar masses 1/n, as `DetonationResult.molecular_weight_ratio`.

Tests are skipped automatically if the `cea` package is not installed.
"""
import os
from functools import lru_cache

import numpy as np
import pytest

cea = pytest.importorskip("cea")

import goddard
from goddard import DetonationBranch, GasChemistry
from conftest import find_data_dir, assert_close_rel
from test_fac_cea import run_isolated
from test_shock_cea import (
    MOLE_FRACTIONS as SHOCK_TUBE_MOLE_FRACTIONS,
    PRESSURE as SHOCK_TUBE_PRESSURE,
    RATIO_RTOL as SHOCK_RATIO_RTOL,
    SHOCK_SPEEDS,
    TEMPERATURE as SHOCK_TUBE_TEMPERATURE,
    _cea_shock_tube,
    _product_species as _shock_tube_product_species,
)

DATA_DIR = find_data_dir()
NASA9_GAS = os.path.join(DATA_DIR, "nasa9_gas.yaml")

BAR = 1e5  # Pa

# Same database, and both codes converge the CJ conditions to ~1e-5 or better. Observed worst
# error over the CJ cases 3.0e-5 (CJ velocity, shock-tube mixture). The overdriven comparison
# reuses the shock-tube cases and their tolerance, SHOCK_RATIO_RTOL.
RATIO_RTOL = 1e-4

# (id, reactant mole fractions, T1 [K], P1 [bar])
CJ_CASES = [
    ("ex6-298K-1bar", {"H2": 2.0, "O2": 1.0}, 298.15, 1.0),
    ("ex6-298K-20bar", {"H2": 2.0, "O2": 1.0}, 298.15, 20.0),
    ("ex6-500K-1bar", {"H2": 2.0, "O2": 1.0}, 500.0, 1.0),
    ("ex6-500K-20bar", {"H2": 2.0, "O2": 1.0}, 500.0, 20.0),
    ("h2o2-lean", {"H2": 1.0, "O2": 1.0}, 298.15, 1.0),
    ("h2o2-rich", {"H2": 4.0, "O2": 1.0}, 298.15, 1.0),
    ("h2-air", {"H2": 2.0, "O2": 1.0, "N2": 3.76}, 298.15, 1.01325),
    ("shock-tube-argon", SHOCK_TUBE_MOLE_FRACTIONS, SHOCK_TUBE_TEMPERATURE, SHOCK_TUBE_PRESSURE / BAR),
]


@lru_cache(maxsize=1)
def _nasa9_species():
    """Gas species of nasa9_gas.yaml, by name."""
    import cantera as ct
    return {sp.name: sp for sp in ct.Species.list_from_file(NASA9_GAS)}


def _elements(mole_fractions):
    species = _nasa9_species()
    return frozenset(e for name in mole_fractions for e in species[name].composition)


def _product_species(elements):
    """Every gas species of nasa9_gas.yaml whose elements are among `elements`."""
    return tuple(sorted(name for name, sp in _nasa9_species().items() if set(sp.composition) <= elements))


def _cea_chapman_jouguet(product_names, mole_fractions, T1, P1_bar):
    """Solve a CJ detonation with `cea.DetonationSolver`, inside the forked child."""
    names = list(mole_fractions)
    reactants = cea.Mixture(names)
    products = cea.Mixture(list(product_names))
    solver = cea.DetonationSolver(products, reactants=reactants)
    solution = cea.DetonationSolution(solver)
    weights = reactants.moles_to_weights(np.array([mole_fractions[n] for n in names]))
    solver.solve(solution, weights, T1, P1_bar)
    return {
        "converged": bool(solution.converged),
        "velocity": float(solution.velocity),
        "Mach": float(solution.Mach),
        "P_P1": float(solution.P_P1),
        "T_T1": float(solution.T_T1),
        "rho_rho1": float(solution.rho_rho1),
        "M_M1": float(solution.M_M1),
        "gamma_s": float(solution.gamma_s),
        "sonic_velocity": float(solution.sonic_velocity),
    }


def _goddard_solver(product_names, mole_fractions, T1, P1):
    gas = goddard.Gas(NASA9_GAS, "gas", GasChemistry.FROZEN, set(product_names))
    gas.set_state_TPX(T1, P1, mole_fractions)
    return goddard.DetonationSolver(gas)


@pytest.mark.parametrize("mole_fractions,T1,P1_bar", [c[1:] for c in CJ_CASES],
                         ids=[c[0] for c in CJ_CASES])
def test_chapman_jouguet_matches_cea(mole_fractions, T1, P1_bar):
    products = _product_species(_elements(mole_fractions))
    expected = run_isolated(_cea_chapman_jouguet, products, mole_fractions, T1, P1_bar)
    assert expected["converged"], "CEA did not converge"

    solver = _goddard_solver(products, mole_fractions, T1, P1_bar * BAR)
    cj = solver.chapman_jouguet()
    assert cj.valid
    assert_close_rel(cj.velocity, expected["velocity"], RATIO_RTOL, "CJ velocity")
    assert_close_rel(cj.mach_in, expected["Mach"], RATIO_RTOL, "CJ Mach")
    assert_close_rel(cj.static_pressure_ratio, expected["P_P1"], RATIO_RTOL, "P/P1")
    assert_close_rel(cj.static_temperature_ratio, expected["T_T1"], RATIO_RTOL, "T/T1")
    assert_close_rel(cj.density_ratio, expected["rho_rho1"], RATIO_RTOL, "rho/rho1")
    assert_close_rel(cj.molecular_weight_ratio, expected["M_M1"], RATIO_RTOL, "M/M1")

    products_state = solver.post_detonation_state()
    assert_close_rel(products_state.gamma_s, expected["gamma_s"], RATIO_RTOL, "gamma_s")
    assert_close_rel(products_state.speed_of_sound, expected["sonic_velocity"], RATIO_RTOL,
                     "sound speed")


def test_overdriven_matches_cea_equilibrium_shock():
    products = _shock_tube_product_species()
    reference = run_isolated(_cea_shock_tube, products, False, False)
    solver = _goddard_solver(products, SHOCK_TUBE_MOLE_FRACTIONS, SHOCK_TUBE_TEMPERATURE,
                             SHOCK_TUBE_PRESSURE)
    assert solver.chapman_jouguet().velocity < min(SHOCK_SPEEDS)

    for u1, expected in zip(SHOCK_SPEEDS, reference):
        assert expected["converged"], f"CEA did not converge at u1={u1}"
        label = f"u1={u1:g}"
        r = solver.detonation_from_velocity(u1, DetonationBranch.OVERDRIVEN)
        assert r.valid, label
        assert_close_rel(r.mach_in, expected["Mach"][0], SHOCK_RATIO_RTOL, f"{label} Mach1")
        assert_close_rel(r.mach_out, expected["Mach"][1], SHOCK_RATIO_RTOL, f"{label} Mach2")
        assert_close_rel(r.static_pressure_ratio, expected["P21"], SHOCK_RATIO_RTOL, f"{label} P2/P1")
        assert_close_rel(r.static_temperature_ratio, expected["T21"], SHOCK_RATIO_RTOL, f"{label} T2/T1")
        assert_close_rel(1.0 / r.density_ratio, expected["rho12"], SHOCK_RATIO_RTOL, f"{label} rho1/rho2")
