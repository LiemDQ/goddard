"""Shared test infrastructure for Goddard Python tests.

Provides dataclasses for test case specification, builder functions for both
Goddard and CEA problems, and comparison utilities with configurable tolerances.
"""
import os
from dataclasses import dataclass

import numpy as np
import pytest
import goddard

# ---------------------------------------------------------------------------
# Data directory
# ---------------------------------------------------------------------------

def find_data_dir() -> str:
    """Find the Goddard data/ directory by walking up from this test file."""
    here = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(here))
    data_dir = os.path.join(project_root, "data")
    if not os.path.isdir(data_dir):
        pytest.skip(f"Data directory not found at {data_dir}")
    return data_dir


# ---------------------------------------------------------------------------
# Test case specification
# ---------------------------------------------------------------------------

@dataclass
class ReactantSpec:
    """Specifies a single reactant (fuel or oxidizer) for both Goddard and CEA."""
    cea_name: str                # CEA species name, e.g. "H2" or "O2(L)"
    cantera_composition: str     # Cantera composition string, e.g. "H2:1"
    temperature: float           # Reactant temperature in K


@dataclass
class RocketTestCase:
    """Declarative specification of a rocket test case.

    Contains all information needed to set up equivalent problems in both
    Goddard (via Cantera) and CEA.
    """
    name: str
    fuel: ReactantSpec
    oxidizer: ReactantSpec
    of_ratio: float
    chamber_pressure_bar: float
    area_ratios: list[float]
    thermo_file: str             # e.g. "h2o2.yaml"
    phase_name: str              # e.g. "ohmech"
    species: set[str]
    nozzle_chemistry: str = "equilibrium"   # "equilibrium" or "frozen"
    frozen_NFZ: int = 1


@dataclass
class ComparisonTolerances:
    """Tolerances for comparing Goddard vs CEA results.

    Values from cea_comparer.hpp with slight relaxation for cross-code comparison.
    """
    temperature_rel: float = 1e-2       # 1%
    pressure_rel: float = 1e-2          # 1%
    density_rel: float = 1e-2           # 1%
    enthalpy_rel: float = 1e-2          # 1%
    entropy_rel: float = 1e-2           # 1%
    molecular_weight_rel: float = 1e-3  # 0.1%
    gamma_rel: float = 1e-2             # 1%
    sound_speed_rel: float = 1e-2       # 1%
    thermo_deriv_rel: float = 5e-2      # 5%
    mass_fraction_abs: float = 5e-3     # absolute (relaxed: different thermo databases)
    cstar_rel: float = 1e-2             # 1%
    isp_rel: float = 1e-2               # 1%
    ivac_rel: float = 1e-2              # 1%
    cf_rel: float = 2e-2                # 2%


# ---------------------------------------------------------------------------
# Goddard problem builder
# ---------------------------------------------------------------------------

def cantera_save_state(yaml_path: str, phase_name: str,
                       composition: str, temperature: float,
                       pressure: float = 101325.0) -> list[float]:
    """Create a Cantera state vector matching C++ ThermoPhase::saveState() format.

    The format is [T, density, Y_0, Y_1, ..., Y_n] as a flat list.
    """
    import cantera as ct
    sol = ct.Solution(yaml_path, phase_name)
    sol.TPX = temperature, pressure, composition
    return [sol.T, sol.density] + list(sol.Y)


def build_goddard_problem(case: RocketTestCase):
    """Build a Goddard RocketProblem from a RocketTestCase specification."""
    import goddard

    data_dir = find_data_dir()
    yaml_path = os.path.join(data_dir, case.thermo_file)

    chem_params = goddard.ChemicalParameters()
    chem_params.thermo_file = yaml_path
    chem_params.species = case.species
    chem_params.cantera_fuel_state = goddard.ThermodynamicState(
        case.fuel.temperature, 101325.0, case.fuel.cantera_composition) 
    chem_params.cantera_oxidizer_state = goddard.ThermodynamicState(
        case.oxidizer.temperature, 101325.0, case.oxidizer.cantera_composition)
    chem_params.OF_ratios = [case.of_ratio]

    case_params = goddard.RocketCaseParameters()
    case_params.name = case.name
    case_params.problem_type = "rocket"
    case_params.combustor_options = goddard.infinite_area_combustor(
        [case.chamber_pressure_bar * 1e5])  # bar -> Pa

    if case.nozzle_chemistry == "equilibrium":
        case_params.nozzle_options = goddard.equilibrium_nozzle(*case.area_ratios)
    else:
        case_params.nozzle_options = goddard.frozen_nozzle(
            *case.area_ratios, frozen_NFZ=case.frozen_NFZ)

    return goddard.RocketProblem(chem_params, [case_params], case.phase_name)


# ---------------------------------------------------------------------------
# CEA problem builder
# ---------------------------------------------------------------------------

def solve_cea_problem(case: RocketTestCase):
    """Solve a CEA rocket problem from a RocketTestCase specification.

    Returns the CEA RocketSolution object. Caller must have `cea` installed.
    """
    import cea

    reac_names = [case.fuel.cea_name, case.oxidizer.cea_name]
    T_reactant = np.array([case.fuel.temperature, case.oxidizer.temperature])
    fuel_weights = np.array([1.0, 0.0])
    oxidant_weights = np.array([0.0, 1.0])

    reac = cea.Mixture(reac_names)
    prod = cea.Mixture(reac_names, products_from_reactants=True)

    solver = cea.RocketSolver(prod, reactants=reac)
    solution = cea.RocketSolution(solver)

    weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, case.of_ratio)
    hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant) / cea.R

    solve_kwargs = dict(iac=True, hc=hc)
    if case.area_ratios:
        solve_kwargs["supar"] = case.area_ratios
    if case.nozzle_chemistry == "frozen":
        # CEA station numbering: 1=combustor, 2=throat
        # frozen_NFZ=1 in Goddard (freeze at throat) maps to n_frz=2 in CEA
        solve_kwargs["n_frz"] = case.frozen_NFZ + 1

    # CEA requires at least one pressure ratio; use a large one to get all stations
    pi_p = [1000.0]
    solver.solve(solution, weights, case.chamber_pressure_bar, pi_p, **solve_kwargs)
    return solution


# ---------------------------------------------------------------------------
# Comparison utilities
# ---------------------------------------------------------------------------

def assert_close_rel(actual: float, expected: float, rtol: float, name: str):
    """Assert two values are within relative tolerance."""
    if abs(expected) < 1e-12:
        assert abs(actual - expected) < rtol, (
            f"{name}: actual={actual}, expected={expected}, "
            f"abs_diff={abs(actual - expected)}, tol={rtol}")
    else:
        rel_err = abs(actual - expected) / abs(expected)
        assert rel_err <= rtol, (
            f"{name}: actual={actual:.6g}, expected={expected:.6g}, "
            f"rel_err={rel_err:.4e}, tol={rtol}")


def assert_close_abs(actual: float, expected: float, atol: float, name: str):
    """Assert two values are within absolute tolerance."""
    assert abs(actual - expected) <= atol, (
        f"{name}: actual={actual:.6g}, expected={expected:.6g}, "
        f"abs_diff={abs(actual - expected):.4e}, tol={atol}")


def compare_thermo_states(goddard_state, cea_solution, cea_station_idx: int,
                          tol: ComparisonTolerances = None, label: str = ""):
    """Compare a Goddard ThermoStateInfo against a CEA solution at a station.

    Checks temperature, pressure, molecular weight, gamma_s, and speed of sound.
    """
    if tol is None:
        tol = ComparisonTolerances()
    prefix = f"{label} " if label else ""

    # Temperature (both in K)
    assert_close_rel(goddard_state.temperature,
                     cea_solution.T[cea_station_idx],
                     tol.temperature_rel, f"{prefix}temperature")

    # Pressure (Goddard Pa, CEA bar)
    assert_close_rel(goddard_state.pressure / 1e5,
                     cea_solution.P[cea_station_idx],
                     tol.pressure_rel, f"{prefix}pressure")

    # Molecular weight (both g/mol)
    assert_close_rel(goddard_state.molecular_weight,
                     cea_solution.M[cea_station_idx],
                     tol.molecular_weight_rel, f"{prefix}molecular_weight")

    # Isentropic exponent
    assert_close_rel(goddard_state.gamma_s,
                     cea_solution.gamma_s[cea_station_idx],
                     tol.gamma_rel, f"{prefix}gamma_s")

    # Speed of sound (both m/s)
    assert_close_rel(goddard_state.speed_of_sound,
                     cea_solution.sonic_velocity[cea_station_idx],
                     tol.sound_speed_rel, f"{prefix}speed_of_sound")


# ---------------------------------------------------------------------------
# Predefined test cases
# ---------------------------------------------------------------------------

H2O2_SPECIES = {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"}

H2_O2_GAS_EQUILIBRIUM = RocketTestCase(
    name="h2_o2_gas_eq",
    fuel=ReactantSpec(cea_name="H2", cantera_composition="H2:1", temperature=300.0),
    oxidizer=ReactantSpec(cea_name="O2", cantera_composition="O2:1", temperature=300.0),
    of_ratio=6.0,
    chamber_pressure_bar=206.84,
    area_ratios=[15.0, 35.0],
    thermo_file="h2o2.yaml",
    phase_name="ohmech",
    species=H2O2_SPECIES,
    nozzle_chemistry="equilibrium",
)

H2_O2_GAS_FROZEN = RocketTestCase(
    name="h2_o2_gas_frz",
    fuel=ReactantSpec(cea_name="H2", cantera_composition="H2:1", temperature=300.0),
    oxidizer=ReactantSpec(cea_name="O2", cantera_composition="O2:1", temperature=300.0),
    of_ratio=6.0,
    chamber_pressure_bar=206.84,
    area_ratios=[15.0, 35.0],
    thermo_file="h2o2.yaml",
    phase_name="ohmech",
    species=H2O2_SPECIES,
    nozzle_chemistry="frozen",
    frozen_NFZ=1,
)

ALL_TEST_CASES = [H2_O2_GAS_EQUILIBRIUM, H2_O2_GAS_FROZEN]
EQUILIBRIUM_CASES = [H2_O2_GAS_EQUILIBRIUM]
FROZEN_CASES = [H2_O2_GAS_FROZEN]
