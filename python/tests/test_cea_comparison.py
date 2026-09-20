"""Compare Goddard rocket simulation results against NASA CEA Python module.

Tests are skipped automatically if the `cea` package is not installed.
Run with: pytest python/tests/test_cea_comparison.py -v
"""
import os

import numpy as np
import pytest

cea = pytest.importorskip("cea")

import goddard
from conftest import (
    find_data_dir,
    RocketTestCase,
    ComparisonTolerances,
    build_goddard_problem,
    solve_cea_problem,
    compare_thermo_states,
    assert_close_rel,
    assert_close_abs,
    ReactantSpec
)


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

DEFAULT_CASES = ALL_TEST_CASES
EXPANSION_CASES = DEFAULT_CASES


# ---------------------------------------------------------------------------
# CEA station index helpers
# ---------------------------------------------------------------------------

def discover_cea_stations(cea_sol, area_ratios=None):
    """Discover CEA station mapping from solution arrays.

    CEA station layout for IAC with pi_p and supar:
      station 0 = chamber (Mach=0)
      station 1 = throat (Mach~1.0, ae_at=1.0)
      stations 2..N = exits (pressure ratio exits first, then area ratio exits)

    When area_ratios is provided, exit indices are matched by ae_at value.

    Returns a dict with keys 'chamber', 'throat', and 'exits' (list of indices).
    """
    mach = cea_sol.Mach
    ae_at = cea_sol.ae_at
    n = cea_sol.num_pts

    # Chamber is always first
    chamber_idx = 0

    # Find throat (Mach closest to 1.0)
    throat_idx = None
    for i in range(n):
        if 0.99 < mach[i] < 1.01:
            throat_idx = i
            break

    # Match exits by area ratio
    exit_indices = []
    if area_ratios is not None:
        for target_ar in area_ratios:
            best_idx = None
            best_err = float("inf")
            for i in range(n):
                if mach[i] > 1.0 and i != throat_idx:
                    err = abs(ae_at[i] - target_ar) / target_ar
                    if err < best_err:
                        best_err = err
                        best_idx = i
            if best_idx is not None and best_err < 0.01:
                exit_indices.append(best_idx)
    else:
        # Fallback: all supersonic stations after throat
        if throat_idx is not None:
            for i in range(throat_idx + 1, n):
                if mach[i] > 1.0:
                    exit_indices.append(i)

    return {
        "chamber": chamber_idx,
        "throat": throat_idx,
        "exits": exit_indices,
    }


# ---------------------------------------------------------------------------
# Thermodynamic state comparison tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("case", DEFAULT_CASES, ids=lambda c: c.name)
def test_chamber_state(case: RocketTestCase):
    """Compare chamber thermodynamic state between Goddard and CEA."""
    problem = build_goddard_problem(case)
    results = problem.solve()
    chamber = results.chamber(0, case.name).thermo

    cea_sol = solve_cea_problem(case)
    stations = discover_cea_stations(cea_sol, case.area_ratios)

    compare_thermo_states(chamber, cea_sol, stations["chamber"], label="chamber")


@pytest.mark.parametrize("case", EXPANSION_CASES, ids=lambda c: c.name)
def test_throat_state(case: RocketTestCase):
    """Compare throat thermodynamic state between Goddard and CEA."""
    problem = build_goddard_problem(case)
    results = problem.solve()
    throat = results.throat(0, case.name).thermo

    cea_sol = solve_cea_problem(case)
    stations = discover_cea_stations(cea_sol, case.area_ratios)
    assert stations["throat"] is not None, "CEA solution has no throat station"

    compare_thermo_states(throat, cea_sol, stations["throat"], label="throat")


@pytest.mark.parametrize("case", EXPANSION_CASES, ids=lambda c: c.name)
def test_exit_states(case: RocketTestCase):
    """Compare exit/expansion thermodynamic states between Goddard and CEA."""
    problem = build_goddard_problem(case)
    results = problem.solve()
    exits = [s.thermo for s in results.exits(0, case.name)]
    assert len(exits) == len(case.area_ratios), (
        f"Expected {len(case.area_ratios)} exit states, got {len(exits)}")

    cea_sol = solve_cea_problem(case)
    stations = discover_cea_stations(cea_sol, case.area_ratios)
    assert len(stations["exits"]) >= len(case.area_ratios), (
        f"CEA has {len(stations['exits'])} exits, expected {len(case.area_ratios)}")

    for i, exit_state in enumerate(exits):
        cea_idx = stations["exits"][i]
        compare_thermo_states(
            exit_state, cea_sol, cea_idx,
            label=f"exit[{i}] AR={case.area_ratios[i]}")


# ---------------------------------------------------------------------------
# Performance comparison tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("case", EXPANSION_CASES, ids=lambda c: c.name)
def test_performance(case: RocketTestCase):
    """Compare every field of RocketProblemResults.performance() against CEA.

    Both specific impulses are velocities [m/s] in both codes. CEA's
    `coefficient_of_thrust` is exactly `Isp / c_star` at every station (checked below), the
    same matched-nozzle definition Goddard uses, so CF compares directly with no pressure term
    to account for.
    """
    problem = build_goddard_problem(case)
    results = problem.solve()
    assert len(results.exits(0, case.name)) == len(case.area_ratios)

    cea_sol = solve_cea_problem(case)
    stations = discover_cea_stations(cea_sol, case.area_ratios)
    tol = ComparisonTolerances()

    chamber_pressure_bar = cea_sol.P[stations["chamber"]]

    for i, requested_area_ratio in enumerate(case.area_ratios):
        cea_idx = stations["exits"][i]
        label = f"exit[{i}] AR={requested_area_ratio}"
        performance = results.performance(0, i, case.name)

        # CEA's own CF definition, confirmed station by station.
        cea_cf = cea_sol.coefficient_of_thrust[cea_idx]
        assert cea_cf == pytest.approx(
            cea_sol.Isp[cea_idx] / cea_sol.c_star[cea_idx], rel=1e-6)

        assert_close_rel(performance.area_ratio, requested_area_ratio,
                         tol.area_ratio_rel, f"{label} area_ratio (requested)")
        assert_close_rel(performance.area_ratio, cea_sol.ae_at[cea_idx],
                         tol.area_ratio_rel, f"{label} area_ratio (CEA ae_at)")
        assert_close_rel(performance.pressure_ratio,
                         chamber_pressure_bar / cea_sol.P[cea_idx],
                         tol.pressure_ratio_rel, f"{label} pressure_ratio")
        assert_close_rel(performance.cstar, cea_sol.c_star[cea_idx],
                         tol.cstar_rel, f"{label} c_star")
        assert_close_rel(performance.isp, cea_sol.Isp[cea_idx],
                         tol.isp_rel, f"{label} Isp")
        assert_close_rel(performance.ivac, cea_sol.Isp_vacuum[cea_idx],
                         tol.ivac_rel, f"{label} Ivac")
        assert_close_rel(performance.mach_number, cea_sol.Mach[cea_idx],
                         tol.mach_rel, f"{label} Mach")
        assert_close_rel(performance.CF, cea_cf, tol.cf_rel, f"{label} CF")

        # The static entry point, fed the station states directly, must give the same numbers.
        direct = goddard.RocketProblemResults.calculate_performance(
            results.chamber(0, case.name).thermo,
            results.throat(0, case.name).thermo,
            results.exits(0, case.name)[i].thermo)
        for field in ("pressure_ratio", "area_ratio", "mach_number", "cstar", "CF",
                      "isp", "ivac"):
            assert getattr(direct, field) == pytest.approx(
                getattr(performance, field), rel=1e-12), f"{label} {field}"


# ---------------------------------------------------------------------------
# Species composition comparison tests
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("case", DEFAULT_CASES, ids=lambda c: c.name)
def test_chamber_composition(case: RocketTestCase):
    """Compare chamber species mass fractions between Goddard and CEA."""
    problem = build_goddard_problem(case)
    results = problem.solve()
    chamber = results.chamber(0, case.name).thermo

    cea_sol = solve_cea_problem(case)
    stations = discover_cea_stations(cea_sol, case.area_ratios)
    tol = ComparisonTolerances()

    goddard_comp = chamber.composition
    cea_mass_fracs = cea_sol.mass_fractions

    for species, goddard_frac in goddard_comp.items():
        if goddard_frac < 1e-6:
            continue
        # CEA may use different species names; skip if not found
        if species not in cea_mass_fracs:
            continue
        cea_frac = cea_mass_fracs[species][stations["chamber"]]
        assert_close_abs(goddard_frac, cea_frac, tol.mass_fraction_abs,
                         f"chamber mass_frac[{species}]")


# ---------------------------------------------------------------------------
# Isochoric (constant-volume) combustion
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("of_ratio", [4.0, 6.0, 8.0])
@pytest.mark.parametrize("initial_pressure_pa", [1e5, 1e6])
def test_isochoric_combustion_matches_cea(of_ratio, initial_pressure_pa):
    """Compare UV equilibrium of gaseous H2/O2 against cea.EqSolver."""
    from goddard import Combustor, CombustorOptions, CombustionProcess, Gas, GasChemistry

    T_reactant = 300.0
    yaml_path = os.path.join(find_data_dir(), "h2o2.yaml")

    gas = Gas(yaml_path, "ohmech", species=H2O2_SPECIES)
    combustor = Combustor(gas, {"H2": 1.0}, {"O2": 1.0})
    options = CombustorOptions(process=CombustionProcess.ISOCHORIC)
    burnt = combustor.solve(np.array([T_reactant]), np.array([initial_pressure_pa]),
                            np.array([of_ratio]), options)

    # Specific volume of the unburnt reactants sets the constant-volume constraint in CEA.
    reactants = Gas(yaml_path, "ohmech", species=H2O2_SPECIES)
    reactants.set_state_TPX(T_reactant, initial_pressure_pa,
                            f"H2:{1.0 / 2.01588}, O2:{of_ratio / 31.9988}")
    specific_volume = 1.0 / reactants.density

    reac = cea.Mixture(["H2", "O2"])
    prod = cea.Mixture(["H2", "O2"], products_from_reactants=True)
    solver = cea.EqSolver(prod, reactants=reac)
    solution = cea.EqSolution(solver)
    weights = reac.of_ratio_to_weights(np.array([0.0, 1.0]), np.array([1.0, 0.0]), of_ratio)
    u_reactants = reac.calc_property(cea.ENERGY, weights, np.array([T_reactant, T_reactant]))
    solver.solve(solution, cea.UV, u_reactants / cea.R, specific_volume, weights)
    assert solution.converged

    gas.restore_state(burnt.get_state(0))
    gas.chemistry = GasChemistry.EQUILIBRIUM
    # Tighter than ComparisonTolerances: gas-phase H2/O2 needs no cross-database allowance.
    assert_close_rel(gas.temperature, solution.T, 3e-3, "UV temperature")
    assert_close_rel(gas.pressure, solution.P * 1e5, 3e-3, "UV pressure")
    assert_close_rel(gas.mean_molecular_weight, solution.M, 1e-3, "UV molecular weight")
    assert_close_rel(gas.gamma_s, solution.gamma_s, 3e-3, "UV gamma_s")
