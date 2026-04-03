"""Compare Goddard rocket simulation results against NASA CEA Python module.

Tests are skipped automatically if the `cea` package is not installed.
Run with: pytest python/tests/test_cea_comparison.py -v
"""
import pytest

cea = pytest.importorskip("cea")

import goddard
from conftest import (
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

# Frozen cases excluded from default runs pending gamma_s fix.
# Use FROZEN_CASES when ready to enable them.
DEFAULT_CASES = ALL_TEST_CASES


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


@pytest.mark.parametrize("case", DEFAULT_CASES, ids=lambda c: c.name)
def test_throat_state(case: RocketTestCase):
    """Compare throat thermodynamic state between Goddard and CEA."""
    problem = build_goddard_problem(case)
    results = problem.solve()
    throat = results.throat(0, case.name).thermo

    cea_sol = solve_cea_problem(case)
    stations = discover_cea_stations(cea_sol, case.area_ratios)
    assert stations["throat"] is not None, "CEA solution has no throat station"

    compare_thermo_states(throat, cea_sol, stations["throat"], label="throat")


@pytest.mark.parametrize("case", DEFAULT_CASES, ids=lambda c: c.name)
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

@pytest.mark.parametrize("case", DEFAULT_CASES, ids=lambda c: c.name)
def test_performance(case: RocketTestCase):
    """Compare rocket performance metrics between Goddard and CEA.

    Note: calculate_performance() uses an ideal-gas c* approximation.
    We compute c* = Pc / (rho_t * a_t) directly for a fair comparison.
    Isp is computed from enthalpy difference and does not depend on c*.
    """
    import math

    problem = build_goddard_problem(case)
    results = problem.solve()

    chamber = results.chamber(0, case.name).thermo
    throat = results.throat(0, case.name).thermo
    exits = [s.thermo for s in results.exits(0, case.name)]
    assert len(exits) > 0

    cea_sol = solve_cea_problem(case)
    stations = discover_cea_stations(cea_sol, case.area_ratios)
    tol = ComparisonTolerances()

    # c* from first principles: c* = Pc / (rho_t * a_t)
    goddard_cstar = chamber.pressure / (throat.density * throat.speed_of_sound)
    g0 = 9.80665

    for i, exit_state in enumerate(exits):
        cea_idx = stations["exits"][i]
        label = f"exit[{i}] AR={case.area_ratios[i]}"

        # c* (both in m/s)
        assert_close_rel(goddard_cstar, cea_sol.c_star[cea_idx],
                         tol.cstar_rel, f"{label} c_star")

        # Isp from enthalpy difference: v_e = sqrt(2*(h_c - h_e)), Isp = v_e / g0
        exit_velocity = math.sqrt(2.0 * (chamber.enthalpy - exit_state.enthalpy))
        goddard_isp_m_s = exit_velocity
        assert_close_rel(goddard_isp_m_s, cea_sol.Isp[cea_idx],
                         tol.isp_rel, f"{label} Isp")

        # Ivac = v_e + Pe*Ae/(mdot) = v_e + Pe/(rho_t*a_t)*(Ae/At)
        # Using c* = Pc/(rho_t*a_t) and Ae/At from throat/exit density/velocity:
        area_ratio = (throat.density * throat.speed_of_sound) / \
                     (exit_state.density * exit_velocity)
        goddard_ivac_m_s = exit_velocity + \
            (exit_state.pressure * area_ratio) / \
            (throat.density * throat.speed_of_sound)
        assert_close_rel(goddard_ivac_m_s, cea_sol.Isp_vacuum[cea_idx],
                         tol.ivac_rel, f"{label} Ivac")


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
