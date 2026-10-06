"""Solvers work on their own copy of the Gas passed to them, and state accessors return copies."""
import os

import numpy as np
import pytest

import goddard
from goddard import ExpansionType, Gas, GasChemistry
from conftest import find_data_dir


@pytest.fixture
def h2o2_yaml():
    return os.path.join(find_data_dir(), "h2o2.yaml")


def snapshot(gas):
    return gas.temperature, gas.pressure, np.array(gas.mass_fractions), gas.chemistry


def assert_unchanged(gas, before):
    T, P, Y, chemistry = before
    assert gas.temperature == T
    assert gas.pressure == P
    np.testing.assert_array_equal(np.array(gas.mass_fractions), Y)
    assert gas.chemistry == chemistry


def test_nozzle_leaves_gas_unchanged(h2o2_yaml):
    gas = Gas(h2o2_yaml, phase_name="ohmech", chemistry=GasChemistry.EQUILIBRIUM)
    gas.set_state_TPX(3000.0, 2e6, "H2:2, O2:1")
    gas.equilibrate("HP")
    before = snapshot(gas)

    goddard.Nozzle(gas).solve(ExpansionType.SUPERSONIC_AREA_RATIO, 10.0)
    assert_unchanged(gas, before)


def test_combustor_leaves_gas_unchanged(h2o2_yaml):
    gas = Gas(h2o2_yaml, phase_name="ohmech", chemistry=GasChemistry.EQUILIBRIUM)
    gas.set_state_TPX(300.0, 1e6, "H2:1")
    before = snapshot(gas)

    combustor = goddard.Combustor(gas, {"H2": 1.0}, {"O2": 1.0})
    combustor.solve_adiabatic(300.0, 300.0, np.array([1e6]), np.array([6.0]),
                              goddard.infinite_area_combustor([1e6]))
    assert_unchanged(gas, before)


def test_shock_solver_leaves_gas_and_earlier_states_unchanged(h2o2_yaml):
    gas = Gas(h2o2_yaml, phase_name="ohmech", chemistry=GasChemistry.FROZEN)
    gas.set_state_TPX(300.0, 1e5, "H2:2, O2:1, AR:7")
    before = snapshot(gas)

    solver = goddard.ShockSolver(gas)
    solver.normal_shock(3.0)
    assert_unchanged(gas, before)

    first = solver.post_shock_state()
    T_first = first.temperature
    pre = solver.pre_shock_state()
    solver.normal_shock(5.0)
    second = solver.post_shock_state()
    assert first.temperature == T_first
    assert second.temperature > T_first
    assert pre.temperature == pytest.approx(300.0)


def test_detonation_solver_leaves_gas_unchanged(h2o2_yaml):
    gas = Gas(h2o2_yaml, phase_name="ohmech", chemistry=GasChemistry.FROZEN)
    gas.set_state_TPX(300.0, 1e5, "H2:2, O2:1")
    before = snapshot(gas)

    solver = goddard.DetonationSolver(gas)
    products = solver.post_detonation_state()
    von_neumann = solver.von_neumann_state()
    assert_unchanged(gas, before)
    # Each accessor returns its own copy.
    assert products.temperature != von_neumann.temperature
