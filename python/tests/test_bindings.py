"""Behaviour of the Python bindings themselves: errors, read-only arrays, reprs, Gas methods."""
import os

import numpy as np
import pytest

import goddard
from goddard import Gas, GasChemistry
from conftest import find_data_dir


@pytest.fixture
def h2o2_yaml():
    return os.path.join(find_data_dir(), "h2o2.yaml")


@pytest.fixture
def burnt(h2o2_yaml):
    gas = Gas(h2o2_yaml, phase_name="ohmech", chemistry=GasChemistry.EQUILIBRIUM)
    gas.set_state_TPX(3000.0, 2e6, "H2:2, O2:1")
    gas.equilibrate("HP")
    return gas


def test_convergence_error_carries_solver_details(h2o2_yaml):
    """An equilibrium shock below the CJ speed of an exothermic mixture does not converge."""
    gas = Gas(h2o2_yaml, phase_name="ohmech", chemistry=GasChemistry.EQUILIBRIUM)
    gas.set_state_TPX(300.0, 1e5, "H2:2, O2:1")
    solver = goddard.ShockSolver(gas)
    with pytest.raises(goddard.ConvergenceError) as info:
        solver.normal_shock(1.5)
    error = info.value
    assert isinstance(error, RuntimeError)
    assert error.iterations > 0
    assert error.tolerance > 0.0
    assert np.isfinite(error.residual)


def test_profile_arrays_are_read_only():
    profile = goddard.conical_nozzle(4.0)
    with pytest.raises(ValueError):
        profile.x[0] = 1.0
    x = np.array(profile.x)
    x[0] += 0.0
    profile.x = x  # whole-array assignment still works
    np.testing.assert_array_equal(profile.x, x)


def test_gas_state_setters_round_trip(burnt):
    T, rho = burnt.temperature, burnt.density
    copy = burnt.clone()
    copy.set_state_TP(1000.0, 1e5)
    copy.set_state_TD(T, rho)
    assert copy.pressure == pytest.approx(burnt.pressure, rel=1e-10)

    u = burnt.enthalpy_mass - burnt.pressure / burnt.density
    copy.set_state_TP(1000.0, 1e5)
    copy.set_state_UV(u, 1.0 / rho)
    assert copy.temperature == pytest.approx(T, rel=1e-8)


def test_gas_flow_helpers(burnt):
    velocity = 1000.0
    expected = burnt.temperature * 8314.46261815324 / (
        burnt.pressure * velocity * burnt.mean_molecular_weight)
    assert burnt.area_per_mdot(velocity) == pytest.approx(expected, rel=1e-6)
    assert 1000.0 < burnt.cstar() < 3000.0


def test_gas_mixture_ratios(h2o2_yaml):
    gas = Gas(h2o2_yaml, phase_name="ohmech")
    gas.set_state_TP(300.0, 1e5)
    gas.set_OF_ratio(8.0, "H2:1", "O2:1", basis="mass")
    assert gas.stoich_OF_ratio("H2:1", "O2:1", basis="mass") == pytest.approx(7.94, rel=1e-2)
    assert gas.equivalence_ratio() == pytest.approx(7.94 / 8.0, rel=1e-2)
    assert gas.fuel_fraction("H2:1", "O2:1", basis="mass") == pytest.approx(1.0 / 9.0, rel=1e-8)
    gas.set_equivalence_ratio(1.0, "H2:1", "O2:1")
    assert gas.equivalence_ratio("H2:1", "O2:1") == pytest.approx(1.0, rel=1e-8)
    with pytest.raises(ValueError):
        gas.fuel_fraction("H2:1", "O2:1", basis="volume")
    assert "temperature" in gas.report()


def test_reprs(burnt):
    assert repr(burnt).startswith("<Gas ")
    results = goddard.Nozzle(burnt).solve(goddard.ExpansionType.SUPERSONIC_AREA_RATIO, 5.0)
    assert repr(results).startswith("<NozzleResults ")
    assert repr(results.throat).startswith("<ThroatCondition ")
    assert repr(results.expansions[0]).startswith("<NozzleStation A/At=")
    assert repr(results.expansions[0].thermo).startswith("<ThermodynamicState T=")
    shock = goddard.ShockSolver(burnt.clone()).normal_shock(2.0)
    assert repr(shock).startswith("<ShockResult valid=True")


def test_thermo_array_shape_properties(h2o2_yaml):
    gas = Gas(h2o2_yaml, phase_name="ohmech", chemistry=GasChemistry.EQUILIBRIUM)
    combustor = goddard.Combustor(gas, {"H2": 1.0}, {"O2": 1.0})
    states = combustor.solve(300.0, 300.0, np.array([1e6, 2e6]), np.array([6.0]),
                                       goddard.infinite_area_combustor([1e6, 2e6]))
    assert states.size == 2
    assert states.ndim == 2
    assert list(states.shape) == [1, 2]
