"""Tests for the Python ThermoArray: array properties, indexing and compositions."""
import os

import numpy as np
import pytest

import goddard
from conftest import find_data_dir

BAR = 1e5  # Pa


@pytest.fixture(scope="module")
def states():
    """H2/O2 combustion products on a (mixture ratio, pressure) grid of shape (3, 2)."""
    gas = goddard.Gas(os.path.join(find_data_dir(), "h2o2.yaml"), "ohmech")
    combustor = goddard.Combustor(gas, {"H2": 1.0}, {"O2": 1.0})
    return combustor.solve(300.0, 300.0, np.array([10 * BAR, 50 * BAR]),
                           np.array([4.0, 5.0, 6.0]), goddard.infinite_area_combustor())


SCALAR_PROPERTIES = ["temperature", "pressure", "density", "enthalpy_mass", "entropy_mass",
                     "cp_mass", "cv_mass", "mean_molecular_weight", "mixture_molecular_weight",
                     "gas_mass_fraction", "gamma_s", "equilibrium_gamma_s", "speed_of_sound",
                     "equilibrium_speed_of_sound"]


def test_array_properties_are_the_gas_properties_of_each_entry(states):
    """states.<name>[i, j] is states[i, j].<name> for every scalar Gas property."""
    assert states.shape == (3, 2)
    for name in SCALAR_PROPERTIES:
        values = getattr(states, name)
        assert values.shape == states.shape, name
        for i in range(3):
            for j in range(2):
                assert values[i, j] == pytest.approx(getattr(states[i, j], name), rel=1e-12), \
                    f"{name} at ({i}, {j})"


def test_grid_axes_follow_the_solve_arguments(states):
    """Axis 0 is the mixture ratio, axis 1 the pressure: pressure is constant along axis 0."""
    assert states.pressure[:, 0] == pytest.approx(np.full(3, 10 * BAR), rel=1e-9)
    assert states.pressure[:, 1] == pytest.approx(np.full(3, 50 * BAR), rel=1e-9)
    # More oxidizer leaves more O2 in the fuel-rich products.
    o2 = states.species_names.index("O2")
    assert np.all(np.diff(states.mole_fractions[:, 0, o2]) > 0.0)


def test_internal_energy_is_enthalpy_minus_flow_work(states):
    expected = states.enthalpy_mass - states.pressure / states.density
    assert states.internal_energy_mass == pytest.approx(expected, rel=1e-12)


def test_compositions_have_a_species_axis_and_sum_to_one(states):
    n_species = len(states.species_names)
    for name in ("mole_fractions", "mass_fractions", "mixture_mass_fractions"):
        fractions = getattr(states, name)
        assert fractions.shape == states.shape + (n_species,), name
        assert fractions.sum(axis=-1) == pytest.approx(np.ones(states.shape), rel=1e-12), name
    assert states.mass_fractions[2, 1] == pytest.approx(states[2, 1].mass_fractions, rel=1e-12)


def test_indexing_returns_an_independent_gas(states):
    gas = states[1, 0]
    temperature = states.temperature[1, 0]
    gas.set_state_TP(500.0, 1e5)
    assert states.temperature[1, 0] == temperature
    assert states[-1, -1].temperature == states.temperature[2, 1]
    with pytest.raises(IndexError):
        states[3, 0]
    with pytest.raises(IndexError):
        states[0]


def test_snapshot_and_evaluate(states):
    snapshot = states.snapshot((2, 1), goddard.CompositionModel.EQUILIBRIUM)
    assert snapshot.composition_model == goddard.CompositionModel.EQUILIBRIUM
    assert snapshot.gamma_s == pytest.approx(states.equilibrium_gamma_s[2, 1], rel=1e-12)
    cstar = states.evaluate(lambda gas: gas.cstar(goddard.CompositionModel.EQUILIBRIUM))
    assert cstar.shape == states.shape
    assert cstar[0, 1] == pytest.approx(
        states[0, 1].cstar(goddard.CompositionModel.EQUILIBRIUM), rel=1e-12)
