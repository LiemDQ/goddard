"""Tests for the Gas Python binding."""
import os
import math

import pytest

from conftest import find_data_dir


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def h2o2_yaml():
    """Path to h2o2.yaml data file."""
    return os.path.join(find_data_dir(), "h2o2.yaml")


@pytest.fixture
def gas(h2o2_yaml):
    """A frozen Gas from h2o2.yaml at default state (300 K, 1 atm)."""
    from goddard import Gas, GasChemistry
    return Gas(h2o2_yaml, phase_name="ohmech", chemistry=GasChemistry.FROZEN)


@pytest.fixture
def eq_gas(h2o2_yaml):
    """An equilibrium Gas from h2o2.yaml."""
    from goddard import Gas, GasChemistry
    return Gas(h2o2_yaml, phase_name="ohmech", chemistry=GasChemistry.EQUILIBRIUM)


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

def test_yaml_construction(gas):
    assert gas.temperature > 0
    assert gas.pressure > 0


def test_yaml_construction_with_species(h2o2_yaml):
    from goddard import Gas, GasChemistry
    # Must include all species referenced by the phase definition
    species = {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"}
    g = Gas(h2o2_yaml, phase_name="ohmech", species=species,
            chemistry=GasChemistry.FROZEN)
    assert g.temperature > 0


def test_solution_handle_construction(h2o2_yaml):
    from goddard import Gas, GasChemistry, create_solution
    species = {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"}
    sol = create_solution(h2o2_yaml, "ohmech", species)
    g = Gas(sol, chemistry=GasChemistry.FROZEN)
    assert g.temperature > 0


# ---------------------------------------------------------------------------
# Property access (should be attributes, not method calls)
# ---------------------------------------------------------------------------

def test_properties_are_attributes(gas):
    # These should all be accessible as attributes, not methods
    assert isinstance(gas.temperature, float)
    assert isinstance(gas.pressure, float)
    assert isinstance(gas.density, float)
    assert isinstance(gas.enthalpy_mass, float)
    assert isinstance(gas.entropy_mass, float)
    assert isinstance(gas.cp_mass, float)
    assert isinstance(gas.cv_mass, float)
    assert isinstance(gas.mean_molecular_weight, float)
    assert isinstance(gas.gamma_s, float)
    assert isinstance(gas.speed_of_sound, float)


def test_properties_positive(gas):
    assert gas.temperature > 0
    assert gas.pressure > 0
    assert gas.density > 0
    assert gas.cp_mass > 0
    assert gas.cv_mass > 0
    assert gas.mean_molecular_weight > 0
    assert gas.gamma_s > 1.0
    assert gas.speed_of_sound > 0


# ---------------------------------------------------------------------------
# State setters
# ---------------------------------------------------------------------------

def test_set_state_TP(gas):
    gas.set_state_TP(1000.0, 5e5)
    assert abs(gas.temperature - 1000.0) < 0.1
    assert abs(gas.pressure - 5e5) < 1.0


def test_set_state_TPX(gas):
    gas.set_state_TPX(2000.0, 1e6, "H2:2, O2:1")
    assert abs(gas.temperature - 2000.0) < 0.1
    assert abs(gas.pressure - 1e6) < 1.0


def test_set_state_HP(gas):
    gas.set_state_TP(1500.0, 1e5)
    h = gas.enthalpy_mass
    p = gas.pressure
    gas.set_state_TP(300.0, 1e5)  # change state
    gas.set_state_HP(h, p)
    assert abs(gas.temperature - 1500.0) < 1.0


def test_set_state_SP(gas):
    gas.set_state_TP(1500.0, 1e5)
    s = gas.entropy_mass
    p = gas.pressure
    gas.set_state_TP(300.0, 1e5)  # change state
    gas.set_state_SP(s, p)
    assert abs(gas.temperature - 1500.0) < 1.0


# ---------------------------------------------------------------------------
# State save/restore
# ---------------------------------------------------------------------------

def test_save_restore_roundtrip(gas):
    gas.set_state_TP(1500.0, 2e6)
    T_orig = gas.temperature
    P_orig = gas.pressure

    state = gas.save_state()
    assert isinstance(state, list)
    assert len(state) > 0

    gas.set_state_TP(300.0, 1e5)
    assert abs(gas.temperature - 300.0) < 0.1

    gas.restore_state(state)
    assert abs(gas.temperature - T_orig) < 0.1
    assert abs(gas.pressure - P_orig) < 1.0


# ---------------------------------------------------------------------------
# Derived calculations
# ---------------------------------------------------------------------------

def test_stagnation_enthalpy(gas):
    gas.set_state_TP(300.0, 1e5)
    h = gas.enthalpy_mass
    v = 100.0
    h0 = gas.stagnation_enthalpy(v)
    assert abs(h0 - (h + 0.5 * v**2)) < 1.0


def test_stagnation_pressure(gas):
    gas.set_state_TP(300.0, 1e5)
    p0 = gas.stagnation_pressure(100.0)
    assert p0 > gas.pressure


def test_mach(gas):
    gas.set_state_TP(300.0, 1e5)
    a = gas.speed_of_sound
    m = gas.mach(a)
    assert abs(m - 1.0) < 1e-6


def test_isenthalpic_velocity_explicit(gas):
    gas.set_state_TP(1000.0, 1e5)
    h = gas.enthalpy_mass
    H0 = h + 0.5 * 500.0**2
    v = gas.isenthalpic_velocity(H0)
    assert abs(v - 500.0) < 1.0


def test_isenthalpic_velocity_stored(gas):
    gas.set_state_TP(1000.0, 1e5)
    h = gas.enthalpy_mass
    H0 = h + 0.5 * 500.0**2
    gas.stagnation_enthalpy_ref = H0
    v = gas.isenthalpic_velocity()
    assert abs(v - 500.0) < 1.0


# ---------------------------------------------------------------------------
# Expansion properties and equilibrium
# ---------------------------------------------------------------------------

def test_expansion_properties(gas):
    from goddard import ExpansionProperties
    gas.set_state_TP(1000.0, 1e5)
    ep = gas.expansion_properties()
    assert isinstance(ep, ExpansionProperties)
    assert ep.gamma_s > 1.0
    assert ep.spec_heat_p > 0


def test_equilibrate(eq_gas):
    eq_gas.set_state_TPX(3000.0, 1e6, "H2:2, O2:1")
    T_before = eq_gas.temperature
    eq_gas.equilibrate("HP")
    # Temperature should change after equilibration
    assert eq_gas.temperature != pytest.approx(T_before, rel=1e-3)


# ---------------------------------------------------------------------------
# Snapshot
# ---------------------------------------------------------------------------

def test_snapshot(gas):
    from goddard import ThermoStateInfo
    gas.set_state_TP(1000.0, 1e5)
    info = gas.snapshot()
    assert isinstance(info, ThermoStateInfo)
    assert abs(info.temperature - 1000.0) < 0.1
    assert abs(info.pressure - 1e5) < 1.0


# ---------------------------------------------------------------------------
# Solution accessor
# ---------------------------------------------------------------------------

def test_solution_property(gas):
    from goddard import SolutionHandle
    sol = gas.solution
    assert sol is not None


# ---------------------------------------------------------------------------
# Chemistry
# ---------------------------------------------------------------------------

def test_chemistry_read_write(gas):
    from goddard import GasChemistry
    assert gas.chemistry == GasChemistry.FROZEN
    gas.chemistry = GasChemistry.EQUILIBRIUM
    assert gas.chemistry == GasChemistry.EQUILIBRIUM


# ---------------------------------------------------------------------------
# Reference state properties
# ---------------------------------------------------------------------------

def test_stagnation_enthalpy_ref(gas):
    gas.stagnation_enthalpy_ref = 1e6
    assert abs(gas.stagnation_enthalpy_ref - 1e6) < 1.0


def test_reference_entropy(gas):
    gas.set_state_TP(300.0, 1e5)
    s = gas.entropy_mass
    gas.reference_entropy = s
    assert abs(gas.reference_entropy - s) < 1.0


# ---------------------------------------------------------------------------
# Convenience functions
# ---------------------------------------------------------------------------

def test_gas_from_yaml(h2o2_yaml):
    from goddard import gas_from_yaml, GasChemistry
    g = gas_from_yaml(h2o2_yaml, phase_name="ohmech",
                      chemistry=GasChemistry.FROZEN)
    assert g.temperature > 0
    assert g.pressure > 0


# ---------------------------------------------------------------------------
# Nozzle from Gas
# ---------------------------------------------------------------------------

def test_nozzle_from_gas(h2o2_yaml):
    from goddard import Gas, GasChemistry, Nozzle, ExpansionType

    g = Gas(h2o2_yaml, phase_name="ohmech", chemistry=GasChemistry.EQUILIBRIUM)
    g.set_state_TPX(3000.0, 2e6, "H2:2, O2:1")
    g.equilibrate("HP")

    nozzle = Nozzle(g)
    result = nozzle.solve(ExpansionType.SUPERSONIC_AREA_RATIO, 5.0)
    assert result.throat.converged
    assert len(result.expansions) == 1


# ---------------------------------------------------------------------------
# from_cantera round-trip
# ---------------------------------------------------------------------------

def test_from_cantera():
    ct = pytest.importorskip("cantera")
    from goddard import from_cantera, GasChemistry

    ct_sol = ct.Solution("h2o2.yaml", "ohmech")
    ct_sol.TPX = 2500.0, 3e6, "H2:2, O2:1"

    gas = from_cantera(ct_sol, chemistry=GasChemistry.FROZEN)

    assert abs(gas.temperature - 2500.0) < 1.0
    assert abs(gas.pressure - 3e6) < 100.0
