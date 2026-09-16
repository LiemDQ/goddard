"""Tests for the Gas Python binding."""
import os
import math

import numpy as np
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
    from goddard import ThermodynamicState
    gas.set_state_TP(1000.0, 1e5)
    info = gas.snapshot()
    assert isinstance(info, ThermodynamicState)
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
    # A gas-only mixture can never sit at a condensed phase transition.
    assert result.throat.pinned_transition is False
    assert result.expansions[0].pinned_transition is False


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


# ---------------------------------------------------------------------------
# Condensed species
# ---------------------------------------------------------------------------

@pytest.fixture
def nasa9_files():
    """Paths to the NASA9 gas, condensed and reactant data files."""
    data_dir = find_data_dir()
    return (os.path.join(data_dir, "nasa9_gas.yaml"),
            os.path.join(data_dir, "nasa9_condensed.yaml"),
            os.path.join(data_dir, "nasa9_reactants.yaml"))


@pytest.fixture
def water_gas(nasa9_files):
    """H2/O2 products carrying the two water candidates."""
    from goddard import Gas, GasChemistry

    gas_file, condensed_file, _ = nasa9_files
    return Gas(gas_file, phase_name="products",
               species={"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2"},
               condensed_file=condensed_file,
               condensed_species={"H2O(L)", "H2O(cr)"},
               chemistry=GasChemistry.EQUILIBRIUM)


def test_condensed_candidates(water_gas):
    assert water_gas.has_condensed_candidates
    assert not water_gas.has_condensed_phases
    assert set(water_gas.condensed_species_names) == {"H2O(L)", "H2O(cr)"}
    assert water_gas.condensed_moles == [0.0, 0.0]
    assert water_gas.gas_mass_fraction == 1.0
    assert not water_gas.at_phase_transition
    assert water_gas.pinned_polymorphs == (-1, -1)


def test_condensed_state_vector_is_extended(water_gas, h2o2_yaml):
    from goddard import Gas

    gas_only = Gas(h2o2_yaml, phase_name="ohmech")
    assert len(water_gas.save_state()) == water_gas.num_species + 2 + 2
    assert len(gas_only.save_state()) == gas_only.num_species + 2


def test_no_condensed_candidates_by_default(gas):
    assert not gas.has_condensed_candidates
    assert gas.condensed_species_names == []
    assert gas.gas_mass_fraction == 1.0
    assert gas.mixture_molecular_weight == pytest.approx(gas.mean_molecular_weight)


def test_equilibrate_tp_condenses_water(water_gas, nasa9_files):
    import goddard
    from goddard import reactant_gas

    _, _, reactant_file = nasa9_files
    stream = reactant_gas(reactant_file, {"H2": 100.0, "O2": 60.0}, 298.15, 101325.0)
    amounts = dict(zip(stream.element_names, stream.element_moles))
    elements = np.array([amounts[element] for element in water_gas.element_names])

    pressure = 0.05 * 101325.0
    water_gas.set_element_moles(elements, 300.0, pressure)
    water_gas.equilibrate_TP(300.0, pressure)

    liquid = water_gas.condensed_moles[water_gas.condensed_species_names.index("H2O(L)")]
    assert liquid > 0.0
    assert water_gas.has_condensed_phases
    assert water_gas.gas_mass_fraction < 1.0
    assert water_gas.mixture_molecular_weight < water_gas.mean_molecular_weight
    assert water_gas.last_equilibrium_solve_count() > 0
    assert sum(water_gas.mixture_mass_fractions) == pytest.approx(1.0)

    # One entry per *present* condensed species, in candidate order.
    assert len(water_gas.condensed_enthalpy_RT) == 1
    assert len(water_gas.condensed_cp_R) == 1
    assert water_gas.condensed_molar_masses[0] == pytest.approx(18.0153, rel=1e-3)
    assert water_gas.condensed_stoich_coeffs.shape == (1, len(water_gas.element_names))

    derivatives = goddard.equilibrium_derivatives(water_gas)
    assert len(derivatives.dn_condensed_dlogT_P) == 1
    assert len(derivatives.dn_condensed_dlogP_T) == 1
    assert not derivatives.pinned_transition


def test_equilibrium_options_are_writable(water_gas):
    from goddard import EquilibriumOptions

    assert water_gas.equilibrium_options.max_steps == 20000
    options = EquilibriumOptions()
    options.max_steps = 5000
    water_gas.equilibrium_options = options
    assert water_gas.equilibrium_options.max_steps == 5000


def test_equilibrium_properties_of_a_gas(eq_gas):
    import goddard

    eq_gas.set_state_TPX(3000.0, 2e6, "H2:2, O2:1")
    eq_gas.equilibrate("HP")

    properties = goddard.equilibrium_properties(eq_gas)
    frozen = goddard.frozen_properties(eq_gas)
    derivatives = goddard.equilibrium_derivatives(eq_gas)

    assert properties.gamma_s > 1.0
    assert properties.spec_heat_p > frozen.spec_heat_p
    assert not properties.pinned_transition
    assert len(derivatives.dn_condensed_dlogT_P) == 0
    assert properties.speed_of_sound == pytest.approx(eq_gas.speed_of_sound, rel=1e-9)
