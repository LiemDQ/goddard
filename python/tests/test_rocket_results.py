"""Tests for RocketProblemResults accessors and the RocketProblem convenience helpers."""
import os

import pytest

import goddard
from conftest import find_data_dir


H2O2_SPECIES = {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"}
BAR = 1e5  # Pa


def build_h2o2_problem(name, pressures, area_ratios, of_ratios=(6.0,), nozzle_options=None):
    """H2/O2 RocketProblem with gaseous reactants at 300 K on h2o2.yaml (composition path)."""
    chem_params = goddard.ChemicalParameters()
    chem_params.thermo_file = os.path.join(find_data_dir(), "h2o2.yaml")
    chem_params.species = H2O2_SPECIES
    chem_params.fuel_state = goddard.PhaseSpecification(300.0, pressures[0], "H2:1")
    chem_params.oxidizer_state = goddard.PhaseSpecification(300.0, pressures[0], "O2:1")
    chem_params.OF_ratios = list(of_ratios)

    case_params = goddard.RocketCaseParameters()
    case_params.name = name
    case_params.problem_type = "rocket"
    case_params.combustor_options = goddard.infinite_area_combustor(pressures)
    case_params.nozzle_options = (nozzle_options if nozzle_options is not None
                                  else goddard.equilibrium_nozzle(*area_ratios))
    return goddard.RocketProblem(chem_params, [case_params], "ohmech")


# ---------------------------------------------------------------------------
# Several chamber pressures
# ---------------------------------------------------------------------------

PRESSURES = [20.0 * BAR, 50.0 * BAR]
AREA_RATIOS = [5.0, 20.0]


@pytest.fixture(scope="module")
def two_pressure_results():
    return build_h2o2_problem("two_pressures", PRESSURES, AREA_RATIOS).solve()


def test_chamber_accessor_selects_pressure(two_pressure_results):
    """chamber() returns the station of the requested chamber pressure."""
    results = two_pressure_results
    for p_idx, pressure in enumerate(PRESSURES):
        chamber = results.chamber(0, pressure_index=p_idx, case_name="two_pressures")
        assert chamber.pressure_index == p_idx
        assert chamber.thermo.pressure == pytest.approx(pressure, rel=1e-9)
    assert (results.chamber(0, pressure_index=0).thermo.pressure
            != results.chamber(0, pressure_index=1).thermo.pressure)


def test_performance_area_ratio_matches_requested_at_every_pressure(two_pressure_results):
    """performance() combines throat and exit of one pressure.

    The area ratio comes from mass flux continuity, rho_t a_t / (rho_e v_e). Taking the throat
    of one pressure and the exit of another would give about AR * P0/P1 instead of AR.
    Tolerance: ComparisonTolerances.area_ratio_rel (0.1%).
    """
    results = two_pressure_results
    for p_idx in range(len(PRESSURES)):
        exits = results.exits(0, pressure_index=p_idx)
        assert len(exits) == len(AREA_RATIOS)
        for exit_index, area_ratio in enumerate(AREA_RATIOS):
            assert exits[exit_index].pressure_index == p_idx
            performance = results.performance(0, pressure_index=p_idx, exit_index=exit_index)
            assert performance.area_ratio == pytest.approx(area_ratio, rel=1e-3), (
                f"pressure_index={p_idx}, exit_index={exit_index}")


def test_accessors_reject_positional_arguments_after_of_index(two_pressure_results):
    """pressure_index was inserted before exit_index and case_name, so they are keyword-only.

    An old call such as performance(0, 1) must fail rather than read another operating point.
    """
    results = two_pressure_results
    with pytest.raises(TypeError):
        results.performance(0, 1)
    with pytest.raises(TypeError):
        results.chamber(0, "two_pressures")


# ---------------------------------------------------------------------------
# Nozzle option defaults
# ---------------------------------------------------------------------------

def test_frozen_nfz_defaults_to_the_chamber_everywhere():
    """C++, the keyword constructor and frozen_nozzle() agree: frozen_NFZ = 0, the chamber.

    0 is CEA's default freezing point (nfz = frozen_NFZ + 1 = 1).
    """
    assert goddard.NozzleOptions().frozen_NFZ == 0
    assert goddard.NozzleOptions(chemistry=goddard.GasChemistry.FROZEN).frozen_NFZ == 0
    assert goddard.frozen_nozzle(10.0).frozen_NFZ == 0


# ---------------------------------------------------------------------------
# Combustor helpers accept numpy arrays
# ---------------------------------------------------------------------------

def test_combustor_helpers_accept_numpy_pressures():
    """A numpy array of pressures must not be tested for truth value (ambiguous for size > 1)."""
    np = pytest.importorskip("numpy")
    pressures = np.array([20.0 * BAR, 50.0 * BAR])

    for opts in (goddard.infinite_area_combustor(pressures),
                 goddard.finite_mass_flux_combustor(1000.0, pressures),
                 goddard.finite_contraction_ratio_combustor(2.0, pressures)):
        assert list(opts.pressures) == pytest.approx(list(pressures))

    assert list(goddard.infinite_area_combustor().pressures) == []
