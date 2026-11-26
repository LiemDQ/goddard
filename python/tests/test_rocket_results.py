"""Tests for the RocketSolver front end: RocketPoint, RocketResult and RocketResults."""
import math
import os
from collections.abc import Sequence

import numpy as np
import pytest

import goddard
from conftest import find_data_dir


H2O2_SPECIES = {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"}
BAR = 1e5  # Pa


def h2o2_solver():
    """RocketSolver for H2/O2 on the h2o2.yaml mechanism (composition path)."""
    products = goddard.Gas(os.path.join(find_data_dir(), "h2o2.yaml"), "ohmech",
                           species=H2O2_SPECIES)
    return goddard.RocketSolver(products)


def h2o2_point(chamber_pressure, expansion_ratios, mixture_ratio=6.0, label="",
               expansion_type=goddard.ExpansionType.SUPERSONIC_AREA_RATIO,
               combustor=None, nozzle=None):
    """H2/O2 RocketPoint with gaseous reactants at 300 K."""
    return goddard.RocketPoint(
        fuel=goddard.PhaseSpecification(300.0, 101325.0, "H2:1"),
        oxidizer=goddard.PhaseSpecification(300.0, 101325.0, "O2:1"),
        mixture_ratio=mixture_ratio,
        chamber_pressure=chamber_pressure,
        combustor=combustor if combustor is not None else goddard.infinite_area_combustor(),
        nozzle=nozzle if nozzle is not None else goddard.equilibrium_nozzle(),
        expansion_type=expansion_type,
        expansion_ratios=list(expansion_ratios),
        label=label)


# ---------------------------------------------------------------------------
# Several chamber pressures
# ---------------------------------------------------------------------------

PRESSURES = [20.0 * BAR, 50.0 * BAR]
AREA_RATIOS = [5.0, 20.0]


@pytest.fixture(scope="module")
def two_pressure_results():
    points = [h2o2_point(pressure, AREA_RATIOS, label=f"pc_{index}")
              for index, pressure in enumerate(PRESSURES)]
    return h2o2_solver().solve(points)


def test_results_follow_the_order_of_the_points(two_pressure_results):
    """Result k belongs to point k, with its own chamber pressure and label."""
    results = two_pressure_results
    assert isinstance(results, goddard.RocketResults)
    assert isinstance(results, Sequence)
    assert len(results) == len(PRESSURES)
    for index, pressure in enumerate(PRESSURES):
        assert results[index].point.label == f"pc_{index}"
        assert results[index].chamber.thermo.pressure == pytest.approx(pressure, rel=1e-9)
    assert (results[0].chamber.thermo.pressure != results[1].chamber.thermo.pressure)
    assert [result.point.chamber_pressure for result in results] == PRESSURES
    assert list(results.point("chamber_pressure")) == PRESSURES


def test_performance_area_ratio_matches_requested_at_every_pressure(two_pressure_results):
    """performance[e] combines the throat and exit e of one point.

    The area ratio comes from mass flux continuity, rho_t u_t / (rho_e u_e). Taking the throat
    of one pressure and the exit of another would give about AR * P0/P1 instead of AR.
    Tolerance: ComparisonTolerances.area_ratio_rel (0.1%).
    """
    for result in two_pressure_results:
        assert len(result.exits) == len(AREA_RATIOS)
        assert len(result.performance) == len(AREA_RATIOS)
        for exit_index, area_ratio in enumerate(AREA_RATIOS):
            assert result.performance[exit_index].area_ratio == pytest.approx(
                area_ratio, rel=1e-3), f"label={result.point.label}, exit_index={exit_index}"


def test_infinite_area_chamber_stations_coincide(two_pressure_results):
    for result in two_pressure_results:
        assert result.chamber.thermo.temperature == result.stagnation.thermo.temperature
        assert result.chamber.thermo.temperature == result.combustion_end.thermo.temperature
        assert math.isinf(result.contraction_ratio)
        assert result.mass_flux == 0.0


def test_results_arrays_have_one_row_per_point(two_pressure_results):
    results = two_pressure_results
    ivac = results.performance("ivac")
    assert ivac.shape == (len(PRESSURES), len(AREA_RATIOS))
    for row, result in enumerate(results):
        for column, performance in enumerate(result.performance):
            assert ivac[row, column] == performance.ivac

    exit_areas = results.station("exits", "area_ratio")
    assert exit_areas.shape == (len(PRESSURES), len(AREA_RATIOS))
    assert exit_areas == pytest.approx(np.array([AREA_RATIOS] * len(PRESSURES)), rel=1e-3)

    chamber_temperature = results.station("chamber", "temperature")
    assert chamber_temperature.shape == (len(PRESSURES),)
    assert chamber_temperature[0] == results[0].chamber.thermo.temperature
    assert results.station("throat", "mach") == pytest.approx(np.ones(len(PRESSURES)),
                                                              rel=1e-6)

    with pytest.raises(ValueError):
        results.performance("not_a_field")
    with pytest.raises(ValueError):
        results.station("not_a_station", "temperature")


def test_results_arrays_pad_points_with_fewer_exits_with_nan():
    results = h2o2_solver().solve([
        h2o2_point(20.0 * BAR, [5.0, 10.0, 20.0]),
        h2o2_point(20.0 * BAR, [5.0]),
    ])
    ivac = results.performance("ivac")
    assert ivac.shape == (2, 3)
    assert np.isfinite(ivac[0]).all()
    assert np.isfinite(ivac[1, 0])
    assert np.isnan(ivac[1, 1:]).all()
    exit_pressure = results.station("exits", "pressure")
    assert exit_pressure.shape == (2, 3)
    assert np.isnan(exit_pressure[1, 1:]).all()


def test_results_slicing_and_dataframe(two_pressure_results):
    assert isinstance(two_pressure_results[:1], goddard.RocketResults)
    assert len(two_pressure_results[:1]) == 1
    pandas = pytest.importorskip("pandas")
    frame = two_pressure_results.to_dataframe()
    assert isinstance(frame, pandas.DataFrame)
    assert len(frame) == len(PRESSURES) * len(AREA_RATIOS)
    assert frame["chamber_pressure"].tolist() == [p for p in PRESSURES for _ in AREA_RATIOS]


def test_solve_accepts_a_generator_and_a_single_point():
    solver = h2o2_solver()
    generated = solver.solve(h2o2_point(pressure, [10.0]) for pressure in PRESSURES)
    assert isinstance(generated, goddard.RocketResults)
    assert len(generated) == len(PRESSURES)
    assert generated[1].chamber.thermo.pressure == pytest.approx(PRESSURES[1], rel=1e-9)

    single = solver.solve(h2o2_point(PRESSURES[0], [10.0]))
    assert isinstance(single, goddard.RocketResult)
    assert single.chamber.thermo.temperature == pytest.approx(
        generated[0].chamber.thermo.temperature, rel=1e-12)


def test_specific_impulse_in_seconds(two_pressure_results):
    performance = two_pressure_results[0].performance[0]
    assert performance.isp_s == pytest.approx(performance.isp / 9.80665, rel=1e-15)
    assert performance.ivac_s == pytest.approx(performance.ivac / 9.80665, rel=1e-15)


def test_rocket_performance_matches_the_solver(two_pressure_results):
    result = two_pressure_results[0]
    for index, performance in enumerate(result.performance):
        direct = goddard.rocket_performance(
            result.stagnation.thermo, result.throat.thermo, result.exits[index].thermo)
        for field in ("pressure_ratio", "area_ratio", "mach_number", "cstar", "CF",
                      "isp", "ivac"):
            assert getattr(direct, field) == pytest.approx(
                getattr(performance, field), rel=1e-12), field


# ---------------------------------------------------------------------------
# Pressure-ratio exits
# ---------------------------------------------------------------------------

def test_pressure_ratio_exit_reports_area_ratio():
    """A pressure-ratio exit reports Ae/At in area_ratio, matching its performance entry."""
    point = h2o2_point(20.0 * BAR, [10.0, 100.0],
                       expansion_type=goddard.ExpansionType.PRESSURE_RATIO)
    result = h2o2_solver().solve(point)
    assert len(result.exits) == 2
    for index, station in enumerate(result.exits):
        assert station.area_ratio > 1.0
        assert station.area_ratio == pytest.approx(result.performance[index].area_ratio,
                                                   rel=1e-6)
    assert result.exits[1].area_ratio > result.exits[0].area_ratio


# ---------------------------------------------------------------------------
# Convenience builders
# ---------------------------------------------------------------------------

def test_frozen_nfz_defaults_to_the_chamber_everywhere():
    """C++, the keyword constructor and frozen_nozzle() agree: frozen_NFZ = 0, the chamber.

    0 is CEA's default freezing point (nfz = frozen_NFZ + 1 = 1).
    """
    assert goddard.NozzleOptions().frozen_NFZ == 0
    assert goddard.NozzleOptions(chemistry=goddard.GasChemistry.FROZEN).frozen_NFZ == 0
    assert goddard.frozen_nozzle().frozen_NFZ == 0
    assert goddard.frozen_nozzle(frozen_NFZ=1).frozen_NFZ == 1
    assert goddard.frozen_nozzle().chemistry == goddard.GasChemistry.FROZEN


def test_convenience_builders_take_no_ratios_or_pressures():
    """Ratios and pressures belong to the RocketPoint; the old positional arguments fail."""
    assert goddard.equilibrium_nozzle().chemistry == goddard.GasChemistry.EQUILIBRIUM
    with pytest.raises(TypeError):
        goddard.equilibrium_nozzle(10.0)
    with pytest.raises(TypeError):
        goddard.frozen_nozzle(10.0)
    with pytest.raises(TypeError):
        goddard.infinite_area_combustor([20.0 * BAR])
    with pytest.raises(TypeError):
        goddard.finite_mass_flux_combustor(1000.0, [20.0 * BAR])
    with pytest.raises(TypeError):
        goddard.finite_contraction_ratio_combustor(2.0, [20.0 * BAR])

    assert goddard.infinite_area_combustor().type == goddard.CombustorType.INFINITE_AREA
    assert goddard.finite_mass_flux_combustor(1000.0).mass_flux == 1000.0
    assert goddard.finite_contraction_ratio_combustor(2.0).contraction_ratio == 2.0
    assert not hasattr(goddard.CombustorOptions(), "pressures")
    assert not hasattr(goddard.NozzleOptions(), "expansion_ratios")


# ---------------------------------------------------------------------------
# Options for planned features raise instead of being ignored
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("attribute, value", [
    ("include_transport", True),
    ("include_ionized_species", True),
    ("trace_cutoff", 1e-5),
])
def test_unimplemented_solver_options_raise(attribute, value):
    solver = h2o2_solver()
    setattr(solver, attribute, value)
    with pytest.raises(NotImplementedError):
        solver.solve(h2o2_point(20.0 * BAR, [5.0]))


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def test_report_labels_mixture_ratio_by_type(two_pressure_results):
    assert "O/F=" in two_pressure_results[0].report()
    assert "O/F=" in two_pressure_results.report()

    combustor = goddard.infinite_area_combustor()
    combustor.mixture_type = goddard.MixtureRatioType.FUEL_FRAC
    point = h2o2_point(20.0 * BAR, [5.0], mixture_ratio=0.15, combustor=combustor)
    report = h2o2_solver().solve(point).report()
    assert "FUEL FRACTION=" in report
    assert "O/F=" not in report


def test_report_prints_the_point_label(two_pressure_results):
    report = two_pressure_results.report()
    assert "CASE = pc_0" in report
    assert "CASE = pc_1" in report
    assert goddard.rocket_report(two_pressure_results[0]).count("CASE = pc_0") == 1
    assert "CASE = pc_1" not in two_pressure_results[0].report()
