"""Failure reporting: SolveStatus, FailureReason, the exceptions, and `.status` on every result."""
import math
import os

import numpy as np
import pytest

import goddard
from goddard import (
    ConvergenceError, FailureReason, Gas, GasChemistry, SolveError, SolveStatus,
)
from conftest import find_data_dir
from test_rocket_results import h2o2_point, h2o2_solver

BAR = 1e5  # Pa


@pytest.fixture(scope="module")
def h2o2_yaml():
    return os.path.join(find_data_dir(), "h2o2.yaml")


def h2o2_gas(h2o2_yaml, composition="H2:2, O2:1", T=300.0, P=1e5):
    gas = Gas(h2o2_yaml, phase_name="ohmech")
    gas.set_state_TPX(T, P, composition)
    return gas


# ---------------------------------------------------------------------------
# Types
# ---------------------------------------------------------------------------

def test_failure_reason_members():
    names = {"NONE", "NOT_CONVERGED", "NO_SOLUTION", "OUT_OF_DATA_RANGE", "THERMO_FAILURE",
             "INCOMPLETE", "NOT_COMPUTED"}
    assert names <= set(dir(FailureReason))
    assert len({getattr(FailureReason, name) for name in names}) == len(names)


def test_exception_hierarchy():
    assert issubclass(SolveError, RuntimeError)
    assert issubclass(ConvergenceError, SolveError)
    assert not issubclass(SolveError, ConvergenceError)


def failed_status():
    """A NO_SOLUTION status: a subsonic normal shock."""
    return goddard.normal_shock(0.5, 1.4).status


def test_status_attributes_of_a_failed_solve():
    status = failed_status()
    assert isinstance(status, SolveStatus)
    assert status.ok is False
    assert status.reason == FailureReason.NO_SOLUTION
    assert isinstance(status.message, str) and status.message
    assert isinstance(status.stage, str)
    assert status.iterations == -1
    assert math.isnan(status.tolerance)
    assert math.isnan(status.residual)
    assert list(status.warnings) == []
    assert "NO_SOLUTION" in status.describe()
    assert status.message in status.describe()


def test_status_attributes_of_a_successful_solve():
    status = goddard.normal_shock(2.0, 1.4).status
    assert status.ok is True
    assert status.reason == FailureReason.NONE
    assert status.message == ""
    assert status.describe() == "ok"
    status.raise_if_failed()  # does nothing


def test_status_is_read_only():
    status = failed_status()
    with pytest.raises(AttributeError):
        status.reason = FailureReason.NONE
    with pytest.raises(AttributeError):
        status.ok = True


def test_raise_if_failed_raises_solve_error_for_no_solution():
    status = failed_status()
    with pytest.raises(SolveError) as info:
        status.raise_if_failed()
    assert not isinstance(info.value, ConvergenceError)
    assert info.value.status.reason == FailureReason.NO_SOLUTION
    assert info.value.status.message == status.message


def test_raise_if_failed_raises_convergence_error_with_iteration_data(h2o2_yaml):
    """An equilibrium shock below the CJ speed of an exothermic mixture does not converge."""
    gas = h2o2_gas(h2o2_yaml)
    status = goddard.ShockSolver(gas, GasChemistry.EQUILIBRIUM).normal_shock(2.0).status
    assert status.reason == FailureReason.NOT_CONVERGED
    with pytest.raises(ConvergenceError) as info:
        status.raise_if_failed()
    error = info.value
    assert error.status.reason == FailureReason.NOT_CONVERGED
    assert error.iterations == status.iterations > 0
    assert error.tolerance == status.tolerance > 0.0
    assert error.residual == status.residual
    with pytest.raises(SolveError):  # a ConvergenceError is also a SolveError
        status.raise_if_failed()


# ---------------------------------------------------------------------------
# Exceptions that leave the library
# ---------------------------------------------------------------------------

def test_unknown_species_and_files_are_value_errors(h2o2_yaml):
    gas = h2o2_gas(h2o2_yaml)
    with pytest.raises(ValueError, match="NOT_A_SPECIES"):
        gas.set_state_TPX(300.0, 1e5, "H2:1, NOT_A_SPECIES:1")
    with pytest.raises(ValueError):
        Gas(os.path.join(find_data_dir(), "no_such_file.yaml"), phase_name="ohmech")
    with pytest.raises(ValueError):
        Gas(h2o2_yaml, phase_name="no_such_phase")


def test_misspelt_equilibrium_request_is_a_value_error(h2o2_yaml):
    gas = h2o2_gas(h2o2_yaml)
    with pytest.raises(ValueError, match="PT"):
        gas.equilibrate("PT")
    with pytest.raises(ValueError, match="gibs"):
        gas.equilibrate("TP", "gibs")


def test_unknown_species_in_a_combustor_is_a_value_error(h2o2_yaml):
    gas = h2o2_gas(h2o2_yaml)
    combustor = goddard.Combustor(gas, {"H2": 1.0}, {"XX": 1.0})
    with pytest.raises(ValueError):
        combustor.solve(300.0, 300.0, np.array([10 * BAR]), np.array([6.0]),
                        goddard.infinite_area_combustor())


def test_rejected_state_is_a_solve_error_not_a_backend_error(h2o2_yaml):
    gas = h2o2_gas(h2o2_yaml)
    with pytest.raises(SolveError) as info:
        gas.set_state_TP(0.0, 1e5)
    error = info.value
    assert type(error).__module__.startswith("goddard")
    assert not isinstance(error, ConvergenceError)
    assert error.status.reason == FailureReason.THERMO_FAILURE
    assert "set_state_TP" in error.status.message


def test_scalar_shock_functions_raise_instead_of_returning_nan():
    with pytest.raises(SolveError) as info:
        goddard.oblique_shock_max_deflection(0.8, 1.4)
    assert info.value.status.reason == FailureReason.NO_SOLUTION
    with pytest.raises(ValueError):
        goddard.oblique_shock_max_deflection(2.0, 0.9)  # gamma <= 1 is an invalid argument


def test_no_backend_exception_type_is_exposed():
    assert not hasattr(goddard, "CanteraError")
    assert not hasattr(goddard, "FmtError")


# ---------------------------------------------------------------------------
# `.status` on the result classes
# ---------------------------------------------------------------------------

def test_shock_results_carry_status_and_nan_fields():
    ok = goddard.reflected_shock(2.0, 1.4)
    assert ok.status.ok and ok.incident.status.ok and ok.reflected.status.ok
    assert not hasattr(ok, "valid")

    bad = goddard.reflected_shock(0.5, 1.4)
    assert bad.status.reason == FailureReason.NO_SOLUTION
    assert bad.incident.status.reason == FailureReason.NO_SOLUTION
    assert bad.reflected.status.reason == FailureReason.NOT_COMPUTED
    assert math.isnan(bad.incident.mach_out)
    assert math.isnan(bad.reflected.static_pressure_ratio)

    oblique = goddard.oblique_shock_from_deflection(
        3.0, goddard.oblique_shock_max_deflection(3.0, 1.4) + 1e-3, 1.4)
    assert oblique.status.reason == FailureReason.NO_SOLUTION
    assert math.isnan(oblique.beta) and math.isnan(oblique.mach_out)
    assert goddard.oblique_shock_from_deflection(3.0, 0.2, 1.4).status.ok


def test_shock_solver_post_state_follows_the_most_recent_solve(h2o2_yaml):
    solver = goddard.ShockSolver(h2o2_gas(h2o2_yaml, "N2:0.79, O2:0.21"), GasChemistry.FROZEN)
    with pytest.raises(SolveError) as info:
        solver.post_shock_state()
    assert info.value.status.reason == FailureReason.NOT_COMPUTED

    assert solver.normal_shock(2.0).status.ok
    assert solver.post_shock_state().temperature > 300.0

    assert solver.normal_shock(0.5).status.reason == FailureReason.NO_SOLUTION
    with pytest.raises(SolveError) as info:
        solver.post_shock_state()
    assert info.value.status.reason == FailureReason.NO_SOLUTION


def test_detonation_results_carry_status(h2o2_yaml):
    gamma, heat_release = 1.4, 10.0
    assert goddard.chapman_jouguet_detonation(gamma, heat_release).status.ok
    slow = goddard.detonation(0.9, gamma, heat_release)
    assert slow.status.reason == FailureReason.NO_SOLUTION
    assert math.isnan(slow.mach_in) and math.isnan(slow.static_pressure_ratio)
    with pytest.raises(ValueError):
        goddard.chapman_jouguet_detonation(1.0, heat_release)

    solver = goddard.DetonationSolver(h2o2_gas(h2o2_yaml))
    reflected = solver.reflected_detonation(0.99)
    assert reflected.status.reason == FailureReason.NO_SOLUTION
    assert reflected.reflected.status.reason == FailureReason.NOT_COMPUTED
    assert solver.reflected_detonation().status.ok


def test_detonation_solver_is_lazy_and_failures_are_in_the_result(h2o2_yaml):
    # Argon releases no heat. Construction does not solve, so it cannot fail.
    inert = goddard.DetonationSolver(h2o2_gas(h2o2_yaml, "AR:1"))
    with pytest.raises(SolveError) as info:
        inert.post_detonation_state()  # nothing has been solved yet
    assert info.value.status.reason == FailureReason.NOT_COMPUTED
    result = inert.chapman_jouguet()
    assert result.status.reason == FailureReason.NO_SOLUTION
    assert math.isnan(result.velocity)
    with pytest.raises(SolveError) as info:
        inert.post_detonation_state()
    assert info.value.status.reason == FailureReason.NO_SOLUTION

    solver = goddard.DetonationSolver(h2o2_gas(h2o2_yaml))
    with pytest.raises(SolveError):
        solver.von_neumann_state()
    cj = solver.chapman_jouguet()
    assert cj.status.ok
    assert solver.chapman_jouguet().velocity == cj.velocity  # cached
    assert solver.post_detonation_state().temperature > 2000.0


def test_nozzle_results_carry_status_per_station(h2o2_yaml):
    nozzle = goddard.Nozzle(h2o2_gas(h2o2_yaml, "H2O:1", 3500.0, 70 * BAR))
    results = nozzle.solve_ratios(goddard.ExpansionType.SUPERSONIC_AREA_RATIO, [5.0, 10.0])
    assert results.status.ok and results.throat.status.ok
    assert all(station.status.ok for station in results.expansions)
    for name in ("valid", "converged"):
        assert not hasattr(results, name)

    with pytest.raises(ValueError):
        nozzle.solve(goddard.ExpansionType.SUPERSONIC_AREA_RATIO, 1.0)
    with pytest.raises(ValueError):
        nozzle.solve(goddard.ExpansionType.PRESSURE_RATIO, 0.5)


def test_failed_throat_makes_every_station_not_computed(h2o2_yaml):
    """Cold H2/O2 is not at equilibrium, which the equilibrium throat iteration does not tolerate."""
    gas = h2o2_gas(h2o2_yaml, "H2:2, O2:1", 300.0, 70 * BAR)
    results = goddard.Nozzle(gas).solve_ratios(
        goddard.ExpansionType.SUPERSONIC_AREA_RATIO, [5.0, 10.0])
    assert results.throat.status.reason == FailureReason.NOT_CONVERGED
    assert math.isnan(results.throat.velocity)
    assert [station.status.reason for station in results.expansions] == [FailureReason.NOT_COMPUTED] * 2
    assert all(math.isnan(station.mach) for station in results.expansions)
    assert results.status.reason == FailureReason.NOT_CONVERGED


def test_finite_area_chamber_status_and_choked_flux(h2o2_yaml):
    gas = h2o2_gas(h2o2_yaml, "H2O:1", 3500.0, 70 * BAR)
    gas.equilibrate("HP")
    nozzle = goddard.Nozzle(gas)
    injector_state = gas.save_state() if hasattr(gas, "save_state") else None
    if injector_state is None:
        pytest.skip("Gas.save_state is not exposed")
    chamber = nozzle.solve_finite_area_chamber(
        injector_state, goddard.CombustorType.FINITE_CONTRACTION_RATIO, 3.0)
    assert chamber.status.ok
    choked = nozzle.solve_finite_area_chamber(
        injector_state, goddard.CombustorType.FINITE_MASS_FLUX, 1e6)
    assert choked.status.reason == FailureReason.NO_SOLUTION
    assert choked.injector.status.ok


def test_kinetic_nozzle_status_incomplete_when_steps_run_out(h2o2_yaml):
    gas = h2o2_gas(h2o2_yaml, "H2O:0.9, H2:0.05, O2:0.03, OH:0.02", 3500.0, 70 * BAR)
    gas.equilibrate("HP")
    profile = goddard.conical_nozzle(4.0)
    nozzle = goddard.KineticNozzle(gas, profile, 0.01, 1.0)
    truncated = nozzle.solve(dt_max=1e-6, dx_max=0.1, max_steps=3)
    assert truncated.status.reason == FailureReason.INCOMPLETE
    assert len(truncated.stations) == 3
    assert all(math.isfinite(station.velocity) for station in truncated.stations)
    assert not hasattr(truncated, "reached_exit")
    with pytest.raises(ValueError):
        nozzle.solve(dt_max=0.0)
    with pytest.raises(ValueError):
        nozzle.solve(max_steps=0)


def test_moc_result_status_agrees_with_failure_code():
    ok = goddard.moc_analysis(goddard.conical_nozzle(4.0))
    assert isinstance(ok.status, SolveStatus)
    assert ok.status.ok and ok.failure.code == goddard.MocErrorCode.NONE
    assert not hasattr(ok, "converged")
    assert hasattr(ok, "reached_exit_plane")

    # A forced Kliegel-Levine start line that misses the wall angle fails at initialization
    options = goddard.MocOptions(
        start_line=goddard.MocStartLine.KLIEGEL_LEVINE, kl_max_wall_angle_error=0.01,
        flow_type=goddard.MocFlowKind.AXISYMMETRIC, chemistry=GasChemistry.PERFECT_GAS,
        mode=goddard.MocMode.ANALYSIS, gamma=1.4, num_characteristics=8,
        nozzle_profile=goddard.conical_nozzle(8.0))
    options.geometry.throat_radius = 1.0
    failed = goddard.MocNozzle(options).solve()
    assert failed.failure.code == goddard.MocErrorCode.INITIALIZATION_FAILED
    assert failed.status.reason == FailureReason.NOT_CONVERGED
    assert failed.status.stage == "MoC initialization"
    assert "INITIALIZATION_FAILED" in failed.status.message


# ---------------------------------------------------------------------------
# ThermoArray
# ---------------------------------------------------------------------------

def combustion_states(h2o2_yaml, fuel_temperature):
    gas = Gas(h2o2_yaml, phase_name="ohmech")
    combustor = goddard.Combustor(gas, {"H2": 1.0}, {"O2": 1.0})
    pressures = np.array([10 * BAR, 20 * BAR])
    return combustor.solve(fuel_temperature, 300.0, pressures, np.array([4.0, 6.0]),
                           goddard.infinite_area_combustor())


def test_thermoarray_without_failures(h2o2_yaml):
    states = combustion_states(h2o2_yaml, 300.0)
    assert states.num_failed == 0
    assert states.status.ok
    assert states.ok.dtype == bool and states.ok.shape == states.shape and states.ok.all()
    assert states.status_at((0, 0)).ok
    assert np.isfinite(states.temperature).all()


def test_thermoarray_reports_failed_operating_points(h2o2_yaml):
    """A 1e5 K fuel stream has an enthalpy no combustion temperature can match."""
    states = combustion_states(h2o2_yaml, 1e5)
    assert states.num_failed == states.size == 4
    assert not states.ok.any()
    assert states.status.reason == FailureReason.THERMO_FAILURE
    assert "combustion at O/F" in states.status.stage
    assert states.status_at((0, 1)).reason == FailureReason.THERMO_FAILURE
    assert np.isnan(states.temperature).all()
    assert repr(states) == "<ThermoArray shape=(2, 2) 4 of 4 entries failed>"
    # A failed entry holds the state it had before the failed solve, not a solution.
    with pytest.raises(SolveError) as caught:
        states.get_state((0, 0))
    assert caught.value.status.reason == FailureReason.THERMO_FAILURE
    with pytest.raises(SolveError):
        states[0, 0]
    with pytest.raises(IndexError):
        states.status_at((2, 0))


def test_combustor_rejects_nonpositive_temperature_and_pressure(h2o2_yaml):
    gas = Gas(h2o2_yaml, phase_name="ohmech")
    combustor = goddard.Combustor(gas, {"H2": 1.0}, {"O2": 1.0})
    options = goddard.infinite_area_combustor()
    with pytest.raises(ValueError):
        combustor.solve(0.0, 300.0, np.array([10 * BAR]), np.array([6.0]), options)
    with pytest.raises(ValueError):
        combustor.solve(300.0, 300.0, np.array([-1.0]), np.array([6.0]), options)


# ---------------------------------------------------------------------------
# RocketSolver
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def partial_failure_results():
    """At an area ratio of 1e5 the equilibrium exit cannot be solved; the exit at 10 can."""
    points = [h2o2_point(70 * BAR, [10.0, 1e5], mixture_ratio=mixture_ratio, label="sweep")
              for mixture_ratio in (4.0, 6.0)]
    return h2o2_solver().solve(points)


def test_failed_exit_does_not_abort_the_sweep(partial_failure_results):
    results = partial_failure_results
    assert results.failed() == [0, 1]
    assert len(results) == 2

    assert results.status.reason == FailureReason.THERMO_FAILURE
    assert not results.status.ok
    for result in results:
        assert not result.status.ok
        assert result.chamber.status.ok
        assert result.throat.status.ok
        assert len(result.exits) == 2
        assert result.exits[0].status.ok and not result.exits[1].status.ok
        assert result.exits[1].status.reason == FailureReason.THERMO_FAILURE
        assert math.isnan(result.exits[1].thermo.temperature)
        assert math.isnan(result.exits[1].area_ratio)


def test_performance_of_failed_exit_is_nan(partial_failure_results):
    for result in partial_failure_results:
        good = result.performance[0]
        assert good.status.ok and good.isp > 2000.0
        bad = result.performance[1]
        assert bad.status.reason == FailureReason.THERMO_FAILURE
        for name in ("isp", "ivac", "CF", "cstar", "mach_number", "pressure_ratio", "area_ratio"):
            assert math.isnan(getattr(bad, name)), name

    ivac = partial_failure_results.performance("ivac")
    assert np.isfinite(ivac[:, 0]).all() and np.isnan(ivac[:, 1]).all()


def test_report_marks_failed_stations(partial_failure_results):
    report = partial_failure_results[0].report()
    assert "STATIONS NOT SOLVED (NaN ABOVE)" in report
    assert "OPERATING POINT NOT SOLVED" not in report


def test_rocket_solver_argument_errors():
    solver = h2o2_solver()
    good = h2o2_point(70 * BAR, [10.0])

    unset = goddard.RocketPoint(mixture_ratio=6.0, chamber_pressure=70 * BAR,
                                expansion_ratios=[10.0])
    with pytest.raises(ValueError):
        solver.solve(unset)
    with pytest.raises(ValueError):
        solver.validate(unset)

    # Nothing is solved when a later point of a list is invalid; the message names it.
    with pytest.raises(ValueError, match="point 1: "):
        solver.solve([good, unset])

    # A reactant species that is not a product species, with no reactant file.
    missing_species = h2o2_point(70 * BAR, [10.0])
    missing_species.oxidizer = goddard.PhaseSpecification(300.0, 70 * BAR, "NOT_A_SPECIES:1")
    with pytest.raises(ValueError):
        solver.solve(missing_species)

    solver.validate(good)  # does nothing
