"""Smoke tests to verify all bindings are importable and basic types work."""
import pytest


def test_import_module():
    import goddard
    assert hasattr(goddard, "__version__")


def test_import_enums():
    from goddard import CombustorType, GasChemistry, ExpansionType

    assert CombustorType.INFINITE_AREA is not None
    assert CombustorType.FINITE_MASS_FLUX is not None
    assert CombustorType.FINITE_CONTRACTION_RATIO is not None
    assert CombustorType.NONE is not None

    assert GasChemistry.FROZEN is not None
    assert GasChemistry.EQUILIBRIUM is not None
    assert GasChemistry.KINETIC is not None

    assert ExpansionType.SUPERSONIC_AREA_RATIO is not None
    assert ExpansionType.SUBSONIC_AREA_RATIO is not None
    assert ExpansionType.PRESSURE_RATIO is not None


def test_import_structs():
    from goddard import (
        CombustorOptions, NozzleOptions, ThroatCondition,
        NozzleResult, NozzleResults, RocketCaseParameters,
        ChemicalParameters, ThermoStateInfo, RocketPerformance,
        RocketState, ExpansionProperties, EquilibriumDerivatives,
    )


def test_import_classes():
    from goddard import (
        RocketProblem, RocketProblemResults, ThermoArray,
        Combustor, Nozzle,
        MixtureRatio, MixtureRatios, SolutionHandle,
    )


def test_import_kinetic_nozzle():
    from goddard import KineticNozzle, KineticNozzleStation, KineticNozzleResults
    assert KineticNozzle is not None
    assert KineticNozzleStation is not None
    assert KineticNozzleResults is not None


def test_import_gas():
    from goddard import Gas
    assert Gas is not None


def test_import_shocks():
    from goddard import (
        ShockResult, ObliqueShockResult, ShockSolver,
        normal_shock, reflected_shock,
        oblique_shock_wave_angle, oblique_shock_deflection_angle,
        oblique_shock_from_wave_angle, oblique_shock_from_deflection,
    )
    assert ShockResult is not None
    assert ShockSolver is not None


def test_perfect_gas_normal_shock():
    from goddard import normal_shock
    r = normal_shock(2.0, 1.4)
    assert r.valid
    assert r.mach_out < 1.0
    assert r.static_pressure_ratio > 1.0


def test_nozzle_options_new_fields():
    from goddard import NozzleOptions, SolverOptions
    opts = NozzleOptions()
    assert opts.dt_max == 1e-6
    assert opts.dx_max == 1e-3
    assert opts.max_steps == 100000
    opts.solver = SolverOptions(abstol=1e-8)
    assert opts.solver.abstol == 1e-8


def test_import_errors():
    from goddard import ConvergenceError
    assert issubclass(ConvergenceError, Exception)


def test_import_convenience():
    from goddard import (
        OF_ratio, supersonic_ratio, subsonic_ratio, pressure_ratio,
        infinite_area_combustor, finite_mass_flux_combustor,
        finite_contraction_ratio_combustor,
        equilibrium_nozzle, frozen_nozzle, from_cantera,
    )


def test_combustor_options_construction():
    from goddard import CombustorOptions, CombustorType

    opts = CombustorOptions()
    opts.type = CombustorType.INFINITE_AREA
    opts.pressures = [1e6, 2e6]
    assert opts.type == CombustorType.INFINITE_AREA
    assert opts.pressures == [1e6, 2e6]


def test_nozzle_options_construction():
    from goddard import NozzleOptions, GasChemistry, ExpansionType

    opts = NozzleOptions()
    opts.chemistry = GasChemistry.EQUILIBRIUM
    opts.expansion_type = ExpansionType.SUPERSONIC_AREA_RATIO
    opts.expansion_ratios = [3.0, 5.0, 10.0]
    assert opts.chemistry == GasChemistry.EQUILIBRIUM
    assert opts.expansion_ratios == [3.0, 5.0, 10.0]


def test_convenience_of_ratio():
    from goddard import OF_ratio
    assert OF_ratio(2.0, 3.0, 4.0) == [2.0, 3.0, 4.0]


def test_convenience_supersonic_ratio():
    from goddard import supersonic_ratio, ExpansionType, GasChemistry

    opts = supersonic_ratio(3.0, 5.0, 10.0)
    assert opts.expansion_type == ExpansionType.SUPERSONIC_AREA_RATIO
    assert opts.chemistry == GasChemistry.EQUILIBRIUM
    assert opts.expansion_ratios == [3.0, 5.0, 10.0]


def test_convenience_combustor():
    from goddard import infinite_area_combustor, CombustorType

    opts = infinite_area_combustor([1e6])
    assert opts.type == CombustorType.INFINITE_AREA
    assert opts.pressures == [1e6]
