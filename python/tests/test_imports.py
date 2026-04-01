"""Smoke tests to verify all bindings are importable and basic types work."""
import pytest


def test_import_module():
    import goddard
    assert hasattr(goddard, "__version__")


def test_import_enums():
    from goddard import CombustorType, NozzleChemistryType, ExpansionType

    assert CombustorType.INFINITE_AREA is not None
    assert CombustorType.FINITE_MASS_FLUX is not None
    assert CombustorType.FINITE_CONTRACTION_RATIO is not None
    assert CombustorType.NONE is not None

    assert NozzleChemistryType.FROZEN is not None
    assert NozzleChemistryType.EQUILIBRIUM is not None
    assert NozzleChemistryType.KINETIC is not None

    assert ExpansionType.SUPERSONIC_AREA_RATIO is not None
    assert ExpansionType.SUBSONIC_AREA_RATIO is not None
    assert ExpansionType.PRESSURE_RATIO is not None


def test_import_structs():
    from goddard import (
        CombustorOptions, NozzleOptions, ThroatCondition,
        NozzleResult, NozzleResults, RocketCaseParameters,
        ChemicalParameters, ThermoStateInfo, RocketPerformance,
        RocketState, EquilibriumProperties, EquilibriumDerivatives,
    )


def test_import_classes():
    from goddard import (
        RocketProblem, RocketProblemResults, ThermoArray,
        Combustor, NozzleBase, EquilibriumNozzle, FrozenNozzle,
        MixtureRatio, MixtureRatios, SolutionHandle,
    )


def test_import_kinetic_nozzle():
    from goddard import KineticNozzle, KineticNozzleStation, KineticNozzleResults
    assert KineticNozzle is not None
    assert KineticNozzleStation is not None
    assert KineticNozzleResults is not None


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
    from goddard import NozzleOptions, NozzleChemistryType, ExpansionType

    opts = NozzleOptions()
    opts.chemistry = NozzleChemistryType.EQUILIBRIUM
    opts.expansion_type = ExpansionType.SUPERSONIC_AREA_RATIO
    opts.expansion_ratios = [3.0, 5.0, 10.0]
    assert opts.chemistry == NozzleChemistryType.EQUILIBRIUM
    assert opts.expansion_ratios == [3.0, 5.0, 10.0]


def test_convenience_of_ratio():
    from goddard import OF_ratio
    assert OF_ratio(2.0, 3.0, 4.0) == [2.0, 3.0, 4.0]


def test_convenience_supersonic_ratio():
    from goddard import supersonic_ratio, ExpansionType, NozzleChemistryType

    opts = supersonic_ratio(3.0, 5.0, 10.0)
    assert opts.expansion_type == ExpansionType.SUPERSONIC_AREA_RATIO
    assert opts.chemistry == NozzleChemistryType.EQUILIBRIUM
    assert opts.expansion_ratios == [3.0, 5.0, 10.0]


def test_convenience_combustor():
    from goddard import infinite_area_combustor, CombustorType

    opts = infinite_area_combustor([1e6])
    assert opts.type == CombustorType.INFINITE_AREA
    assert opts.pressures == [1e6]


def test_mixture_ratio():
    from goddard import MixtureRatio

    mr = MixtureRatio(3.0, 2.016, 32.0)
    assert mr.OF_ratio == 3.0
    assert mr.M_fuel == 2.016


def test_mixture_ratios():
    import numpy as np
    from goddard import MixtureRatios

    mrs = MixtureRatios(np.array([2.0, 3.0, 4.0]), 2.016, 32.0)
    assert mrs.size() == 3
