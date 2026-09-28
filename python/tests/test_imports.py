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


def test_import_station_type():
    from goddard import StationType

    assert StationType.CHAMBER is not None
    assert StationType.STAGNATION is not None
    assert StationType.COMBUSTION_END is not None
    assert StationType.THROAT is not None
    assert StationType.EXIT is not None


def test_import_structs():
    from goddard import (
        CombustorOptions, NozzleOptions, ThroatCondition,
        NozzleResult, NozzleResults, FiniteAreaChamber, RocketCaseParameters,
        ChemicalParameters, ThermodynamicState, RocketPerformance,
        RocketState, ExpansionProperties, EquilibriumDerivatives,
    )


def test_import_classes():
    from goddard import (
        RocketProblem, RocketProblemResults, ThermoArray,
        Combustor, Nozzle, Gas, MocNozzle, ShockSolver
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
        ShockResult, ReflectedShockResult, ObliqueShockResult, ShockSolver,
        normal_shock, reflected_shock,
        oblique_shock_wave_angle, oblique_shock_deflection_angle,
        oblique_shock_max_deflection, oblique_shock_max_deflection_wave_angle,
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


def test_perfect_gas_reflected_shock():
    from goddard import reflected_shock, normal_shock
    r = reflected_shock(2.0, 1.4)
    assert r.valid
    assert r.incident.static_pressure_ratio == normal_shock(2.0, 1.4).static_pressure_ratio
    # The reflected shock is weaker than the incident one but still a shock
    assert 1.0 < r.reflected.mach_in < 2.0


def test_perfect_gas_max_deflection():
    import math
    from goddard import oblique_shock_max_deflection, oblique_shock_from_deflection
    # NACA 1135 chart 2: theta_max = 22.97 deg at M = 2, gamma = 1.4
    theta_max = oblique_shock_max_deflection(2.0, 1.4)
    assert abs(math.degrees(theta_max) - 22.97) < 0.01
    assert not oblique_shock_from_deflection(2.0, theta_max + 1e-3, 1.4).valid


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


def test_combustor_options_process():
    from goddard import CombustorOptions, CombustionProcess

    opts = CombustorOptions()
    assert opts.process == CombustionProcess.ISOBARIC
    opts.process = CombustionProcess.ISOCHORIC
    assert opts.process == CombustionProcess.ISOCHORIC
    assert CombustorOptions(process=CombustionProcess.ISOCHORIC).process == CombustionProcess.ISOCHORIC


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


def test_import_condensed_api():
    from goddard import (
        EquilibriumOptions, EquilibriumProperty, MixtureRatioType,
        equilibrium_derivatives, equilibrium_properties, frozen_properties,
        reactant_gas, condensed_species,
    )

    assert EquilibriumProperty.ENTHALPY is not None
    assert EquilibriumProperty.ENTROPY is not None
    assert MixtureRatioType.OF_RATIO is not None
    assert EquilibriumOptions().max_steps > 0
    assert condensed_species("f.yaml", ["C(gr)"]) == {
        "condensed_file": "f.yaml", "condensed_species": {"C(gr)"}}
    assert condensed_species("f.yaml") == {"condensed_file": "f.yaml", "all_condensed": True}


def test_import_detonations():
    from goddard import (
        DetonationBranch, DetonationResult, ReflectedDetonationResult, DetonationSolver,
        chapman_jouguet_detonation, detonation, reflected_detonation,
    )
    assert DetonationBranch.OVERDRIVEN is not None
    assert DetonationBranch.UNDERDRIVEN is not None
    assert DetonationSolver is not None


def test_perfect_gas_detonation():
    from goddard import chapman_jouguet_detonation, detonation, DetonationBranch
    # M_CJ = sqrt(H) + sqrt(H + 1), H = (gamma^2 - 1) Q/(2 gamma); Q is chosen so that H = 4.
    gamma = 1.4
    Q = 2.0 * gamma * 4.0 / (gamma * gamma - 1.0)
    cj = chapman_jouguet_detonation(gamma, Q)
    assert cj.valid
    assert abs(cj.mach_in - (2.0 + 5.0**0.5)) < 1e-12
    assert abs(cj.mach_out - 1.0) < 1e-12
    strong = detonation(1.5, gamma, Q, DetonationBranch.OVERDRIVEN)
    weak = detonation(1.5, gamma, Q, DetonationBranch.UNDERDRIVEN)
    assert strong.mach_out < 1.0 < weak.mach_out
    assert not detonation(0.9, gamma, Q).valid
