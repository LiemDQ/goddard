from goddard._core import (
    # Enums
    CombustorType,
    NozzleChemistryType,
    ExpansionType,
    # MoC enums
    MocFlowKind,
    MocChemistry,
    MocMode,
    MocInitialization,
    # Structs
    CombustorOptions,
    NozzleOptions,
    ThroatCondition,
    NozzleResult,
    NozzleResults,
    RocketCaseParameters,
    ChemicalParameters,
    ThermoStateInfo,
    RocketPerformance,
    RocketState,
    RocketProblemCaseResult,
    EquilibriumProperties,
    EquilibriumDerivatives,
    # MoC structs
    ThroatGeometry,
    NozzleProfile,
    MocOptions,
    CharacteristicPoint,
    CharacteristicNet,
    ExitPlane,
    MocResult,
    ThrustCoefficient,
    # Classes
    SolutionHandle,
    RocketProblem,
    RocketProblemResults,
    ThermoArray,
    Combustor,
    NozzleBase,
    EquilibriumNozzle,
    FrozenNozzle,
    MixtureRatio,
    MixtureRatios,
    ThermodynamicState,
    # MoC class + function
    MocNozzle,
    compute_thrust_coefficient,
    # Errors
    ConvergenceError,
    # Functions
    create_solution,
    # Kinetic nozzle
    KineticNozzle,
    KineticNozzleStation,
    KineticNozzleResults
)

from goddard.convenience import (
    OF_ratio,
    supersonic_ratio,
    subsonic_ratio,
    pressure_ratio,
    infinite_area_combustor,
    finite_mass_flux_combustor,
    finite_contraction_ratio_combustor,
    equilibrium_nozzle,
    frozen_nozzle,
    from_cantera,
    moc_design,
    moc_analysis,
)

__version__ = "0.1.0"
