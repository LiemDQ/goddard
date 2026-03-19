from goddard._core import (
    # Enums
    CombustorType,
    NozzleChemistryType,
    ExpansionType,
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
    # Errors
    ConvergenceError,
    # Functions
    create_solution,
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
)

__version__ = "0.1.0"
