import pathlib

from goddard._core import (
    # Enums
    CombustorType,
    CombustionProcess,
    MixtureRatioType,
    GasChemistry,
    ExpansionType,
    StationType,
    # MoC enums
    MocFlowKind,
    MocMode,
    MocLogLevel,
    MocErrorCode,
    MocStepLimiter,
    MocStartLine,
    CharacteristicFamily,
    ChainTermination,
    SolverOptions,
    # Structs
    CombustorOptions,
    NozzleOptions,
    ThroatCondition,
    NozzleResult,
    NozzleResults,
    FiniteAreaChamber,
    RocketCaseParameters,
    ChemicalParameters,
    ThermodynamicState,
    PhaseSpecification,
    RocketStation,
    RocketPerformance,
    RocketState,
    ExpansionProperties,
    EquilibriumDerivatives,
    EquilibriumOptions,
    EquilibriumProperty,
    # MoC structs
    NozzleGeometry,
    NozzleProfile,
    MocOptions,
    CharacteristicPoint,
    CharacteristicNet,
    ChainMetadata,
    PointMembership,
    ExitPlane,
    MocResult,
    MocFailure,
    MocPassDiagnostics,
    MocInitDiagnostics,
    MocFrontShear,
    MocCrossings,
    ThrustCoefficient,
    # Classes
    SolutionHandle,
    RocketProblem,
    RocketProblemResults,
    ThermoArray,
    Combustor,
    Nozzle,
    # MoC class + functions
    MocNozzle,
    compute_thrust_coefficient,
    find_like_characteristic_crossings,
    validate_moc_options,
    summarize_front_shear,
    # Errors
    ConvergenceError,
    # Functions
    add_data_directory,
    create_solution,
    equilibrium_derivatives,
    equilibrium_properties,
    frozen_properties,
    # Kinetic nozzle
    KineticNozzle,
    KineticNozzleStation,
    KineticNozzleResults,
    # Gas wrapper
    Gas,
    # Shocks
    ShockResult,
    ReflectedShockResult,
    ObliqueShockResult,
    ShockSolver,
    normal_shock,
    reflected_shock,
    oblique_shock_wave_angle,
    oblique_shock_deflection_angle,
    oblique_shock_max_deflection_wave_angle,
    oblique_shock_max_deflection,
    oblique_shock_from_wave_angle,
    oblique_shock_from_deflection,
    # Detonations
    DetonationBranch,
    DetonationResult,
    ReflectedDetonationResult,
    DetonationSolver,
    chapman_jouguet_detonation,
    detonation,
    reflected_detonation,
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
    gas_from_yaml,
    reactant_gas,
    condensed_species,
    moc_design,
    moc_rao_design,
    moc_analysis,
    conical_nozzle,
    rao_nozzle,
    bezier_nozzle,
    pass_diagnostics_table,
)

from goddard import plotting

# Thermodynamic data files shipped with the package (absent in a source checkout, where the
# library registers the repository's data/ directory itself).
data_dir = pathlib.Path(__file__).parent / "data"
if data_dir.is_dir():
    add_data_directory(str(data_dir))

__version__ = "0.1.0"
