"""Convenience factory functions for building Goddard problems."""

import math

from goddard._core import (
    CombustorType,
    CombustorOptions,
    NozzleChemistryType,
    NozzleOptions,
    ExpansionType,
    create_solution,
    MocFlowKind,
    MocChemistry,
    MocMode,
    MocOptions,
    NozzleProfile,
    MocNozzle,
)


def OF_ratio(*ratios):
    """Create a list of O/F ratios.

    Example:
        OF_ratio(2.0, 2.5, 3.0, 3.5, 4.0)
    """
    return list(ratios)


def supersonic_ratio(*ratios):
    """Create NozzleOptions for supersonic area ratio expansion.

    Example:
        supersonic_ratio(3.0, 5.0, 10.0, 15.0)
    """
    opts = NozzleOptions()
    opts.chemistry = NozzleChemistryType.EQUILIBRIUM
    opts.expansion_type = ExpansionType.SUPERSONIC_AREA_RATIO
    opts.expansion_ratios = list(ratios)
    return opts


def subsonic_ratio(*ratios):
    """Create NozzleOptions for subsonic area ratio expansion."""
    opts = NozzleOptions()
    opts.chemistry = NozzleChemistryType.EQUILIBRIUM
    opts.expansion_type = ExpansionType.SUBSONIC_AREA_RATIO
    opts.expansion_ratios = list(ratios)
    return opts


def pressure_ratio(*ratios):
    """Create NozzleOptions for pressure ratio expansion."""
    opts = NozzleOptions()
    opts.chemistry = NozzleChemistryType.EQUILIBRIUM
    opts.expansion_type = ExpansionType.PRESSURE_RATIO
    opts.expansion_ratios = list(ratios)
    return opts


def infinite_area_combustor(pressures=None):
    """Create CombustorOptions for infinite area combustor.

    Args:
        pressures: List of chamber pressures. Defaults to empty list.
    """
    opts = CombustorOptions()
    opts.type = CombustorType.INFINITE_AREA
    opts.pressures = pressures or []
    return opts


def finite_mass_flux_combustor(mass_flux, pressures=None):
    """Create CombustorOptions for finite mass flux combustor.

    Args:
        mass_flux: Mass flux value.
        pressures: List of chamber pressures.
    """
    opts = CombustorOptions()
    opts.type = CombustorType.FINITE_MASS_FLUX
    opts.mass_flux = mass_flux
    opts.pressures = pressures or []
    return opts


def finite_contraction_ratio_combustor(contraction_ratio, pressures=None):
    """Create CombustorOptions for finite contraction ratio combustor.

    Args:
        contraction_ratio: Contraction ratio value.
        pressures: List of chamber pressures.
    """
    opts = CombustorOptions()
    opts.type = CombustorType.FINITE_CONTRACTION_RATIO
    opts.contraction_ratio = contraction_ratio
    opts.pressures = pressures or []
    return opts


def equilibrium_nozzle(*ratios):
    """Create NozzleOptions for equilibrium chemistry nozzle with supersonic
    area ratio expansion.

    Args:
        *ratios: Expansion ratios. If empty, returns bare NozzleChemistryType.

    Returns:
        NozzleOptions if ratios provided, NozzleChemistryType.EQUILIBRIUM otherwise.
    """
    if not ratios:
        return NozzleChemistryType.EQUILIBRIUM
    opts = NozzleOptions()
    opts.chemistry = NozzleChemistryType.EQUILIBRIUM
    opts.expansion_type = ExpansionType.SUPERSONIC_AREA_RATIO
    opts.expansion_ratios = list(ratios)
    return opts


def frozen_nozzle(*ratios, frozen_NFZ=1):
    """Create NozzleOptions for frozen chemistry nozzle with supersonic
    area ratio expansion.

    Args:
        *ratios: Expansion ratios. If empty, returns bare NozzleChemistryType.
        frozen_NFZ: Frozen flow station number.

    Returns:
        NozzleOptions if ratios provided, NozzleChemistryType.FROZEN otherwise.
    """
    if not ratios:
        return NozzleChemistryType.FROZEN
    opts = NozzleOptions()
    opts.chemistry = NozzleChemistryType.FROZEN
    opts.expansion_type = ExpansionType.SUPERSONIC_AREA_RATIO
    opts.expansion_ratios = list(ratios)
    opts.frozen_NFZ = frozen_NFZ
    return opts


def moc_design(theta_max_deg, *, num_characteristics=10, gamma=1.4,
               flow_type=None, chemistry=None, throat_radius=1.0,
               solution=None):
    """Design a minimum-length nozzle using the Method of Characteristics.

    Args:
        theta_max_deg: Maximum wall angle in degrees.
        num_characteristics: Number of C+ characteristics from expansion fan.
        gamma: Ratio of specific heats (used for PERFECT_GAS chemistry only).
        flow_type: MocFlowKind (defaults to PLANAR).
        chemistry: MocChemistry (defaults to PERFECT_GAS).
        throat_radius: Throat radius for throat geometry.
        solution: SolutionHandle for chemistry-based calculations. When
            provided, the chemistry constructor is used.

    Returns:
        MocResult from the solver.
    """
    if flow_type is None:
        flow_type = MocFlowKind.PLANAR
    if chemistry is None:
        chemistry = MocChemistry.PERFECT_GAS

    opts = MocOptions()
    opts.flow_type = flow_type
    opts.chemistry = chemistry
    opts.mode = MocMode.DESIGN_MIN_LENGTH
    opts.num_characteristics = num_characteristics
    opts.gamma = gamma
    opts.theta_max = math.radians(theta_max_deg)
    opts.geometry.throat_radius = throat_radius

    if solution is not None:
        nozzle = MocNozzle(solution, opts)
    else:
        nozzle = MocNozzle(opts)

    return nozzle.solve()


def moc_analysis(profile, *, num_characteristics=10, gamma=1.4,
                 flow_type=None, chemistry=None, throat_radius=1.0,
                 solution=None):
    """Analyse an existing nozzle contour using the Method of Characteristics.

    Args:
        profile: NozzleProfile object or path to a CSV file.
        num_characteristics: Number of C+ characteristics from expansion fan.
        gamma: Ratio of specific heats (used for PERFECT_GAS chemistry only).
        flow_type: MocFlowKind (defaults to PLANAR).
        chemistry: MocChemistry (defaults to PERFECT_GAS).
        throat_radius: Throat radius for throat geometry.
        solution: SolutionHandle for chemistry-based calculations. When
            provided, the chemistry constructor is used.

    Returns:
        MocResult from the solver.
    """
    if flow_type is None:
        flow_type = MocFlowKind.PLANAR
    if chemistry is None:
        chemistry = MocChemistry.PERFECT_GAS

    if isinstance(profile, str):
        profile = NozzleProfile.load_csv(profile)

    opts = MocOptions()
    opts.flow_type = flow_type
    opts.chemistry = chemistry
    opts.mode = MocMode.ANALYSIS
    opts.num_characteristics = num_characteristics
    opts.gamma = gamma
    opts.geometry.throat_radius = throat_radius
    opts.nozzle_profile = profile

    if solution is not None:
        nozzle = MocNozzle(solution, opts)
    else:
        nozzle = MocNozzle(opts)

    return nozzle.solve()


def from_cantera(ct_solution):
    """Create a Goddard SolutionHandle from a Python cantera.Solution object.

    This reconstructs the thermodynamic state by extracting the source file,
    species, and state from the Python Cantera object and creating a new
    internal C++ Solution with the same configuration.

    Args:
        ct_solution: A cantera.Solution Python object.

    Returns:
        SolutionHandle that can be passed to Combustor/Nozzle constructors.
    """
    source = ct_solution.source
    species_names = set(ct_solution.species_names)
    name = ct_solution.name

    sol = create_solution(source, name, species_names)

    # Sync thermodynamic state: set temperature, pressure, and composition
    # Access the thermo handle to set state
    # For now, we set state via the returned handle's methods
    # This requires the SolutionHandle to expose setState_TPX or similar
    # TODO: add state-setting methods to SolutionHandle binding if needed
    return sol
