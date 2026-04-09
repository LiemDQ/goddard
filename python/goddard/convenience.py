"""Convenience factory functions for building Goddard problems."""

import math

from goddard._core import (
    CombustorType,
    CombustorOptions,
    Gas,
    GasChemistry,
    NozzleOptions,
    ExpansionType,
    create_solution,
    MocFlowKind,
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
    opts.chemistry = GasChemistry.EQUILIBRIUM
    opts.expansion_type = ExpansionType.SUPERSONIC_AREA_RATIO
    opts.expansion_ratios = list(ratios)
    return opts


def subsonic_ratio(*ratios):
    """Create NozzleOptions for subsonic area ratio expansion."""
    opts = NozzleOptions()
    opts.chemistry = GasChemistry.EQUILIBRIUM
    opts.expansion_type = ExpansionType.SUBSONIC_AREA_RATIO
    opts.expansion_ratios = list(ratios)
    return opts


def pressure_ratio(*ratios):
    """Create NozzleOptions for pressure ratio expansion."""
    opts = NozzleOptions()
    opts.chemistry = GasChemistry.EQUILIBRIUM
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
        *ratios: Expansion ratios. If empty, returns bare GasChemistry.

    Returns:
        NozzleOptions if ratios provided, GasChemistry.EQUILIBRIUM otherwise.
    """
    if not ratios:
        return GasChemistry.EQUILIBRIUM
    opts = NozzleOptions()
    opts.chemistry = GasChemistry.EQUILIBRIUM
    opts.expansion_type = ExpansionType.SUPERSONIC_AREA_RATIO
    opts.expansion_ratios = list(ratios)
    return opts


def frozen_nozzle(*ratios, frozen_NFZ=1):
    """Create NozzleOptions for frozen chemistry nozzle with supersonic
    area ratio expansion.

    Args:
        *ratios: Expansion ratios. If empty, returns bare GasChemistry.
        frozen_NFZ: Frozen flow station number.

    Returns:
        NozzleOptions if ratios provided, GasChemistry.FROZEN otherwise.
    """
    if not ratios:
        return GasChemistry.FROZEN
    opts = NozzleOptions()
    opts.chemistry = GasChemistry.FROZEN
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
        chemistry: GasChemistry (defaults to PERFECT_GAS).
        throat_radius: Throat radius for throat geometry.
        solution: SolutionHandle for chemistry-based calculations. When
            provided, the chemistry constructor is used.

    Returns:
        MocResult from the solver.
    """
    if flow_type is None:
        flow_type = MocFlowKind.PLANAR
    if chemistry is None:
        chemistry = GasChemistry.PERFECT_GAS

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
        chemistry: GasChemistry (defaults to PERFECT_GAS).
        throat_radius: Throat radius for throat geometry.
        solution: SolutionHandle for chemistry-based calculations. When
            provided, the chemistry constructor is used.

    Returns:
        MocResult from the solver.
    """
    if flow_type is None:
        flow_type = MocFlowKind.PLANAR
    if chemistry is None:
        chemistry = GasChemistry.PERFECT_GAS

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


def from_cantera(ct_solution, chemistry=GasChemistry.FROZEN):
    """Create a Gas from a Python cantera.Solution object.

    Reconstructs the thermodynamic state by extracting the source file,
    species, and state from the Python Cantera object and creating a new
    internal C++ Solution with the same configuration.

    Args:
        ct_solution: A cantera.Solution Python object.
        chemistry: GasChemistry mode for the Gas (default: FROZEN).

    Returns:
        Gas object with state synchronized from the input Cantera solution.
    """
    source = ct_solution.source
    species_names = set(ct_solution.species_names)
    name = ct_solution.name

    gas = Gas(source, phase_name=name, species=species_names,
              chemistry=chemistry)

    T = ct_solution.T
    P = ct_solution.P
    X = ct_solution.X
    comp_str = ", ".join(
        f"{sp}:{x}" for sp, x in zip(ct_solution.species_names, X) if x > 0
    )
    gas.set_state_TPX(T, P, comp_str)

    return gas


def gas_from_yaml(yaml_file, phase_name="", species=None,
                  chemistry=GasChemistry.FROZEN):
    """Create a Gas from a YAML thermodynamic data file.

    Args:
        yaml_file: Path to a Cantera YAML thermodynamic data file.
        phase_name: Name of the phase in the YAML file (default: "").
        species: Set of species names to include (default: all).
        chemistry: GasChemistry mode (default: FROZEN).

    Returns:
        Gas object.
    """
    return Gas(yaml_file, phase_name=phase_name,
               species=species or set(), chemistry=chemistry)
