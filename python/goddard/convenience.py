"""Convenience factory functions for building Goddard problems."""

import math
import os

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


def conical_nozzle(area_ratio, *, r_expansion_curve=0.382, r_throat=1.0,
                   angle_deg=15.0, n_points=50):
    """Generate a conical nozzle contour.

    Args:
        area_ratio: Exit area divided by throat area.
        r_expansion_curve: Radius of curvature of the throat expansion arc, in the
            same length units as ``r_throat``. Zero gives a sharp throat corner.
        r_throat: Throat radius, in length units.
        angle_deg: Conical expansion half-angle, in degrees.
        n_points: Number of points in the contour.

    Returns:
        NozzleProfile.
    """
    return NozzleProfile.conical(
        area_ratio=area_ratio, r_expansion_curve=r_expansion_curve,
        r_throat=r_throat, angle_deg=angle_deg, n_points=n_points)


def rao_nozzle(area_ratio, *, r_throat=1.0, length_frac=0.8, n_points=50):
    """Generate a thrust-optimized parabolic (Rao TOP) nozzle contour.

    Args:
        area_ratio: Exit area divided by throat area. Must be at least 3.55, the
            lower bound of Rao's tabulated data.
        r_throat: Throat radius, in length units.
        length_frac: Length as a fraction of a comparable 15-degree conical nozzle.
            Must be between 0.6 and 1.0.
        n_points: Number of points in the contour.

    Returns:
        NozzleProfile.
    """
    return NozzleProfile.rao_top(
        area_ratio=area_ratio, r_throat=r_throat,
        length_frac=length_frac, n_points=n_points)


def bezier_nozzle(area_ratio, theta_n_deg, theta_e_deg, *,
                  r_expansion_curve=0.382, r_throat=1.0, length_frac=0.8,
                  n_points=50):
    """Generate a parabolic nozzle contour parametrized by Bezier curves.

    Args:
        area_ratio: Exit area divided by throat area.
        theta_n_deg: Maximum wall angle from the centerline, in degrees.
        theta_e_deg: Wall angle at the nozzle exit, in degrees.
        r_expansion_curve: Radius of curvature of the expansion region, as a
            fraction of the throat radius.
        r_throat: Throat radius, in length units.
        length_frac: Length as a fraction of a comparable 15-degree conical nozzle.
        n_points: Number of points in the contour.

    Returns:
        NozzleProfile.
    """
    return NozzleProfile.bezier(
        area_ratio=area_ratio, theta_n_deg=theta_n_deg, theta_e_deg=theta_e_deg,
        r_expansion_curve=r_expansion_curve, r_throat=r_throat,
        length_frac=length_frac, n_points=n_points)


def _apply_mesh_control(opts, max_front_spacing_factor, min_front_spacing_factor,
                        front_spacing_growth, max_cell_aspect_ratio,
                        max_front_points):
    """Set the marching-front mesh controls that were explicitly requested."""
    if max_front_spacing_factor is not None:
        opts.max_front_spacing_factor = max_front_spacing_factor
    if min_front_spacing_factor is not None:
        opts.min_front_spacing_factor = min_front_spacing_factor
    if front_spacing_growth is not None:
        opts.front_spacing_growth = front_spacing_growth
    if max_cell_aspect_ratio is not None:
        opts.max_cell_aspect_ratio = max_cell_aspect_ratio
    if max_front_points is not None:
        opts.max_front_points = max_front_points


def _build_nozzle(opts, gas, solution):
    """Pick the MocNozzle constructor matching the supplied gas description."""
    if gas is not None:
        return MocNozzle(gas, opts)
    if solution is not None:
        return MocNozzle(solution, opts)
    return MocNozzle(opts)


def moc_design(theta_max_deg, *, num_characteristics=10, gamma=1.4,
               flow_type=None, chemistry=None, throat_radius=1.0,
               exit_mach=None, geometry=None, log_level=None,
               max_front_spacing_factor=None, min_front_spacing_factor=None,
               front_spacing_growth=None, max_cell_aspect_ratio=None,
               max_front_points=None, gas=None, solution=None):
    """Design a minimum-length nozzle using the Method of Characteristics.

    Args:
        theta_max_deg: Maximum wall angle in degrees.
        num_characteristics: Number of C+ characteristics from expansion fan.
        gamma: Ratio of specific heats (used for PERFECT_GAS chemistry only).
        flow_type: MocFlowKind (defaults to PLANAR).
        chemistry: GasChemistry (defaults to PERFECT_GAS).
        throat_radius: Throat radius, used when ``geometry`` is not given.
        exit_mach: Target exit Mach number, if the design is Mach-driven.
        geometry: A NozzleGeometry, overriding ``throat_radius``.
        log_level: MocLogLevel; pass ``MocLogLevel.DEBUG`` for a kernel trace.
        max_front_spacing_factor: Upper bound on marching-front point spacing.
        min_front_spacing_factor: Lower bound on marching-front point spacing.
        front_spacing_growth: How much target spacing grows with nozzle radius.
        max_cell_aspect_ratio: Largest tolerated mesh-cell side ratio.
        max_front_points: Safety cap on marching-front size (0 selects the default).
        gas: A Gas for chemistry-based calculations.
        solution: A SolutionHandle, as an alternative to ``gas``.

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
    if geometry is not None:
        opts.geometry = geometry
    else:
        opts.geometry.throat_radius = throat_radius
    if exit_mach is not None:
        opts.exit_mach = exit_mach
    if log_level is not None:
        opts.log_level = log_level
    _apply_mesh_control(opts, max_front_spacing_factor, min_front_spacing_factor,
                        front_spacing_growth, max_cell_aspect_ratio, max_front_points)

    return _build_nozzle(opts, gas, solution).solve()


def moc_rao_design(expansion_ratio, *, length_frac=0.8, num_characteristics=10,
                   gamma=1.4, flow_type=None, chemistry=None, throat_radius=1.0,
                   geometry=None, log_level=None,
                   max_front_spacing_factor=None, min_front_spacing_factor=None,
                   front_spacing_growth=None, max_cell_aspect_ratio=None,
                   max_front_points=None, gas=None, solution=None):
    """Design a Rao thrust-optimized nozzle using the Method of Characteristics.

    The solver generates the Rao contour from ``expansion_ratio`` and
    ``length_frac``, then marches the flow field through it.

    Args:
        expansion_ratio: Exit area divided by throat area.
        length_frac: Length as a fraction of a comparable 15-degree conical nozzle.
        num_characteristics: Number of C+ characteristics from expansion fan.
        gamma: Ratio of specific heats (used for PERFECT_GAS chemistry only).
        flow_type: MocFlowKind (defaults to AXISYMMETRIC, the mode Rao contours are for).
        chemistry: GasChemistry (defaults to PERFECT_GAS).
        throat_radius: Throat radius, used when ``geometry`` is not given.
        geometry: A NozzleGeometry, overriding ``throat_radius``, ``expansion_ratio``
            and ``length_frac``.
        log_level: MocLogLevel; pass ``MocLogLevel.DEBUG`` for a kernel trace.
        max_front_spacing_factor: Upper bound on marching-front point spacing.
        min_front_spacing_factor: Lower bound on marching-front point spacing.
        front_spacing_growth: How much target spacing grows with nozzle radius.
        max_cell_aspect_ratio: Largest tolerated mesh-cell side ratio.
        max_front_points: Safety cap on marching-front size (0 selects the default).
        gas: A Gas for chemistry-based calculations.
        solution: A SolutionHandle, as an alternative to ``gas``.

    Returns:
        MocResult from the solver.
    """
    if flow_type is None:
        flow_type = MocFlowKind.AXISYMMETRIC
    if chemistry is None:
        chemistry = GasChemistry.PERFECT_GAS

    opts = MocOptions()
    opts.flow_type = flow_type
    opts.chemistry = chemistry
    opts.mode = MocMode.DESIGN_RAO
    opts.num_characteristics = num_characteristics
    opts.gamma = gamma
    if geometry is not None:
        opts.geometry = geometry
    else:
        opts.geometry.throat_radius = throat_radius
        opts.geometry.expansion_ratio = expansion_ratio
        opts.geometry.length_fraction = length_frac
    if log_level is not None:
        opts.log_level = log_level
    _apply_mesh_control(opts, max_front_spacing_factor, min_front_spacing_factor,
                        front_spacing_growth, max_cell_aspect_ratio, max_front_points)

    return _build_nozzle(opts, gas, solution).solve()


def moc_analysis(profile, *, num_characteristics=10, gamma=1.4,
                 flow_type=None, chemistry=None, throat_radius=1.0,
                 geometry=None, log_level=None,
                 max_front_spacing_factor=None, min_front_spacing_factor=None,
                 front_spacing_growth=None, max_cell_aspect_ratio=None,
                 max_front_points=None, gas=None, solution=None):
    """Analyse an existing nozzle contour using the Method of Characteristics.

    Args:
        profile: NozzleProfile object, or a path (str or os.PathLike) to a CSV file.
        num_characteristics: Number of C+ characteristics from expansion fan.
        gamma: Ratio of specific heats (used for PERFECT_GAS chemistry only).
        flow_type: MocFlowKind (defaults to PLANAR).
        chemistry: GasChemistry (defaults to PERFECT_GAS).
        throat_radius: Throat radius, used when ``geometry`` is not given.
        geometry: A NozzleGeometry, overriding ``throat_radius``.
        log_level: MocLogLevel; pass ``MocLogLevel.DEBUG`` for a kernel trace.
        max_front_spacing_factor: Upper bound on marching-front point spacing.
        min_front_spacing_factor: Lower bound on marching-front point spacing.
        front_spacing_growth: How much target spacing grows with nozzle radius.
        max_cell_aspect_ratio: Largest tolerated mesh-cell side ratio.
        max_front_points: Safety cap on marching-front size (0 selects the default).
        gas: A Gas for chemistry-based calculations.
        solution: A SolutionHandle, as an alternative to ``gas``.

    Returns:
        MocResult from the solver.
    """
    if flow_type is None:
        flow_type = MocFlowKind.PLANAR
    if chemistry is None:
        chemistry = GasChemistry.PERFECT_GAS

    if isinstance(profile, (str, os.PathLike)):
        profile = NozzleProfile.load_csv(os.fspath(profile))

    opts = MocOptions()
    opts.flow_type = flow_type
    opts.chemistry = chemistry
    opts.mode = MocMode.ANALYSIS
    opts.num_characteristics = num_characteristics
    opts.gamma = gamma
    if geometry is not None:
        opts.geometry = geometry
    else:
        opts.geometry.throat_radius = throat_radius
    opts.nozzle_profile = profile
    if log_level is not None:
        opts.log_level = log_level
    _apply_mesh_control(opts, max_front_spacing_factor, min_front_spacing_factor,
                        front_spacing_growth, max_cell_aspect_ratio, max_front_points)

    return _build_nozzle(opts, gas, solution).solve()


def pass_diagnostics_table(result):
    """Columnarize a result's per-pass marching-front diagnostics.

    Turns ``result.pass_diagnostics`` -- a list of MocPassDiagnostics -- into a dict
    of numpy arrays, which is the shape a grid-convergence study or a plot wants.

    Args:
        result: A MocResult.

    Returns:
        dict mapping field name to a numpy array with one entry per kernel pass.
    """
    import numpy as np

    fields = ("pass_index", "front_points", "min_spacing", "max_spacing",
              "mean_spacing", "target_spacing", "max_cell_aspect",
              "min_spacelike_margin", "inserted", "retired")
    diagnostics = result.pass_diagnostics
    return {
        name: np.array([getattr(d, name) for d in diagnostics])
        for name in fields
    }


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
