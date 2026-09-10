"""Smoke tests for MoC bindings — no solver calls, just imports and construction."""
import inspect

import numpy as np
import pytest


def test_import_moc_enums():
    from goddard import (
        MocFlowKind, GasChemistry, MocMode, MocLogLevel, MocErrorCode,
        CharacteristicFamily, ChainTermination,
    )

    assert MocFlowKind.PLANAR is not None
    assert MocFlowKind.AXISYMMETRIC is not None

    assert GasChemistry.PERFECT_GAS is not None
    assert GasChemistry.FROZEN is not None
    assert GasChemistry.EQUILIBRIUM is not None

    assert MocMode.DESIGN_MIN_LENGTH is not None
    assert MocMode.DESIGN_RAO is not None
    assert MocMode.DESIGN_CENTERLINE is not None
    assert MocMode.ANALYSIS is not None

    assert MocLogLevel.NORMAL is not None
    assert MocLogLevel.DEBUG is not None

    assert MocErrorCode.NONE is not None
    assert MocErrorCode.INCOMPLETE_MARCH is not None

    assert CharacteristicFamily.PLUS is not None
    assert CharacteristicFamily.MINUS is not None
    assert ChainTermination.MERGED is not None


def test_import_moc_marching_enums():
    from goddard import MocMarchScheme, MocStepLimiter, MocStartLine

    assert MocMarchScheme.AUTO is not None
    assert MocMarchScheme.DIRECT is not None
    assert MocMarchScheme.INVERSE is not None

    assert MocStepLimiter.NONE is not None
    assert MocStepLimiter.CFL is not None
    assert MocStepLimiter.WALL_FOOT is not None
    assert MocStepLimiter.WALL_TURN is not None
    assert MocStepLimiter.EXIT is not None

    assert MocStartLine.AUTO is not None
    assert MocStartLine.KLIEGEL_LEVINE is not None
    assert MocStartLine.CENTERED_FAN is not None


def test_import_moc_structs():
    from goddard import (
        NozzleGeometry, NozzleProfile, MocOptions,
        CharacteristicPoint, CharacteristicNet, ChainMetadata, PointMembership,
        ExitPlane, MocResult, MocFailure, MocPassDiagnostics, MocCrossings,
        ThrustCoefficient,
    )


def test_import_moc_class():
    from goddard import (
        MocNozzle, compute_thrust_coefficient,
        find_like_characteristic_crossings, validate_moc_options,
    )


def test_import_convenience_moc():
    from goddard import (
        moc_design, moc_rao_design, moc_analysis,
        conical_nozzle, rao_nozzle, bezier_nozzle, pass_diagnostics_table,
    )


def test_import_plotting_without_matplotlib_at_import_time():
    """`import goddard` must not pull in matplotlib; plotting imports it lazily."""
    from goddard import plotting

    assert hasattr(plotting, "plot_characteristic_net")
    assert hasattr(plotting, "plot_field")


def test_nozzle_geometry_construction():
    from goddard import NozzleGeometry

    # Default construction
    g = NozzleGeometry()
    assert isinstance(g.throat_radius, float)
    assert g.throat_radius == pytest.approx(1.0)
    assert g.length_fraction == pytest.approx(0.8)
    assert g.expansion_ratio == pytest.approx(5.0)

    # Kwargs construction
    g2 = NozzleGeometry(throat_radius=2.0,
                        upstream_wall_curvature_radius=3.0,
                        downstream_wall_curvature_radius=0.5,
                        length_fraction=0.9,
                        expansion_ratio=12.0)
    assert g2.throat_radius == 2.0
    assert g2.upstream_wall_curvature_radius == 3.0
    assert g2.downstream_wall_curvature_radius == 0.5
    assert g2.length_fraction == pytest.approx(0.9)
    assert g2.expansion_ratio == pytest.approx(12.0)

    # Field write
    g2.throat_radius = 5.0
    assert g2.throat_radius == 5.0


def test_nozzle_profile_construction():
    from goddard import NozzleProfile

    p = NozzleProfile()
    assert len(p) == 0

    # From coordinate sequences
    p2 = NozzleProfile([0.0, 1.0, 2.0], [1.0, 1.5, 2.0])
    assert len(p2) == 3
    np.testing.assert_allclose(p2.x, [0.0, 1.0, 2.0])

    with pytest.raises(ValueError):
        NozzleProfile([0.0, 1.0], [1.0])


def test_nozzle_profile_push_back():
    from goddard import NozzleProfile

    p = NozzleProfile()
    p.push_back(0.0, 1.0)
    p.push_back(1.0, 1.5)
    p.push_back(2.0, 2.0)

    assert len(p) == 3
    # Coordinates come back as numpy arrays.
    assert isinstance(p.x, np.ndarray)
    np.testing.assert_allclose(p.x, [0.0, 1.0, 2.0])
    np.testing.assert_allclose(p.y, [1.0, 1.5, 2.0])


def test_nozzle_profile_indexing_and_iteration():
    from goddard import NozzleProfile

    p = NozzleProfile([0.0, 1.0, 2.0], [1.0, 1.5, 2.0])
    assert p[0] == pytest.approx((0.0, 1.0))
    assert p[-1] == pytest.approx((2.0, 2.0))
    assert len(list(p)) == 3

    with pytest.raises(IndexError):
        p[3]


def test_moc_options_defaults():
    from goddard import MocOptions, MocFlowKind, GasChemistry, MocMode

    opts = MocOptions()
    assert opts.flow_type == MocFlowKind.PLANAR
    assert opts.chemistry == GasChemistry.PERFECT_GAS
    assert opts.mode == MocMode.DESIGN_MIN_LENGTH
    assert opts.num_characteristics == 10
    assert opts.gamma == pytest.approx(1.4)
    assert opts.solver_options.reltol == pytest.approx(1e-5)
    assert opts.solver_options.abstol == pytest.approx(1e-10)
    assert opts.theta_max == pytest.approx(0.0)
    assert opts.exit_mach == pytest.approx(0.0)
    assert len(opts.theta_schedule) == 0


def test_moc_options_mesh_control_defaults():
    from goddard import MocOptions

    opts = MocOptions()
    assert opts.max_front_spacing_factor == pytest.approx(1.5)
    assert opts.min_front_spacing_factor == pytest.approx(0.35)
    assert opts.front_spacing_growth == pytest.approx(1.0)
    assert opts.max_cell_aspect_ratio == pytest.approx(6.0)
    assert opts.max_front_points == 0
    assert opts.initial_line_axial_shift == pytest.approx(0.1)


def test_moc_options_kwargs():
    from goddard import MocOptions, MocFlowKind, GasChemistry, MocMode

    opts = MocOptions(
        flow_type=MocFlowKind.AXISYMMETRIC,
        chemistry=GasChemistry.FROZEN,
        mode=MocMode.ANALYSIS,
        num_characteristics=20,
        gamma=1.3,
        theta_max=0.3,
        max_front_spacing_factor=1.8,
        max_front_points=64,
    )
    assert opts.flow_type == MocFlowKind.AXISYMMETRIC
    assert opts.chemistry == GasChemistry.FROZEN
    assert opts.mode == MocMode.ANALYSIS
    assert opts.num_characteristics == 20
    assert opts.gamma == pytest.approx(1.3)
    assert opts.theta_max == pytest.approx(0.3)
    assert opts.max_front_spacing_factor == pytest.approx(1.8)
    assert opts.max_front_points == 64


def test_moc_options_mutability():
    from goddard import MocOptions, MocFlowKind

    opts = MocOptions()
    opts.flow_type = MocFlowKind.AXISYMMETRIC
    assert opts.flow_type == MocFlowKind.AXISYMMETRIC

    opts.geometry.throat_radius = 2.5
    assert opts.geometry.throat_radius == pytest.approx(2.5)

    opts.theta_schedule = [0.1, 0.2, 0.3]
    np.testing.assert_allclose(opts.theta_schedule, [0.1, 0.2, 0.3])


def test_validate_moc_options_rejects_bad_input():
    from goddard import MocOptions, validate_moc_options

    validate_moc_options(MocOptions())

    with pytest.raises(ValueError):
        validate_moc_options(MocOptions(num_characteristics=-3))


def test_moc_options_marching_defaults():
    from goddard import MocOptions, MocMarchScheme, MocStartLine

    opts = MocOptions()
    assert opts.march_scheme == MocMarchScheme.AUTO
    assert opts.start_line == MocStartLine.AUTO
    assert opts.inverse_cfl == pytest.approx(0.8)
    assert opts.max_wall_turn_per_step == pytest.approx(0.0175)
    assert opts.front_tilt_decay == pytest.approx(0.9)
    assert opts.kl_max_wall_angle_error == pytest.approx(0.25)


def test_moc_options_marching_kwargs():
    from goddard import MocOptions, MocMarchScheme, MocStartLine

    opts = MocOptions(
        march_scheme=MocMarchScheme.INVERSE,
        start_line=MocStartLine.KLIEGEL_LEVINE,
        inverse_cfl=0.5,
        max_wall_turn_per_step=0.01,
        front_tilt_decay=0.7,
        kl_max_wall_angle_error=0.1,
    )
    assert opts.march_scheme == MocMarchScheme.INVERSE
    assert opts.start_line == MocStartLine.KLIEGEL_LEVINE
    assert opts.inverse_cfl == pytest.approx(0.5)
    assert opts.max_wall_turn_per_step == pytest.approx(0.01)
    assert opts.front_tilt_decay == pytest.approx(0.7)
    assert opts.kl_max_wall_angle_error == pytest.approx(0.1)

    # Also settable after construction, like every other MocOptions field.
    opts.march_scheme = MocMarchScheme.DIRECT
    assert opts.march_scheme == MocMarchScheme.DIRECT


def test_characteristic_net_fronts_default_empty():
    from goddard import CharacteristicNet

    net = CharacteristicNet()
    assert net.fronts == []


def test_moc_pass_diagnostics_step_fields_default():
    from goddard import MocPassDiagnostics, MocStepLimiter

    diag = MocPassDiagnostics()
    assert diag.step_dx == pytest.approx(0.0)
    assert diag.step_limiter == MocStepLimiter.NONE


def test_moc_result_min_theta_and_init_diagnostics_defaults():
    from goddard import MocResult, MocStartLine

    result = MocResult()
    assert result.min_theta == pytest.approx(0.0)
    assert result.min_theta_x == pytest.approx(0.0)
    assert result.min_theta_y == pytest.approx(0.0)

    # MocInitDiagnostics is reached through MocResult.init_diagnostics rather than
    # imported directly (it is bound, but not re-exported at package scope).
    diagnostics = result.init_diagnostics
    assert diagnostics.wall_bc_residual == pytest.approx(0.0)
    assert diagnostics.kplus_wall_end == pytest.approx(0.0)
    assert diagnostics.start_line_used == MocStartLine.AUTO


def test_convenience_moc_functions_accept_march_scheme_and_start_line():
    from goddard.convenience import moc_analysis, moc_rao_design

    for func in (moc_analysis, moc_rao_design):
        params = inspect.signature(func).parameters
        assert "march_scheme" in params
        assert params["march_scheme"].default is None
        assert "start_line" in params
        assert params["start_line"].default is None


def test_import_plot_fronts():
    from goddard import plotting

    assert hasattr(plotting, "plot_fronts")
    assert "plot_fronts" in plotting.__all__
