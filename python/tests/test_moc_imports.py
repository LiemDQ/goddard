"""Smoke tests for MoC bindings — no solver calls, just imports and construction."""
import pytest


def test_import_moc_enums():
    from goddard import MocFlowKind, GasChemistry, MocMode

    assert MocFlowKind.PLANAR is not None
    assert MocFlowKind.AXISYMMETRIC is not None

    assert GasChemistry.PERFECT_GAS is not None
    assert GasChemistry.FROZEN is not None
    assert GasChemistry.EQUILIBRIUM is not None

    assert MocMode.DESIGN_MIN_LENGTH is not None
    assert MocMode.DESIGN_RAO is not None
    assert MocMode.ANALYSIS is not None


def test_import_moc_structs():
    from goddard import (
        ThroatGeometry, NozzleProfile, MocOptions,
        CharacteristicPoint, CharacteristicNet, ExitPlane,
        MocResult, ThrustCoefficient,
    )


def test_import_moc_class():
    from goddard import MocNozzle, compute_thrust_coefficient


def test_import_convenience_moc():
    from goddard import moc_design, moc_analysis


def test_throat_geometry_construction():
    from goddard import ThroatGeometry

    # Default construction
    g = ThroatGeometry()
    assert isinstance(g.throat_radius, float)

    # Kwargs construction
    g2 = ThroatGeometry(throat_radius=2.0,
                        upstream_wall_curvature_radius=3.0,
                        downstream_wall_curvature_radius=0.5)
    assert g2.throat_radius == 2.0
    assert g2.upstream_wall_curvature_radius == 3.0
    assert g2.downstream_wall_curvature_radius == 0.5

    # Field write
    g2.throat_radius = 5.0
    assert g2.throat_radius == 5.0


def test_nozzle_profile_construction():
    from goddard import NozzleProfile

    p = NozzleProfile()
    assert len(p) == 0


def test_nozzle_profile_push_back():
    from goddard import NozzleProfile

    p = NozzleProfile()
    p.push_back(0.0, 1.0)
    p.push_back(1.0, 1.5)
    p.push_back(2.0, 2.0)

    assert len(p) == 3
    assert p.x == [0.0, 1.0, 2.0]
    assert p.y == [1.0, 1.5, 2.0]


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
    assert opts.theta_schedule == []


def test_moc_options_kwargs():
    from goddard import MocOptions, MocFlowKind, GasChemistry, MocMode

    opts = MocOptions(
        flow_type=MocFlowKind.AXISYMMETRIC,
        chemistry=GasChemistry.FROZEN,
        mode=MocMode.ANALYSIS,
        num_characteristics=20,
        gamma=1.3,
        theta_max=0.3,
    )
    assert opts.flow_type == MocFlowKind.AXISYMMETRIC
    assert opts.chemistry == GasChemistry.FROZEN
    assert opts.mode == MocMode.ANALYSIS
    assert opts.num_characteristics == 20
    assert opts.gamma == pytest.approx(1.3)
    assert opts.theta_max == pytest.approx(0.3)


def test_moc_options_mutability():
    from goddard import MocOptions, MocFlowKind

    opts = MocOptions()
    opts.flow_type = MocFlowKind.AXISYMMETRIC
    assert opts.flow_type == MocFlowKind.AXISYMMETRIC

    opts.geometry.throat_radius = 2.5
    assert opts.geometry.throat_radius == pytest.approx(2.5)
