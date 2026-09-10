"""Integration tests for MoC bindings — runs solver, checks structure not accuracy."""
import math

import numpy as np
import pytest


@pytest.fixture
def planar_result():
    """A basic planar perfect-gas design solve (theta_max=15 deg, N=10)."""
    from goddard import MocOptions, MocFlowKind, GasChemistry, MocMode, MocNozzle
    opts = MocOptions()
    opts.flow_type = MocFlowKind.PLANAR
    opts.chemistry = GasChemistry.PERFECT_GAS
    opts.mode = MocMode.DESIGN_MIN_LENGTH
    opts.num_characteristics = 10
    opts.gamma = 1.4
    opts.theta_max = math.radians(15.0)
    return MocNozzle(opts).solve()


def test_planar_perfect_gas_solve(planar_result):
    assert planar_result.converged
    assert planar_result.exit_mach > 1.0
    assert planar_result.nozzle_length > 0.0
    assert planar_result.area_ratio > 1.0
    assert planar_result.messages == []


def test_result_types(planar_result):
    from goddard import CharacteristicNet, NozzleProfile, ExitPlane
    assert isinstance(planar_result.converged, bool)
    assert isinstance(planar_result.exit_mach, float)
    assert isinstance(planar_result.nozzle_length, float)
    assert isinstance(planar_result.area_ratio, float)
    assert isinstance(planar_result.net, CharacteristicNet)
    assert isinstance(planar_result.profile, NozzleProfile)
    assert isinstance(planar_result.exit_plane, ExitPlane)


def test_characteristic_net_structure(planar_result):
    from goddard import CharacteristicFamily

    net = planar_result.net
    assert len(net) > 0
    assert bool(net) is True
    assert not net.empty()
    # A finished solve has marched every chain to termination.
    assert not net.has_active_chains()

    # The topology lives in the chains: each is an index list into net.points.
    assert len(net.chains) > 0
    assert len(net.chains) == len(net.chain_metadata)
    assert all(chain.dtype == np.int64 for chain in net.chains)
    assert max(int(i) for chain in net.chains for i in chain) < len(net.points)

    families = {meta.family for meta in net.chain_metadata}
    assert CharacteristicFamily.PLUS in families
    assert CharacteristicFamily.MINUS in families

    # membership maps point index -> the chains through it
    assert len(net.membership) == len(net.points)

    assert net.leading_point(0) is not None


def test_characteristic_net_columnar_arrays(planar_result):
    net = planar_result.net
    n = len(net.points)
    for name in ("x", "y", "theta", "nu", "mach", "mu",
                 "pressure", "temperature", "gamma_s", "V"):
        column = getattr(net, name)
        assert isinstance(column, np.ndarray)
        assert column.dtype == np.float64
        assert column.shape == (n,)


def test_characteristic_net_boundaries(planar_result):
    net = planar_result.net
    assert len(net.wall_points()) > 0
    assert len(net.axis_points()) > 0
    # outflow_points may legitimately be empty for a min-length design.
    assert isinstance(net.outflow_points(), list)
    assert net.wall_point_indices.dtype == np.int64
    assert net.axis_point_indices.dtype == np.int64
    assert isinstance(net.wall_x, np.ndarray)


def test_result_diagnostic_fields_readable(planar_result):
    from goddard import MocCrossings, MocFailure, MocErrorCode

    assert isinstance(planar_result.failure, MocFailure)
    assert planar_result.failure.code == MocErrorCode.NONE
    assert isinstance(planar_result.crossings, MocCrossings)
    assert planar_result.crossings.count == 0
    assert isinstance(planar_result.exit_coverage, float)
    assert isinstance(planar_result.reached_exit_plane, bool)
    assert planar_result.inserted_characteristics >= 0
    assert planar_result.retired_characteristics >= 0

    assert len(planar_result.pass_diagnostics) > 0
    diag = planar_result.pass_diagnostics[0]
    # `pass` is a Python keyword, so the field is exposed as pass_index.
    assert isinstance(diag.pass_index, int)
    assert isinstance(diag.max_spacing, float)


def test_find_like_characteristic_crossings(planar_result):
    from goddard import find_like_characteristic_crossings

    crossings = find_like_characteristic_crossings(planar_result.net)
    assert crossings.count == 0


def test_exit_plane_vectors_same_length(planar_result):
    ep = planar_result.exit_plane
    n = len(ep.y)
    assert n > 0
    assert len(ep.mach) == n
    assert len(ep.theta) == n
    assert len(ep.pressure) == n
    assert len(ep.temperature) == n
    assert len(ep.gamma_s) == n
    assert len(ep.velocity) == n
    # All Mach numbers should be supersonic
    assert all(m > 1.0 for m in ep.mach)


def test_wall_points_populated(planar_result):
    net = planar_result.net
    assert len(net.wall_x) > 0
    assert len(net.wall_y) > 0
    # wall_points() is a method on the net, not a field.
    assert len(net.wall_points()) > 0
    # wall_x should be monotonically increasing
    assert np.all(np.diff(net.wall_x) > 0)


def test_profile_matches_wall(planar_result):
    profile = planar_result.profile
    net = planar_result.net
    assert len(profile.x) == len(net.wall_x)
    # Profile starts at or near x=0 (throat)
    assert profile.x[0] == pytest.approx(0.0, abs=1e-6)


def test_nozzle_profile_csv_round_trip(tmp_path):
    from goddard import NozzleProfile
    p = NozzleProfile()
    coords = [(0.0, 1.0), (0.5, 1.2), (1.0, 1.5), (1.5, 1.8), (2.0, 2.0)]
    for x, y in coords:
        p.push_back(x, y)

    csv_path = str(tmp_path / "profile.csv")
    p.save_csv(csv_path)

    p2 = NozzleProfile.load_csv(csv_path)
    assert len(p2) == len(p)
    for i, (x, y) in enumerate(coords):
        xi, yi = p2.at(i)
        assert xi == pytest.approx(x, rel=1e-9)
        assert yi == pytest.approx(y, rel=1e-9)


def test_nozzle_profile_methods(planar_result):
    profile = planar_result.profile
    # Query at a mid-point x
    mid_x = profile.x[len(profile.x) // 2]
    slope = profile.slope_at(mid_x)
    theta = profile.theta_at(mid_x)
    max_t = profile.max_theta()

    assert isinstance(slope, float)
    assert isinstance(theta, float)
    assert isinstance(max_t, float)
    # max_theta should be positive for a diverging nozzle
    assert max_t > 0.0


def test_analysis_round_trip(planar_result):
    """Design a nozzle, then analyse the resulting profile — exit Mach should match."""
    from goddard import MocOptions, MocFlowKind, GasChemistry, MocMode, MocNozzle

    design_mach = planar_result.exit_mach
    profile = planar_result.profile

    opts = MocOptions()
    opts.flow_type = MocFlowKind.PLANAR
    opts.chemistry = GasChemistry.PERFECT_GAS
    opts.mode = MocMode.ANALYSIS
    opts.num_characteristics = 10
    opts.gamma = 1.4
    opts.nozzle_profile = profile

    result = MocNozzle(opts).solve()
    assert result.converged
    # Analysis should reproduce design Mach within 5%
    assert result.exit_mach == pytest.approx(design_mach, rel=0.05)


def test_thrust_coefficient(planar_result):
    from goddard import compute_thrust_coefficient, MocFlowKind

    cf = compute_thrust_coefficient(planar_result, MocFlowKind.PLANAR)
    assert cf.Cf_vacuum > 0.0
    # momentum + pressure ≈ Cf_vacuum (vacuum means zero ambient)
    assert cf.momentum_thrust + cf.pressure_thrust == pytest.approx(cf.Cf_vacuum, rel=1e-6)


def test_thrust_coefficient_ambient(planar_result):
    from goddard import compute_thrust_coefficient, MocFlowKind

    cf_vac = compute_thrust_coefficient(planar_result, MocFlowKind.PLANAR, 0.0)
    cf_amb = compute_thrust_coefficient(planar_result, MocFlowKind.PLANAR, 0.1)
    # Vacuum Cf should be greater than Cf at positive back pressure
    assert cf_vac.Cf_vacuum > cf_amb.Cf


def test_convenience_moc_design():
    from goddard import moc_design

    result = moc_design(15.0)
    assert result.converged
    assert result.exit_mach > 1.0


def test_convenience_moc_analysis():
    from goddard import moc_design, moc_analysis

    design = moc_design(15.0)
    assert design.converged

    result = moc_analysis(design.profile)
    assert result.converged


def test_convenience_analysis_from_csv(tmp_path):
    from goddard import moc_design, moc_analysis

    design = moc_design(15.0)
    csv_path = str(tmp_path / "wall.csv")
    design.profile.save_csv(csv_path)

    result = moc_analysis(csv_path)
    assert result.converged


def test_convenience_axisymmetric():
    from goddard import moc_design, MocFlowKind

    result = moc_design(15.0, flow_type=MocFlowKind.AXISYMMETRIC)
    assert result.converged
    assert result.exit_mach > 1.0


def test_convenience_analysis_from_pathlib(tmp_path):
    from goddard import moc_design, moc_analysis

    design = moc_design(15.0)
    csv_path = tmp_path / "wall.csv"
    design.profile.save_csv(str(csv_path))

    # Only path handling is under test: a pathlib.Path must be accepted where a
    # str is, and give the same answer. Whether the solve converges is a solver
    # question, not a bindings one.
    from_path = moc_analysis(csv_path)
    from_str = moc_analysis(str(csv_path))
    assert from_path.exit_mach == pytest.approx(from_str.exit_mach)
    assert len(from_path.net) == len(from_str.net)


# ---------------------------------------------------------------------------
# Contour generators — smoke level: shape and basic geometry only.
# ---------------------------------------------------------------------------

AREA_RATIO = 8.0


@pytest.mark.parametrize("factory,args", [
    ("conical_nozzle", (AREA_RATIO,)),
    ("rao_nozzle", (AREA_RATIO,)),
    ("bezier_nozzle", (AREA_RATIO, 30.0, 8.0)),
])
def test_generated_contours(factory, args):
    import goddard

    profile = getattr(goddard, factory)(*args)
    assert len(profile) > 0
    assert np.all(np.diff(profile.x) > 0), "x must be strictly increasing"
    # Every generator is defined to reach the requested area ratio.
    assert profile.y[-1] == pytest.approx(math.sqrt(AREA_RATIO), rel=1e-9)
    assert profile.length() > 0.0
    assert profile.max_theta() > 0.0


def test_throat_expansion_curve():
    from goddard import NozzleProfile

    arc = NozzleProfile.throat_expansion_curve(
        theta_n_deg=30.0, r_expansion_curve=0.382, r_throat=1.0)
    assert len(arc) > 0
    assert np.all(np.diff(arc.x) > 0)
    assert arc.y[0] == pytest.approx(1.0)


def test_generator_rejects_bad_input():
    from goddard import NozzleProfile

    with pytest.raises(ValueError):
        NozzleProfile.conical(area_ratio=0.5, r_expansion_curve=0.382)
    with pytest.raises(ValueError):
        NozzleProfile.rao_top(area_ratio=AREA_RATIO, length_frac=0.1)


def test_generated_contour_is_analysable():
    from goddard import conical_nozzle, moc_analysis

    # A generated contour can be fed straight to the solver. Whether it converges
    # is a solver question, not a bindings one; only reachability is asserted here.
    result = moc_analysis(conical_nozzle(AREA_RATIO))
    assert result is not None
    assert len(result.net) > 0


def test_conical_axisymmetric_analysis_converges():
    """D.md item 7 (Package D validation suite): a conical AR=4 axisymmetric analysis
    through goddard.moc_analysis converges, exit_plane.y increases strictly from the axis
    to the exit radius, and the inverse-march front record (net.fronts) is non-empty.

    flow_type=AXISYMMETRIC + mode=ANALYSIS (moc_analysis's default) resolves
    MocMarchScheme.AUTO to INVERSE (see MocMarchScheme in moc.hpp) -- the same default
    path MocDefaultOptionsConvergence.ConicalDefaultThroatConvergesAcrossN exercises in
    C++ (test_moc_convergence.cpp) -- so this is a Python-side check that the same fix
    is reachable through the bindings, not a new numerical claim.
    """
    from goddard import conical_nozzle, moc_analysis, MocFlowKind

    profile = conical_nozzle(4.0)
    result = moc_analysis(profile, flow_type=MocFlowKind.AXISYMMETRIC, num_characteristics=15)

    assert result.converged, result.failure.message

    y = np.asarray(result.exit_plane.y)
    assert len(y) > 1
    assert y[0] == pytest.approx(0.0, abs=1e-6)
    assert np.all(np.diff(y) > 0), "exit_plane.y must increase strictly from the axis outward"
    # generate_conical_nozzle sets the exit radius to sqrt(area_ratio); the inverse march's
    # exit plane lands exactly on the requested station (exit_coverage measured at 1.0 for
    # this configuration in the C++ sweep, tools/moc_sweep.cpp), so 1% covers roundoff.
    assert y[-1] == pytest.approx(math.sqrt(4.0), rel=1e-2)

    # TODO(Package E): CharacteristicNet.fronts is a C++ field (Package B, the inverse
    # march's per-pass front record) not yet exposed in python/src/bind_moc.cpp as of this
    # test. Skip cleanly, keyed on the attribute's absence, so this starts enforcing the
    # instant the binding lands rather than silently staying green forever.
    if not hasattr(result.net, "fronts"):
        pytest.skip("CharacteristicNet.fronts not yet bound (Package E); "
                     "remove this skip once bind_moc.cpp exposes net.fronts")
    assert len(result.net.fronts) > 0


def test_rao_design_mode_runs():
    from goddard import moc_rao_design

    # DESIGN_RAO needs geometry.expansion_ratio / length_fraction, so this whole
    # mode was unreachable from Python before. Convergence is not asserted.
    result = moc_rao_design(AREA_RATIO)
    assert result is not None
    assert result.failure is not None
    assert len(result.pass_diagnostics) > 0


def test_pass_diagnostics_table():
    from goddard import moc_design, pass_diagnostics_table

    table = pass_diagnostics_table(moc_design(15.0))
    assert "pass_index" in table
    assert "max_spacing" in table
    lengths = {len(column) for column in table.values()}
    assert len(lengths) == 1


# ---------------------------------------------------------------------------
# Plotting — smoke level: the helpers run and hand back Axes.
# ---------------------------------------------------------------------------

def test_plotting_helpers(planar_result):
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib.axes import Axes

    from goddard import plotting as gp

    assert isinstance(gp.plot_characteristic_net(planar_result), Axes)
    assert isinstance(gp.plot_field(planar_result, field="mach"), Axes)
    assert isinstance(gp.plot_profile(planar_result.profile), Axes)
    assert isinstance(gp.plot_exit_plane(planar_result, field="mach"), Axes)

    axes = gp.plot_front_diagnostics(planar_result)
    assert len(axes) == 3
    assert all(isinstance(ax, Axes) for ax in axes)


def test_mesh_node_mask_excludes_bootstrap_points(planar_result):
    pytest.importorskip("matplotlib")
    from goddard.plotting import mesh_node_mask

    net = planar_result.net
    mask = mesh_node_mask(net)
    assert mask.shape == (len(net.points),)
    # The seeded throat-lip wall point belongs to no chain and carries placeholder
    # zeros; masking it is what keeps the Mach contour from starting at 0.
    assert not mask.all()
    assert np.asarray(net.mach)[mask].min() >= 1.0
