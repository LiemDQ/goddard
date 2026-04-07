"""Integration tests for MoC bindings — runs solver, checks structure not accuracy."""
import math
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
    net = planar_result.net
    N = 10
    # num_c_plus should equal N (one per fan characteristic)
    assert net.num_c_plus == N
    # Each wavefront has one fewer point than the previous (triangular)
    assert len(net.wavefronts) > 0
    # First wavefront should have N points; last should have 1
    assert len(net.wavefronts[0]) == N
    assert len(net.wavefronts[-1]) == 1


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
    assert len(net.wall_points) > 0
    # wall_x should be monotonically increasing
    assert all(net.wall_x[i] < net.wall_x[i+1] for i in range(len(net.wall_x)-1))


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
