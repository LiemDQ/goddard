// Acceptance tests for the inverse (reference-plane) marching kernel (Package B).
// See instructions/moc_fix/B.md for the full spec these implement.
#include "goddard/moc.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/prandtlmeyer.hpp"
#include "goddard/nozzle.hpp"
#include "cantera/core.h"
#include <cmath>
#include <vector>
#include <string>
#include "gtest/gtest.h"

using namespace Goddard;

static constexpr double DEG = M_PI / 180.0;

namespace {

// Exposes the protected test-only front-override hook (see MocNozzle::m_inverse_front_override)
// so a hand-built front with a known exact solution can be marched without a throat/KL/fan
// construction consistent with it.
class InverseMarchTestNozzle : public MocNozzle {
public:
    using MocNozzle::MocNozzle;
    void set_front_override(std::vector<CharacteristicPoint> front) {
        m_inverse_front_override = std::move(front);
    }
};

// Builds a manufactured front point given (x, y, theta, mach, gamma) directly, bypassing
// ThermodynamicContext (perfect-gas only, matching every test in this file): nu and mu follow
// from the perfect-gas Prandtl-Meyer relations, and pressure/temperature are set to harmless
// placeholders since check_point_validity does not examine them.
CharacteristicPoint make_point(double x, double y, double theta, double mach, double gamma) {
    CharacteristicPoint pt{};
    pt.x = x;
    pt.y = y;
    pt.theta = theta;
    pt.mach = mach;
    pt.nu = prandtl_meyer(mach, gamma);
    pt.mu = mach_to_mu(mach);
    pt.gamma_s = gamma;
    pt.V = mach;
    pt.pressure = 1.0;
    pt.temperature = 1.0;
    pt.K_plus = pt.theta - pt.nu;
    pt.K_minus = pt.theta + pt.nu;
    return pt;
}

// Isentropic 1D area-Mach relation A/A* (quasi-1D mass conservation), matching
// test_moc_convergence.cpp's helper of the same name.
double area_ratio_1d(double mach, double gamma) {
    double t = (2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * mach * mach);
    return std::pow(t, (gamma + 1.0) / (2.0 * (gamma - 1.0))) / mach;
}

// Supersonic Mach from A/A* (bisection; relation is monotone for M > 1), matching
// test_moc_convergence.cpp's helper of the same name.
double mach_from_area_ratio_1d(double area_ratio, double gamma) {
    double lo = 1.0 + 1e-9, hi = 50.0;
    for (int i = 0; i < 200; i++) {
        double mid = 0.5 * (lo + hi);
        if (area_ratio_1d(mid, gamma) > area_ratio) hi = mid;
        else lo = mid;
    }
    return 0.5 * (lo + hi);
}

MocOptions make_inverse_options(MocFlowKind flow, double gamma, int n) {
    MocOptions opts;
    opts.flow_type = flow;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::ANALYSIS;
    opts.gamma = gamma;
    opts.num_characteristics = n;
    opts.geometry.throat_radius = 1.0;
    opts.march_scheme = MocMarchScheme::INVERSE;
    return opts;
}

} // namespace

// ============================================================
// Test 1: uniform flow (theta=0, M=const) is exact in a straight duct for both flow kinds.
// ============================================================

class InverseMarchUniformFlow : public ::testing::TestWithParam<MocFlowKind> {};

TEST_P(InverseMarchUniformFlow, StaysUniform) {
    const MocFlowKind flow = GetParam();
    const double gamma = 1.4;
    const double mach = 2.0;
    const int n = 21;

    MocOptions opts = make_inverse_options(flow, gamma, n);
    opts.nozzle_profile.push_back({1.0, 1.0});
    opts.nozzle_profile.push_back({5.0, 1.0});

    std::vector<CharacteristicPoint> front0(n);
    for (int i = 0; i < n; i++) {
        const double y = static_cast<double>(i) / static_cast<double>(n - 1);
        front0[i] = make_point(1.0, y, 0.0, mach, gamma);
    }

    InverseMarchTestNozzle solver(opts);
    solver.set_front_override(front0);
    MocResult result = solver.solve();

    ASSERT_TRUE(result.converged) << to_string(result.failure.code) << ": " << result.failure.message;
    ASSERT_GE(result.net.fronts.size(), 2u) << "march must have advanced past F_0";

    for (size_t f = 0; f < result.net.fronts.size(); f++) {
        for (size_t idx : result.net.fronts[f]) {
            const CharacteristicPoint& pt = result.net.points[idx];
            EXPECT_NEAR(pt.mach, mach, 1e-10) << "front " << f << " point at y=" << pt.y;
            EXPECT_NEAR(pt.theta, 0.0, 1e-10) << "front " << f << " point at y=" << pt.y;
        }
    }
    EXPECT_NEAR(result.net.fronts.back().size(), static_cast<size_t>(n), 0u);
    EXPECT_NEAR(result.net.wall_x.back(), 5.0, 1e-9);
}

INSTANTIATE_TEST_SUITE_P(PlanarAndAxisymmetric, InverseMarchUniformFlow,
    ::testing::Values(MocFlowKind::PLANAR, MocFlowKind::AXISYMMETRIC),
    [](const ::testing::TestParamInfo<MocFlowKind>& info) {
        return info.param == MocFlowKind::PLANAR ? "Planar" : "Axisymmetric";
    });

// ============================================================
// Test 2: radial source flow from a virtual apex is an exact solution of the full (not just
// 1D) supersonic potential-flow equations, so it can be evaluated pointwise anywhere -- not
// only on the spherical cap the flow is naturally uniform on -- and used as a manufactured
// solution for a genuine grid-convergence check of the inverse kernel's order of accuracy.
// A straight cone through the apex is a streamline of this flow (theta = atan2(y, x - x0) is
// constant along any ray from the apex), so it is also a valid, exact wall.
// ============================================================

namespace {

struct SourceFlow {
    double x0;     // apex x
    double r_star; // critical radius (A/A* = 1 there)
    double gamma;

    double mach_at(double x, double y) const {
        const double r = std::hypot(x - x0, y);
        const double area_ratio = std::pow(r / r_star, 2); // axisymmetric (spherical cap)
        return mach_from_area_ratio_1d(area_ratio, gamma);
    }
    double mach_at_planar(double x, double y) const {
        const double r = std::hypot(x - x0, y);
        const double area_ratio = r / r_star; // planar (cylindrical cap)
        return mach_from_area_ratio_1d(area_ratio, gamma);
    }
    double theta_at(double x, double y) const {
        return std::atan2(y, x - x0);
    }
};

std::vector<CharacteristicPoint> build_source_flow_front(
    const SourceFlow& flow, MocFlowKind flow_kind, double x_station, double y_wall, int n)
{
    std::vector<CharacteristicPoint> front(n);
    for (int i = 0; i < n; i++) {
        const double y = y_wall * static_cast<double>(i) / static_cast<double>(n - 1);
        const double mach = (flow_kind == MocFlowKind::AXISYMMETRIC)
            ? flow.mach_at(x_station, y) : flow.mach_at_planar(x_station, y);
        const double theta = flow.theta_at(x_station, y);
        front[i] = make_point(x_station, y, theta, mach, flow.gamma);
    }
    return front;
}

} // namespace

class InverseMarchSourceFlow : public ::testing::TestWithParam<MocFlowKind> {};

TEST_P(InverseMarchSourceFlow, SecondOrderConvergence) {
    const MocFlowKind flow_kind = GetParam();
    const double gamma = 1.4;
    const double half_angle = 15.0 * DEG;

    // Apex placed so the flow is comfortably supersonic (M ~ 2.7-3) across x in [1, 4] and
    // the cone (a streamline of this exact flow) is a valid wall there.
    SourceFlow flow{-0.3, 0.5, gamma};

    const double x_start = 1.0;
    const double x_end = 4.0;
    const double y_wall_start = (x_start - flow.x0) * std::tan(half_angle);

    NozzleProfile profile;
    profile.push_back({0.9, (0.9 - flow.x0) * std::tan(half_angle)});
    profile.push_back({x_end, (x_end - flow.x0) * std::tan(half_angle)});

    std::vector<int> levels = {11, 21, 41};
    std::vector<double> linf_errors;

    for (int n : levels) {
        MocOptions opts = make_inverse_options(flow_kind, gamma, n);
        opts.nozzle_profile = profile;

        std::vector<CharacteristicPoint> front0 =
            build_source_flow_front(flow, flow_kind, x_start, y_wall_start, n);

        InverseMarchTestNozzle solver(opts);
        solver.set_front_override(front0);
        MocResult result = solver.solve();

        ASSERT_TRUE(result.converged) << "N=" << n << ": " << to_string(result.failure.code)
            << " -- " << result.failure.message;
        ASSERT_FALSE(result.net.fronts.empty());
        EXPECT_NEAR(result.net.wall_x.back(), x_end, 1e-6) << "N=" << n;

        double linf = 0.0;
        for (size_t idx : result.net.fronts.back()) {
            const CharacteristicPoint& pt = result.net.points[idx];
            const double mach_exact = (flow_kind == MocFlowKind::AXISYMMETRIC)
                ? flow.mach_at(pt.x, pt.y) : flow.mach_at_planar(pt.x, pt.y);
            linf = std::max(linf, std::abs(pt.mach - mach_exact));
        }
        linf_errors.push_back(linf);
        RecordProperty("linf_error_N" + std::to_string(n), std::to_string(linf));
    }

    for (size_t i = 0; i + 1 < levels.size(); i++) {
        EXPECT_LT(linf_errors[i + 1], linf_errors[i] / 3.0)
            << "L-inf Mach error must drop by >= 3x per doubling of N (second order): "
            << "N=" << levels[i] << " err=" << linf_errors[i]
            << ", N=" << levels[i + 1] << " err=" << linf_errors[i + 1];
    }
}

INSTANTIATE_TEST_SUITE_P(PlanarAndAxisymmetric, InverseMarchSourceFlow,
    ::testing::Values(MocFlowKind::PLANAR, MocFlowKind::AXISYMMETRIC),
    [](const ::testing::TestParamInfo<MocFlowKind>& info) {
        return info.param == MocFlowKind::PLANAR ? "Planar" : "Axisymmetric";
    });

// ============================================================
// Test 3: planar DESIGN_MIN_LENGTH -> ANALYSIS (fan start line) round trip.
// ============================================================

namespace {

MocResult solve_planar_min_length_design(double gamma, double theta_max, int n) {
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = gamma;
    opts.theta_max = theta_max;
    opts.num_characteristics = n;
    opts.geometry.throat_radius = 1.0;
    MocNozzle solver(opts);
    return solver.solve();
}

} // namespace

TEST(InverseMarch, PlanarDesignRoundTrip) {
    const double gamma = 1.4;
    const double theta_max = 15.0 * DEG;
    const double expected_exit_mach = mach_from_prandtl_meyer(2.0 * theta_max, gamma);

    double err[2] = {0.0, 0.0};
    bool ok[2] = {false, false};
    int levels[2] = {16, 32};
    for (int i = 0; i < 2; i++) {
        MocResult design = solve_planar_min_length_design(gamma, theta_max, levels[i]);
        ASSERT_TRUE(design.converged) << "design N=" << levels[i];

        MocOptions opts = make_inverse_options(MocFlowKind::PLANAR, gamma, levels[i]);
        opts.geometry.downstream_wall_curvature_radius = -1.0; // fan start line
        opts.nozzle_profile = design.profile;

        MocNozzle solver(opts);
        MocResult analysis = solver.solve();
        ok[i] = analysis.converged;
        EXPECT_TRUE(analysis.converged) << "N=" << levels[i] << ": "
            << to_string(analysis.failure.code) << " -- " << analysis.failure.message;
        if (!analysis.converged) continue;

        err[i] = std::abs(analysis.exit_mach - expected_exit_mach);
        RecordProperty("exit_mach_error_N" + std::to_string(levels[i]), std::to_string(err[i]));
    }

    if (ok[1]) {
        EXPECT_LT(err[1], 5e-3) << "N=32 exit Mach must match PM^-1(2*theta_max) within 5e-3 "
            << "(err=" << err[1] << ")";
    }
    if (ok[0] && ok[1]) {
        EXPECT_LT(err[1], 0.6 * err[0])
            << "Round-trip error must shrink with N (N=16 err=" << err[0]
            << ", N=32 err=" << err[1] << ")";
    }
}

// ============================================================
// Test 4: conical nozzle, clean (large-radius) throat, KL start line. Exit-plane mass flow
// vs. the 1-D critical value, mirroring MocNozzle::start_line_mass_flow_error's construction
// (moc_nozzle.cpp) but applied to the exit plane instead of the start line, since that
// helper is protected and not accessible from a test.
// ============================================================

namespace {

double perfect_gas_mass_flux(double mach, double gamma) {
    const double exponent = -(gamma + 1.0) / (2.0 * (gamma - 1.0));
    return mach * std::pow(1.0 + 0.5 * (gamma - 1.0) * mach * mach, exponent);
}

double exit_plane_mass_flow_error(const MocResult& result, double gamma, double r_throat) {
    const ExitPlane& ep = result.exit_plane;
    if (ep.y.size() < 2) return std::numeric_limits<double>::quiet_NaN();
    double mdot = 0.0;
    for (size_t i = 1; i < ep.y.size(); i++) {
        const double dy = ep.y[i] - ep.y[i - 1];
        const double theta_avg = 0.5 * (ep.theta[i] + ep.theta[i - 1]);
        const double y_avg = 0.5 * (ep.y[i] + ep.y[i - 1]);
        const double flux = 0.5 * (perfect_gas_mass_flux(ep.mach[i], gamma)
                                   + perfect_gas_mass_flux(ep.mach[i - 1], gamma));
        // The exit plane is a vertical cut (constant x): dx=0, so the downstream normal flux
        // reduces to cos(theta_avg)*dy (see start_line_mass_flow_error for the general form).
        const double normal_flux = std::cos(theta_avg) * dy;
        mdot += flux * normal_flux * (2.0 * M_PI * y_avg);
    }
    const double mdot_reference = perfect_gas_mass_flux(1.0, gamma) * (M_PI * r_throat * r_throat);
    return (mdot - mdot_reference) / mdot_reference;
}

} // namespace

// Both marching kernels place a weak compression converging on the axis near x = 3.9-4.3
// in this geometry, steepening under refinement (instructions/moc_fix/diagnosis.md A8):
// the internal shock known for conical nozzles with circular-arc throats. The march passes
// through it isentropically, so the exit Mach right at it (AR = 4, exit at x = 4.0) is not a
// convergence metric; AR = 8's exit lies well beyond it and is. Mass conservation across the
// exit plane is the check that holds for both.
TEST(InverseMarch, ConicalConvergesAtCleanThroat) {
    const double gamma = 1.4;
    const int levels[3] = {15, 31, 61};

    for (double ar : {4.0, 8.0}) {
        double exit_mach[3] = {0, 0, 0};
        bool ok[3] = {false, false, false};

        for (int li = 0; li < 3; li++) {
            const int n = levels[li];
            MocOptions opts = make_inverse_options(MocFlowKind::AXISYMMETRIC, gamma, n);
            opts.geometry.downstream_wall_curvature_radius = 2.0; // KL start line, clean throat
            opts.nozzle_profile = NozzleProfile::generate_conical_nozzle(ar, 2.0, 1.0, 15.0, 120);

            MocNozzle solver(opts);
            MocResult result = solver.solve();

            ok[li] = result.converged;
            EXPECT_TRUE(result.converged) << "AR=" << ar << " N=" << n << ": "
                << to_string(result.failure.code) << " -- " << result.failure.message;
            if (!result.converged) continue;
            EXPECT_TRUE(result.reached_exit_plane) << "AR=" << ar << " N=" << n;

            // Wall Mach must not decrease along the contour. The first wall point is the
            // start line's own (series) value; comparisons start at the first solved one.
            // 1e-3 absorbs the facet-to-facet wall-angle interpolation on a 120-point arc.
            std::vector<double> wall_mach;
            for (size_t idx : result.net.wall_point_indices) {
                wall_mach.push_back(result.net.points[idx].mach);
            }
            for (size_t i = 2; i < wall_mach.size(); i++) {
                EXPECT_GE(wall_mach[i], wall_mach[i - 1] - 1e-3)
                    << "AR=" << ar << " N=" << n << ": wall Mach dropped at wall point " << i;
            }

            // The compression on the axis: recorded, and bounded so a folded solution
            // (a runaway to large negative angles) cannot pass as a weak feature.
            const std::string tag = "_AR" + std::to_string(static_cast<int>(ar)) + "_N" + std::to_string(n);
            RecordProperty("min_theta_deg" + tag, std::to_string(result.min_theta / DEG));
            RecordProperty("min_theta_x" + tag, std::to_string(result.min_theta_x));
            EXPECT_GT(result.min_theta, -2.0 * DEG) << "AR=" << ar << " N=" << n
                << ": flow angle dipped to " << result.min_theta / DEG << " deg at x=" << result.min_theta_x;

            const double mdot_err = exit_plane_mass_flow_error(result, gamma, 1.0);
            EXPECT_LT(std::abs(mdot_err), 0.01)
                << "AR=" << ar << " N=" << n << ": exit-plane mass flow must be within 1% of the "
                << "throat value (got " << mdot_err << ")";

            exit_mach[li] = result.exit_mach;
            RecordProperty("exit_mach" + tag, std::to_string(exit_mach[li]));
            RecordProperty("mdot_err" + tag, std::to_string(mdot_err));
            RecordProperty("passes" + tag, std::to_string(result.pass_diagnostics.size()));
        }

        if (ar > 4.0 && ok[0] && ok[1] && ok[2]) {
            // Measured 3.4990 / 3.4989 / 3.4985 at N = 15 / 31 / 61: the increments are already
            // at the 1e-4 level, where a ratio test says nothing, so the assertion is on the
            // spread. A regression that moved the exit Mach by more than this would show.
            EXPECT_LT(std::abs(exit_mach[2] - exit_mach[1]), 2e-3)
                << "AR=" << ar << ": exit Mach at N=61 should agree with N=31 to 2e-3 (M15="
                << exit_mach[0] << ", M31=" << exit_mach[1] << ", M61=" << exit_mach[2] << ")";
            EXPECT_LT(std::abs(exit_mach[1] - exit_mach[0]), 5e-3)
                << "AR=" << ar << ": exit Mach at N=31 should agree with N=15 to 5e-3";
        }
    }
}

// The default throat (r_arc = 0.382), where the KL line is wall-corrected (Package A). The
// compression on the axis is stronger here (about -3 deg at N=31 and -5 deg at N=61 near
// x = 3.4) and steepens with N: a forming shock. The march still reaches the exit plane and
// conserves mass to about 1%, which is what this asserts; the flow downstream of the
// compression is approximate and is not compared across N.
TEST(InverseMarch, ConicalDefaultThroatReachesExit) {
    const double gamma = 1.4;
    for (double ar : {4.0, 8.0}) {
        for (int n : {31, 61}) {
            MocOptions opts = make_inverse_options(MocFlowKind::AXISYMMETRIC, gamma, n);
            opts.geometry.downstream_wall_curvature_radius = 0.382;
            opts.nozzle_profile = NozzleProfile::generate_conical_nozzle(ar, 0.382, 1.0, 15.0, 60);
            MocResult result = MocNozzle(opts).solve();
            EXPECT_TRUE(result.converged) << "AR=" << ar << " N=" << n << ": "
                << to_string(result.failure.code) << " -- " << result.failure.message;
            if (!result.converged) continue;
            EXPECT_TRUE(result.reached_exit_plane) << "AR=" << ar << " N=" << n;
            const double mdot_err = exit_plane_mass_flow_error(result, gamma, 1.0);
            EXPECT_LT(std::abs(mdot_err), 0.02) << "AR=" << ar << " N=" << n << ": mass-flow error " << mdot_err;
            EXPECT_LT(result.min_theta, -1.0 * DEG)
                << "AR=" << ar << " N=" << n << ": the axis compression is expected here; "
                << "if it has gone, the geometry or the physics changed -- update this test";
            EXPECT_GT(result.min_theta, -8.0 * DEG) << "AR=" << ar << " N=" << n;
            const std::string tag = "_AR" + std::to_string(static_cast<int>(ar)) + "_N" + std::to_string(n);
            RecordProperty("min_theta_deg" + tag, std::to_string(result.min_theta / DEG));
            RecordProperty("exit_mach" + tag, std::to_string(result.exit_mach));
            RecordProperty("mdot_err" + tag, std::to_string(mdot_err));
        }
    }
}

// Axisymmetric design -> analysis round trip. The DIRECT minimum-length design carries a
// known +0.08 exit-Mach bias against the 1-D area-Mach relation (it does not conserve mass
// exactly; see MocConvergence.AxiDesign1DConsistencyBounded), so its exit Mach is not the
// reference here. The analysis of the designed contour is judged on what must hold for any
// correct isentropic solution: the exit-plane mass flow matches the throat's, increasingly
// so with N, and the area-averaged exit Mach matches the 1-D value for the contour's
// area ratio.
TEST(InverseMarch, AxiDesignRoundTripConservesMass) {
    const double gamma = 1.4;
    const double theta_max = 12.0 * DEG;
    double mdot_err[3] = {0, 0, 0};
    const int levels[3] = {8, 16, 32};
    for (int i = 0; i < 3; i++) {
        MocOptions d;
        d.flow_type = MocFlowKind::AXISYMMETRIC;
        d.chemistry = GasChemistry::PERFECT_GAS;
        d.mode = MocMode::DESIGN_MIN_LENGTH;
        d.gamma = gamma;
        d.theta_max = theta_max;
        d.num_characteristics = levels[i];
        d.geometry.throat_radius = 1.0;
        MocResult design = MocNozzle(d).solve();
        ASSERT_TRUE(design.converged) << "design N=" << levels[i];

        MocOptions a = make_inverse_options(MocFlowKind::AXISYMMETRIC, gamma, levels[i]);
        a.geometry.downstream_wall_curvature_radius = -1.0; // fan start line
        a.nozzle_profile = design.profile;
        MocResult analysis = MocNozzle(a).solve();
        ASSERT_TRUE(analysis.converged) << "N=" << levels[i] << ": "
            << to_string(analysis.failure.code) << " -- " << analysis.failure.message;

        mdot_err[i] = exit_plane_mass_flow_error(analysis, gamma, 1.0);
        EXPECT_LT(std::abs(mdot_err[i]), 0.02) << "N=" << levels[i] << ": mass-flow error " << mdot_err[i];

        // Area-averaged exit Mach against the 1-D value for the achieved area ratio.
        const ExitPlane& ep = analysis.exit_plane;
        double mach_area = 0.0, area = 0.0;
        for (size_t k = 1; k < ep.y.size(); k++) {
            const double dA = M_PI * (ep.y[k] * ep.y[k] - ep.y[k - 1] * ep.y[k - 1]);
            mach_area += 0.5 * (ep.mach[k] + ep.mach[k - 1]) * dA;
            area += dA;
        }
        const double mach_mean = mach_area / area;
        const double mach_1d = mach_from_area_ratio_1d(analysis.area_ratio, gamma);
        EXPECT_NEAR(mach_mean, mach_1d, 0.04 * mach_1d) << "N=" << levels[i]
            << ": area-mean exit Mach " << mach_mean << " vs 1-D " << mach_1d;
        RecordProperty("exit_mach_axis_N" + std::to_string(levels[i]), std::to_string(analysis.exit_mach));
        RecordProperty("design_exit_mach_N" + std::to_string(levels[i]), std::to_string(design.exit_mach));
        RecordProperty("mdot_err_N" + std::to_string(levels[i]), std::to_string(mdot_err[i]));
    }
    EXPECT_LT(std::abs(mdot_err[2]), std::abs(mdot_err[1]))
        << "mass-flow error must shrink from N=16 (" << mdot_err[1] << ") to N=32 (" << mdot_err[2] << ")";
    EXPECT_LT(std::abs(mdot_err[1]), std::abs(mdot_err[0]));
}

// ============================================================
// Test 5: frozen/equilibrium chemistry smoke test on the AR=4 conical contour, using the
// MocFrozenTest fixture's gas (test_moc_phases.cpp).
// ============================================================

class InverseMarchChemistrySmoke : public ::testing::Test {
protected:
    void SetUp() override {
        gas = Cantera::newSolution("h2o2.yaml", "ohmech");
        gas->thermo()->setState_TPX(3000.0, 3e6, "H2O:0.8, OH:0.1, H2:0.05, O2:0.05");
    }
    std::shared_ptr<Cantera::Solution> gas;
};

TEST_F(InverseMarchChemistrySmoke, FrozenAndEquilibriumSmoke) {
    NozzleProfile profile = NozzleProfile::generate_conical_nozzle(4.0, 2.0, 1.0, 15.0, 120);

    auto solve_with = [&](GasChemistry chemistry) {
        MocOptions opts;
        opts.flow_type = MocFlowKind::AXISYMMETRIC;
        opts.chemistry = chemistry;
        opts.mode = MocMode::ANALYSIS;
        // N=15: N=31/61 fail to converge at r_arc=2.0 with the default front_tilt_decay --
        // see InverseMarch.ConicalConvergesAtCleanThroat and the report ("changes to the
        // algorithm"). N=15 is the largest of the three N tested there that converges.
        opts.num_characteristics = 15;
        opts.geometry.throat_radius = 1.0;
        opts.geometry.downstream_wall_curvature_radius = 2.0;
        opts.march_scheme = MocMarchScheme::INVERSE;
        opts.nozzle_profile = profile;

        gas->thermo()->setState_TPX(3000.0, 3e6, "H2O:0.8, OH:0.1, H2:0.05, O2:0.05");
        MocNozzle nozzle(Gas(gas, chemistry), opts);
        return nozzle.solve();
    };

    MocResult frozen_result = solve_with(GasChemistry::FROZEN);
    MocResult equil_result = solve_with(GasChemistry::EQUILIBRIUM);

    ASSERT_TRUE(frozen_result.converged) << to_string(frozen_result.failure.code)
        << " -- " << frozen_result.failure.message;
    ASSERT_TRUE(equil_result.converged) << to_string(equil_result.failure.code)
        << " -- " << equil_result.failure.message;

    RecordProperty("frozen_exit_mach", std::to_string(frozen_result.exit_mach));
    RecordProperty("equilibrium_exit_mach", std::to_string(equil_result.exit_mach));

    // 1-D reference at the same area ratio, same ordering as MocFrozenTest.AxiFrozenVs1D /
    // AxiEquilibriumVs1D / AxiFrozenVsEquilibriumDiffers (test_moc_phases.cpp).
    gas->thermo()->setState_TPX(3000.0, 3e6, "H2O:0.8, OH:0.1, H2:0.05, O2:0.05");
    NozzleOptions nozzle_opts_frozen;
    nozzle_opts_frozen.chemistry = GasChemistry::FROZEN;
    Nozzle nozzle_1d_frozen(Gas(gas, GasChemistry::FROZEN), nozzle_opts_frozen);
    auto result_1d_frozen = nozzle_1d_frozen.solve(
        ExpansionType::SUPERSONIC_AREA_RATIO, frozen_result.area_ratio);
    ASSERT_TRUE(result_1d_frozen.throat.converged);
    ASSERT_FALSE(result_1d_frozen.expansions.empty());
    ASSERT_TRUE(result_1d_frozen.expansions[0].converged);
    gas->thermo()->restoreState(result_1d_frozen.expansions[0].state);
    const double mach_1d_frozen = Goddard::mach(
        *gas->thermo(), result_1d_frozen.throat.H_stagnation, result_1d_frozen.expansions[0].gamma_s);

    gas->thermo()->setState_TPX(3000.0, 3e6, "H2O:0.8, OH:0.1, H2:0.05, O2:0.05");
    NozzleOptions nozzle_opts_equil;
    Nozzle nozzle_1d_equil(Gas(gas, GasChemistry::EQUILIBRIUM), nozzle_opts_equil);
    auto result_1d_equil = nozzle_1d_equil.solve(
        ExpansionType::SUPERSONIC_AREA_RATIO, equil_result.area_ratio);
    ASSERT_TRUE(result_1d_equil.throat.converged);
    ASSERT_FALSE(result_1d_equil.expansions.empty());
    ASSERT_TRUE(result_1d_equil.expansions[0].converged);
    gas->thermo()->restoreState(result_1d_equil.expansions[0].state);
    const double mach_1d_equil = Goddard::mach(
        *gas->thermo(), result_1d_equil.throat.H_stagnation, result_1d_equil.expansions[0].gamma_s);

    RecordProperty("frozen_1d_mach", std::to_string(mach_1d_frozen));
    RecordProperty("equilibrium_1d_mach", std::to_string(mach_1d_equil));

    EXPECT_NEAR(frozen_result.exit_mach, mach_1d_frozen, 0.5)
        << "Inverse-march frozen exit Mach should be in the same ballpark as 1D";
    EXPECT_NEAR(equil_result.exit_mach, mach_1d_equil, 0.5)
        << "Inverse-march equilibrium exit Mach should be in the same ballpark as 1D";

    // Same ordering as MocFrozenTest.AxiFrozenVsEquilibriumDiffers: frozen and equilibrium
    // must give different exit Mach, and (since neither 1D reference is degenerate here) the
    // inverse-march results must differ in the same direction as their 1D counterparts.
    EXPECT_GT(std::abs(frozen_result.exit_mach - equil_result.exit_mach), 1e-6)
        << "Frozen and equilibrium should produce different exit Mach";
    if (std::abs(mach_1d_frozen - mach_1d_equil) > 1e-6) {
        EXPECT_EQ(frozen_result.exit_mach > equil_result.exit_mach,
                  mach_1d_frozen > mach_1d_equil)
            << "Inverse-march frozen/equilibrium ordering should match the 1D ordering";
    }
}
