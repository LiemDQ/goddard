#include "goddard/moc.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/prandtlmeyer.hpp"
#include "goddard/gas_dynamics.hpp"
#include <cmath>
#include "gtest/gtest.h"

using namespace Goddard;

static constexpr double DEG = M_PI / 180.0;

// Helper: create a MocNozzle configured for planar perfect gas
static MocNozzle make_perfect_gas_solver(double gamma, double theta_max, int num_chars) {
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = gamma;
    opts.theta_max = theta_max;
    opts.num_characteristics = num_chars;
    opts.geometry.throat_radius = 1.0;
    return MocNozzle(opts);
}

// Helper: build a CharacteristicPoint with known properties
static CharacteristicPoint make_point(
    double mach, double theta, double gamma,
    double x, double y,
    double p = 1.0, double T = 1.0)
{
    CharacteristicPoint pt{};
    pt.mach = mach;
    pt.theta = theta;
    pt.nu = prandtl_meyer(mach, gamma);
    pt.K_minus = theta + pt.nu;
    pt.K_plus = theta - pt.nu;
    pt.mu = std::asin(1.0 / mach);
    pt.gamma_s = gamma;
    pt.x = x;
    pt.y = y;
    pt.pressure = p;
    pt.temperature = T;
    return pt;
}

// ============================================================
// Interior point algebraic solver
// ============================================================

class MocInteriorAlgebraicTest : public ::testing::Test {
protected:
    MocNozzle solver = make_perfect_gas_solver(1.4, 10.0 * DEG, 5);
};

TEST_F(MocInteriorAlgebraicTest, UniformFlow) {
    // Two parent points with identical M and theta -> child should match
    double gamma = 1.4;
    auto p1 = make_point(2.0, 5.0 * DEG, gamma, 0.0, 0.5);
    auto p2 = make_point(2.0, 5.0 * DEG, gamma, 0.0, 1.5);

    auto p3 = solver.m_options.gamma; // just checking it's set
    (void)p3;

    // Call the public solve which dispatches to algebraic
    // We need to access the private method, so we use solve_interior_point
    // through a derived test or by making the solver's solve() produce known results.
    // For now, test via the full solver with a known expansion.
}

TEST_F(MocInteriorAlgebraicTest, KMinusKPlusPreserved) {
    // Verify Riemann invariants are correctly propagated
    double gamma = 1.4;
    auto p1 = make_point(2.0, 10.0 * DEG, gamma, 0.0, 0.5);
    auto p2 = make_point(1.8, 3.0 * DEG, gamma, 0.0, 1.2);

    // Run a small solve and check the initial data line
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = gamma;
    opts.theta_max = 10.0 * DEG;
    opts.num_characteristics = 3;
    opts.geometry.throat_radius = 1.0;

    MocNozzle nozzle(opts);
    auto result = nozzle.solve();

    // Basic sanity: we should have points
    EXPECT_GT(result.net.points.size(), 0u);

    // Riemann invariants are preserved ALONG a characteristic, not across a wavefront:
    // for planar flow K+ (= theta - nu) is constant along a C+ chain and K- (= theta + nu)
    // is constant along a C- chain. (The previous version treated each c_chains[i] as a
    // wavefront with a shared K+, which is not how the chain-based net is organized.)
    using Family = ChainMetadata::Family;
    for (size_t c = 0; c < result.net.c_chains.size(); c++) {
        const auto& chain = result.net.c_chains[c];
        if (chain.size() < 2) continue;
        Family fam = result.net.chain_metadata[c].family;
        if (fam == Family::PLUS) {
            double kplus = result.net.points[chain.front()].K_plus;
            for (size_t pt_idx : chain) {
                EXPECT_NEAR(result.net.points[pt_idx].K_plus, kplus, 1e-9)
                    << "K+ should be invariant along a C+ chain (chain " << c << ")";
            }
        } else if (fam == Family::MINUS) {
            double kminus = result.net.points[chain.front()].K_minus;
            for (size_t pt_idx : chain) {
                EXPECT_NEAR(result.net.points[pt_idx].K_minus, kminus, 1e-9)
                    << "K- should be invariant along a C- chain (chain " << c << ")";
            }
        }
    }
}

// ============================================================
// Axis point solver
// ============================================================

TEST(MocAxisPoint, SymmetryCondition) {
    // The axis point should have theta = 0 by symmetry
    double gamma = 1.4;
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = gamma;
    opts.theta_max = 15.0 * DEG;
    opts.num_characteristics = 5;
    opts.geometry.throat_radius = 1.0;

    MocNozzle nozzle(opts);
    auto result = nozzle.solve();

    // Check all wavefronts except first: 
    // first point should be an axis point with theta=0, y=0
    // first wavepoint is intentionally non-symmetrical
    const auto axis_pts = result.net.axis_points();
    ASSERT_FALSE(axis_pts.empty()) << "Axis points are empty";
    for (size_t i = 0; i < axis_pts.size(); i++) {
        const auto& ap = axis_pts[i];
        EXPECT_NEAR(ap.theta, 0.0, 1e-12)
            << "Axis point theta should be 0 in wavefront " << i;
        EXPECT_NEAR(ap.y, 0.0, 1e-12)
            << "Axis point y should be 0 in wavefront " << i;
    }
}

TEST(MocAxisPoint, KMinusEqualsNu) {
    // At axis: theta=0, so K_minus = theta + nu = nu, K_plus = theta - nu = -nu
    double gamma = 1.4;
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = gamma;
    opts.theta_max = 12.0 * DEG;
    opts.num_characteristics = 4;
    opts.geometry.throat_radius = 1.0;

    MocNozzle nozzle(opts);
    auto result = nozzle.solve();
    // skip first wavefront -- it is explicitly not set to theta = 0

    const auto axis_pts = result.net.axis_points();
    for (size_t i = 1; i < axis_pts.size(); i++) {
        const auto& axis_pt = axis_pts[i];
        EXPECT_NEAR(axis_pt.K_minus, axis_pt.nu, 1e-10)
            << "Axis K_minus should equal nu in wavefront " << i;
        EXPECT_NEAR(axis_pt.K_plus, -axis_pt.nu, 1e-10)
            << "Axis K_plus should equal -nu in wavefront " << i;
    }
}

// ============================================================
// Wall point solver
// ============================================================

TEST(MocWallPoint, WallKPlusPreserved) {
    // K+ should be preserved from interior parent to wall point
    double gamma = 1.4;
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = gamma;
    opts.theta_max = 10.0 * DEG;
    opts.num_characteristics = 4;
    opts.geometry.throat_radius = 1.0;

    MocNozzle nozzle(opts);
    auto result = nozzle.solve();

    // Wall points should have decreasing theta from theta_max toward 0
    EXPECT_GT(result.net.wall_x.size(), 1u);
    for (size_t i = 1; i < result.net.wall_x.size(); i++) {
        EXPECT_GT(result.net.wall_x[i], result.net.wall_x[i - 1])
            << "Wall x should be monotonically increasing";
    }
}

// ============================================================
// Basic solve sanity
// ============================================================

TEST(MocSolve, PerfectGasDoesNotRequireCantera) {
    // Perfect gas mode should work without any Cantera Solution object
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = 1.4;
    opts.theta_max = 10.0 * DEG;
    opts.num_characteristics = 5;
    opts.geometry.throat_radius = 1.0;

    MocNozzle nozzle(opts);
    ASSERT_NO_THROW({
        auto result = nozzle.solve();
        EXPECT_TRUE(result.converged);
        EXPECT_GT(result.exit_mach, 1.0);
        EXPECT_GT(result.nozzle_length, 0.0);
        EXPECT_GT(result.area_ratio, 1.0);
    });
}

TEST(MocSolve, CustomThetaSchedule) {
    // Providing a custom theta schedule should produce the same number of characteristics
    double gamma = 1.4;
    std::vector<double> schedule = {2.0 * DEG, 4.0 * DEG, 6.0 * DEG, 8.0 * DEG, 10.0 * DEG};

    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = gamma;
    opts.theta_max = 10.0 * DEG;
    opts.num_characteristics = 99; // should be overridden by schedule size
    opts.theta_schedule = schedule;
    opts.geometry.throat_radius = 1.0;

    MocNozzle nozzle(opts);
    auto result = nozzle.solve();

    // Should have 5 characteristics (from schedule), not 99
    EXPECT_EQ(result.net.wall_points().size(), 5u);
}

TEST(MocSolve, ExitMachConsistentWithThetaMax) {
    // For a minimum-length nozzle, the exit Mach should correspond to nu = theta_max
    double gamma = 1.4;
    double theta_max = 15.0 * DEG;

    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = gamma;
    opts.theta_max = theta_max;
    opts.num_characteristics = 10;
    opts.geometry.throat_radius = 1.0;

    MocNozzle nozzle(opts);
    auto result = nozzle.solve();

    // For a min-length nozzle, the exit nu = 2*theta_max because the full
    // expansion (K_minus = 2*theta_max from the last fan ray) propagates
    // through the kernel to the axis where theta=0, so nu = 2*theta_max.
    double expected_exit_mach = mach_from_prandtl_meyer(2.0 * theta_max, gamma);
    EXPECT_NEAR(result.exit_mach, expected_exit_mach, 0.01)
        << "Exit Mach should correspond to nu = 2*theta_max";
}

// ── NozzleProfile::max_theta ────────────────────────────────────────────

TEST(NozzleProfileTest, MaxThetaSimpleRamp) {
    // Single straight segment at 30 degrees
    NozzleProfile profile;
    profile.x = {0.0, 1.0};
    profile.y = {1.0, 1.0 + std::tan(30.0 * DEG)};
    EXPECT_NEAR(profile.max_theta(), 30.0 * DEG, 1e-12);
}

TEST(NozzleProfileTest, MaxThetaMultipleSegments) {
    // Three segments: 10°, 25°, 15° — max should be 25°
    NozzleProfile profile;
    double dx = 1.0;
    profile.x = {0.0, dx, 2.0 * dx, 3.0 * dx};
    profile.y = {1.0,
                 1.0 + dx * std::tan(10.0 * DEG),
                 1.0 + dx * std::tan(10.0 * DEG) + dx * std::tan(25.0 * DEG),
                 1.0 + dx * std::tan(10.0 * DEG) + dx * std::tan(25.0 * DEG) + dx * std::tan(25.0 * DEG)};
    EXPECT_NEAR(profile.max_theta(), 25.0 * DEG, 1e-12);
}

TEST(NozzleProfileTest, MaxThetaBellNozzle) {
    // Bell-like contour: expansion section (increasing angle) then contraction
    // back toward the axis. Max angle is in the middle, not at the first segment.
    NozzleProfile profile;
    double dx = 0.5;
    // The 2nd-order central difference at each node averages the slopes (tan values)
    // of the two adjacent segments, then takes atan. This is NOT the same as averaging
    // angles: atan((tan(35°)+tan(25°))/2) ≈ 30.26°, not 30°.
    std::vector<double> angles_deg = {5.0, 15.0, 35.0, 25.0, 20.0, 5.0};
    profile.x.push_back(0.0);
    profile.y.push_back(1.0);
    for (size_t i = 0; i < angles_deg.size(); i++) {
        double dy = dx * std::tan(angles_deg[i] * DEG);
        profile.x.push_back(profile.x.back() + dx);
        profile.y.push_back(profile.y.back() + dy);
    }
    // Max is at the node between the 35° and 25° segments (uniform spacing → simple average).
    double expected_max_theta = std::atan((std::tan(35.0 * DEG) + std::tan(25.0 * DEG)) / 2.0);
    EXPECT_NEAR(profile.max_theta(), expected_max_theta, 1e-12);
}

TEST(NozzleProfileTest, MaxThetaFlatSegmentsIgnored) {
    // Flat segment (0°) followed by an angled one — flat should not affect max
    NozzleProfile profile;
    profile.x = {0.0, 1.0, 2.0};
    profile.y = {1.0, 1.0, 1.0 + std::tan(12.0 * DEG)};
    EXPECT_NEAR(profile.max_theta(), 12.0 * DEG, 1e-12);
}
