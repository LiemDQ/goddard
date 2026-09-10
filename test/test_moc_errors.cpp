#include "goddard/moc.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/characteristics.hpp"
#include "goddard/profile.hpp"
#include "goddard/gas.hpp"
#include "cantera/core.h"
#include <cmath>
#include <limits>
#include <string>
#include "gtest/gtest.h"

using namespace Goddard;

static constexpr double DEG = M_PI / 180.0;

namespace {

CharacteristicPoint make_valid_point() {
    CharacteristicPoint pt{};
    pt.mach = 2.0;
    pt.theta = 10.0 * DEG;
    pt.nu = 20.0 * DEG;
    pt.mu = std::asin(1.0 / pt.mach);
    pt.x = 1.0;
    pt.y = 1.0;
    return pt;
}

// Every point that actually entered the net via a unit process must
// independently pass the validity check -- this is the invariant that makes
// converged == true trustworthy. The one exception is CharacteristicNet's
// throat-lip seed point (CharacteristicNet::seed_wall_point, used to anchor
// the wall march for a centered-fan initial data line): it owns no
// characteristic chain (both membership optionals empty) and is not the
// output of any unit process, so it intentionally carries incomplete flow
// state (e.g. mach == 0) by design and is exempt.
void expect_all_points_valid(const MocResult& result, double tol) {
    ASSERT_EQ(result.net.points.size(), result.net.membership.size());
    for (size_t i = 0; i < result.net.points.size(); i++) {
        const PointMembership& mem = result.net.membership[i];
        if (!mem.c_plus_chain_idx.has_value() && !mem.c_minus_chain_idx.has_value()) {
            continue;
        }
        const CharacteristicPoint& pt = result.net.points[i];
        EXPECT_EQ(check_point_validity(pt, tol), MocErrorCode::NONE)
            << "point at (" << pt.x << ", " << pt.y << ") failed validity check";
    }
}

} // namespace

// ============================================================
// check_point_validity: unit tests for each violated condition
// ============================================================

TEST(CheckPointValidity, ValidPointReturnsNone) {
    auto pt = make_valid_point();
    EXPECT_EQ(check_point_validity(pt, 1e-10), MocErrorCode::NONE);
}

TEST(CheckPointValidity, NegativeNuDetected) {
    auto pt = make_valid_point();
    pt.nu = -0.01;
    EXPECT_EQ(check_point_validity(pt, 1e-10), MocErrorCode::NEGATIVE_NU);
}

TEST(CheckPointValidity, NegativeThetaBeyondToleranceDetected) {
    auto pt = make_valid_point();
    pt.theta = -0.01;
    EXPECT_EQ(check_point_validity(pt, 1e-10), MocErrorCode::NEGATIVE_THETA);
}

TEST(CheckPointValidity, SmallNegativeThetaWithinToleranceIsValid) {
    // Axis points legitimately pin theta to exactly 0.0, and interior points
    // can carry small negative roundoff; neither should trip NEGATIVE_THETA.
    auto pt = make_valid_point();
    const double tol = 1e-8;
    pt.theta = -0.5 * tol;
    EXPECT_EQ(check_point_validity(pt, tol), MocErrorCode::NONE);
}

TEST(CheckPointValidity, ExactlyZeroThetaIsValid) {
    // Axis points pin theta to exactly 0.0.
    auto pt = make_valid_point();
    pt.theta = 0.0;
    EXPECT_EQ(check_point_validity(pt, 1e-10), MocErrorCode::NONE);
}

TEST(CheckPointValidity, SubsonicMachDetected) {
    auto pt = make_valid_point();
    pt.mach = 0.8;
    EXPECT_EQ(check_point_validity(pt, 1e-10), MocErrorCode::SUBSONIC_MACH);
}

TEST(CheckPointValidity, NonFiniteYDetected) {
    // Keep nu/theta/mach independently valid so the earlier checks (which
    // precede finiteness in priority order) do not mask the NaN.
    auto pt = make_valid_point();
    pt.y = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(check_point_validity(pt, 1e-10), MocErrorCode::NONFINITE_VALUE);
}

TEST(CheckPointValidity, InfiniteXDetected) {
    auto pt = make_valid_point();
    pt.x = std::numeric_limits<double>::infinity();
    EXPECT_EQ(check_point_validity(pt, 1e-10), MocErrorCode::NONFINITE_VALUE);
}

TEST(CheckPointValidity, NonFiniteMuDetected) {
    auto pt = make_valid_point();
    pt.mu = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(check_point_validity(pt, 1e-10), MocErrorCode::NONFINITE_VALUE);
}

// ============================================================
// Regression: representative previously-passing configurations must still
// report converged == true with an all-NONE failure record, and every point
// that entered the net must independently pass check_point_validity.
// ============================================================

TEST(MocErrorsRegression, AndersonMinLengthConverges) {
    // Anderson, "Modern Compressible Flow" Ch. 11, Table 11.1.
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.gamma = 1.40;
    opts.theta_max = 18.375 * DEG;
    opts.theta_schedule = {
        0.375 * DEG, 3.375 * DEG, 6.375 * DEG, 9.375 * DEG,
        12.375 * DEG, 15.375 * DEG, 18.375 * DEG};
    opts.num_characteristics = static_cast<int>(opts.theta_schedule.size());
    opts.geometry.throat_radius = 1.0;

    MocNozzle nozzle(opts);
    MocResult result;
    ASSERT_NO_THROW({ result = nozzle.solve(); });

    ASSERT_TRUE(result.converged);
    EXPECT_EQ(result.failure.code, MocErrorCode::NONE);
    EXPECT_NEAR(result.exit_mach, 2.4, 0.01)
        << "Exit Mach should match Anderson Table 11.1 (nu_exit = 2*theta_max)";

    expect_all_points_valid(result, opts.solver_options.abstol);
}

class MocErrorsFrozenTest : public ::testing::Test {
protected:
    void SetUp() override {
        gas = Cantera::newSolution("h2o2.yaml", "ohmech");
        gas->thermo()->setState_TPX(3000.0, 3e6, "H2O:0.8, OH:0.1, H2:0.05, O2:0.05");
    }
    std::shared_ptr<Cantera::Solution> gas;
};

TEST_F(MocErrorsFrozenTest, BasicFrozenCaseConverges) {
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::FROZEN;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.theta_max = 15.0 * DEG;
    opts.num_characteristics = 7;
    opts.geometry.throat_radius = 0.05;

    MocNozzle nozzle(Gas(gas, opts.chemistry), opts);
    MocResult result;
    ASSERT_NO_THROW({ result = nozzle.solve(); });

    ASSERT_TRUE(result.converged);
    EXPECT_EQ(result.failure.code, MocErrorCode::NONE);

    expect_all_points_valid(result, opts.solver_options.abstol);
}

TEST_F(MocErrorsFrozenTest, BasicEquilibriumCaseConverges) {
    MocOptions opts;
    opts.flow_type = MocFlowKind::PLANAR;
    opts.chemistry = GasChemistry::EQUILIBRIUM;
    opts.mode = MocMode::DESIGN_MIN_LENGTH;
    opts.theta_max = 15.0 * DEG;
    opts.num_characteristics = 7;
    opts.geometry.throat_radius = 0.05;

    MocNozzle nozzle(Gas(gas, opts.chemistry), opts);
    MocResult result;
    ASSERT_NO_THROW({ result = nozzle.solve(); });

    ASSERT_TRUE(result.converged);
    EXPECT_EQ(result.failure.code, MocErrorCode::NONE);

    expect_all_points_valid(result, opts.solver_options.abstol);
}

TEST(MocErrorsRegression, DesignRaoLeavesNoInvalidPoints) {
    // Mirrors MocDesignRao.BasicSolveReachesRecordedCoverage (test_moc_rao.cpp), which
    // explains why this configuration is asserted rather than skipped. The load-bearing
    // assertion here is different from that one's: whatever the march achieves before it
    // stops, every point it did commit to the net must be valid. A solver that fails is
    // acceptable; a solver that leaves NaNs or subsonic points behind is not, and that
    // check must not be skipped away just because the march did not finish.
    MocOptions opts;
    opts.flow_type = MocFlowKind::AXISYMMETRIC;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_RAO;
    opts.gamma = 1.23;
    opts.num_characteristics = 15;
    opts.geometry.throat_radius = 1.0;
    opts.geometry.expansion_ratio = 5.0;
    opts.geometry.length_fraction = 0.8;

    MocNozzle nozzle(opts);
    MocResult result;
    ASSERT_NO_THROW({ result = nozzle.solve(); });

    RecordProperty("exit_coverage", std::to_string(result.exit_coverage));
    RecordProperty("failure_code", std::string(to_string(result.failure.code)));

    EXPECT_EQ(result.converged, result.failure.code == MocErrorCode::NONE)
        << "converged and failure.code disagree: " << result.failure.message;
    EXPECT_GT(result.exit_mach, 1.0);
    EXPECT_GT(result.area_ratio, 1.0);

    expect_all_points_valid(result, opts.solver_options.abstol);
}

// ============================================================
// Failure honesty: a configuration known to trigger huge marching steps (and,
// before this phase, a silently-laundered Prandtl-Meyer inversion failure)
// must report converged == false with a specific error code, must never
// throw, and -- the load-bearing assertion -- must never leave a NaN in the
// net. See instructions/moc_convergence_roadmap.md Sec. 1/4 for the original
// diagnosis of this laundering bug.
// ============================================================

TEST(MocErrorsFailureHonesty, AxisymmetricConicalAnalysisReportsFailureHonestly) {
    // AR=4/N=8 was the original known-broken configuration this test targeted; the
    // Phase 2 near-axis-void fix (dual-family KL-line seeding via a
    // downstream-shifted start line) now makes it converge. AR=8 remains a reliably
    // non-converging configuration for the DIRECT kernel: with the wall-consistent
    // Kliegel-Levine line it marches to the compression that converges on the axis near
    // x = 2.7 (instructions/moc_fix/diagnosis.md A8) and stops there on the flow-angle
    // check. The INVERSE kernel, now the default for axisymmetric analysis, passes through
    // that compression, so the scheme is pinned to DIRECT here to keep exercising the
    // failure-honesty machinery.
    MocOptions opts;
    opts.march_scheme = MocMarchScheme::DIRECT;
    opts.flow_type = MocFlowKind::AXISYMMETRIC;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::ANALYSIS;
    opts.gamma = 1.4;
    opts.num_characteristics = 8;
    opts.geometry.throat_radius = 1.0;
    // downstream_wall_curvature_radius left at its positive default (0.382):
    // generate_initial_data_line() takes the Kliegel-Levine transonic path.
    opts.nozzle_profile = NozzleProfile::generate_conical_nozzle(8.0, 0.382, 1.0, 15.0, 60);

    MocNozzle nozzle(opts);
    MocResult result;
    ASSERT_NO_THROW({ result = nozzle.solve(); })
        << "solve() must never throw for a numerical/convergence failure";

    EXPECT_FALSE(result.converged);
    EXPECT_NE(result.failure.code, MocErrorCode::NONE)
        << "A non-converging solve must carry a specific, non-NONE error code";

    // The laundering path this phase closes: a PM-inversion or table-lookup
    // failure must never flow into mach_to_mu()/corrector-step averaging and
    // produce a NaN that silently enters the net. No point that was actually
    // inserted into the net -- not even ones computed before the abort -- may
    // be non-finite or otherwise invalid.
    for (const auto& pt : result.net.points) {
        EXPECT_TRUE(std::isfinite(pt.x)) << "x is non-finite";
        EXPECT_TRUE(std::isfinite(pt.y)) << "y is non-finite";
        EXPECT_TRUE(std::isfinite(pt.theta)) << "theta is non-finite";
        EXPECT_TRUE(std::isfinite(pt.nu)) << "nu is non-finite";
        EXPECT_TRUE(std::isfinite(pt.mach)) << "mach is non-finite";
        EXPECT_TRUE(std::isfinite(pt.mu)) << "mu is non-finite";
    }
    expect_all_points_valid(result, opts.solver_options.abstol);
}
