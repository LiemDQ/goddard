#include "goddard/profile.hpp"
#include "goddard/numerics.hpp"
#include <cmath>
#include <stdexcept>
#include "gtest/gtest.h"

using namespace Goddard;

namespace {
constexpr double DEG = M_PI / 180.0;

// area_ratio^2 recovered from radius_at(x_max()) should match the requested
// area ratio (r_exit/r_throat)^2, since generation always sets r_throat=1.0
// in these tests.
double recovered_area_ratio(const NozzleProfile& profile, double r_throat) {
    double r_exit = profile.radius_at(profile.x_max());
    return (r_exit * r_exit) / (r_throat * r_throat);
}

bool all_finite_over_profile(const NozzleProfile& profile, int n_samples) {
    for (int i = 0; i < n_samples; i++) {
        double frac = static_cast<double>(i) / (n_samples - 1);
        double x_query = profile.x_min() + frac * (profile.x_max() - profile.x_min());
        if (!std::isfinite(profile.slope_at(x_query)) || !std::isfinite(profile.theta_at(x_query))) {
            return false;
        }
    }
    return true;
}
} // namespace

// ============================================================
// Input validation
// ============================================================

TEST(ConicalNozzleValidationTest, InvalidAngleThrows) {
    EXPECT_THROW(NozzleProfile::generate_conical_nozzle(20.0, 0.0, 1.0, 0.0, 50), std::invalid_argument);
    EXPECT_THROW(NozzleProfile::generate_conical_nozzle(20.0, 0.0, 1.0, 90.0, 50), std::invalid_argument);
}

TEST(ConicalNozzleValidationTest, InvalidAreaRatioThrows) {
    EXPECT_THROW(NozzleProfile::generate_conical_nozzle(1.0, 0.0, 1.0, 15.0, 50), std::invalid_argument);
}

TEST(ConicalNozzleValidationTest, InvalidThroatRadiusThrows) {
    EXPECT_THROW(NozzleProfile::generate_conical_nozzle(20.0, 0.0, 0.0, 15.0, 50), std::invalid_argument);
}

TEST(ConicalNozzleValidationTest, InvalidExpansionCurveRadiusThrows) {
    // r_expansion_curve is a radius (in throat radii) and must be non-negative;
    // 0 is the valid sharp-corner case, so only negative values should throw.
    EXPECT_THROW(NozzleProfile::generate_conical_nozzle(20.0, -1.0, 1.0, 15.0, 50), std::invalid_argument);
}

TEST(BezierNozzleValidationTest, InvalidLengthFracThrows) {
    EXPECT_THROW(NozzleProfile::generate_bezier_nozzle(20.0, 30.0, 8.0, 0.382, 1.0, 0.0, 50), std::invalid_argument);
    EXPECT_THROW(NozzleProfile::generate_bezier_nozzle(20.0, 30.0, 8.0, 0.382, 1.0, 1.01, 50), std::invalid_argument);
}

TEST(BezierNozzleValidationTest, LengthFracOfExactlyOneIsValid) {
    // length_frac == 1.0 (a "100% bell") is a normal point on Rao's chart and must
    // not throw -- regression for a boundary mismatch with generate_Rao_TOP_nozzle,
    // which itself allows length_frac up to and including 1.0.
    NozzleProfile profile;
    EXPECT_NO_THROW(profile = NozzleProfile::generate_bezier_nozzle(20.0, 30.0, 8.0, 0.382, 1.0, 1.0, 50));
    EXPECT_GT(profile.size(), 0u);
}

TEST(BezierNozzleValidationTest, InvalidAreaRatioThrows) {
    EXPECT_THROW(NozzleProfile::generate_bezier_nozzle(1.0, 30.0, 8.0, 0.382, 1.0, 0.8, 50), std::invalid_argument);
}

TEST(BezierNozzleValidationTest, InvalidThroatRadiusThrows) {
    EXPECT_THROW(NozzleProfile::generate_bezier_nozzle(20.0, 30.0, 8.0, 0.382, 0.0, 0.8, 50), std::invalid_argument);
}

TEST(RaoTopNozzleValidationTest, InvalidLengthFracThrows) {
    // The Rao chart only covers length_frac in [0.6, 1.0).
    EXPECT_THROW(NozzleProfile::generate_Rao_TOP_nozzle(20.0, 1.0, 0.3, 50), std::invalid_argument);
}

TEST(RaoTopNozzleValidationTest, InvalidAreaRatioThrows) {
    // The Rao chart's minimum area ratio is 10^0.550258197 =~ 3.55.
    EXPECT_THROW(NozzleProfile::generate_Rao_TOP_nozzle(2.0, 1.0, 0.8, 50), std::invalid_argument);
}

TEST(RaoTopNozzleValidationTest, InvalidThroatRadiusThrows) {
    EXPECT_THROW(NozzleProfile::generate_Rao_TOP_nozzle(20.0, 0.0, 0.8, 50), std::invalid_argument);
}

// ============================================================
// Conical nozzle (control group; also the Bug 6 / find_index regression)
// ============================================================

TEST(ConicalNozzleTest, XCoordinatesMonotonicIncreasing) {
    NozzleProfile profile = NozzleProfile::generate_conical_nozzle(20.0, 0.0, 1.0, 15.0, 50);
    for (size_t i = 1; i < profile.size(); i++) {
        EXPECT_GT(profile.at(i).first, profile.at(i - 1).first) << "at index " << i;
    }
}

TEST(ConicalNozzleTest, RadiusNonDecreasing) {
    NozzleProfile profile = NozzleProfile::generate_conical_nozzle(20.0, 0.0, 1.0, 15.0, 50);
    for (size_t i = 1; i < profile.size(); i++) {
        EXPECT_GE(profile.at(i).second, profile.at(i - 1).second) << "at index " << i;
    }
}

TEST(ConicalNozzleTest, AreaRatioMatchesRequested) {
    // Regression for Bug 6: NozzleProfile::find_index used to always throw for a
    // query at or beyond the last point, breaking radius_at(x_max()) exactly.
    double area_ratio = 20.0;
    NozzleProfile profile = NozzleProfile::generate_conical_nozzle(area_ratio, 0.0, 1.0, 15.0, 50);
    EXPECT_NEAR(recovered_area_ratio(profile, 1.0), area_ratio, max_fp_error(area_ratio, 1e-9, 1e-9));
}

TEST(ConicalNozzleTest, HalfAngleMatchesSpecifiedAngle) {
    NozzleProfile profile = NozzleProfile::generate_conical_nozzle(20.0, 0.0, 1.0, 15.0, 50);
    EXPECT_NEAR(profile.max_theta(), 15.0 * DEG, 1e-9);
}

// ============================================================
// Conical nozzle throat-arc geometry (r_expansion_curve > 0)
// ============================================================

TEST(ConicalNozzleTest, ThroatArcLeavesTangentToConicalSection) {
    // area_ratio=4.0, r_expansion_curve=0.382 (the Kliegel-Levine default
    // curvature radius), r_throat=1.0, angle=15deg, n_points=10.
    // Geometry independently verified by hand and by debug_probes/profile_check.cpp:
    // the arc's last point (index n_points/2 - 1 = 4) is (0.098869, 1.013016),
    // and the exit (index 9) is (3.782342, 2.000000).
    NozzleProfile profile = NozzleProfile::generate_conical_nozzle(4.0, 0.382, 1.0, 15.0, 10);

    ASSERT_EQ(profile.size(), 10u);

    for (size_t i = 1; i < profile.size(); i++) {
        EXPECT_GT(profile.at(i).first, profile.at(i - 1).first) << "at index " << i;
    }

    double y_exit = profile.at(profile.size() - 1).second;
    EXPECT_NEAR(y_exit, std::sqrt(4.0) * 1.0, 1e-9);

    // n_arc = n_points/2 = 5 points (indices 0..4) form the throat arc; the cone
    // proper starts at index 5. The cone section is parametrized as a straight
    // line anchored at the arc's endpoint, so this departure segment's slope must
    // equal tan(angle) exactly -- a mismatch would mean the cone continues at the
    // wrong angle (e.g. a degrees/radians mixup) rather than tangentially from
    // where the arc left off.
    auto [x_arc_end, y_arc_end] = profile.at(4);
    auto [x_cone_next, y_cone_next] = profile.at(5);
    double departure_slope = (y_cone_next - y_arc_end) / (x_cone_next - x_arc_end);
    EXPECT_NEAR(departure_slope, std::tan(15.0 * DEG), 1e-9);
}

// ============================================================
// Bezier nozzle (Bugs 1, 3, 6)
// ============================================================

TEST(BezierNozzleTest, XCoordinatesMonotonicIncreasing) {
    // Regression for Bug 1 (tan(15) treated as radians): pre-fix, `length` came
    // out negative and the whole contour ran backwards.
    NozzleProfile profile = NozzleProfile::generate_bezier_nozzle(20.0, 30.0, 8.0, 0.382, 1.0, 0.8, 50);
    for (size_t i = 1; i < profile.size(); i++) {
        EXPECT_GT(profile.at(i).first, profile.at(i - 1).first) << "at index " << i;
    }
}

TEST(BezierNozzleTest, RadiusNonDecreasing) {
    NozzleProfile profile = NozzleProfile::generate_bezier_nozzle(20.0, 30.0, 8.0, 0.382, 1.0, 0.8, 50);
    for (size_t i = 1; i < profile.size(); i++) {
        EXPECT_GE(profile.at(i).second, profile.at(i - 1).second) << "at index " << i;
    }
}

TEST(BezierNozzleTest, AreaRatioMatchesRequested) {
    double area_ratio = 20.0;
    NozzleProfile profile = NozzleProfile::generate_bezier_nozzle(area_ratio, 30.0, 8.0, 0.382, 1.0, 0.8, 50);
    EXPECT_NEAR(recovered_area_ratio(profile, 1.0), area_ratio, max_fp_error(area_ratio, 1e-9, 1e-9));
}

TEST(BezierNozzleTest, NoDuplicateOrZeroLengthSegments) {
    // Regression for Bug 3: the arc segment's last point and the Bezier segment's
    // first point used to both evaluate to exactly the same (x, r).
    NozzleProfile profile = NozzleProfile::generate_bezier_nozzle(20.0, 30.0, 8.0, 0.382, 1.0, 0.8, 50);
    for (size_t i = 1; i < profile.size(); i++) {
        auto [x1, y1] = profile.at(i - 1);
        auto [x2, y2] = profile.at(i);
        bool is_duplicate = std::abs(x2 - x1) < 1e-12 && std::abs(y2 - y1) < 1e-12;
        EXPECT_FALSE(is_duplicate) << "duplicate point at index " << i;
    }
}

TEST(BezierNozzleTest, SlopeAndThetaAreFiniteEverywhere) {
    // Regression for Bug 3's actual symptom: the duplicate junction point caused
    // slope_at_idx to divide 0/0, producing NaN near that x.
    NozzleProfile profile = NozzleProfile::generate_bezier_nozzle(20.0, 30.0, 8.0, 0.382, 1.0, 0.8, 50);
    EXPECT_TRUE(all_finite_over_profile(profile, 200));
}

TEST(BezierNozzleTest, WallAngleWithinPhysicalBounds) {
    NozzleProfile profile = NozzleProfile::generate_bezier_nozzle(20.0, 30.0, 8.0, 0.382, 1.0, 0.8, 50);
    double max_theta_deg = profile.max_theta() / DEG;
    EXPECT_GT(max_theta_deg, 0.0);
    // The Rao chart's largest tabulated theta_n is 40.31 degrees, so 45 degrees
    // is a safe physical bound for any theta_n in its covered range.
    EXPECT_LT(max_theta_deg, 45.0);
}

TEST(BezierNozzleTest, MaxThetaNearSpecifiedThetaN) {
    // At n_points=50, the second-order finite-difference smoothing at the
    // arc/Bezier junction node measurably underestimates theta_n (~0.24 degrees
    // for this configuration) -- expected discretization behavior, not a bug.
    double theta_n = 30.0;
    NozzleProfile profile = NozzleProfile::generate_bezier_nozzle(20.0, theta_n, 8.0, 0.382, 1.0, 0.8, 50);
    EXPECT_NEAR(profile.max_theta() / DEG, theta_n, 0.5);
}

// ============================================================
// Rao TOP nozzle (Bug 2 -- the critical regression)
// ============================================================

// Grid points from RAO_PARABOLIC_NOZZLE_THETA_N/_THETA_E in src/profile.cpp,
// chosen so (area_ratio, length_frac) land exactly on a table node. At an exact
// node the bicubic spline reproduces the tabulated value exactly (t=u=0), so
// generate_Rao_TOP_nozzle must produce the same contour as calling
// generate_bezier_nozzle directly with these theta_n/theta_e values -- if the
// interpolator's axes are transposed, this fails dramatically (wrong angles
// looked up).
struct RaoGridPoint {
    double area_ratio;
    double length_frac;
    double theta_n;
    double theta_e;
};

class RaoTopNozzleGridPointTest : public ::testing::TestWithParam<RaoGridPoint> {};

TEST_P(RaoTopNozzleGridPointTest, MatchesDirectBezierAtExactGridPoint) {
    RaoGridPoint p = GetParam();
    NozzleProfile rao_profile = NozzleProfile::generate_Rao_TOP_nozzle(p.area_ratio, 1.0, p.length_frac, 50);
    NozzleProfile direct_profile = NozzleProfile::generate_bezier_nozzle(
        p.area_ratio, p.theta_n, p.theta_e, 0.382, 1.0, p.length_frac, 50);

    ASSERT_EQ(rao_profile.size(), direct_profile.size());
    for (size_t i = 0; i < rao_profile.size(); i++) {
        auto [rx, ry] = rao_profile.at(i);
        auto [dx, dy] = direct_profile.at(i);
        EXPECT_NEAR(rx, dx, max_fp_error(dx, 1e-9, 1e-9)) << "x mismatch at index " << i;
        EXPECT_NEAR(ry, dy, max_fp_error(dy, 1e-9, 1e-9)) << "r mismatch at index " << i;
    }
}

INSTANTIATE_TEST_SUITE_P(
    GridNodes, RaoTopNozzleGridPointTest,
    ::testing::Values(
        // row 0 (60%), col 0
        RaoGridPoint{3.550243958072517, 0.60, 25.5617, 21.7742},
        // row 2 (80%), col 13
        RaoGridPoint{20.143187577910176, 0.80, 28.8714, 9.0158},
        // row 3 (90%), col 25 (last column)
        RaoGridPoint{100.0018874468128, 0.90, 32.5528, 5.9733}));

TEST(RaoTopNozzleTest, LargeAreaRatioClampsInsteadOfThrowing) {
    // Regression for Bug 4 (x_max()/y_max() off-by-one) combined with Bug 2
    // (transposed axes): area ratios beyond the chart's max (100) should clamp
    // to the 100:1 column rather than throwing or reading garbage.
    NozzleProfile clamped;
    EXPECT_NO_THROW(clamped = NozzleProfile::generate_Rao_TOP_nozzle(1000.0, 1.0, 0.8, 50));
    NozzleProfile at_chart_max = NozzleProfile::generate_Rao_TOP_nozzle(100.0, 1.0, 0.8, 50);
    EXPECT_NEAR(clamped.max_theta() / DEG, at_chart_max.max_theta() / DEG, 0.05)
        << "Clamped result should closely match the chart's 100:1 column";
}

TEST(RaoTopNozzleTest, TooSmallAreaRatioThrows) {
    // Above generate_Rao_TOP_nozzle's own area_ratio > 1.0 style check, but below
    // the Rao chart's tabulated minimum (~3.55).
    EXPECT_THROW(NozzleProfile::generate_Rao_TOP_nozzle(2.0, 1.0, 0.8, 50), std::invalid_argument);
}

// ============================================================
// NozzleProfile boundary queries (isolated Bug 6 regression)
// ============================================================

TEST(NozzleProfileBoundaryQueryTest, RadiusAtExactXMaxDoesNotThrow) {
    NozzleProfile profile;
    profile.x = {0.0, 1.0, 2.0};
    profile.y = {1.0, 1.5, 2.0};
    double radius = 0.0;
    EXPECT_NO_THROW(radius = profile.radius_at(2.0));
    EXPECT_NEAR(radius, 2.0, 1e-12);
}

TEST(NozzleProfileBoundaryQueryTest, RadiusBeyondXMaxThrows) {
    NozzleProfile profile;
    profile.x = {0.0, 1.0, 2.0};
    profile.y = {1.0, 1.5, 2.0};
    EXPECT_THROW(profile.radius_at(2.5), std::runtime_error);
}
