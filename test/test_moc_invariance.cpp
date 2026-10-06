#include "goddard/moc.hpp"
#include "goddard/moc_nozzle.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>
#include "gtest/gtest.h"

using namespace Goddard;

// ============================================================
// Bit-level invariance of the minimum-length design path.
//
// DESIGN_MIN_LENGTH planar and axisymmetric output must stay unchanged by
// any work on the Kliegel-Levine start line: min-length always initializes
// with a centered expansion fan (build_start_line, moc_initialization.cpp)
// and uses the chain ladder (DirectMarch) only, unaffected by anything on
// the inverse-march path. That invariance was not enforced by any test, so
// a KL-path change could silently move the one configuration validated
// against Anderson Table 11.1.
//
// The references were captured from the tree at the baseline commit with
// %.17g and are exact to the last bit there. The comparison tolerance is
// 1e-12 relative rather than exact equality so the test survives a compiler
// or optimization-flag change without going quiet about a real regression;
// any actual change to the march moves these numbers far more than that.
// ============================================================

namespace {

struct WallPoint { double x, y; };

struct GoldenCase {
    const char* name;
    MocFlowKind flow;
    double gamma;
    double theta_max_deg;
    int num_characteristics;
    double exit_mach;
    double area_ratio;
    double nozzle_length;
    std::vector<WallPoint> wall;
};

MocResult solve_golden(const GoldenCase& c) {
    MocOptions o;
    o.flow_type = c.flow;
    o.chemistry = GasChemistry::PERFECT_GAS;
    o.mode = MocMode::DESIGN_MIN_LENGTH;
    o.gamma = c.gamma;
    o.theta_max = c.theta_max_deg * M_PI / 180.0;
    o.num_characteristics = c.num_characteristics;
    o.geometry.throat_radius = 1.0;
    MocNozzle solver(o);
    return solver.solve();
}

// Relative comparison, falling back to absolute near zero (wall[0] is exactly
// the throat lip at x = 0).
void expect_matches(double actual, double expected, const char* what, size_t idx) {
    const double scale = std::max(1.0, std::abs(expected));
    EXPECT_NEAR(actual, expected, 1e-12 * scale)
        << what << " at index " << idx << " moved: got " << actual
        << ", golden " << expected;
}

std::vector<GoldenCase> golden_cases() {
    return {
        GoldenCase{
            "planar", MocFlowKind::PLANAR, 1.4, 15.0, 20,
            2.1339050332318474, 1.8894505839966758, 5.7525398692310965,
            {
                {0.0,                  1.0},
                {0.84848684094848181,  1.2267562636210831},
                {1.4393456899580401,   1.3799174783727617},
                {1.7157072809023792,   1.4475259131248843},
                {1.9457014930234722,   1.5004601690089807},
                {2.1581412210345166,   1.5462967986351306},
                {2.3633269123960465,   1.5876327152960039},
                {2.5663458100988668,   1.62564355910336},
                {2.7702127376501844,   1.6609274072012852},
                {2.97695386701223,     1.6937961022812793},
                {3.1880680299060695,   1.7243987739468074},
                {3.4047521753479462,   1.7527815177819939},
                {3.6280231724077194,   1.7789188686949713},
                {3.8587889513427034,   1.8027313780510814},
                {4.097892821716087,    1.8240957364637855},
                {4.3461426830770007,   1.8428506051261013},
                {4.6043313015788492,   1.8587998167277315},
                {4.873251102039065,    1.8717138660588886},
                {5.1537054986910533,   1.8813302209203042},
                {5.4465179990578525,   1.8873527680692519},
                {5.7525398692310965,   1.8894505839966758},
            }},
        GoldenCase{
            "axisymmetric", MocFlowKind::AXISYMMETRIC, 1.4, 12.0, 20,
            3.0322104501599916, 4.0581414538115386, 8.3750502491154535,
            {
                {0.0,                  1.0},
                {0.91502911251641772,  1.1939947432958415},
                {1.8244493139662288,   1.3810985288530089},
                {2.2817368434035665,   1.4699644237231251},
                {2.642309415529887,    1.5359398920539205},
                {2.9762877751102863,   1.5932710189157973},
                {3.2955438466260936,   1.6444767407150247},
                {3.611497178272292,    1.6916043046251925},
                {3.92810533556885,     1.7352853749458101},
                {4.2488171773160017,   1.7759533575284077},
                {4.5758241760473801,   1.8137802259827738},
                {4.9108522294067667,   1.8488157449385167},
                {5.2552234457721951,   1.8810141710246369},
                {5.6099884914708635,   1.9102631906388341},
                {5.9759709108738948,   1.9363991426644018},
                {6.3538013410468865,   1.9592191070366725},
                {6.7439452489488216,   1.9784907298895873},
                {7.1467530654196851,   1.9939613396819718},
                {7.5625585713003733,   2.0053656020351727},
                {7.9875177651361913,   2.0123576720398484},
                {8.3750502491154535,   2.0144829246760914},
            }},
    };
}

} // namespace

TEST(MocInvariance, DesignMinLengthUnchanged) {
    for (const GoldenCase& c : golden_cases()) {
        SCOPED_TRACE(c.name);
        MocResult r = solve_golden(c);

        ASSERT_TRUE(r.converged) << c.name << ": " << r.failure.message;
        expect_matches(r.exit_mach, c.exit_mach, "exit_mach", 0);
        expect_matches(r.area_ratio, c.area_ratio, "area_ratio", 0);
        expect_matches(r.nozzle_length, c.nozzle_length, "nozzle_length", 0);

        ASSERT_EQ(r.net.wall_x.size(), c.wall.size()) << c.name << ": wall point count changed";
        ASSERT_EQ(r.net.wall_y.size(), c.wall.size()) << c.name;
        for (size_t i = 0; i < c.wall.size(); i++) {
            expect_matches(r.net.wall_x[i], c.wall[i].x, "wall_x", i);
            expect_matches(r.net.wall_y[i], c.wall[i].y, "wall_y", i);
        }
    }
}

// ============================================================
// Throat-radius scaling invariance.
//
// The inviscid compatibility relations contain lengths only as ratios (the axisymmetric
// source term is sin(mu) sin(theta) / (y cos(theta -/+ mu)) times dx), so a nozzle k times
// larger has the same flow field at corresponding points. Solving with
// NozzleGeometry::throat_radius = k (and, in analysis, a contour built with r_throat = k)
// must therefore reproduce the throat_radius = 1 solve exactly up to roundoff: every
// dimensionless output (Mach numbers, area ratio, thrust coefficient, flow angles) is
// unchanged, and every length (net coordinates, contour, nozzle length, exit-plane
// stations) is multiplied by k.
//
// 1e-9 relative: the solver works in throat radii internally, so the only differences
// between the two solves are the roundoff of converting the input contour to throat radii
// and the outputs back (a few ulps). Any unit inconsistency (a start line built in throat
// radii against a wall in physical units, a curvature ratio divided by the throat radius
// twice, an area ratio taken against the wrong reference) moves these numbers by O(1).
// ============================================================

namespace {

constexpr double kThroatScale = 2.5;
constexpr double kScaleRelTol = 1e-9;

/** Largest relative deviation of `scaled` from k * `unit`, with an absolute floor of k. */
double max_length_deviation(const std::vector<double>& scaled, const std::vector<double>& unit,
                            double k)
{
    double worst = 0.0;
    for (size_t i = 0; i < scaled.size() && i < unit.size(); i++) {
        const double deviation = std::abs(scaled[i] - k * unit[i]) / (k * std::max(1.0, std::abs(unit[i])));
        if (std::isnan(deviation)) return std::numeric_limits<double>::infinity();
        worst = std::max(worst, deviation);
    }
    return worst;
}

/** Largest relative deviation between two dimensionless arrays, with an absolute floor of 1. */
double max_invariant_deviation(const std::vector<double>& scaled, const std::vector<double>& unit) {
    double worst = 0.0;
    for (size_t i = 0; i < scaled.size() && i < unit.size(); i++) {
        const double deviation = std::abs(scaled[i] - unit[i]) / std::max(1.0, std::abs(unit[i]));
        if (std::isnan(deviation)) return std::numeric_limits<double>::infinity();
        worst = std::max(worst, deviation);
    }
    return worst;
}

void expect_invariant(double scaled, double unit, const char* what) {
    EXPECT_NEAR(scaled, unit, kScaleRelTol * std::max(1.0, std::abs(unit)))
        << what << " must not depend on the throat radius";
}

void expect_length_scaled(double scaled, double unit, double k, const char* what) {
    EXPECT_NEAR(scaled, k * unit, kScaleRelTol * k * std::max(1.0, std::abs(unit)))
        << what << " must scale with the throat radius";
}

/**
 * Area-weighted mean Mach number over the exit plane: sum(M dA) / sum(dA), with dA = y dy
 * (axisymmetric) or dy (planar), trapezoidal. Dimensionless, so scale invariant.
 */
double exit_plane_mean_mach(const MocResult& result, MocFlowKind flow) {
    const ExitPlane& ep = result.exit_plane;
    double weighted = 0.0;
    double area = 0.0;
    for (size_t i = 1; i < ep.y.size(); i++) {
        const double dy = ep.y[i] - ep.y[i - 1];
        const double weight = (flow == MocFlowKind::AXISYMMETRIC) ? 0.5 * (ep.y[i] + ep.y[i - 1]) : 1.0;
        weighted += 0.5 * (ep.mach[i] + ep.mach[i - 1]) * weight * dy;
        area += weight * dy;
    }
    return (area > 0.0) ? weighted / area : 0.0;
}

/** Copy of `profile` with every coordinate multiplied by k, built point by point. */
NozzleProfile scale_profile(const NozzleProfile& profile, double k) {
    NozzleProfile out;
    for (size_t i = 0; i < profile.size(); i++) {
        out.push_back({k * profile.x[i], k * profile.y[i]});
    }
    out.throat_index = profile.throat_index;
    return out;
}

/**
 * Every invariant and every scaled length of a throat_radius = k solve against the
 * throat_radius = 1 solve of the same nozzle.
 */
void expect_scaled_solution(const MocResult& scaled, const MocResult& unit, double k,
                            MocFlowKind flow)
{
    ASSERT_TRUE(unit.converged) << "reference (throat_radius = 1): "
        << to_string(unit.failure.code) << " -- " << unit.failure.message;
    ASSERT_TRUE(scaled.converged) << "throat_radius = " << k << ": "
        << to_string(scaled.failure.code) << " -- " << scaled.failure.message;

    // Dimensionless performance.
    expect_invariant(scaled.area_ratio, unit.area_ratio, "area_ratio");
    expect_invariant(scaled.exit_mach, unit.exit_mach, "exit_mach");
    expect_invariant(scaled.exit_coverage, unit.exit_coverage, "exit_coverage");
    expect_invariant(scaled.min_theta, unit.min_theta, "min_theta");
    expect_invariant(exit_plane_mean_mach(scaled, flow), exit_plane_mean_mach(unit, flow),
                     "area-weighted exit-plane Mach");

    // Thrust coefficient is force / (p0 * A_throat): the exit-plane integral grows with k
    // (planar) or k^2 (axisymmetric), and so does the throat area it is normalized by.
    const double ambient_pressure_ratio = 0.01;
    const ThrustCoefficient cf_scaled = compute_thrust_coefficient(scaled, flow, ambient_pressure_ratio);
    const ThrustCoefficient cf_unit = compute_thrust_coefficient(unit, flow, ambient_pressure_ratio);
    expect_invariant(cf_scaled.Cf_vacuum, cf_unit.Cf_vacuum, "Cf_vacuum");
    expect_invariant(cf_scaled.Cf, cf_unit.Cf, "Cf");
    expect_invariant(cf_scaled.momentum_thrust, cf_unit.momentum_thrust, "momentum_thrust");

    // Lengths.
    expect_length_scaled(scaled.nozzle_length, unit.nozzle_length, k, "nozzle_length");
    expect_length_scaled(scaled.min_theta_x, unit.min_theta_x, k, "min_theta_x");
    expect_length_scaled(scaled.min_theta_y, unit.min_theta_y, k, "min_theta_y");
    expect_length_scaled(scaled.init_diagnostics.wall_gap, unit.init_diagnostics.wall_gap, k,
                         "init_diagnostics.wall_gap");

    ASSERT_EQ(scaled.net.wall_x.size(), unit.net.wall_x.size()) << "wall point count changed";
    EXPECT_LT(max_length_deviation(scaled.net.wall_x, unit.net.wall_x, k), kScaleRelTol) << "net.wall_x";
    EXPECT_LT(max_length_deviation(scaled.net.wall_y, unit.net.wall_y, k), kScaleRelTol) << "net.wall_y";

    ASSERT_EQ(scaled.profile.size(), unit.profile.size()) << "profile point count changed";
    EXPECT_LT(max_length_deviation(scaled.profile.x, unit.profile.x, k), kScaleRelTol) << "profile.x";
    EXPECT_LT(max_length_deviation(scaled.profile.y, unit.profile.y, k), kScaleRelTol) << "profile.y";

    ASSERT_EQ(scaled.exit_plane.y.size(), unit.exit_plane.y.size()) << "exit-plane point count changed";
    EXPECT_LT(max_length_deviation(scaled.exit_plane.y, unit.exit_plane.y, k), kScaleRelTol) << "exit_plane.y";
    EXPECT_LT(max_invariant_deviation(scaled.exit_plane.mach, unit.exit_plane.mach), kScaleRelTol)
        << "exit_plane.mach";
    EXPECT_LT(max_invariant_deviation(scaled.exit_plane.theta, unit.exit_plane.theta), kScaleRelTol)
        << "exit_plane.theta";

    // The whole net, point by point.
    ASSERT_EQ(scaled.net.points.size(), unit.net.points.size()) << "net point count changed";
    std::vector<double> xs, xu, ys, yu, mach_s, mach_u;
    for (size_t i = 0; i < unit.net.points.size(); i++) {
        xs.push_back(scaled.net.points[i].x);
        xu.push_back(unit.net.points[i].x);
        ys.push_back(scaled.net.points[i].y);
        yu.push_back(unit.net.points[i].y);
        mach_s.push_back(scaled.net.points[i].mach);
        mach_u.push_back(unit.net.points[i].mach);
    }
    EXPECT_LT(max_length_deviation(xs, xu, k), kScaleRelTol) << "net point x";
    EXPECT_LT(max_length_deviation(ys, yu, k), kScaleRelTol) << "net point y";
    EXPECT_LT(max_invariant_deviation(mach_s, mach_u), kScaleRelTol) << "net point Mach";

    // Per-pass front geometry (inverse march only; empty for minimum-length design).
    ASSERT_EQ(scaled.pass_diagnostics.size(), unit.pass_diagnostics.size()) << "pass count changed";
    std::vector<double> step_s, step_u, spacing_s, spacing_u;
    for (size_t i = 0; i < unit.pass_diagnostics.size(); i++) {
        step_s.push_back(scaled.pass_diagnostics[i].step_dx);
        step_u.push_back(unit.pass_diagnostics[i].step_dx);
        spacing_s.push_back(scaled.pass_diagnostics[i].max_spacing);
        spacing_u.push_back(unit.pass_diagnostics[i].max_spacing);
    }
    EXPECT_LT(max_length_deviation(step_s, step_u, k), kScaleRelTol) << "pass_diagnostics.step_dx";
    EXPECT_LT(max_length_deviation(spacing_s, spacing_u, k), kScaleRelTol) << "pass_diagnostics.max_spacing";
}

/**
 * Wall points the march placed on a prescribed contour must lie on that contour as given,
 * in its own units. The first wall point (the start line's) is skipped: its station comes
 * from a fixed-point iteration onto the contour, not from radius_at() directly.
 */
void expect_wall_on_contour(const MocResult& result, const NozzleProfile& contour, double k) {
    double worst = 0.0;
    for (size_t i = 1; i < result.net.wall_x.size(); i++) {
        const double x = result.net.wall_x[i];
        ASSERT_GE(x, contour.x_min());
        ASSERT_LE(x, contour.x_max() * (1.0 + 1e-12));
        const double x_clamped = std::min(x, contour.x_max());
        worst = std::max(worst, std::abs(result.net.wall_y[i] - contour.radius_at(x_clamped)) / k);
    }
    EXPECT_LT(worst, kScaleRelTol) << "wall points must lie on the contour given in physical units";
}

MocOptions analysis_options(MocFlowKind flow, int n, double throat_radius) {
    MocOptions opts;
    opts.flow_type = flow;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::ANALYSIS;
    opts.gamma = 1.4;
    opts.num_characteristics = n;
    opts.geometry.throat_radius = throat_radius;
    return opts;
}

} // namespace

// Axisymmetric conical analysis, as MocDefaultOptionsConvergence.ConicalDefaultThroatConvergesAcrossN
// (test_moc_convergence.cpp) runs it at AR = 4, N = 15: AUTO start line (Kliegel-Levine with
// the wall-consistency correction at r_arc = 0.382), the inverse march, and a contour whose
// arc radius is 0.382 throat radii in both solves.
TEST(MocThroatScaling, ConicalAxisymmetricAnalysisIsScaleInvariant) {
    const double k = kThroatScale;
    const MocFlowKind flow = MocFlowKind::AXISYMMETRIC;

    MocOptions unit_opts = analysis_options(flow, 15, 1.0);
    unit_opts.geometry.downstream_wall_curvature_radius = 0.382;
    unit_opts.nozzle_profile = NozzleProfile::generate_conical_nozzle(4.0, 0.382, 1.0, 15.0, 60);

    MocOptions scaled_opts = analysis_options(flow, 15, k);
    scaled_opts.geometry.downstream_wall_curvature_radius = 0.382; // a multiple of throat_radius
    const NozzleProfile scaled_contour = NozzleProfile::generate_conical_nozzle(4.0, 0.382, k, 15.0, 60);
    scaled_opts.nozzle_profile = scaled_contour;

    const MocResult unit = MocNozzle(unit_opts).solve();
    const MocResult scaled = MocNozzle(scaled_opts).solve();

    EXPECT_EQ(scaled.init_diagnostics.start_line_used, unit.init_diagnostics.start_line_used)
        << "the start line chosen must not depend on the throat radius";
    expect_scaled_solution(scaled, unit, k, flow);
    expect_wall_on_contour(scaled, scaled_contour, k);
    // The contour's exit radius is sqrt(4) throat radii.
    EXPECT_NEAR(scaled.net.wall_y.back(), 2.0 * k, kScaleRelTol * k);
}

// Planar design -> analysis round trip with a centered-fan start line, the path whose
// initial front is built from the throat lip: design a minimum-length contour at
// throat_radius = k, analyse it at throat_radius = k, and compare with the same chain at
// throat_radius = 1. Also checks that the design's contour comes back in the same units the
// analysis reads.
TEST(MocThroatScaling, PlanarFanDesignAnalysisRoundTripIsScaleInvariant) {
    const double k = kThroatScale;
    const MocFlowKind flow = MocFlowKind::PLANAR;

    auto design_at = [&](double throat_radius) {
        MocOptions opts;
        opts.flow_type = flow;
        opts.chemistry = GasChemistry::PERFECT_GAS;
        opts.mode = MocMode::DESIGN_MIN_LENGTH;
        opts.gamma = 1.4;
        opts.theta_max = 15.0 * M_PI / 180.0;
        opts.num_characteristics = 16;
        opts.geometry.throat_radius = throat_radius;
        return MocNozzle(opts).solve();
    };
    auto analyse = [&](const NozzleProfile& contour, double throat_radius) {
        MocOptions opts = analysis_options(flow, 16, throat_radius);
        opts.geometry.downstream_wall_curvature_radius = -1.0; // centered-fan start line
        opts.nozzle_profile = contour;
        return MocNozzle(opts).solve();
    };

    const MocResult unit_design = design_at(1.0);
    const MocResult scaled_design = design_at(k);
    ASSERT_TRUE(unit_design.converged);
    ASSERT_TRUE(scaled_design.converged);
    {
        SCOPED_TRACE("design");
        expect_scaled_solution(scaled_design, unit_design, k, flow);
    }
    // The designed contour starts at the throat lip, (0, throat_radius).
    ASSERT_FALSE(scaled_design.profile.y.empty());
    EXPECT_NEAR(scaled_design.profile.y.front(), k, kScaleRelTol * k);

    const MocResult unit = analyse(unit_design.profile, 1.0);
    const MocResult scaled = analyse(scaled_design.profile, k);
    {
        SCOPED_TRACE("analysis");
        expect_scaled_solution(scaled, unit, k, flow);
        expect_wall_on_contour(scaled, scale_profile(unit_design.profile, k), k);
    }
}

// Rao design generates its own contour from expansion_ratio and length_fraction; with
// throat_radius = k it must be the Rao contour of a throat of radius k, marched to the same
// nondimensional solution. Configuration of MocDesignRao.BasicSolveReachesRecordedCoverage
// (test_moc_rao.cpp).
TEST(MocThroatScaling, RaoDesignIsScaleInvariant) {
    const double k = kThroatScale;
    const MocFlowKind flow = MocFlowKind::AXISYMMETRIC;

    auto design_at = [&](double throat_radius) {
        MocOptions opts;
        opts.flow_type = flow;
        opts.chemistry = GasChemistry::PERFECT_GAS;
        opts.mode = MocMode::DESIGN_RAO;
        opts.gamma = 1.23;
        opts.num_characteristics = 15;
        opts.geometry.throat_radius = throat_radius;
        opts.geometry.expansion_ratio = 5.0;
        opts.geometry.length_fraction = 0.8;
        return MocNozzle(opts).solve();
    };

    const MocResult unit = design_at(1.0);
    const MocResult scaled = design_at(k);
    expect_scaled_solution(scaled, unit, k, flow);
    expect_wall_on_contour(scaled, NozzleProfile::generate_Rao_TOP_nozzle(5.0, k, 0.8), k);
}

// Minimum-length design against the golden contours above: with throat_radius = k the wall
// is the golden wall times k, and the area ratio and exit Mach are the golden values.
TEST(MocThroatScaling, DesignMinLengthContourScalesWithThroatRadius) {
    const double k = kThroatScale;
    for (const GoldenCase& c : golden_cases()) {
        SCOPED_TRACE(c.name);
        MocOptions o;
        o.flow_type = c.flow;
        o.chemistry = GasChemistry::PERFECT_GAS;
        o.mode = MocMode::DESIGN_MIN_LENGTH;
        o.gamma = c.gamma;
        o.theta_max = c.theta_max_deg * M_PI / 180.0;
        o.num_characteristics = c.num_characteristics;
        o.geometry.throat_radius = k;
        const MocResult scaled = MocNozzle(o).solve();
        const MocResult unit = solve_golden(c);

        ASSERT_TRUE(scaled.converged) << scaled.failure.message;
        expect_invariant(scaled.exit_mach, c.exit_mach, "exit_mach");
        expect_invariant(scaled.area_ratio, c.area_ratio, "area_ratio");
        expect_length_scaled(scaled.nozzle_length, c.nozzle_length, k, "nozzle_length");

        ASSERT_EQ(scaled.net.wall_x.size(), c.wall.size());
        for (size_t i = 0; i < c.wall.size(); i++) {
            expect_length_scaled(scaled.net.wall_x[i], c.wall[i].x, k, "wall_x");
            expect_length_scaled(scaled.net.wall_y[i], c.wall[i].y, k, "wall_y");
        }
        expect_scaled_solution(scaled, unit, k, c.flow);
    }
}

// ============================================================
// Throat-radius validation.
// ============================================================

TEST(MocThroatScaling, NonPositiveThroatRadiusIsRejected) {
    for (double throat_radius : {0.0, -1.0}) {
        MocOptions o;
        o.flow_type = MocFlowKind::PLANAR;
        o.chemistry = GasChemistry::PERFECT_GAS;
        o.mode = MocMode::DESIGN_MIN_LENGTH;
        o.gamma = 1.4;
        o.theta_max = 15.0 * M_PI / 180.0;
        o.num_characteristics = 10;
        o.geometry.throat_radius = throat_radius;
        EXPECT_THROW({ MocNozzle nozzle(o); }, std::invalid_argument) << "throat_radius = " << throat_radius;
    }
}

// An analysis contour whose narrowest radius is not the throat radius is in different units
// from the geometry (here a contour in throat radii with a physical throat radius of 2.5): the
// start line would be built at a throat that is not on the wall.
TEST(MocThroatScaling, AnalysisContourMustMatchThroatRadius) {
    MocOptions o = analysis_options(MocFlowKind::AXISYMMETRIC, 15, kThroatScale);
    o.nozzle_profile = NozzleProfile::generate_conical_nozzle(4.0, 0.382, 1.0, 15.0, 60);
    EXPECT_THROW({ MocNozzle nozzle(o); }, std::invalid_argument);

    // The same contour in the geometry's units is accepted.
    o.nozzle_profile = NozzleProfile::generate_conical_nozzle(4.0, 0.382, kThroatScale, 15.0, 60);
    EXPECT_NO_THROW({ MocNozzle nozzle(o); });
}
