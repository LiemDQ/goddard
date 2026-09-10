#include "goddard/moc.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/profile.hpp"
#include <cmath>
#include <string>
#include "gtest/gtest.h"

using namespace Goddard;

static constexpr double DEG = M_PI / 180.0;

// Helper: create a MocNozzle configured for axisymmetric perfect-gas Rao design.
// downstream_wall_curvature_radius is left at its default (0.382, positive), so
// generate_initial_data_line() takes the Kliegel-Levine transonic path -- the
// combination documented as recommended for DESIGN_RAO/ANALYSIS.
static MocNozzle make_rao_solver(
    double gamma, double expansion_ratio, double length_fraction, int num_chars)
{
    MocOptions opts;
    opts.flow_type = MocFlowKind::AXISYMMETRIC;
    opts.chemistry = GasChemistry::PERFECT_GAS;
    opts.mode = MocMode::DESIGN_RAO;
    opts.gamma = gamma;
    opts.num_characteristics = num_chars;
    opts.geometry.throat_radius = 1.0;
    opts.geometry.expansion_ratio = expansion_ratio;
    opts.geometry.length_fraction = length_fraction;
    return MocNozzle(opts);
}

// Flipped 2026-09-09 (D.md item 2 / Addendum): DESIGN_RAO is seeded by the
// Kliegel-Levine start line and, before Package A/B, shared the unresolved axisymmetric
// marching failure that blocked conical ANALYSIS at AR >= 4 (see
// instructions/moc_convergence_roadmap.md) -- this configuration did not reach the exit
// plane and was asserted on a recorded coverage floor (0.8195) rather than `converged`,
// "for the same reason... a boolean that flips... tells you less than the number that
// drifted." Package A's wall-consistent KL line and Package B's inverse march (the
// default for DESIGN_RAO) together remove that structural failure, so the coverage-floor
// pattern -- which "exists only because nothing converged" (D.md item 2) -- is replaced
// with an assertion of convergence.
TEST(MocDesignRao, BasicSolveReachesRecordedCoverage) {
    auto solver = make_rao_solver(1.23, 5.0, 0.8, 15);
    auto result = solver.solve();

    RecordProperty("exit_coverage", std::to_string(result.exit_coverage));
    RecordProperty("exit_mach", std::to_string(result.exit_mach));
    RecordProperty("failure_code", std::string(to_string(result.failure.code)));

    // Honesty: converged and "no failure code" must agree, whichever way they go.
    EXPECT_EQ(result.converged, result.failure.code == MocErrorCode::NONE)
        << "converged and failure.code disagree: " << result.failure.message;

    ASSERT_TRUE(result.converged) << to_string(result.failure.code)
        << " -- " << result.failure.message;
    EXPECT_TRUE(result.reached_exit_plane);
    EXPECT_GT(result.exit_mach, 1.0);
    EXPECT_GT(result.area_ratio, 1.0);
    EXPECT_GE(result.exit_coverage, 0.99)
        << "DESIGN_RAO AR=5/Lf=0.8 should now reach (near-)complete exit-radius coverage";
}

// ============================================================
// D.md item 5: Rao contours must outperform a conical nozzle of the same length and
// area ratio -- the whole point of the Rao/thrust-optimized-parabola shape (Rao 1958;
// see generate_Rao_TOP_nozzle's references). This is "the originally-planned... test"
// noted in the project memory as dropped during the pre-fix debugging sessions because
// no (AR, Lf) combination reliably converged; Packages A and B are what make it possible
// now.
// ============================================================

namespace {

// A conical nozzle (throat arc r_expansion_curve=0.382, matching DESIGN_RAO's default
// throat) with the same area ratio as `target` but half-angle theta_n chosen so its
// length matches target_length -- the fair-comparison contour item 5 asks for ("same
// length and area ratio"), rather than the fixed-15-degree cone
// generate_conical_nozzle's theta_n defaults to. Nozzle length decreases monotonically
// with theta_n at fixed area ratio (a steeper cone reaches the same exit area sooner),
// so bisection on theta_n over (1, 45) degrees is well posed for the (AR, Lf) grid this
// test uses, where Lf <= 1.0 means the Rao contour is never longer than the reference
// 15-degree cone.
NozzleProfile conical_nozzle_matching_length(double area_ratio, double target_length, size_t n_points) {
    double lo_deg = 1.0, hi_deg = 45.0;
    for (int i = 0; i < 60; i++) {
        const double mid_deg = 0.5 * (lo_deg + hi_deg);
        NozzleProfile p = NozzleProfile::generate_conical_nozzle(area_ratio, 0.382, 1.0, mid_deg, n_points);
        if (p.length() > target_length) lo_deg = mid_deg; else hi_deg = mid_deg;
    }
    return NozzleProfile::generate_conical_nozzle(area_ratio, 0.382, 1.0, 0.5 * (lo_deg + hi_deg), n_points);
}

} // namespace

TEST(MocDesignRao, CfExceedsConicalOfSameLengthAndAreaRatio) {
    for (double gamma : {1.23, 1.4}) {
        for (double ar : {5.0, 10.0, 20.0}) {
            for (double lf : {0.6, 0.8, 1.0}) {
                const std::string tag = "_AR" + std::to_string(static_cast<int>(ar))
                    + "_Lf" + std::to_string(lf) + "_g" + std::to_string(gamma);

                MocResult rao = make_rao_solver(gamma, ar, lf, 31).solve();
                ASSERT_TRUE(rao.converged) << tag << ": " << to_string(rao.failure.code)
                    << " -- " << rao.failure.message;
                EXPECT_TRUE(rao.reached_exit_plane) << tag;

                ThrustCoefficient rao_cf = compute_thrust_coefficient(rao, MocFlowKind::AXISYMMETRIC);
                EXPECT_GT(rao_cf.Cf_vacuum, 0.0) << tag << ": Rao Cf_vacuum must be physical";
                EXPECT_NEAR(rao_cf.momentum_thrust + rao_cf.pressure_thrust, rao_cf.Cf_vacuum,
                            std::max(1e-8, 1e-10 * rao_cf.Cf_vacuum))
                    << tag << ": momentum + pressure thrust must sum to the vacuum coefficient";

                // The comparison cone: same throat, same area ratio, same length as the
                // Rao contour actually achieved (rao.nozzle_length/rao.area_ratio, not the
                // requested targets -- the achieved area ratio can differ slightly from
                // the request, see MocKlInitConvergence's area_ratio tolerances).
                NozzleProfile cone = conical_nozzle_matching_length(rao.area_ratio, rao.nozzle_length, 200);
                ASSERT_NEAR(cone.length(), rao.nozzle_length, 1e-3 * rao.nozzle_length)
                    << tag << ": comparison-cone length bisection did not converge";

                MocOptions cone_opts;
                cone_opts.flow_type = MocFlowKind::AXISYMMETRIC;
                cone_opts.chemistry = GasChemistry::PERFECT_GAS;
                cone_opts.mode = MocMode::ANALYSIS;
                cone_opts.gamma = gamma;
                cone_opts.num_characteristics = 31;
                cone_opts.geometry.throat_radius = 1.0;
                cone_opts.geometry.downstream_wall_curvature_radius = 0.382;
                cone_opts.nozzle_profile = cone;
                // march_scheme, start_line left at their AUTO defaults, as in
                // MocDefaultOptionsConvergence.ConicalDefaultThroatConvergesAcrossN
                // (test_moc_convergence.cpp): AUTO -> INVERSE for axisymmetric ANALYSIS.

                MocResult cone_result = MocNozzle(cone_opts).solve();
                ASSERT_TRUE(cone_result.converged) << tag << " (comparison cone): "
                    << to_string(cone_result.failure.code) << " -- " << cone_result.failure.message;
                EXPECT_TRUE(cone_result.reached_exit_plane) << tag << " (comparison cone)";

                ThrustCoefficient cone_cf = compute_thrust_coefficient(cone_result, MocFlowKind::AXISYMMETRIC);
                EXPECT_GT(cone_cf.Cf_vacuum, 0.0) << tag << ": conical Cf_vacuum must be physical";

                RecordProperty("rao_Cf_vacuum" + tag, std::to_string(rao_cf.Cf_vacuum));
                RecordProperty("conical_Cf_vacuum" + tag, std::to_string(cone_cf.Cf_vacuum));
                RecordProperty("cone_theta_n_deg" + tag, std::to_string(cone.max_theta() / DEG));

                EXPECT_GT(rao_cf.Cf_vacuum, cone_cf.Cf_vacuum)
                    << tag << ": Rao contour should out-thrust a conical nozzle of the same "
                    << "length and area ratio (Rao " << rao_cf.Cf_vacuum << " vs conical "
                    << cone_cf.Cf_vacuum << ")";
            }
        }
    }
}
