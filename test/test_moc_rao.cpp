#include "goddard/moc.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/profile.hpp"
#include <string>
#include "gtest/gtest.h"

using namespace Goddard;

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

// DESIGN_RAO is seeded by the Kliegel-Levine start line and shares the unresolved
// axisymmetric marching failure that blocks conical ANALYSIS at AR >= 4 (see
// instructions/moc_convergence_roadmap.md). This configuration does not currently reach
// the exit plane.
//
// It is asserted rather than skipped. A GTEST_SKIP on non-convergence makes the case
// disappear from the suite exactly when it regresses, which is the moment its signal is
// worth most -- a solver that gets worse and a solver that stops being tested look
// identical in the summary line. Instead the test pins what the solve actually achieves
// today: it must not throw, it must report its failure honestly, and it must still cover
// at least the fraction of the exit radius recorded below. Raise the floor when the
// coverage improves; a drop fails here rather than going quiet.
//
// Measured 2026-09-02 (after the Kliegel-Levine critical-velocity-ratio and wall-anchoring
// fixes): exit_coverage 0.8195, exit_mach 2.977, NON_DOWNSTREAM_POINT.
TEST(MocDesignRao, BasicSolveReachesRecordedCoverage) {
    auto solver = make_rao_solver(1.23, 5.0, 0.8, 15);
    auto result = solver.solve();

    RecordProperty("exit_coverage", std::to_string(result.exit_coverage));
    RecordProperty("exit_mach", std::to_string(result.exit_mach));
    RecordProperty("failure_code", std::string(to_string(result.failure.code)));

    // Honesty: converged and "no failure code" must agree, whichever way they go.
    EXPECT_EQ(result.converged, result.failure.code == MocErrorCode::NONE)
        << "converged and failure.code disagree: " << result.failure.message;

    // Whatever the march achieves, the part it did solve has to be physical.
    EXPECT_GT(result.exit_mach, 1.0);
    EXPECT_GT(result.area_ratio, 1.0);

    EXPECT_GE(result.exit_coverage, 0.80)
        << "DESIGN_RAO AR=5/Lf=0.8 coverage regressed below its recorded floor (was 0.8195): "
        << result.failure.message;
}
