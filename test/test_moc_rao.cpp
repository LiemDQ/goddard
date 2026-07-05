#include "goddard/moc.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/profile.hpp"
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

// NOTE: DESIGN_RAO's Kliegel-Levine-seeded kernel still has a known, unresolved
// geometric-consistency issue beyond the bootstrap-ordering bug fixed alongside
// this test (see CharacteristicNet::add_initial_data_line and
// MocNozzle::sort_plus_edges_by_proximity): many (area_ratio, length_fraction)
// combinations either fail to converge or converge with a wall contour that
// doesn't actually reach the target area ratio. Until that's fixed, this file
// only smoke-tests the one configuration confirmed to converge to a
// self-consistent result; it does not yet attempt the Rao-vs-conical thrust
// comparison.
TEST(MocDesignRao, BasicSolveConverges) {
    auto solver = make_rao_solver(1.23, 5.0, 0.8, 15);
    auto result = solver.solve();

    EXPECT_TRUE(result.converged);
    EXPECT_GT(result.exit_mach, 1.0);
    EXPECT_GT(result.area_ratio, 1.0);
}
