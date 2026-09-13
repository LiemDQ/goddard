// Regression tests for the free-function unit processes in moc_unit_processes.hpp/.cpp,
// focused on the generic iterative interior-point solver (solve_interior_point_iterative)
// and intersect_ray_with_wall -- neither is called by either kernel today, so this file is
// their only exercise.
#include "goddard/moc.hpp"
#include "goddard/moc_context.hpp"
#include "goddard/moc_thermo.hpp"
#include "goddard/moc_unit_processes.hpp"
#include "goddard/characteristics.hpp"
#include "goddard/profile.hpp"
#include <cmath>
#include "gtest/gtest.h"

using namespace Goddard;

static constexpr double DEG = M_PI / 180.0;

namespace {

// Builds a MocSolveContext for a perfect-gas solve of the given flow type. The wall profile
// is empty: none of the tests here need a contour (the interior-point solvers never consult
// ctx.wall, and the wall test builds its own profile directly).
MocOptions make_options(MocFlowKind flow_type) {
    MocOptions options;
    options.flow_type = flow_type;
    options.chemistry = GasChemistry::PERFECT_GAS;
    options.gamma = 1.4;
    options.solver_options.abstol = 1e-10;
    return options;
}

// Builds a point with the given Mach number's full thermodynamic state (via MocThermo, as
// the brief prescribes), then the caller sets x/y/theta and calls update_Ks().
CharacteristicPoint make_point(const MocThermo& thermo, double mach, double x, double y, double theta) {
    CharacteristicPoint pt{};
    MocErrorCode err = thermo.set_state_from_mach(pt, mach);
    EXPECT_EQ(err, MocErrorCode::NONE);
    pt.x = x;
    pt.y = y;
    pt.theta = theta;
    pt.update_Ks();
    return pt;
}

// Zero-source lambdas, matching the SourceTermFn signature documented on
// solve_interior_point_iterative: double(const CharacteristicPoint& parent, double new_y).
double zero_source(const CharacteristicPoint&, double) { return 0.0; }

} // namespace

// ── Planar: the iterative solver must reproduce the planar algebraic solver exactly ────────
//
// With zero-source lambdas the iterative solver's compatibility relations are the planar
// ones: theta + nu = K- of the C- parent and theta - nu = K+ of the C+ parent. Its Newton
// residual, (nu - nu_1) + (nu - nu_2) - (theta_1 - theta_2) - (S1 - S2), is exactly those
// two relations with theta eliminated, so for ANY pair of parents -- different flow angles
// included -- the root is the algebraic solver's nu, and the two solvers agree to solver
// tolerance, not merely to discretization order.
TEST(MocUnitProcesses, IterativeSolverReproducesPlanarAlgebraic) {
    MocOptions options = make_options(MocFlowKind::PLANAR);
    NozzleProfile wall;
    MocThermo thermo = MocThermo::perfect_gas(options.gamma);
    MocLog log;
    MocSolveContext ctx{options, wall, thermo, log};

    // Different flow angles at the two parents: theta_1 - theta_2 enters the residual and
    // must be carried with the right sign.
    CharacteristicPoint c_minus_parent = make_point(thermo, 2.0, 0.0, 0.5, 0.10);
    CharacteristicPoint c_plus_parent = make_point(thermo, 2.2, 0.1, 0.2, 0.05);

    PointResult algebraic = solve_interior_point(ctx, c_minus_parent, c_plus_parent);
    ASSERT_EQ(algebraic.error, MocErrorCode::NONE);

    PointResult iterative = solve_interior_point_iterative(
        ctx, c_minus_parent, c_plus_parent, zero_source, zero_source);
    ASSERT_EQ(iterative.error, MocErrorCode::NONE);

    EXPECT_NEAR(iterative.point.x, algebraic.point.x, 1e-12);
    EXPECT_NEAR(iterative.point.y, algebraic.point.y, 1e-12);
    EXPECT_NEAR(iterative.point.theta, algebraic.point.theta, 1e-12);
    EXPECT_NEAR(iterative.point.nu, algebraic.point.nu, 1e-12);
    EXPECT_NEAR(iterative.point.mach, algebraic.point.mach, 1e-12);
}

// ── Axisymmetric: the iterative solver (fed the axisymmetric source terms) must agree ─────
// with the axisymmetric algebraic solver to second order
//
// The two solvers are not the same discretization here. The algebraic solver
// (solve_interior_point) integrates the axisymmetric source term in dx-form, averaged
// between each parent and the intersection by a predictor-corrector; the iterative solver's
// cminus_source_term/cplus_source_term integrate it in dy-form from each parent to the shared
// new_y, inside a Newton iteration on Mach. Both are second-order discretizations of the
// same compatibility relations (Zucrow & Hoffman ch. 17), so on a smooth flow field their
// disagreement must shrink like h^2 as the parent separation h shrinks: by a factor of ~4
// per halving of h. The cell below has both theta and Mach varying with y, so the
// theta_1 - theta_2 term and both source terms are all exercised with nonzero values.
namespace {

// A small "cell": two parents straddling y = 0.5 by +/- h/2, sampled from a smooth flow field
// (theta and Mach linear in y) so that halving h is a meaningful refinement of the same
// underlying continuous problem, not two unrelated point pairs.
struct CellDifference {
    double dtheta, dnu, dmach;
};

CellDifference axisymmetric_cell_difference(const MocSolveContext& ctx, const MocThermo& thermo, double h) {
    const double y_center = 0.5;
    const double theta_center = 0.08, theta_slope = 0.25;  // rad, rad per unit y
    const double mach_center = 2.10, mach_slope = 0.6;     // per unit y

    const double y_minus = y_center + 0.5 * h;  // C- parent: larger y (see find_pair_partner)
    const double y_plus = y_center - 0.5 * h;   // C+ parent: smaller y

    CharacteristicPoint c_minus_parent = make_point(
        thermo, mach_center + mach_slope * (y_minus - y_center),
        0.0, y_minus, theta_center + theta_slope * (y_minus - y_center));
    CharacteristicPoint c_plus_parent = make_point(
        thermo, mach_center + mach_slope * (y_plus - y_center),
        0.0, y_plus, theta_center + theta_slope * (y_plus - y_center));

    // cminus_source_term/cplus_source_term are overloaded (a two-point dy-form also exists,
    // used by the wall solvers); bind explicitly to the (parent, new_y) overload the template
    // expects, since a template argument cannot deduce a type from an overload set directly.
    auto cminus = [](const CharacteristicPoint& p, double new_y) { return cminus_source_term(p, new_y); };
    auto cplus = [](const CharacteristicPoint& p, double new_y) { return cplus_source_term(p, new_y); };

    PointResult algebraic = solve_interior_point(ctx, c_minus_parent, c_plus_parent);
    PointResult iterative = solve_interior_point_iterative(
        ctx, c_minus_parent, c_plus_parent, cminus, cplus);

    EXPECT_EQ(algebraic.error, MocErrorCode::NONE);
    EXPECT_EQ(iterative.error, MocErrorCode::NONE);

    return CellDifference{
        std::abs(iterative.point.theta - algebraic.point.theta),
        std::abs(iterative.point.nu - algebraic.point.nu),
        std::abs(iterative.point.mach - algebraic.point.mach)
    };
}

} // namespace

TEST(MocUnitProcesses, IterativeSolverAgreesWithAxisymmetricAlgebraicToSecondOrder) {
    MocOptions options = make_options(MocFlowKind::AXISYMMETRIC);
    NozzleProfile wall;
    MocThermo thermo = MocThermo::perfect_gas(options.gamma);
    MocLog log;
    MocSolveContext ctx{options, wall, thermo, log};

    const double h = 0.01;
    CellDifference diff_h = axisymmetric_cell_difference(ctx, thermo, h);
    CellDifference diff_h_half = axisymmetric_cell_difference(ctx, thermo, h / 2.0);

    RecordProperty("dtheta_h", diff_h.dtheta);
    RecordProperty("dnu_h", diff_h.dnu);
    RecordProperty("dmach_h", diff_h.dmach);
    RecordProperty("dtheta_h_half", diff_h_half.dtheta);
    RecordProperty("dnu_h_half", diff_h_half.dnu);

    // Measured at h = 0.01 (perfect gas, gamma = 1.4, the flow field above): dtheta 5.5e-6,
    // dnu 6.3e-6, dmach 1.3e-5; at h/2 they are 1.37e-6, 1.56e-6 (ratios 4.02, 4.03). The
    // bound is a few times the largest measured value.
    constexpr double bound = 5e-5;
    EXPECT_LT(diff_h.dtheta, bound);
    EXPECT_LT(diff_h.dnu, bound);
    EXPECT_LT(diff_h.dmach, bound);

    // Convergence order: halving h shrinks the disagreement by ~4x. Allow a band around 4,
    // since this is an asymptotic rate measured at one finite h, and a floor under the h/2
    // residual so a chance near-zero difference cannot blow up the ratio.
    constexpr double floor = 1e-13;
    double ratio_theta = diff_h.dtheta / std::max(diff_h_half.dtheta, floor);
    double ratio_nu = diff_h.dnu / std::max(diff_h_half.dnu, floor);
    EXPECT_GT(ratio_theta, 3.0) << "dtheta(h)=" << diff_h.dtheta << " dtheta(h/2)=" << diff_h_half.dtheta;
    EXPECT_LT(ratio_theta, 5.5) << "dtheta(h)=" << diff_h.dtheta << " dtheta(h/2)=" << diff_h_half.dtheta;
    EXPECT_GT(ratio_nu, 3.0) << "dnu(h)=" << diff_h.dnu << " dnu(h/2)=" << diff_h_half.dnu;
    EXPECT_LT(ratio_nu, 5.5) << "dnu(h)=" << diff_h.dnu << " dnu(h/2)=" << diff_h_half.dnu;
}

// ── intersect_ray_with_wall ──────────────────────────────────────────────────────────────

TEST(MocUnitProcesses, IntersectRayWithWallHitsConicalContour) {
    NozzleProfile wall = NozzleProfile::generate_conical_nozzle(
        /*area_ratio=*/4.0, /*r_expansion_curve=*/0.0, /*r_throat=*/1.0, /*theta_n=*/15.0, 50);

    // From near the axis, well inside the nozzle, a steeply rising ray (45 deg, steeper than
    // the wall's 15 deg half-angle) is guaranteed to catch the diverging wall within its span.
    const double x0 = 0.5, y0 = 0.3;
    const double angle = 45.0 * DEG;

    std::optional<std::pair<double, double>> hit = intersect_ray_with_wall(x0, y0, angle, wall);
    ASSERT_TRUE(hit.has_value());
    auto [x_hit, y_hit] = *hit;
    ASSERT_GE(x_hit, wall.x_min());
    ASSERT_LE(x_hit, wall.x_max());
    EXPECT_NEAR(y_hit, wall.radius_at(x_hit), 1e-9);
}

TEST(MocUnitProcesses, IntersectRayWithWallMissesWhenPointingAway) {
    NozzleProfile wall = NozzleProfile::generate_conical_nozzle(4.0, 0.0, 1.0, 15.0, 50);

    // A ray starting well inside the nozzle (below the wall everywhere) angled downward: it
    // descends while the wall only ever rises, so it can never catch the contour.
    const double x0 = 0.5, y0 = 0.3;
    const double angle = -30.0 * DEG;

    EXPECT_EQ(intersect_ray_with_wall(x0, y0, angle, wall), std::nullopt);
}

TEST(MocUnitProcesses, IntersectRayWithWallThrowsOnDegenerateProfile) {
    NozzleProfile empty_wall;
    EXPECT_THROW(intersect_ray_with_wall(0.0, 0.5, 0.0, empty_wall), std::runtime_error);

    NozzleProfile single_point_wall;
    single_point_wall.x = {0.0};
    single_point_wall.y = {1.0};
    EXPECT_THROW(intersect_ray_with_wall(0.0, 0.5, 0.0, single_point_wall), std::runtime_error);
}
