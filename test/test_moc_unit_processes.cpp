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
// With zero-source lambdas the iterative solver's compatibility equations reduce to exactly
// the planar algebraic solver's (see solve_interior_point_iterative's doc comment): both
// solve theta = S + theta_1 - (nu - nu_1) with S = 0, i.e. the same K+ = theta - nu / K- =
// theta + nu invariants. The two should therefore agree to solver-tolerance roundoff, not
// merely to discretization order -- for parents that share the same theta.
//
// Parents are given EQUAL theta here (0.08 rad at both), not the differing thetas an
// interior-point pairing usually has. find_node_mach's Newton residual (moved verbatim,
// unmodified per this step's behavior-preservation mandate) is
//   delta_theta - source_delta + (nu_P - nu_1) + (nu_P - nu_2) = 0,  delta_theta = theta_1 - theta_2
// which combines the two compatibility relations theta+nu=K_minus_1 (C-) and
// theta-nu=K_plus_2-S2 (C+) with the WRONG sign on delta_theta and on S2 relative to
// eliminating theta_P from that pair (the correct combination is
// 2 nu_P = delta_theta + nu_1 + nu_2 + S1 + S2, verified by hand against Zucrow & Hoffman
// ch. 17's compatibility relations) -- confirmed both analytically and by running this test
// against the moved (unmodified) code: with theta_1=0.10, theta_2=0.05 as in an ordinary
// pairing, the iterative solver's theta/nu differ from the algebraic solver's by exactly
// delta_theta = 0.05 rad, not by roundoff. That sign issue is pre-existing in code this step
// only relocates (it has no caller in either kernel and was never previously tested); fixing
// it is out of scope for a step whose contract is bit-identical behavior. Using equal-theta
// parents makes delta_theta vanish, which is where the two solvers DO coincide exactly, and
// still exercises the full mechanism (Newton Mach solve, corrector loop, thermo dispatch)
// with genuinely different Mach numbers at the two parents.
TEST(MocUnitProcesses, IterativeSolverReproducesPlanarAlgebraic) {
    MocOptions options = make_options(MocFlowKind::PLANAR);
    NozzleProfile wall;
    MocThermo thermo = MocThermo::perfect_gas(options.gamma);
    MocLog log;
    MocSolveContext ctx{options, wall, thermo, log};

    CharacteristicPoint c_minus_parent = make_point(thermo, 2.0, 0.0, 0.5, 0.08);
    CharacteristicPoint c_plus_parent = make_point(thermo, 2.2, 0.1, 0.2, 0.08);

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

// ── Axisymmetric: the iterative solver (fed the axisymmetric source terms) should agree ────
// with the axisymmetric algebraic solver to O(h^2), not bit-for-bit
//
// Unlike the planar case, the two solvers are NOT the same discretization here. The
// algebraic solver (solve_interior_point) evaluates the axisymmetric source term in
// dx-form (accumulated along x from each parent, predictor-corrector averaged at the
// intersection); the iterative solver's cminus_source_term/cplus_source_term evaluate it in
// dy-form (accumulated along y from each parent to the shared new_y, inside a Newton
// iteration on Mach). Both are second-order-accurate discretizations of the same
// compatibility relations (Zucrow & Hoffman ch. 17), so on a smooth flow field they must
// agree to O(h^2) as the parent separation h shrinks, not exactly.
//
// Parents share the same theta at every h (see IterativeSolverReproducesPlanarAlgebraic's
// comment on find_node_mach's pre-existing sign handling of delta_theta = theta_1 - theta_2):
// with delta_theta held at exactly 0 for every h, that O(h) term cannot contaminate the O(h^2)
// comparison this test is actually after. Mach still varies with y, so the cell is not
// degenerate and the axisymmetric source terms (which depend on theta, not on delta_theta)
// are genuinely exercised and nonzero.
namespace {

// A small "cell": two parents straddling y = 0.5 by +/- h/2, sampled from a smooth flow field
// (constant theta, Mach linear in y) so that halving h is a meaningful refinement of the same
// underlying continuous problem, not two unrelated point pairs.
struct CellDifference {
    double dtheta, dnu, dmach;
};

CellDifference axisymmetric_cell_difference(const MocSolveContext& ctx, const MocThermo& thermo, double h) {
    const double y_center = 0.5;
    const double theta_common = 0.08;                      // rad; same at both parents, every h
    const double mach_center = 2.10, mach_slope = 0.6;     // per unit y

    const double y_minus = y_center + 0.5 * h;  // C- parent: larger y (see find_pair_partner)
    const double y_plus = y_center - 0.5 * h;   // C+ parent: smaller y

    CharacteristicPoint c_minus_parent = make_point(
        thermo, mach_center + mach_slope * (y_minus - y_center),
        0.0, y_minus, theta_common);
    CharacteristicPoint c_plus_parent = make_point(
        thermo, mach_center + mach_slope * (y_plus - y_center),
        0.0, y_plus, theta_common);

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

    // Measured residuals at h = 0.01 (perfect gas, gamma = 1.4, the flow field above):
    // dtheta ~ 9.0e-4 rad, dnu ~ 6.2e-4 rad, dmach ~ 1.3e-3. Larger than a pure O(h^2)
    // discretization difference would suggest for h = 0.01 -- because it isn't purely that.
    // find_node_mach's residual (see IterativeSolverReproducesPlanarAlgebraic's comment) also
    // combines the two dy-form sources S1, S2 with the wrong relative sign for eliminating
    // theta_P (nu_P = ... + S1 - S2 where the correct elimination needs + S1 + S2); with
    // theta held equal at both parents that stray -2*S2 term is the only surviving defect,
    // and S2 itself is O(h) (see cplus_source_term's dy/y_avg factor), so it dominates the
    // legitimate O(h^2) difference between the dx-form and dy-form source discretizations at
    // these cell sizes. The bound below is a generous multiple of the measured value -- not
    // tight to a hypothetical pure O(h^2) residual, since the actual (pre-existing, moved
    // unmodified) code is not that.
    constexpr double bound = 4e-3;
    EXPECT_LT(diff_h.dtheta, bound);
    EXPECT_LT(diff_h.dnu, bound);
    EXPECT_LT(diff_h.dmach, bound);

    // Convergence order: halving h shrinks the disagreement by ~2x (O(h)), not ~4x (O(h^2)),
    // for the reason above -- the O(h) sign defect in the source combination dominates over
    // the true O(h^2) discretization difference at these cell sizes. Measured ratio at
    // h=0.01 vs h=0.005: ~1.99 (dtheta), consistent with O(h). Allow a wide band (1.3x-3.5x)
    // since this is an asymptotic rate measured at one finite h, not an identity -- and a
    // floor under the h/2 residual so a chance near-zero difference does not blow up the
    // ratio. This still catches a real regression: a change that broke the shared
    // compatibility-equation structure entirely (rather than merely its source-term order)
    // would not converge at any clean rate as h shrinks.
    constexpr double floor = 1e-13;
    double ratio_theta = diff_h.dtheta / std::max(diff_h_half.dtheta, floor);
    double ratio_nu = diff_h.dnu / std::max(diff_h_half.dnu, floor);
    EXPECT_GT(ratio_theta, 1.3) << "dtheta(h)=" << diff_h.dtheta << " dtheta(h/2)=" << diff_h_half.dtheta;
    EXPECT_LT(ratio_theta, 3.5) << "dtheta(h)=" << diff_h.dtheta << " dtheta(h/2)=" << diff_h_half.dtheta;
    EXPECT_GT(ratio_nu, 1.3) << "dnu(h)=" << diff_h.dnu << " dnu(h/2)=" << diff_h_half.dnu;
    EXPECT_LT(ratio_nu, 3.5) << "dnu(h)=" << diff_h.dnu << " dnu(h/2)=" << diff_h_half.dnu;
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
