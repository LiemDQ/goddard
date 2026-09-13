#pragma once
#include <algorithm>
#include <cmath>
#include <optional>
#include <string_view>
#include <utility>
#include "goddard/moc.hpp"
#include "goddard/moc_context.hpp"
#include "goddard/moc_thermo.hpp"
#include "goddard/characteristics.hpp"
#include "goddard/profile.hpp"

namespace Goddard {

/**
 * A CharacteristicPoint together with the error code (if any) produced while
 * computing it. Unit-process solvers return this instead of a bare
 * CharacteristicPoint so a numerical failure is reported as data -- instead of
 * being laundered into the point's fields (e.g. a sentinel Mach number) or
 * thrown as an exception.
 */
struct PointResult {
    CharacteristicPoint point;
    MocErrorCode error = MocErrorCode::NONE;
};

/**
 * Validity of a computed point: finite fields, supersonic, non-negative Prandtl-Meyer angle,
 * and (by default) non-negative flow angle.
 *
 * The flow-angle requirement is what the chain ladder (`DirectMarch`, minimum-length design)
 * uses to catch a march that has gone wrong. It is optional because a negative flow angle is
 * physically admissible in the analysis of an arbitrary contour (any locally converging
 * wall) and because, wherever the exact angle is zero -- the uniform exit region of a
 * minimum-length nozzle, the axis neighbourhood -- a computed angle is discretization noise
 * of either sign, which a tolerance of the solver's abstol cannot absorb.
 *
 * @param pt  Point to validate.
 * @param tol Absolute tolerance for the nu/theta non-negativity checks (axis points
 *            legitimately pin theta to exactly 0.0, and interior points can carry
 *            small negative roundoff).
 * @return MocErrorCode::NONE if the point is valid, else the first violated
 *         condition in priority order: nu, theta, mach, finiteness.
 */
MocErrorCode check_point_validity(const CharacteristicPoint& pt, double tol,
                                  bool require_nonnegative_theta = true);

// -- Classical unit processes (pure free functions) --------------------------------------
// Every process below takes the whole solve's fixed inputs (options, wall contour, thermo
// dispatch, log) as a MocSolveContext, and the point(s) it is solving from. These are the
// chain ladder's (`DirectMarch`, minimum-length design) unit processes; the reference-plane
// march's (`InverseMarch`, analysis and Rao design) own interior/axis/wall processes are
// siblings that live with that kernel (moc_inverse_march.cpp) but share the source-term
// helpers and axis corrector below.

/** Interior point: dispatches on ctx.options.flow_type to the planar or axisymmetric solve. */
PointResult solve_interior_point(
    const MocSolveContext& ctx,
    const CharacteristicPoint& c_minus_parent,
    const CharacteristicPoint& c_plus_parent);

/**
 * Compute flow properties at centerline for axisymmetric flow.
 *
 * A special method is needed because the axisymmetric compatibility
 * equations have a singularity on the axis of rotation.
 */
PointResult solve_axis_point(
    const MocSolveContext& ctx,
    const CharacteristicPoint& off_axis_parent);

/**
 * The first point of the initial expansion fan's marched data line is a special case, as it
 * lies on the axis but is assigned a nonzero theta. This is because the calculations are
 * started on the characteristic line along which theta is known.
 *
 * This leads to a small physical inconsistency, but it is necessary to bootstrap the
 * downstream marching. (The planar and axisymmetric bootstraps are identical -- both simply
 * carry the expansion point's theta and nu onto the axis -- so there is one body, not two.)
 */
PointResult solve_fan_axis_point(
    const MocSolveContext& ctx,
    const CharacteristicPoint& expansion_point);

/**
 * Minimum-length design only: compute flow and position at the next wall point from the
 * previous wall point and its interior parent, given the wall angle theta_wall (from the
 * DESIGN_MIN_LENGTH theta schedule). Check PointResult::error for a numerical failure.
 */
PointResult solve_wall_point_design(
    const MocSolveContext& ctx,
    const CharacteristicPoint& interior_parent,
    const CharacteristicPoint& previous_wall_point,
    double theta_wall);

/**
 * Source term for axisymmetric flow along the C+ characteristic, two-point (dy) form: used by
 * the wall solvers, which know both endpoints of the step.
 */
double cplus_source_term(const CharacteristicPoint& p1, const CharacteristicPoint& p3);

/**
 * Source term for axisymmetric flow along the C+ characteristic, evaluated at the parent with
 * a prescribed new y (kept for the generic iterative interior solver below).
 */
double cplus_source_term(const CharacteristicPoint& p1, double new_y);

/**
 * Source term for axisymmetric flow along the C- characteristic, evaluated at the parent with
 * a prescribed new y (kept for the generic iterative interior solver below).
 */
double cminus_source_term(const CharacteristicPoint& p1, double new_y);

/**
 * The dx-form axisymmetric source term shared by the axisymmetric interior unit process and
 * the inverse march's interior process: sin(mu_avg) sin(theta_avg) / (y_avg cos(char_angle_avg)).
 *
 * Use only where the caller's existing expression multiplies and divides in exactly this
 * order; where it does not (e.g. the inverse interior process's mirrored-foot branch, which
 * evaluates sin(theta)/y at the new point instead of averaging), the expression stays inline
 * at the call site -- bit-identity matters more than sharing the helper there.
 */
double axisymmetric_source(double mu_avg, double theta_avg, double y_avg, double char_angle_avg);

/**
 * The axis corrector shared by solve_axis_point and the inverse march's axis process: the
 * dy/y_avg = -2 limit of the axisymmetric C- source term as the far point approaches the axis
 * (an algebraic identity there, not a small-perturbation approximation).
 *
 * The |sin(theta_avg - mu_avg)| > 1e-12 guard is carried over from the inverse-march copy;
 * the chain-ladder copy this now also covers never exercised it (mu sits in 30-60 deg for
 * every reachable state), so adding it here does not change the dump.
 */
double axis_source_correction(const CharacteristicPoint& foot, const CharacteristicPoint& axis_point);

/**
 * Find where a ray from (x, y) travelling at `angle` intersects `wall`.
 *
 * @return The (x, y) of the intersection, or std::nullopt if the ray does not meet the
 *         profile within its span.
 * @throws std::runtime_error if `wall` has fewer than two points.
 */
std::optional<std::pair<double, double>> intersect_ray_with_wall(
    double x, double y, double angle, const NozzleProfile& wall);

/**
 * Unwrap a PointResult, throwing ConvergenceError with a formatted message identifying
 * `context` and the failing point/error code if the unit process that produced it failed.
 *
 * Used wherever a unit-process failure must be communicated by exception rather than as data
 * -- the initial data line (before the kernel starts marching) and the inverse march's
 * initial-front construction, both of which run inside solve()'s existing initialization
 * exception boundary.
 */
CharacteristicPoint unwrap_or_throw(const PointResult& result, std::string_view context);

/**
 * Iteratively find the Mach number of a node from two upstream states and the net source-term
 * difference between them (used by solve_interior_point_iterative below). Not called by
 * either kernel today -- kept as a generic building block for an arbitrary-source interior
 * process (see solve_interior_point_iterative).
 */
double find_node_mach(
    const MocThermo& thermo,
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2,
    double source_delta,
    double mach_guess,
    double abstol,
    MocLog& log);

/**
 * Generic iterative interior-point process for an arbitrary characteristic source term,
 * supplied as two callables (one per family) of the form
 * `double(const CharacteristicPoint& parent, double new_y)`.
 *
 * This is the general compatibility-equation solver the axisymmetric interior process
 * specializes: passing zero-source lambdas reproduces the planar algebraic solve, and passing
 * cminus_source_term/cplus_source_term reproduces (to the discretization's own order, not
 * bit-for-bit -- it is a different iteration scheme) the axisymmetric predictor-corrector.
 * Not called by either kernel; kept as a documented general-purpose building block, with a
 * regression test (test/test_moc_unit_processes.cpp) standing in for the kernel call neither
 * exercises it today.
 */
template <typename CMinusSource, typename CPlusSource>
PointResult solve_interior_point_iterative(
    const MocSolveContext& ctx,
    const CharacteristicPoint& c_minus_parent, const CharacteristicPoint& c_plus_parent,
    CMinusSource&& cminus_source, CPlusSource&& cplus_source)
{
    // General procedure:
    // 1. Receive theta and nu from the two incoming characteristics
    // 2. Estimate characteristic slopes at parent nodes
    // 3. Intersect straight line characteristics
    // 4. Compute source terms (if applicable)
    // 5. Solve nonlinear equation in V using upstream states only
    // 6. Get theta from either compatibility equation
    // 7. Correct characteristic slope by averaging with endpoint
    // 8. Correct V by averaging with endpoint
    // 9. Recompute theta and source terms at corrected location
    // 10. [optional] One additional pass if source term changed significantly
    const CharacteristicPoint& p1 = c_minus_parent;
    const CharacteristicPoint& p2 = c_plus_parent;
    CharacteristicPoint p3{};
    double c_minus_angle = p1.theta - p1.mu;
    double c_plus_angle = p2.theta + p2.mu;

    auto [x, y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
    p3.x = x;
    p3.y = y;
    double S1 = cminus_source(p1, y);
    double S2 = cplus_source(p2, y);

    double mach = find_node_mach(ctx.thermo, p1, p2, S1 - S2, 0.5 * (p1.mach + p2.mach),
                                  ctx.options.solver_options.abstol, ctx.log);
    if (mach < 1.0) {
        return {p3, MocErrorCode::PM_INVERSION_FAILED};
    }
    MocErrorCode thermo_error = ctx.thermo.set_state_from_mach(p3, mach);
    if (thermo_error != MocErrorCode::NONE) {
        return {p3, thermo_error};
    }

    // from compatibility equation.
    p3.theta = S1 + p1.theta - (p3.nu - p1.nu);
    p3.K_plus = p3.theta - p3.nu;
    p3.K_minus = p3.theta + p3.nu;

    // corrector step
    const int max_iters = 3;

    for (int i = 0; i < max_iters; i++) {
        c_minus_angle = average_cminus_angle(p1, p3);
        c_plus_angle = average_cplus_angle(p2, p3);
        auto [new_x, new_y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
        p3.x = new_x;
        p3.y = new_y;
        double S1_new = cminus_source(p1, new_y);
        double S2_new = cplus_source(p2, new_y);
        double mach_new = find_node_mach(ctx.thermo, p1, p2, S1_new - S2_new, p3.mach,
                                          ctx.options.solver_options.abstol, ctx.log);
        if (mach_new < 1.0) {
            return {p3, MocErrorCode::PM_INVERSION_FAILED};
        }
        thermo_error = ctx.thermo.set_state_from_mach(p3, mach_new);
        if (thermo_error != MocErrorCode::NONE) {
            return {p3, thermo_error};
        }

        p3.theta = S1_new + p1.theta - (p3.nu - p1.nu);
        p3.K_plus = p3.theta - p3.nu;
        p3.K_minus = p3.theta + p3.nu;
        double residual = std::max(std::abs(S1_new - S1), std::abs(S2_new - S2));
        if (residual < ctx.options.solver_options.abstol) break;
        S1 = S1_new;
        S2 = S2_new;

        if (i == max_iters - 1 && residual >= ctx.options.solver_options.abstol) {
            ctx.log.warning("Iterative interior solver did not converge."
                "Residual={} after {} iterations at ({}, {}).", residual, max_iters, p3.x, p3.y);
        }
    }

    MocErrorCode validity = check_point_validity(p3, ctx.options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        return {p3, validity};
    }
    return {p3, MocErrorCode::NONE};
}

} // namespace Goddard
