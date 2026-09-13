#include <cmath>
#include <format>
#include <stdexcept>
#include "goddard/moc_unit_processes.hpp"
#include "goddard/error.hpp"
#include "goddard/prandtlmeyer.hpp"

namespace Goddard {

MocErrorCode check_point_validity(const CharacteristicPoint& pt, double tol,
                                  bool require_nonnegative_theta) {
    if (pt.nu < -tol) return MocErrorCode::NEGATIVE_NU;
    if (require_nonnegative_theta && pt.theta < -tol) return MocErrorCode::NEGATIVE_THETA;
    if (pt.mach < 1.0) return MocErrorCode::SUBSONIC_MACH;
    if (!std::isfinite(pt.x) || !std::isfinite(pt.y) || !std::isfinite(pt.theta) ||
        !std::isfinite(pt.nu) || !std::isfinite(pt.mach) || !std::isfinite(pt.mu)) {
        return MocErrorCode::NONFINITE_VALUE;
    }
    return MocErrorCode::NONE;
}

namespace {

// Planar algebraic special case (L = M = 0 in the axisymmetric compatibility equations).
PointResult solve_interior_point_planar(
    const MocSolveContext& ctx,
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2)
{
    CharacteristicPoint p3{};
    p3.theta = 0.5*(p1.K_minus + p2.K_plus);
    double mach_guess = 0.5*(p1.mach + p2.mach);

    // Position (x, y) is computed after the thermo update below, so a thermo
    // failure here leaves p3.x/p3.y unset; fall back to a parent's coordinates so
    // the returned point still carries a meaningful location for diagnostics.
    MocErrorCode thermo_error = ctx.thermo.set_state_from_nu(
        p3, 0.5*(p1.K_minus - p2.K_plus), mach_guess);
    if (thermo_error != MocErrorCode::NONE) {
        p3.x = p1.x;
        p3.y = p1.y;
        return {p3, thermo_error};
    }
    p3.update_Ks();

    // compute intersection point
    // assumimg characteristics are straight lines.
    double c_minus_angle = average_cminus_angle(p1, p3);
    double c_plus_angle = average_cplus_angle(p2, p3);

    auto [x, y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
    p3.x = x;
    p3.y = y;

    // validity checks
    if (p3.x < p1.x || p3.x < p2.x) {
        return {p3, MocErrorCode::NON_DOWNSTREAM_POINT};
    }

    MocErrorCode validity = check_point_validity(p3, ctx.options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        return {p3, validity};
    }

    return {p3, MocErrorCode::NONE};
}

// Axisymmetric flow: the (theta, nu) formulation, consistent with the planar, wall, and axis
// solvers. The Riemann invariants K+ = theta - nu (constant along C+) and K- = theta + nu
// (constant along C-) pick up source terms in axisymmetric flow:
//   along C+:  d(theta - nu) = -L dx,   L = sin(mu) sin(theta) / (y cos(theta + mu))
//   along C-:  d(theta + nu) = +M dx,   M = sin(mu) sin(theta) / (y cos(theta - mu))
// With L = M = 0 this reduces exactly to the planar solver. Working in nu (rather than the
// velocity form) keeps it exact for perfect gas, where nu = nu(M); the velocity form would
// require the *true* velocity for cot(mu) dV/V = dnu to hold, but V is stored as M here.
PointResult solve_interior_point_axisymmetric(
    const MocSolveContext& ctx,
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2)
{
    CharacteristicPoint p3{};

    // Predictor: straight characteristics, source terms evaluated at the parents.
    double c_minus_angle = p1.theta - p1.mu;
    double c_plus_angle = p2.theta + p2.mu;

    auto [x, y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
    p3.x = x;
    p3.y = y;

    // NOTE: watch out for singularities when a parent is on the centerline.
    // The averaging of y should prevent issues.
    double L = axisymmetric_source(p2.mu, p2.theta, 0.5*(p2.y+p3.y), p2.theta + p2.mu);
    double M = axisymmetric_source(p1.mu, p1.theta, 0.5*(p1.y+p3.y), p1.theta - p1.mu);

    double K_plus = (p2.theta - p2.nu) - L*(p3.x - p2.x);
    double K_minus = (p1.theta + p1.nu) + M*(p3.x - p1.x);
    p3.theta = 0.5*(K_minus + K_plus);
    MocErrorCode thermo_error = ctx.thermo.set_state_from_nu(
        p3, 0.5*(K_minus - K_plus), 0.5*(p1.mach + p2.mach));
    if (thermo_error != MocErrorCode::NONE) {
        return {p3, thermo_error};
    }

    // Corrector: average the characteristic angles and source terms across the step.
    double c_minus_angle_corrected = average_cminus_angle(p1, p3);
    double c_plus_angle_corrected = average_cplus_angle(p2, p3);
    auto [x_corrected, y_corrected] = characteristic_intersection_with_angle(
        p1, p2,
        c_minus_angle_corrected,
        c_plus_angle_corrected);
    p3.x = x_corrected;
    p3.y = y_corrected;
    double mu1_avg = 0.5*(p1.mu + p3.mu);
    double mu2_avg = 0.5*(p2.mu + p3.mu);
    double theta1_avg = 0.5*(p1.theta + p3.theta);
    double theta2_avg = 0.5*(p2.theta + p3.theta);
    double M_avg = axisymmetric_source(mu1_avg, theta1_avg, 0.5*(p1.y+p3.y), c_minus_angle_corrected);
    double L_avg = axisymmetric_source(mu2_avg, theta2_avg, 0.5*(p2.y+p3.y), c_plus_angle_corrected);

    K_plus = (p2.theta - p2.nu) - L_avg*(p3.x - p2.x);
    K_minus = (p1.theta + p1.nu) + M_avg*(p3.x - p1.x);
    p3.theta = 0.5*(K_minus + K_plus);
    thermo_error = ctx.thermo.set_state_from_nu(p3, 0.5*(K_minus - K_plus), p3.mach);
    if (thermo_error != MocErrorCode::NONE) {
        return {p3, thermo_error};
    }
    p3.update_Ks();

    // validity checks
    if (p3.x < p1.x || p3.x < p2.x) {
        return {p3, MocErrorCode::NON_DOWNSTREAM_POINT};
    }

    MocErrorCode validity = check_point_validity(p3, ctx.options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        return {p3, validity};
    }

    return {p3, MocErrorCode::NONE};
}

// Predictor: planar K+ preservation. Flow solve: compatibility equation along C+ from the
// interior parent (K+ is preserved: theta - nu = const along C+; this only applies for planar
// flow). For axisymmetric flow, K+ is not preserved along C+: solve_wall_point_design applies
// a predictor-corrector on top of this, using the source term correction which requires y
// (computed only once the wall point's position is known).
PointResult solve_wall_flow(
    const MocSolveContext& ctx,
    const CharacteristicPoint& interior_parent,
    double theta_wall)
{
    CharacteristicPoint wall_point{};
    // geometric constraint: the flow at the wall must follow the wall curvature
    wall_point.theta = theta_wall;

    wall_point.K_plus = interior_parent.K_plus;
    MocErrorCode thermo_error = ctx.thermo.set_state_from_nu(
        wall_point,
        wall_point.theta - wall_point.K_plus,
        interior_parent.mach);
    if (thermo_error != MocErrorCode::NONE) {
        return {wall_point, thermo_error};
    }
    wall_point.K_minus = wall_point.theta + wall_point.nu;

    return {wall_point, MocErrorCode::NONE};
}

} // namespace

PointResult solve_interior_point(
    const MocSolveContext& ctx,
    const CharacteristicPoint& c_minus_parent,
    const CharacteristicPoint& c_plus_parent)
{
    switch (ctx.options.flow_type) {
        case MocFlowKind::PLANAR: {
            return solve_interior_point_planar(ctx, c_minus_parent, c_plus_parent);
        }
        case MocFlowKind::AXISYMMETRIC: {
            return solve_interior_point_axisymmetric(ctx, c_minus_parent, c_plus_parent);
        }
        default:
            throw std::runtime_error("Invalid MoC flow type specified.");
    }
}

PointResult solve_axis_point(const MocSolveContext& ctx, const CharacteristicPoint& off_axis_parent) {
    CharacteristicPoint axis_point{};
    axis_point.y = 0.0;
    axis_point.theta = 0.0; // symmetry condition

    switch (ctx.options.flow_type) {
        case MocFlowKind::PLANAR: {

            axis_point.K_minus = off_axis_parent.K_minus;
            MocErrorCode thermo_error = ctx.thermo.set_state_from_nu(
                axis_point,
                axis_point.K_minus - axis_point.theta,
                off_axis_parent.mach);
            if (thermo_error != MocErrorCode::NONE) {
                axis_point.x = off_axis_parent.x;
                return {axis_point, thermo_error};
            }

            axis_point.K_plus = axis_point.theta - axis_point.nu;

            double c_minus_angle = 0.5 * (off_axis_parent.theta + axis_point.theta)
                - 0.5 * (off_axis_parent.mu + axis_point.mu);

            axis_point.x = off_axis_parent.x - off_axis_parent.y/tan(c_minus_angle);
            break;
        }
        case MocFlowKind::AXISYMMETRIC: {

            const double theta = off_axis_parent.theta;
            const double mu = off_axis_parent.mu;
            // Predictor method recommended by Zucrow & Hoffman, ch. 17
            // Use parent point as initial guess for source term.
            // The C- compatibility per dy is d(theta+nu) = [sin(theta)/(M sin(theta-mu))] (dy/y).
            // The descent from the parent to the axis has dy = -y_parent (signed), which makes
            // the source contribution positive (sin(theta-mu) < 0).
            axis_point.K_minus = off_axis_parent.K_minus;

            const double dy = -off_axis_parent.y; // y_axis - y_parent
            double source_pred = sin(theta) / (off_axis_parent.mach * sin(theta - mu))
                * dy / off_axis_parent.y;

            // dtheta + dnu = S => nu_2 - nu_1 = S + theta_1 = S + K_minus_1
            double nu_pred = axis_point.K_minus + source_pred;

            MocErrorCode thermo_error = ctx.thermo.set_state_from_nu(
                axis_point,
                nu_pred,
                off_axis_parent.mach);
            if (thermo_error != MocErrorCode::NONE) {
                axis_point.x = off_axis_parent.x;
                return {axis_point, thermo_error};
            }

            axis_point.K_plus = axis_point.theta - axis_point.nu;

            double c_minus_angle = average_cminus_angle(off_axis_parent, axis_point);
            axis_point.x = off_axis_parent.x - off_axis_parent.y / tan(c_minus_angle);

            // Corrector: apply the averaged source term with the same signed dy (see
            // axis_source_correction: dy/y_avg = -2 is the finite sin(theta)/y limit at the
            // axis, not a small-perturbation approximation).
            double source_corrected = axis_source_correction(off_axis_parent, axis_point);

            // Corrected nu: K_minus from parent, + source contribution
            double nu_corrected = off_axis_parent.K_minus + source_corrected;
            thermo_error = ctx.thermo.set_state_from_nu(axis_point, nu_corrected, axis_point.mach);
            if (thermo_error != MocErrorCode::NONE) {
                return {axis_point, thermo_error};
            }
            axis_point.K_minus = axis_point.nu; // theta=0
            axis_point.K_plus = -axis_point.nu;

            // Recompute position with corrected slope
            c_minus_angle = average_cminus_angle(off_axis_parent, axis_point);
            axis_point.x = off_axis_parent.x - off_axis_parent.y / tan(c_minus_angle);
            break;
        }
        default:
            throw std::runtime_error("Invalid flow type specified.");
    }

    // solve_axis_point has no non-downstream guard today (a documented gap): the
    // resulting axis point's x must not precede its parent's.
    if (axis_point.x < off_axis_parent.x) {
        return {axis_point, MocErrorCode::NON_DOWNSTREAM_POINT};
    }

    MocErrorCode validity = check_point_validity(axis_point, ctx.options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        return {axis_point, validity};
    }

    return {axis_point, MocErrorCode::NONE};
}

PointResult solve_fan_axis_point(const MocSolveContext& ctx, const CharacteristicPoint& expansion_point) {
    // The planar and axisymmetric bootstraps are identical: both simply carry the expansion
    // point's theta and nu onto the axis (the axisymmetric source term is not applied here --
    // the initial expansion fan is modeled as a centered expansion at a point, and the first
    // axis point is a direct consequence of that geometric construction; source terms enter
    // in subsequent solve_axis_point calls).
    CharacteristicPoint point{};
    point.y = 0.0; //point always lies on axis
    point.theta = expansion_point.theta;
    // nu = expansion point nu follows from geometric analysis
    MocErrorCode thermo_error = ctx.thermo.set_state_from_nu(
        point, expansion_point.nu, expansion_point.mach);
    if (thermo_error != MocErrorCode::NONE) {
        point.x = expansion_point.x;
        return {point, thermo_error};
    }
    point.K_plus = point.theta - point.nu;
    point.K_minus = expansion_point.K_minus;

    double c_minus_angle = average_cminus_angle(expansion_point, point);
    point.x = expansion_point.x - expansion_point.y / tan(c_minus_angle);

    MocErrorCode validity = check_point_validity(point, ctx.options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        return {point, validity};
    }
    return {point, MocErrorCode::NONE};
}

PointResult solve_wall_point_design(
    const MocSolveContext& ctx,
    const CharacteristicPoint& interior_parent,
    const CharacteristicPoint& previous_wall_point,
    double theta_wall)
{
    PointResult flow_result = solve_wall_flow(ctx, interior_parent, theta_wall);
    if (flow_result.error != MocErrorCode::NONE) {
        flow_result.point.x = interior_parent.x;
        flow_result.point.y = interior_parent.y;
        return flow_result;
    }
    CharacteristicPoint wall_point = flow_result.point;

    // the angle between the interior point and the wall point
    // is given by the C+ characteristic.
    double c_plus_angle = average_cplus_angle(interior_parent, wall_point);

    // the angle between the previous wall point and the current wall point
    // is given by wall_theta.
    double prev_wall_angle = 0.5*(previous_wall_point.theta + wall_point.theta);

    auto [x,y] = characteristic_intersection_with_angle(
        interior_parent, previous_wall_point,
        c_plus_angle, prev_wall_angle);

    wall_point.x = x;
    wall_point.y = y;

    switch (ctx.options.flow_type) {
        case MocFlowKind::PLANAR: { break; }
        case MocFlowKind::AXISYMMETRIC: {
            // Wall corrector step
            // Unlike the case of an interior point, we only have one characteristic line
            // and theta is already known (specified by wall).
            // The general compatibility equation for the C+ characteristic is:
            // $theta_P - theta_B - S_{BP} = nu_P - nu_B$
            // where S is the source term, B is the interior point and P is the wall point.

            // The solution procedure is:
            // 1. Calculate source term from estimated y.
            // 2. Calculate nu from axisymmetric compatibility equation.
            // 3. Obtain thermodynamic values for given value of nu.
            // 4. Recompute new position based on new C+ angle
            // 5. Repeat until source residual is satisfactory.
            // In practice only 1-2 iterations should be needed.
            double old_S = 0.0;
            double residual = 0.0;
            const int max_iter = 4;
            for (int i = 0; i < max_iter; i++){

                double S = cplus_source_term(interior_parent, wall_point);
                residual = std::abs(S - old_S);
                old_S = S;
                if (residual < ctx.options.solver_options.abstol) break;

                // C+ compatibility nu_P = theta_P - theta_B + nu_B + S, matching the interior
                // solver's K+ = theta - nu decreasing by the source along C+ (d(theta-nu) = -S).
                // The source enters with a +S sign here, NOT -S.
                MocErrorCode thermo_error = ctx.thermo.set_state_from_nu(
                    wall_point,
                    wall_point.theta-interior_parent.theta + S + interior_parent.nu,
                    interior_parent.mach);
                if (thermo_error != MocErrorCode::NONE) {
                    return {wall_point, thermo_error};
                }

                wall_point.update_Ks();
                double corrected_cplus_angle = average_cplus_angle(interior_parent, wall_point);
                auto [new_x,new_y] = characteristic_intersection_with_angle(
                    interior_parent,
                    previous_wall_point,
                    corrected_cplus_angle,
                    prev_wall_angle);

                wall_point.x = new_x;
                wall_point.y = new_y;

                if (i == max_iter - 1 && residual >= ctx.options.solver_options.abstol) {
                    ctx.log.warning("Wall design source term iteration did not converge."
                        "Residual={} at ({}, {}).", residual, wall_point.x, wall_point.y);
                }
            }
            break;
        }
        default:
            throw std::runtime_error("Invalid flow type specified.");
    }

    MocErrorCode validity = check_point_validity(wall_point, ctx.options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        return {wall_point, validity};
    }

    return {wall_point, MocErrorCode::NONE};
}

double cplus_source_term(const CharacteristicPoint& p, double new_y)
{
    double y_avg = 0.5 * (p.y + new_y);
    double dy = new_y - p.y;
    return -sin(p.theta)/(p.mach * sin(p.theta + p.mu)) * dy/y_avg;
}

double cplus_source_term(const CharacteristicPoint& p1, const CharacteristicPoint& p3)
{
    double y_avg = 0.5 * (p1.y + p3.y);
    double dy = p3.y - p1.y;
    double theta_avg = 0.5 * (p1.theta + p3.theta);
    double theta_plus_mu_avg = 0.5 * ((p1.theta + p1.mu) + (p3.theta + p3.mu));
    double mach_avg = 0.5 * (p1.mach + p3.mach);
    return sin(theta_avg)/(mach_avg * sin(theta_plus_mu_avg)) * dy/y_avg;
}

double cminus_source_term(const CharacteristicPoint& p, double new_y)
{
    double y_avg = 0.5 * (p.y + new_y);
    double dy = new_y - p.y;
    return sin(p.theta)/(p.mach * sin(p.theta - p.mu)) * dy/y_avg;
}

double axisymmetric_source(double mu_avg, double theta_avg, double y_avg, double char_angle_avg) {
    return sin(mu_avg) * sin(theta_avg) / (y_avg * cos(char_angle_avg));
}

double axis_source_correction(const CharacteristicPoint& foot, const CharacteristicPoint& axis_point) {
    double theta_avg = 0.5 * foot.theta;
    double mu_avg = 0.5 * (foot.mu + axis_point.mu);
    double M_avg = 0.5 * (foot.mach + axis_point.mach);
    double denom = sin(theta_avg - mu_avg);
    return std::abs(denom) > 1e-12 ? sin(theta_avg) / (M_avg * denom) * (-2.0) : 0.0;
}

std::optional<std::pair<double, double>> intersect_ray_with_wall(
    double x, double y, double angle, const NozzleProfile& wall)
{
    double c_plus_slope = tan(angle);
    double x_hit = 0.0, y_hit = 0.0;
    bool found = false;

    // A wall hit needs at least one segment; a degenerate profile would underflow the
    // loop bound and the extrapolation indices below.
    if (wall.x.size() < 2) {
        throw std::runtime_error("intersect_ray_with_wall: nozzle profile has fewer than two points.");
    }

    for (size_t k = 0; k < wall.x.size() - 1; k++) {
        double wall_slope = (wall.y[k+1] - wall.y[k])
            / (wall.x[k+1] - wall.x[k]);
        double wall_intercept = wall.y[k] - wall_slope * wall.x[k];
        double char_intercept = y - c_plus_slope * x;

        double x_intercept = (wall_intercept - char_intercept) / (c_plus_slope - wall_slope);

        // check if intersection is within the wall segment
        if (x_intercept >= wall.x[k] && x_intercept <= wall.x[k+1]) {
            x_hit = x_intercept;
            y_hit = y + c_plus_slope * (x_hit - x);
            found = true;
            break;
        }
    }
    if (!found) {
        // if the intersection cannot be found within the nozzle profile, report nothing.
        return std::nullopt;
    }
    return std::make_pair(x_hit, y_hit);
}

CharacteristicPoint unwrap_or_throw(const PointResult& result, std::string_view context) {
    if (result.error != MocErrorCode::NONE) {
        throw ConvergenceError(std::format(
            "{}: unit process failed at ({}, {}) with error {}.",
            context, result.point.x, result.point.y, to_string(result.error)));
    }
    return result.point;
}

double find_node_mach(
    const MocThermo& thermo,
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2,
    double source_delta,
    double mach_guess,
    double abstol,
    MocLog& log)
{
    // initial guess: mach of upstream
    double mach = (mach_guess > 1.0) ? mach_guess : p1.mach;

    double delta_theta = p1.theta - p2.theta;
    const int max_iter = 15;
    // Newton's method to find root
    // Residual function obtained from combining both compatibility equations:
    //
    //$$\epsilon (V_P) = (\theta_{A}- \theta_{B}) + \Delta \nu (V_{A} \rightarrow V_{P})
    // + \Delta \nu (V_{B} \rightarrow V_{P}) - (S_{A}- S_{B})$$
    // where A and B indicate upstream points.
    for (int i = 0; i < max_iter; i++) {
        double nu3, derivative;

        if (thermo.chemistry() == GasChemistry::PERFECT_GAS) {
            nu3 = prandtl_meyer(mach, thermo.gamma());
            // dnu/dM for perfect gas
            derivative = 2.0 * prandtl_meyer_derivative(mach, thermo.gamma());
        } else {
            const PrandtlMeyerTable& table = thermo.table();
            auto [idx, weight] = table.find_mach_index_and_weight(mach);
            nu3 = table.interpolate_at_index(idx, weight, table.nus);
            double V = table.interpolate_at_index(idx, weight, table.velocities);
            derivative = 2.0 * sqrt(mach * mach - 1.0) / V;
        }

        double residual = (nu3 - p1.nu) + (nu3 - p2.nu) - delta_theta - source_delta;
        if (std::abs(residual) < abstol) return mach;

        mach -= residual / derivative;
    }

    // rootfinding has failed
    log.warning("find_node_mach did not converge after {} iterations. "
        "Last Mach={}, delta_theta={}.", max_iter, mach, delta_theta);
    return -1.0;
}

} // namespace Goddard
