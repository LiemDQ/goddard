// Inverse (reference-plane) marching kernel for axisymmetric/planar MoC analysis
// (Package B). See instructions/moc_fix/B.md for the algorithm and
// instructions/moc_fix/diagnosis.md (Summary, A3, A4) for why the DIRECT chain-pairing
// kernel (moc_nozzle.cpp) sheds resolution near the axis in axisymmetric analysis.
//
// The DIRECT kernel and every unit process it uses (solve_interior_point_axisymmetric,
// solve_axis_point, solve_wall_point_analysis, solve_inverse_interior_point, and the
// anonymous trace_back_to_front in moc_nozzle.cpp) are untouched by this file: every
// method here is a new sibling, not a replacement.
#include <cmath>
#include <algorithm>
#include <format>
#include "goddard/moc_nozzle.hpp"
#include "goddard/error.hpp"

namespace Goddard {

std::string_view to_string(MocStepLimiter limiter) {
    switch (limiter) {
        case MocStepLimiter::NONE: return "NONE";
        case MocStepLimiter::CFL: return "CFL";
        case MocStepLimiter::WALL_FOOT: return "WALL_FOOT";
        case MocStepLimiter::WALL_TURN: return "WALL_TURN";
        case MocStepLimiter::EXIT: return "EXIT";
        default: return "UNKNOWN";
    }
}

namespace {

// Ray/polyline intersection: walk upstream from (x, y) along `angle` until the polyline
// `front` (indices into net.points, y strictly increasing) is met. Returns the segment
// index and the parameter along it, or nothing when the ray misses.
//
// Deliberately duplicated from the identically-named helper in moc_nozzle.cpp's anonymous
// namespace (used there by solve_inverse_interior_point for mesh-control rung insertion):
// that helper has internal linkage in a different translation unit, so it cannot be called
// from here. The two are kept byte-for-byte identical in logic.
std::optional<std::pair<size_t, double>> trace_back_to_front(
    double x, double y, double angle,
    const CharacteristicNet& net, const std::vector<size_t>& front)
{
    const double dir_x = -std::cos(angle);
    const double dir_y = -std::sin(angle);

    for (size_t i = 0; i + 1 < front.size(); i++) {
        const CharacteristicPoint& a = net.points[front[i]];
        const CharacteristicPoint& b = net.points[front[i + 1]];
        const double seg_x = b.x - a.x;
        const double seg_y = b.y - a.y;

        const double det = dir_x * (-seg_y) - dir_y * (-seg_x);
        if (std::abs(det) < 1e-14) continue; // ray parallel to this segment

        const double rhs_x = a.x - x;
        const double rhs_y = a.y - y;
        const double s = (rhs_x * (-seg_y) - rhs_y * (-seg_x)) / det;
        const double u = (dir_x * rhs_y - dir_y * rhs_x) / det;

        if (s > 0.0 && u >= 0.0 && u <= 1.0) return std::make_pair(i, u);
    }
    return std::nullopt;
}

/** theta, nu, and a Mach seed at a query y, fit from front data (B.md Sec. "Front
 * interpolation"): quadratic Lagrange through the 3 points bracketing y_query, linear at
 * either end of the front where a third point is not available. */
struct FrontFit { double theta, nu, mach; };

FrontFit lagrange_front_fit(
    const CharacteristicNet& net, const std::vector<size_t>& front,
    size_t seg_idx, double y_query)
{
    const size_t n = front.size();
    size_t i0 = 0, i1 = 0, i2 = 0;
    bool quad = true;
    if (seg_idx > 0) { i0 = seg_idx - 1; i1 = seg_idx; i2 = seg_idx + 1; }
    else if (seg_idx + 2 < n) { i0 = seg_idx; i1 = seg_idx + 1; i2 = seg_idx + 2; }
    else quad = false;

    if (!quad) {
        const CharacteristicPoint& a = net.points[front[seg_idx]];
        const CharacteristicPoint& b = net.points[front[seg_idx + 1]];
        const double t = (b.y > a.y) ? (y_query - a.y) / (b.y - a.y) : 0.0;
        return {a.theta + t * (b.theta - a.theta),
                a.nu + t * (b.nu - a.nu),
                a.mach + t * (b.mach - a.mach)};
    }

    const CharacteristicPoint& p0 = net.points[front[i0]];
    const CharacteristicPoint& p1 = net.points[front[i1]];
    const CharacteristicPoint& p2 = net.points[front[i2]];
    const double y0 = p0.y, y1 = p1.y, y2 = p2.y;
    const double L0 = (y_query - y1) * (y_query - y2) / ((y0 - y1) * (y0 - y2));
    const double L1 = (y_query - y0) * (y_query - y2) / ((y1 - y0) * (y1 - y2));
    const double L2 = (y_query - y0) * (y_query - y1) / ((y2 - y0) * (y2 - y1));
    auto lag = [&](double v0, double v1, double v2) { return v0 * L0 + v1 * L1 + v2 * L2; };
    return {lag(p0.theta, p1.theta, p2.theta), lag(p0.nu, p1.nu, p2.nu), lag(p0.mach, p1.mach, p2.mach)};
}

/** Segment index bracketing y_query along `front` (y increasing), or nullopt if outside. */
std::optional<size_t> bracket_front(
    const CharacteristicNet& net, const std::vector<size_t>& front, double y_query)
{
    for (size_t i = 0; i + 1 < front.size(); i++) {
        const double y0 = net.points[front[i]].y;
        const double y1 = net.points[front[i + 1]].y;
        if (y_query >= y0 - 1e-12 && y_query <= y1 + 1e-12) return i;
    }
    return std::nullopt;
}

/** Geometric position of a ray/front hit, plus whether it was found by the axis mirror. */
struct FootHit { double x, y; size_t seg; bool mirrored; };

/** Resolve a trace_back_to_front hit to its (x, y) position on `front`. */
FootHit resolve_hit(
    const std::pair<size_t, double>& hit, const CharacteristicNet& net,
    const std::vector<size_t>& front, bool mirrored)
{
    const CharacteristicPoint& a = net.points[front[hit.first]];
    const CharacteristicPoint& b = net.points[front[hit.first + 1]];
    const double x = a.x + hit.second * (b.x - a.x);
    const double y = a.y + hit.second * (b.y - a.y);
    return FootHit{x, mirrored ? -y : y, hit.first, mirrored};
}

/** Trace a ray back to `front`, with no axis-mirror fallback: used for the C- family, which
 * (theta - mu, backward direction rising in y) never needs it -- see trace_with_axis_mirror. */
std::optional<FootHit> trace_plain(
    double x0, double y0, double angle,
    const CharacteristicNet& net, const std::vector<size_t>& front)
{
    if (auto hit = trace_back_to_front(x0, y0, angle, net, front)) {
        return resolve_hit(*hit, net, front, false);
    }
    return std::nullopt;
}

/**
 * Trace a ray back to `front`, mirroring across the axis when the direct trace misses.
 *
 * A C+ traced backward from a point close to the axis can cross y=0 before meeting `front`
 * (which only spans y >= 0): by axisymmetric/planar centerline mirror symmetry, tracing the
 * reflected ray from (x0, -y0) at -angle into the (unreflected) front and negating theta on
 * the way back gives the physically correct foot, now reported at y < 0. See B.md Sec.
 * "Interior point", axis mirror. Used only for the C+ family (and the wall's own C+
 * trace, defensively) -- see trace_plain for why C- must not fall back to this.
 */
std::optional<FootHit> trace_with_axis_mirror(
    double x0, double y0, double angle,
    const CharacteristicNet& net, const std::vector<size_t>& front)
{
    if (auto hit = trace_back_to_front(x0, y0, angle, net, front)) {
        return resolve_hit(*hit, net, front, false);
    }
    if (auto hit = trace_back_to_front(x0, -y0, -angle, net, front)) {
        return resolve_hit(*hit, net, front, true);
    }
    return std::nullopt;
}

/** Throw ConvergenceError if a thermodynamic update failed; caught by solve()'s existing
 * initialization exception boundary, matching moc_nozzle.cpp's unwrap_or_throw pattern
 * (not reused directly: that helper has internal linkage in a different translation unit). */
void throw_if_thermo_error(MocErrorCode err, std::string_view context, double x, double y) {
    if (err != MocErrorCode::NONE) {
        throw ConvergenceError(std::format(
            "{}: unit process failed at ({}, {}) with error {}.", context, x, y, to_string(err)));
    }
}

} // namespace

MocMarchScheme MocNozzle::resolve_march_scheme() const {
    if (m_options.march_scheme != MocMarchScheme::AUTO) return m_options.march_scheme;
    if (m_options.mode == MocMode::DESIGN_MIN_LENGTH) return MocMarchScheme::DIRECT;
    if (m_options.flow_type == MocFlowKind::AXISYMMETRIC &&
        (m_options.mode == MocMode::ANALYSIS || m_options.mode == MocMode::DESIGN_RAO)) {
        return MocMarchScheme::INVERSE;
    }
    return MocMarchScheme::DIRECT;
}

std::vector<CharacteristicPoint> MocNozzle::build_inverse_initial_front(
    const std::vector<CharacteristicPoint>& data_line)
{
    if (!m_initial_line_family.has_value()) {
        // Kliegel-Levine topology (or a test-injected override, seeded through the same
        // path -- see m_inverse_front_override): already spans axis to wall.
        return data_line;
    }

    // Centered-fan topology: data_line is the marched C+ through the fan rays, and its top
    // point P sits on the last ray, short of the wall (MocInitialization::
    // initialize_centered_expansion emits the rays; generate_initial_data_line marches them
    // into real positions via solve_initial_axis_point_centered_exp/solve_interior_point).
    // Extend straight up to the wall with the fan's own uniform post-last-ray state -- the
    // region beyond the last ray of a centered fan is uniform -- then correct only the new
    // wall point's nu from the planar C+ compatibility relation with P.
    const CharacteristicPoint& P = data_line.back();
    const NozzleProfile& wall = m_options.nozzle_profile;
    const double y_wall = wall.radius_at(P.x);

    std::vector<CharacteristicPoint> front = data_line;
    if (!(y_wall > P.y)) return front; // degenerate (profile already at/below P): nothing to extend

    const double dy_ref = (data_line.size() >= 2)
        ? (data_line.back().y - data_line[data_line.size() - 2].y)
        : (y_wall - P.y);
    const size_t n_extra = std::max<size_t>(1,
        static_cast<size_t>(std::llround((y_wall - P.y) / std::max(dy_ref, 1e-12))));

    for (size_t k = 1; k < n_extra; k++) {
        const double t = static_cast<double>(k) / static_cast<double>(n_extra);
        CharacteristicPoint pt{};
        pt.x = P.x;
        pt.y = P.y + t * (y_wall - P.y);
        pt.theta = P.theta;
        throw_if_thermo_error(
            update_thermodynamic_state_from_nu(pt, P.nu, P.mach),
            "Inverse march initial front (fan extension)", pt.x, pt.y);
        pt.K_plus = pt.theta - pt.nu;
        pt.K_minus = pt.theta + pt.nu;
        front.push_back(pt);
    }

    CharacteristicPoint wall_pt{};
    wall_pt.x = P.x;
    wall_pt.y = y_wall;
    wall_pt.theta = wall.theta_at(P.x);
    const double nu_w = P.nu + wall_pt.theta - P.theta;
    throw_if_thermo_error(
        update_thermodynamic_state_from_nu(wall_pt, nu_w, P.mach),
        "Inverse march initial front (fan wall point)", wall_pt.x, wall_pt.y);
    wall_pt.K_plus = wall_pt.theta - wall_pt.nu;
    wall_pt.K_minus = wall_pt.theta + wall_pt.nu;
    front.push_back(wall_pt);
    return front;
}

std::vector<size_t> MocNozzle::seed_inverse_front(
    CharacteristicNet& net, const std::vector<CharacteristicPoint>& front_points) const
{
    std::vector<size_t> indices;
    indices.reserve(front_points.size());
    for (const CharacteristicPoint& pt : front_points) {
        net.points.push_back(pt);
        net.membership.push_back(PointMembership{});
        indices.push_back(net.points.size() - 1);
    }
    net.axis_point_indices.push_back(indices.front());
    net.wall_point_indices.push_back(indices.back());
    net.wall_x.push_back(front_points.back().x);
    net.wall_y.push_back(front_points.back().y);
    net.fronts.push_back(indices);
    return indices;
}

PointResult MocNozzle::solve_inverse_march_interior_point(
    double x_new, double y_new,
    const CharacteristicNet& net,
    const std::vector<size_t>& front)
{
    CharacteristicPoint p4{};
    p4.x = x_new;
    p4.y = y_new;

    // Seed from a direct cross-front interpolation at y_new (B.md Sec. 5); only seeds the
    // iteration below, so a linear/quadratic blend across characteristics is fine here even
    // though it would not be an admissible final answer.
    std::optional<size_t> seed_seg = bracket_front(net, front, y_new);
    if (!seed_seg.has_value()) return {p4, MocErrorCode::NON_DOWNSTREAM_POINT};
    FrontFit seed_fit = lagrange_front_fit(net, front, *seed_seg, y_new);
    p4.theta = seed_fit.theta;
    MocErrorCode err = update_thermodynamic_state_from_nu(p4, seed_fit.nu, seed_fit.mach);
    if (err != MocErrorCode::NONE) return {p4, err};

    const bool axisymmetric = (m_options.flow_type == MocFlowKind::AXISYMMETRIC);
    // First iteration: use the seed's own angles for both ends (no foot is known yet).
    double angle_minus = p4.theta - p4.mu;
    double angle_plus = p4.theta + p4.mu;

    constexpr int max_iters = 4;
    for (int iter = 0; iter < max_iters; iter++) {
        std::optional<FootHit> hit_minus = trace_plain(p4.x, p4.y, angle_minus, net, front);
        std::optional<FootHit> hit_plus = trace_with_axis_mirror(p4.x, p4.y, angle_plus, net, front);
        if (!hit_minus.has_value() || !hit_plus.has_value()) {
            return {p4, MocErrorCode::NON_DOWNSTREAM_POINT};
        }

        auto build_foot = [&](const FootHit& h, CharacteristicPoint& out) -> MocErrorCode {
            const double y_for_fit = h.mirrored ? -h.y : h.y;
            FrontFit fit = lagrange_front_fit(net, front, h.seg, y_for_fit);
            out = CharacteristicPoint{};
            out.x = h.x;
            out.y = h.y;
            out.theta = h.mirrored ? -fit.theta : fit.theta;
            return update_thermodynamic_state_from_nu(out, fit.nu, fit.mach);
        };

        CharacteristicPoint foot_minus{}, foot_plus{};
        MocErrorCode err_m = build_foot(*hit_minus, foot_minus);
        if (err_m != MocErrorCode::NONE) return {p4, err_m};
        MocErrorCode err_p = build_foot(*hit_plus, foot_plus);
        if (err_p != MocErrorCode::NONE) return {p4, err_p};

        // Same source terms as solve_interior_point_axisymmetric / solve_inverse_interior_point:
        //   along C+:  d(theta - nu) = -L dx,  L = sin(mu) sin(theta) / (y cos(theta + mu))
        //   along C-:  d(theta + nu) = +M dx,  M = sin(mu) sin(theta) / (y cos(theta - mu))
        double source_minus = 0.0, source_plus = 0.0;
        if (axisymmetric) {
            const double mu_m = 0.5 * (foot_minus.mu + p4.mu);
            const double th_m = 0.5 * (foot_minus.theta + p4.theta);
            const double y_m = 0.5 * (foot_minus.y + p4.y);
            const double ang_m = 0.5 * (angle_minus + (foot_minus.theta - foot_minus.mu));
            if (std::abs(y_m) > 1e-12) source_minus = std::sin(mu_m) * std::sin(th_m) / (y_m * std::cos(ang_m));

            const double mu_p = 0.5 * (foot_plus.mu + p4.mu);
            const double ang_p = 0.5 * (angle_plus + (foot_plus.theta + foot_plus.mu));
            if (hit_plus->mirrored) {
                // Axis-limit ratio evaluated at the new point: the mirrored foot's y is
                // negative while p4.y is small and positive, so their average passes near
                // zero even though neither side is individually degenerate.
                const double ratio = (std::abs(p4.y) > 1e-12) ? std::sin(p4.theta) / p4.y : 0.0;
                source_plus = std::sin(mu_p) * ratio / std::cos(ang_p);
            } else {
                const double th_p = 0.5 * (foot_plus.theta + p4.theta);
                const double y_p = 0.5 * (foot_plus.y + p4.y);
                if (std::abs(y_p) > 1e-12) source_plus = std::sin(mu_p) * std::sin(th_p) / (y_p * std::cos(ang_p));
            }
        }

        const double K_minus = (foot_minus.theta + foot_minus.nu) + source_minus * (p4.x - foot_minus.x);
        const double K_plus = (foot_plus.theta - foot_plus.nu) - source_plus * (p4.x - foot_plus.x);

        const double theta_prev = p4.theta;
        const double nu_prev = p4.nu;
        p4.theta = 0.5 * (K_minus + K_plus);
        err = update_thermodynamic_state_from_nu(p4, 0.5 * (K_minus - K_plus), p4.mach);
        if (err != MocErrorCode::NONE) return {p4, err};

        const bool converged =
            std::abs(p4.theta - theta_prev) < 1e-9 && std::abs(p4.nu - nu_prev) < 1e-9;
        // Predictor-corrector angle for the next trace: average of the refined new point
        // and the foot found this iteration (B.md Sec. "Interior point").
        angle_minus = average_cminus_angle(p4, foot_minus);
        angle_plus = average_cplus_angle(p4, foot_plus);
        if (converged) break;
    }

    p4.K_minus = p4.theta + p4.nu;
    p4.K_plus = p4.theta - p4.nu;

    MocErrorCode validity = check_point_validity(p4, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) return {p4, validity};
    return {p4, MocErrorCode::NONE};
}

PointResult MocNozzle::solve_inverse_march_axis_point(
    double x_new,
    const CharacteristicNet& net,
    const std::vector<size_t>& front)
{
    CharacteristicPoint axis_point{};
    axis_point.x = x_new;
    axis_point.y = 0.0;
    axis_point.theta = 0.0;

    const CharacteristicPoint& seed_src = net.points[front.front()];
    MocErrorCode err = update_thermodynamic_state_from_nu(axis_point, seed_src.nu, seed_src.mach);
    if (err != MocErrorCode::NONE) return {axis_point, err};

    const bool axisymmetric = (m_options.flow_type == MocFlowKind::AXISYMMETRIC);
    double angle_minus = axis_point.theta - axis_point.mu; // seed angle for iteration 1

    constexpr int max_iters = 4;
    for (int iter = 0; iter < max_iters; iter++) {
        std::optional<FootHit> hit = trace_plain(axis_point.x, 0.0, angle_minus, net, front);
        if (!hit.has_value()) return {axis_point, MocErrorCode::NON_DOWNSTREAM_POINT};

        const double y_for_fit = hit->mirrored ? -hit->y : hit->y;
        FrontFit fit = lagrange_front_fit(net, front, hit->seg, y_for_fit);
        CharacteristicPoint foot{};
        foot.x = hit->x;
        foot.y = hit->y;
        foot.theta = hit->mirrored ? -fit.theta : fit.theta;
        err = update_thermodynamic_state_from_nu(foot, fit.nu, fit.mach);
        if (err != MocErrorCode::NONE) return {axis_point, err};

        double nu_new;
        if (axisymmetric) {
            // Reuses solve_axis_point's corrector algebra: y_avg = (foot.y + 0)/2 always,
            // since the axis point sits at y=0, so dy/y_avg = -2 is an algebraic identity
            // here (not a small-perturbation approximation) and generalizes unchanged to a
            // traced-back foot that is not an actual net parent.
            const double theta_avg = 0.5 * foot.theta;
            const double mu_avg = 0.5 * (foot.mu + axis_point.mu);
            const double M_avg = 0.5 * (foot.mach + axis_point.mach);
            const double denom = std::sin(theta_avg - mu_avg);
            const double source = (std::abs(denom) > 1e-12)
                ? std::sin(theta_avg) / (M_avg * denom) * (-2.0) : 0.0;
            nu_new = (foot.theta + foot.nu) + source;
        } else {
            nu_new = foot.theta + foot.nu; // planar: K- exactly preserved, no source
        }

        const double nu_prev = axis_point.nu;
        err = update_thermodynamic_state_from_nu(axis_point, nu_new, axis_point.mach);
        if (err != MocErrorCode::NONE) return {axis_point, err};

        const bool converged = std::abs(axis_point.nu - nu_prev) < 1e-9;
        angle_minus = average_cminus_angle(axis_point, foot);
        if (converged) break;
    }

    axis_point.K_minus = axis_point.nu;
    axis_point.K_plus = -axis_point.nu;

    MocErrorCode validity = check_point_validity(axis_point, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) return {axis_point, validity};
    return {axis_point, MocErrorCode::NONE};
}

PointResult MocNozzle::solve_inverse_march_wall_point(
    double x_new, double y_new,
    const CharacteristicNet& net,
    const std::vector<size_t>& front)
{
    CharacteristicPoint wall_point{};
    wall_point.x = x_new;
    wall_point.y = y_new;
    wall_point.theta = m_options.nozzle_profile.theta_at(x_new); // fixed by the contour

    const CharacteristicPoint& seed_src = net.points[front.back()];
    MocErrorCode err = update_thermodynamic_state_from_nu(wall_point, seed_src.nu, seed_src.mach);
    if (err != MocErrorCode::NONE) return {wall_point, err};

    const bool axisymmetric = (m_options.flow_type == MocFlowKind::AXISYMMETRIC);
    double angle_plus = wall_point.theta + wall_point.mu; // seed angle for iteration 1

    constexpr int max_iters = 4;
    for (int iter = 0; iter < max_iters; iter++) {
        std::optional<FootHit> hit = trace_with_axis_mirror(wall_point.x, wall_point.y, angle_plus, net, front);
        if (!hit.has_value()) return {wall_point, MocErrorCode::NON_DOWNSTREAM_POINT};

        const double y_for_fit = hit->mirrored ? -hit->y : hit->y;
        FrontFit fit = lagrange_front_fit(net, front, hit->seg, y_for_fit);
        CharacteristicPoint foot{};
        foot.x = hit->x;
        foot.y = hit->y;
        foot.theta = hit->mirrored ? -fit.theta : fit.theta;
        err = update_thermodynamic_state_from_nu(foot, fit.nu, fit.mach);
        if (err != MocErrorCode::NONE) return {wall_point, err};

        // K+ transported with cplus_source_term(foot, wall_point), exactly as the DIRECT
        // wall solvers do; here the position is prescribed and the foot (and hence the
        // source term) is what iterates, so convergence is judged on nu rather than on the
        // source-term residual those solvers use.
        const double S = axisymmetric ? cplus_source_term(foot, wall_point) : 0.0;
        const double nu_new = wall_point.theta - foot.theta + S + foot.nu;

        const double nu_prev = wall_point.nu;
        err = update_thermodynamic_state_from_nu(wall_point, nu_new, wall_point.mach);
        if (err != MocErrorCode::NONE) return {wall_point, err};

        const bool converged = std::abs(wall_point.nu - nu_prev) < 1e-9;
        angle_plus = average_cplus_angle(foot, wall_point);
        if (converged) break;
    }

    wall_point.K_plus = wall_point.theta - wall_point.nu;
    wall_point.K_minus = wall_point.theta + wall_point.nu;

    MocErrorCode validity = check_point_validity(wall_point, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) return {wall_point, validity};
    return {wall_point, MocErrorCode::NONE};
}

std::optional<MocFailure> MocNozzle::solve_inverse_characteristic_kernel(CharacteristicNet& net) {
    if (net.fronts.empty() || net.fronts.back().size() < 2) {
        return MocFailure{MocErrorCode::INITIALIZATION_FAILED,
            "Inverse march: initial front has fewer than two points.", 0.0, 0.0, -1};
    }

    const NozzleProfile& wall_profile = m_options.nozzle_profile;
    const double L = wall_profile.x_max();

    std::vector<size_t> front = net.fronts.back();
    const size_t M = front.size();

    // Normalized y-distribution of the initial front, computed once and reused for every
    // subsequent front's interior points (B.md Sec. "New front geometry").
    std::vector<double> s(M);
    {
        const double y_wall0 = net.points[front.back()].y;
        for (size_t j = 0; j < M; j++) {
            s[j] = (y_wall0 > 0.0)
                ? net.points[front[j]].y / y_wall0
                : static_cast<double>(j) / static_cast<double>(M - 1);
        }
    }

    constexpr int max_passes = 20000;
    for (int pass = 0; pass < max_passes; pass++) {
        const CharacteristicPoint& wall_old = net.points[front.back()];
        const CharacteristicPoint& axis_old = net.points[front.front()];

        // --- Step 1: step length ---
        double dy_min = std::numeric_limits<double>::max();
        double slope_max = 0.0;
        for (size_t j = 0; j < M; j++) {
            const CharacteristicPoint& p = net.points[front[j]];
            slope_max = std::max({slope_max, std::tan(p.theta + p.mu), std::tan(p.mu - p.theta)});
            if (j + 1 < M) dy_min = std::min(dy_min, net.points[front[j + 1]].y - p.y);
        }
        const double dx_cfl = (slope_max > 1e-12)
            ? m_options.inverse_cfl * dy_min / slope_max
            : std::numeric_limits<double>::max();

        double dx_foot = std::numeric_limits<double>::max();
        if (M >= 3) {
            const CharacteristicPoint& top = net.points[front[M - 2]];
            const double dy_top = wall_old.y - top.y;
            const double denom = std::tan(top.mu - top.theta) + std::tan(wall_old.theta);
            if (denom > 1e-12) dx_foot = dy_top / denom;
        }

        double dx_turn = std::numeric_limits<double>::max();
        {
            const double probe = std::isfinite(dx_cfl) ? std::max(dx_cfl, 1e-9) : 1e-6;
            const double theta_here = wall_profile.theta_at(wall_old.x);
            const double x_probe = std::min(wall_old.x + probe, L);
            const double dx_probe_actual = x_probe - wall_old.x;
            if (dx_probe_actual > 1e-12) {
                const double theta_probe = wall_profile.theta_at(x_probe);
                const double dtheta_dx = std::abs(theta_probe - theta_here) / dx_probe_actual;
                dx_turn = m_options.max_wall_turn_per_step / std::max(dtheta_dx, 1e-12);
            }
        }

        const double dx_exit = std::max(L - wall_old.x, 0.0);

        double dx = dx_cfl;
        MocStepLimiter limiter = MocStepLimiter::CFL;
        if (dx_foot < dx) { dx = dx_foot; limiter = MocStepLimiter::WALL_FOOT; }
        if (dx_turn < dx) { dx = dx_turn; limiter = MocStepLimiter::WALL_TURN; }
        if (dx_exit < dx) { dx = dx_exit; limiter = MocStepLimiter::EXIT; }

        if (!(dx > 0.0) || !std::isfinite(dx)) {
            return MocFailure{MocErrorCode::NON_DOWNSTREAM_POINT,
                std::format("Inverse march: non-positive step length ({:.3e}) at pass {}.", dx, pass),
                wall_old.x, wall_old.y, pass};
        }

        // --- Step 2: new front geometry ---
        const double x_w = wall_old.x + dx;
        const double y_w = wall_profile.radius_at(std::min(x_w, L));
        const double tilt_raw = m_options.front_tilt_decay * (axis_old.x - wall_old.x);
        double tilt_spacelike = tilt_raw;

        if (y_w > 0.0) {
            double min_cot = std::numeric_limits<double>::max();
            for (size_t j = 0; j < M; j++) {
                const CharacteristicPoint& p = net.points[front[j]];
                const double denom_c = std::tan(p.mu + std::abs(p.theta));
                if (denom_c > 1e-12) min_cot = std::min(min_cot, 1.0 / denom_c);
            }
            if (min_cot < std::numeric_limits<double>::max() && std::abs(tilt_spacelike) / y_w >= min_cot) {
                tilt_spacelike = 0.0;
            }
        }

        // Attempt the whole front (axis, every interior point, wall) at a given tilt
        // candidate. Returns the M solved points, or the failure of whichever point solve
        // failed first.
        struct FrontAttempt {
            bool ok = false;
            std::vector<CharacteristicPoint> points;
            MocFailure failure;
        };
        auto attempt_front = [&](double tilt_try) -> FrontAttempt {
            FrontAttempt out;
            const double x_a_try = x_w + tilt_try;
            std::vector<double> x_new(M), y_new(M);
            x_new[0] = x_a_try; y_new[0] = 0.0;
            x_new[M - 1] = x_w; y_new[M - 1] = y_w;
            for (size_t j = 1; j + 1 < M; j++) {
                y_new[j] = s[j] * y_w;
                x_new[j] = x_a_try + s[j] * (x_w - x_a_try);
            }

            out.points.resize(M);
            PointResult axis_result = solve_inverse_march_axis_point(x_new[0], net, front);
            if (axis_result.error != MocErrorCode::NONE) {
                out.failure = MocFailure{axis_result.error,
                    std::format("Inverse march axis point failed ({}) at pass {}.",
                        to_string(axis_result.error), pass),
                    axis_result.point.x, axis_result.point.y, pass};
                return out;
            }
            out.points[0] = axis_result.point;

            for (size_t j = 1; j + 1 < M; j++) {
                PointResult r = solve_inverse_march_interior_point(x_new[j], y_new[j], net, front);
                if (r.error != MocErrorCode::NONE) {
                    out.failure = MocFailure{r.error,
                        std::format("Inverse march interior point {} failed ({}) at pass {}.",
                            j, to_string(r.error), pass),
                        r.point.x, r.point.y, pass};
                    return out;
                }
                out.points[j] = r.point;
            }

            PointResult wall_result = solve_inverse_march_wall_point(x_new[M - 1], y_new[M - 1], net, front);
            if (wall_result.error != MocErrorCode::NONE) {
                out.failure = MocFailure{wall_result.error,
                    std::format("Inverse march wall point failed ({}) at pass {}.",
                        to_string(wall_result.error), pass),
                    wall_result.point.x, wall_result.point.y, pass};
                return out;
            }
            out.points[M - 1] = wall_result.point;

            for (size_t j = 0; j < M; j++) {
                MocErrorCode v = check_point_validity(out.points[j], m_options.solver_options.abstol);
                if (v != MocErrorCode::NONE) {
                    out.failure = MocFailure{v,
                        std::format("Inverse march point {} invalid ({}) at pass {}.",
                            j, to_string(v), pass),
                        out.points[j].x, out.points[j].y, pass};
                    return out;
                }
            }
            out.ok = true;
            return out;
        };

        FrontAttempt attempt = attempt_front(tilt_spacelike);
        double tilt = tilt_spacelike;

        // Fallback tilt search -- a change from the algorithm as given; see the report.
        // B.md's plain per-pass decay of tilt is computed independently of dx and of the
        // front's own near-axis curvature, so it does not always land in the (sometimes
        // narrow) range of tilt values for which every point's trace back to the previous
        // front actually lands on it, and on it correctly: every old front point near the
        // axis rises near-vertically (mu -> 90 deg there), far steeper than a new point's
        // own trace ray, so a tilt outside the (possibly narrow) window where the ray's
        // shallow climb still meets the old front's steep one either misses outright
        // (NON_DOWNSTREAM_POINT) or lands on the wrong segment, which then fails downstream
        // as an ordinary numerical error (NEGATIVE_THETA, PM_INVERSION_FAILED, ...) rather
        // than a geometric one. The fan-init front's axis end can also sit *upstream* of its
        // wall end (the opposite sign from a Kliegel-Levine line), so the search covers both
        // signs. Coarse scan across the full plausible range first, then bisect around the
        // best bracket found; only engages when the literal formula's own attempt failed.
        if (!attempt.ok) {
            const double tilt_extent = std::max(std::abs(tilt_raw), dx) * 5.0 + dx;
            constexpr int coarse_steps = 4000;
            bool found = false;
            for (int k = -coarse_steps; k <= coarse_steps && !found; k++) {
                const double candidate = tilt_extent * static_cast<double>(k)
                    / static_cast<double>(coarse_steps);
                FrontAttempt trial = attempt_front(candidate);
                if (trial.ok) {
                    // Bisect inward toward the smallest-magnitude working tilt found,
                    // refining against the previous (failing) coarse step for a tighter
                    // bracket.
                    double lo = tilt_extent * static_cast<double>(k - 1)
                        / static_cast<double>(coarse_steps);
                    double hi = candidate;
                    if (k == -coarse_steps) lo = candidate;
                    FrontAttempt best = trial;
                    for (int b = 0; b < 20 && lo != hi; b++) {
                        const double mid = 0.5 * (lo + hi);
                        FrontAttempt mid_trial = attempt_front(mid);
                        if (mid_trial.ok) { hi = mid; best = mid_trial; }
                        else { lo = mid; }
                    }
                    attempt = best;
                    tilt = hi;
                    found = true;
                }
            }
            if (!found) {
                log_debug("FRONT pass={} tilt search exhausted ({} candidates over "
                    "[{:.6f}, {:.6f}]); reporting original failure.",
                    pass, 2 * coarse_steps + 1, -tilt_extent, tilt_extent);
            }
        }

        if (!attempt.ok) {
            return attempt.failure;
        }
        std::vector<CharacteristicPoint>& new_front_points = attempt.points;
        const double x_a = x_w + tilt;

        log_debug("STEP pass={} dx_cfl={:.6f} dx_foot={:.6f} dx_turn={:.6f} dx_exit={:.6f} "
            "dx={:.6f} x_w={:.6f} y_w={:.6f} x_a={:.6f} tilt={:.6f}",
            pass, dx_cfl, dx_foot, dx_turn, dx_exit, dx, x_w, y_w, x_a, tilt);

        // --- Step 7: bookkeeping ---
        std::vector<size_t> new_front = seed_inverse_front(net, new_front_points);

        MocPassDiagnostics diag;
        diag.pass = pass;
        diag.front_points = M;
        diag.step_dx = dx;
        diag.step_limiter = limiter;
        diag.front_axis_x = x_a;
        diag.front_wall_x = x_w;
        {
            double mn = std::numeric_limits<double>::max();
            double mx = 0.0, total = 0.0;
            for (size_t j = 0; j + 1 < M; j++) {
                const double seg_dx = new_front_points[j + 1].x - new_front_points[j].x;
                const double seg_dy = new_front_points[j + 1].y - new_front_points[j].y;
                const double arc = std::hypot(seg_dx, seg_dy);
                mn = std::min(mn, arc);
                mx = std::max(mx, arc);
                total += arc;
            }
            diag.min_spacing = mn;
            diag.max_spacing = mx;
            diag.mean_spacing = total / static_cast<double>(M - 1);
        }
        m_pass_diagnostics.push_back(diag);

        log_debug("FRONT pass={} dx={:.6f} limiter={} x_axis={:.6f} x_wall={:.6f}",
            pass, dx, to_string(limiter), x_a, x_w);

        front = new_front;

        if (limiter == MocStepLimiter::EXIT) {
            return std::nullopt; // the front just built is the exit plane
        }
    }

    return MocFailure{MocErrorCode::MAX_ITERATIONS_REACHED,
        "Inverse march reached its pass safety cap (20000) before the exit plane.",
        0.0, 0.0, max_passes};
}

} // namespace Goddard
