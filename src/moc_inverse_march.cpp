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

/**
 * Riemann invariants (and theta, nu, a Mach seed) at a query y, fit from front data.
 *
 * The invariants K- = theta + nu and K+ = theta - nu are what the march transports, so they
 * are what is interpolated; theta and nu follow. Their profiles along a front are only
 * piecewise smooth: at the edge of a simple-wave region (the last ray of a throat fan, the
 * characteristic from a wall-curvature break) an invariant's slope jumps, typically from a
 * rising ramp to a plateau. A quadratic fit through three points straddling such a kink
 * overshoots the plateau, and in planar flow that overshoot is then carried downstream
 * exactly, along the characteristics, so it accumulates into a ringing of a degree or more
 * around the plateau -- enough to drive the small flow angle next to the axis negative.
 * Measured on the planar minimum-length round trip before this guard: K- of 30.74 deg two
 * points above the axis where the exact value is 30.00 everywhere.
 *
 * Two standard remedies, both applied: the three-point stencil is chosen on the side of the
 * bracketing segment with the smaller second difference (ENO selection, so a stencil is not
 * laid across a kink when a smooth one is available), and the result is clamped to the range
 * of the two bracketing points (so no overshoot survives at all). Smooth regions keep the
 * quadratic's second-order accuracy; a kink degrades locally to linear, which is the right
 * price there.
 */
struct FrontFit { double theta, nu, mach; };

FrontFit lagrange_front_fit(
    const CharacteristicNet& net, const std::vector<size_t>& front,
    size_t seg_idx, double y_query)
{
    const size_t n = front.size();
    const CharacteristicPoint& a = net.points[front[seg_idx]];
    const CharacteristicPoint& b = net.points[front[seg_idx + 1]];

    const double t = (b.y > a.y) ? (y_query - a.y) / (b.y - a.y) : 0.0;
    const double mach_seed = a.mach + t * (b.mach - a.mach);

    auto quadratic = [&](size_t i0, size_t i1, size_t i2, double v0, double v1, double v2) {
        const double y0 = net.points[front[i0]].y, y1 = net.points[front[i1]].y, y2 = net.points[front[i2]].y;
        const double L0 = (y_query - y1) * (y_query - y2) / ((y0 - y1) * (y0 - y2));
        const double L1 = (y_query - y0) * (y_query - y2) / ((y1 - y0) * (y1 - y2));
        const double L2 = (y_query - y0) * (y_query - y1) / ((y2 - y0) * (y2 - y1));
        return v0 * L0 + v1 * L1 + v2 * L2;
    };
    auto second_difference = [&](size_t i0, size_t i1, size_t i2, double CharacteristicPoint::* q) {
        return std::abs(net.points[front[i2]].*q - 2.0 * net.points[front[i1]].*q + net.points[front[i0]].*q);
    };

    const bool has_left = seg_idx >= 1;
    const bool has_right = seg_idx + 2 < n;
    double K_minus, K_plus;
    if (!has_left && !has_right) {
        K_minus = a.K_minus + t * (b.K_minus - a.K_minus);
        K_plus = a.K_plus + t * (b.K_plus - a.K_plus);
    } else {
        size_t i0, i1, i2;
        if (has_left && has_right) {
            const double left = second_difference(seg_idx - 1, seg_idx, seg_idx + 1, &CharacteristicPoint::K_minus)
                              + second_difference(seg_idx - 1, seg_idx, seg_idx + 1, &CharacteristicPoint::K_plus);
            const double right = second_difference(seg_idx, seg_idx + 1, seg_idx + 2, &CharacteristicPoint::K_minus)
                               + second_difference(seg_idx, seg_idx + 1, seg_idx + 2, &CharacteristicPoint::K_plus);
            if (left <= right) { i0 = seg_idx - 1; i1 = seg_idx; i2 = seg_idx + 1; }
            else { i0 = seg_idx; i1 = seg_idx + 1; i2 = seg_idx + 2; }
        } else if (has_left) { i0 = seg_idx - 1; i1 = seg_idx; i2 = seg_idx + 1; }
        else { i0 = seg_idx; i1 = seg_idx + 1; i2 = seg_idx + 2; }
        K_minus = quadratic(i0, i1, i2, net.points[front[i0]].K_minus, net.points[front[i1]].K_minus, net.points[front[i2]].K_minus);
        K_plus = quadratic(i0, i1, i2, net.points[front[i0]].K_plus, net.points[front[i1]].K_plus, net.points[front[i2]].K_plus);
        K_minus = std::clamp(K_minus, std::min(a.K_minus, b.K_minus), std::max(a.K_minus, b.K_minus));
        K_plus = std::clamp(K_plus, std::min(a.K_plus, b.K_plus), std::max(a.K_plus, b.K_plus));
    }
    return {0.5 * (K_minus + K_plus), 0.5 * (K_minus - K_plus), mach_seed};
}

/** Segment index bracketing y_query along `front` (y increasing), or nullopt if outside. */
std::optional<size_t> bracket_front(
    const CharacteristicNet& net, const std::vector<size_t>& front, double y_query)
{
    if (front.size() < 2) return std::nullopt;
    for (size_t i = 0; i + 1 < front.size(); i++) {
        const double y0 = net.points[front[i]].y;
        const double y1 = net.points[front[i + 1]].y;
        if (y_query >= y0 - 1e-12 && y_query <= y1 + 1e-12) return i;
    }
    // Outside the front's span (the new front is slightly taller than the old one at the
    // wall): clamp to the end segment. This only seeds the iteration, so a short
    // extrapolation is harmless.
    return (y_query < net.points[front.front()].y) ? 0 : front.size() - 2;
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

    // Centered-fan topology. data_line is the marched C+ from the fan's first axis point
    // through the rays: a characteristic, which cannot serve as an inverse-march front (a C+
    // traced back from the next front runs parallel to it). Build F_0 instead on the plane
    // x = x0 through that first axis point. Every plane point lies on one of the fan's rays
    // between the throat lip and the ray's marched data-line point, so its invariants are
    // interpolated along that ray between the lip state (K- = 2 theta_i, K+ = 0) and the
    // marched state, blending across the two bracketing rays. In planar flow the invariants
    // are constant along a ray and this is exact; in axisymmetric flow the C- source term
    // accumulated over the ray -- large near the lip, where the rays are near-sonic -- enters
    // through the marched state exactly as the DIRECT kernel's own one-step fan does. Points
    // above the last ray take the last ray's state; the wall point takes the contour's angle
    // with nu from the C+ compatibility relation with the point below it. The lip is (0, 1)
    // in the normalized units initialize_centered_expansion uses.
    const NozzleProfile& wall = m_options.nozzle_profile;
    const size_t n_rays = data_line.size();
    if (n_rays < 2 || m_theta_schedule.size() < n_rays) {
        throw ConvergenceError("Inverse march initial front (fan): data line and theta schedule are inconsistent.");
    }
    const double x0 = data_line.front().x;
    if (!(x0 > 0.0) || x0 < wall.x_min() || x0 > wall.x_max()) {
        throw ConvergenceError(std::format(
            "Inverse march initial front (fan): the fan's first axis point x = {} is not "
            "inside the contour's range [{}, {}].", x0, wall.x_min(), wall.x_max()));
    }
    const double y_wall = wall.radius_at(x0);
    const double lip_x = 0.0;
    const double lip_y = 1.0;

    std::vector<double> ray_angle(n_rays), ray_length(n_rays);
    for (size_t i = 0; i < n_rays; i++) {
        ray_angle[i] = std::atan2(data_line[i].y - lip_y, data_line[i].x - lip_x);
        ray_length[i] = std::hypot(data_line[i].x - lip_x, data_line[i].y - lip_y);
    }
    // Invariants on ray i at a fraction f of the way from the lip to its marched point.
    auto ray_invariants = [&](size_t i, double f, double& K_minus, double& K_plus) {
        const double theta_lip = m_theta_schedule[i];
        const double K_minus_lip = 2.0 * theta_lip;   // nu = theta on a centered fan's rays
        const double K_plus_lip = 0.0;
        K_minus = K_minus_lip + f * (data_line[i].K_minus - K_minus_lip);
        K_plus = K_plus_lip + f * (data_line[i].K_plus - K_plus_lip);
    };

    const size_t n_points = std::max<size_t>(static_cast<size_t>(m_options.num_characteristics), 3);
    std::vector<CharacteristicPoint> front(n_points);
    for (size_t j = 0; j < n_points; j++) {
        const double y = y_wall * static_cast<double>(j) / static_cast<double>(n_points - 1);
        const double phi = std::atan2(y - lip_y, x0 - lip_x);
        const double dist = std::hypot(x0 - lip_x, y - lip_y);
        double K_minus, K_plus;
        if (phi <= ray_angle.front()) {
            ray_invariants(0, std::clamp(dist / ray_length[0], 0.0, 1.0), K_minus, K_plus);
        } else if (phi >= ray_angle.back()) {
            ray_invariants(n_rays - 1, std::clamp(dist / ray_length[n_rays - 1], 0.0, 1.0), K_minus, K_plus);
        } else {
            size_t i = 0;
            while (i + 1 < n_rays && ray_angle[i + 1] < phi) i++;
            const double alpha = (phi - ray_angle[i]) / (ray_angle[i + 1] - ray_angle[i]);
            double Km0, Kp0, Km1, Kp1;
            ray_invariants(i, std::clamp(dist / ray_length[i], 0.0, 1.0), Km0, Kp0);
            ray_invariants(i + 1, std::clamp(dist / ray_length[i + 1], 0.0, 1.0), Km1, Kp1);
            K_minus = Km0 + alpha * (Km1 - Km0);
            K_plus = Kp0 + alpha * (Kp1 - Kp0);
        }
        CharacteristicPoint pt{};
        pt.x = x0;
        pt.y = y;
        pt.theta = 0.5 * (K_minus + K_plus);
        throw_if_thermo_error(
            update_thermodynamic_state_from_nu(pt, 0.5 * (K_minus - K_plus), 1.05),
            "Inverse march initial front (fan plane)", x0, y);
        pt.K_plus = pt.theta - pt.nu;
        pt.K_minus = pt.theta + pt.nu;
        front[j] = pt;
    }
    // The axis point of the plane is the data line's own axis point (it lies on the plane).
    front.front() = data_line.front();

    // Wall point: the contour's own angle, nu from the C+ compatibility relation with the
    // point below (K+ carried up unchanged over the short vertical distance).
    CharacteristicPoint& wall_pt = front.back();
    const CharacteristicPoint below = front[n_points - 2];
    wall_pt = CharacteristicPoint{};
    wall_pt.x = x0;
    wall_pt.y = y_wall;
    wall_pt.theta = wall.theta_at_interpolated(x0);
    throw_if_thermo_error(
        update_thermodynamic_state_from_nu(wall_pt, wall_pt.theta - below.K_plus, below.mach),
        "Inverse march initial front (fan wall point)", wall_pt.x, wall_pt.y);
    wall_pt.K_plus = wall_pt.theta - wall_pt.nu;
    wall_pt.K_minus = wall_pt.theta + wall_pt.nu;
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

    MocErrorCode validity = check_point_validity(p4, m_options.solver_options.abstol, /*require_nonnegative_theta=*/false);
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

    MocErrorCode validity = check_point_validity(axis_point, m_options.solver_options.abstol, false);
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
    wall_point.theta = m_options.nozzle_profile.theta_at_interpolated(x_new); // fixed by the contour

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

    MocErrorCode validity = check_point_validity(wall_point, m_options.solver_options.abstol, false);
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

    // Normalized y-distribution of the initial front, reused for every subsequent front.
    std::vector<double> s(M);
    {
        const double y_wall0 = net.points[front.back()].y;
        for (size_t j = 0; j < M; j++) {
            s[j] = (y_wall0 > 0.0)
                ? net.points[front[j]].y / y_wall0
                : static_cast<double>(j) / static_cast<double>(M - 1);
        }
    }

    // The initial front must be spacelike: a segment leaning upstream (x decreasing with y)
    // must be steeper than the C- through either endpoint, one leaning downstream steeper
    // than the C+. Every later front is this front's shape translated downstream, stretched
    // in y with the wall radius and relaxed toward vertical, none of which can make it less
    // spacelike, so one check here covers the march.
    for (size_t j = 0; j + 1 < M; j++) {
        const CharacteristicPoint& a = net.points[front[j]];
        const CharacteristicPoint& b = net.points[front[j + 1]];
        const double dy = b.y - a.y;
        if (!(dy > 0.0)) {
            return MocFailure{MocErrorCode::INITIALIZATION_FAILED,
                std::format("Inverse march: initial front is not increasing in y at segment {}.", j),
                a.x, a.y, -1};
        }
        // The segment's direction must lie between the C+ direction (theta + mu) and the
        // backward C- direction (theta - mu + pi) of both endpoints.
        const double segment_angle = std::atan2(dy, b.x - a.x);
        const double lower_bound = std::max(a.theta + a.mu, b.theta + b.mu);
        const double upper_bound = std::min(a.theta - a.mu, b.theta - b.mu) + M_PI;
        if (!(segment_angle > lower_bound && segment_angle < upper_bound)) {
            return MocFailure{MocErrorCode::INITIALIZATION_FAILED,
                std::format("Inverse march: initial front segment {} (y {:.4f} to {:.4f}) is not "
                    "spacelike: its direction {:.2f} deg lies outside ({:.2f}, {:.2f}) deg.",
                    j, a.y, b.y, segment_angle * 180.0 / M_PI,
                    lower_bound * 180.0 / M_PI, upper_bound * 180.0 / M_PI),
                a.x, a.y, -1};
        }
    }

    constexpr int max_passes = 20000;
    for (int pass = 0; pass < max_passes; pass++) {
        const CharacteristicPoint& wall_old = net.points[front.back()];

        // The front's shape: every point's axial offset from the wall point.
        std::vector<double> offset(M);
        double max_abs_offset = 0.0;
        for (size_t j = 0; j < M; j++) {
            offset[j] = net.points[front[j]].x - wall_old.x;
            max_abs_offset = std::max(max_abs_offset, std::abs(offset[j]));
        }

        // --- Step 1: step length ---
        // Domain of dependence: a point advanced by dx traces back a height of at most
        // dx * tan(theta +/- mu), which must stay within about one spacing so the feet are
        // interpolated locally.
        double dy_min = std::numeric_limits<double>::max();
        double slope_max = 0.0;
        for (size_t j = 0; j < M; j++) {
            const CharacteristicPoint& p = net.points[front[j]];
            slope_max = std::max({slope_max, std::abs(std::tan(p.theta + p.mu)),
                                  std::abs(std::tan(p.mu - p.theta))});
            if (j + 1 < M) dy_min = std::min(dy_min, net.points[front[j + 1]].y - p.y);
        }
        const double dx_cfl = (slope_max > 1e-12)
            ? m_options.inverse_cfl * dy_min / slope_max
            : std::numeric_limits<double>::max();

        // Wall foot: the top interior point's C- must reach the previous front below its
        // wall point rather than the wall between the two wall points. For a point at
        // offset `off` behind (off < 0) or ahead of (off > 0) the old wall point, with the
        // wall rising at tan(theta_wall):
        //   dx * (1 + tan(theta_wall) cot(mu - theta)) <= dy_top cot(mu - theta) - (1 - lambda) off
        // The relaxation lambda is not known yet, so the bound uses whichever value is
        // conservative for the sign of the offset.
        // Only binding while the C- through the top point rises going backward (mu > theta);
        // once the flow angle exceeds the Mach angle its backward C- descends and cannot
        // reach the wall.
        double dx_foot = std::numeric_limits<double>::max();
        if (M >= 3 && net.points[front[M - 2]].mu - net.points[front[M - 2]].theta > 1e-9) {
            const CharacteristicPoint& top = net.points[front[M - 2]];
            const double dy_top = wall_old.y - top.y;
            const double cot_minus = 1.0 / std::tan(top.mu - top.theta);
            const double keep = (offset[M - 2] < 0.0) ? m_options.front_tilt_decay : 1.0;
            const double numerator = dy_top * cot_minus - keep * offset[M - 2];
            const double denominator = 1.0 + std::tan(wall_old.theta) * cot_minus;
            dx_foot = (numerator > 0.0 && denominator > 0.0) ? numerator / denominator : 0.0;
        }

        // Wall turn: the contour's own angle change ahead of the wall station, read at its
        // vertices (theta_at is piecewise constant per facet, so a one-step difference would
        // see either nothing or a whole facet jump).
        double dx_turn = std::numeric_limits<double>::max();
        {
            const double theta_here = wall_profile.theta_at_interpolated(wall_old.x);
            for (size_t k = 0; k < wall_profile.size(); k++) {
                const double x_v = wall_profile.x[k];
                if (x_v <= wall_old.x + 1e-12) continue;
                if (std::abs(wall_profile.theta_at_interpolated(x_v) - theta_here) > m_options.max_wall_turn_per_step) {
                    dx_turn = x_v - wall_old.x;
                    break;
                }
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
                std::format("Inverse march: non-positive step length ({:.3e}, limiter {}) at pass {}.",
                    dx, to_string(limiter), pass),
                wall_old.x, wall_old.y, pass};
        }

        // --- Step 2: new front geometry ---
        // The previous front's shape, translated by dx, stretched in y with the wall radius,
        // and relaxed toward a vertical plane by removing a fraction lambda of every offset.
        // lambda is capped so the relaxation moves no point by more than half a step: every
        // point advances by between dx/2 and 3dx/2, which is what the step bounds above
        // assume. (A fixed per-pass decay of the tilt, independent of dx, moved the axis end
        // *upstream* whenever the near-axis step was small, which is how the first version
        // of this kernel failed at pass 0 on a Kliegel-Levine front.)
        const double x_w = wall_old.x + dx;
        const double y_w = wall_profile.radius_at(std::min(x_w, L));
        double lambda = 1.0 - m_options.front_tilt_decay;
        if (max_abs_offset > 0.0) lambda = std::min(lambda, 0.5 * dx / max_abs_offset);

        std::vector<double> x_new(M), y_new(M);
        for (size_t j = 0; j < M; j++) {
            x_new[j] = x_w + (1.0 - lambda) * offset[j];
            y_new[j] = s[j] * y_w;
        }
        x_new[M - 1] = x_w;
        y_new[M - 1] = y_w;
        y_new[0] = 0.0;

        // --- Steps 3-6: the unit processes ---
        std::vector<CharacteristicPoint> new_front_points(M);
        {
            PointResult axis_result = solve_inverse_march_axis_point(x_new[0], net, front);
            if (axis_result.error != MocErrorCode::NONE) {
                return MocFailure{axis_result.error,
                    std::format("Inverse march axis point failed ({}) at pass {}.",
                        to_string(axis_result.error), pass),
                    axis_result.point.x, axis_result.point.y, pass};
            }
            new_front_points[0] = axis_result.point;
        }
        for (size_t j = 1; j + 1 < M; j++) {
            PointResult r = solve_inverse_march_interior_point(x_new[j], y_new[j], net, front);
            if (r.error != MocErrorCode::NONE) {
                return MocFailure{r.error,
                    std::format("Inverse march interior point {} failed ({}) at pass {}.",
                        j, to_string(r.error), pass),
                    r.point.x, r.point.y, pass};
            }
            new_front_points[j] = r.point;
        }
        {
            PointResult wall_result = solve_inverse_march_wall_point(x_new[M - 1], y_new[M - 1], net, front);
            if (wall_result.error != MocErrorCode::NONE) {
                return MocFailure{wall_result.error,
                    std::format("Inverse march wall point failed ({}) at pass {}.",
                        to_string(wall_result.error), pass),
                    wall_result.point.x, wall_result.point.y, pass};
            }
            new_front_points[M - 1] = wall_result.point;
        }

        log_debug("STEP pass={} dx_cfl={:.6f} dx_foot={:.6f} dx_turn={:.6f} dx_exit={:.6f} "
            "dx={:.6f} lambda={:.4f} x_w={:.6f} y_w={:.6f} x_a={:.6f}",
            pass, dx_cfl, dx_foot, dx_turn, dx_exit, dx, lambda, x_w, y_w, x_new[0]);

        // --- Step 7: bookkeeping ---
        std::vector<size_t> new_front = seed_inverse_front(net, new_front_points);

        MocPassDiagnostics diag;
        diag.pass = pass;
        diag.front_points = M;
        diag.step_dx = dx;
        diag.step_limiter = limiter;
        diag.front_axis_x = x_new[0];
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
            diag.front_axis_spacing = std::hypot(new_front_points[1].x - new_front_points[0].x,
                                                 new_front_points[1].y - new_front_points[0].y);
            diag.front_wall_spacing = std::hypot(new_front_points[M - 1].x - new_front_points[M - 2].x,
                                                 new_front_points[M - 1].y - new_front_points[M - 2].y);
        }
        m_pass_diagnostics.push_back(diag);

        log_debug("FRONT pass={} dx={:.6f} limiter={} x_axis={:.6f} x_wall={:.6f}",
            pass, dx, to_string(limiter), x_new[0], x_w);

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
