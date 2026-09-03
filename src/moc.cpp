#include <cmath>
#include <algorithm>
#include <format>
#include <vector>
#include <limits>
#include <stdexcept>
#include "goddard/moc.hpp"

namespace Goddard {


namespace {

/** Twice the signed area of triangle (a, b, c): > 0 counter-clockwise, < 0 clockwise. */
double orientation(double ax, double ay, double bx, double by, double cx, double cy) {
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
}

/**
 * True when two segments cross at a point interior to both.
 *
 * Deliberately strict: shared endpoints and collinear overlap do not count. Characteristic
 * chains meet legitimately at wall and axis reflection points, and only a genuine crossing
 * indicates that the family is coalescing.
 */
bool segments_properly_intersect(
    double ax, double ay, double bx, double by,
    double cx, double cy, double dx, double dy)
{
    const double d1 = orientation(ax, ay, bx, by, cx, cy);
    const double d2 = orientation(ax, ay, bx, by, dx, dy);
    const double d3 = orientation(cx, cy, dx, dy, ax, ay);
    const double d4 = orientation(cx, cy, dx, dy, bx, by);
    return ((d1 > 0.0 && d2 < 0.0) || (d1 < 0.0 && d2 > 0.0)) &&
           ((d3 > 0.0 && d4 < 0.0) || (d3 < 0.0 && d4 > 0.0));
}

struct CharacteristicSegment {
    double x0, y0, x1, y1;
    double x_lo, x_hi;
    size_t pt0, pt1;
    size_t chain;
};

} // namespace

MocCrossings find_like_characteristic_crossings(const CharacteristicNet& net) {
    MocCrossings result;

    // A badly broken net can in principle produce enormous numbers of crossings; the count is
    // a diagnostic, not a census, so stop once the answer is unambiguous.
    constexpr size_t max_reported = 10000;
    double first_x = std::numeric_limits<double>::max();

    const ChainMetadata::Family families[2] = {
        ChainMetadata::Family::PLUS, ChainMetadata::Family::MINUS};

    for (ChainMetadata::Family family : families) {
        std::vector<CharacteristicSegment> segments;
        for (size_t c = 0; c < net.c_chains.size(); c++) {
            if (net.chain_metadata[c].family != family) continue;
            const std::vector<size_t>& chain = net.c_chains[c];
            for (size_t k = 0; k + 1 < chain.size(); k++) {
                const CharacteristicPoint& a = net.points[chain[k]];
                const CharacteristicPoint& b = net.points[chain[k + 1]];
                segments.push_back(CharacteristicSegment{
                    a.x, a.y, b.x, b.y,
                    std::min(a.x, b.x), std::max(a.x, b.x),
                    chain[k], chain[k + 1], c});
            }
        }

        // Sorting by the left end lets the inner loop stop as soon as the x-ranges separate,
        // which is what keeps this near-linear on a nozzle: characteristics span a small part
        // of the axial extent each.
        std::sort(segments.begin(), segments.end(),
            [](const CharacteristicSegment& a, const CharacteristicSegment& b) {
                return a.x_lo < b.x_lo;
            });

        for (size_t i = 0; i < segments.size() && result.count < max_reported; i++) {
            const CharacteristicSegment& si = segments[i];
            for (size_t j = i + 1; j < segments.size(); j++) {
                const CharacteristicSegment& sj = segments[j];
                if (sj.x_lo > si.x_hi) break;
                if (si.chain == sj.chain) continue;
                if (si.pt0 == sj.pt0 || si.pt0 == sj.pt1 ||
                    si.pt1 == sj.pt0 || si.pt1 == sj.pt1) continue;

                if (segments_properly_intersect(si.x0, si.y0, si.x1, si.y1,
                                                sj.x0, sj.y0, sj.x1, sj.y1)) {
                    result.count++;
                    const double cross_x = 0.25 * (si.x0 + si.x1 + sj.x0 + sj.x1);
                    if (cross_x < first_x) {
                        first_x = cross_x;
                        result.first_x = cross_x;
                        result.first_y = 0.25 * (si.y0 + si.y1 + sj.y0 + sj.y1);
                        result.first_family = family;
                    }
                }
            }
        }
    }
    return result;
}

void validate_moc_options(const MocOptions& options) {
    if (!(options.max_front_spacing_factor > 1.0)) {
        throw std::invalid_argument(std::format(
            "max_front_spacing_factor must be greater than 1; got {}. At or below 1 nearly "
            "every front segment triggers refinement and the march does not terminate.",
            options.max_front_spacing_factor));
    }
    if (!(options.min_front_spacing_factor > 0.0)) {
        throw std::invalid_argument(std::format(
            "min_front_spacing_factor must be positive; got {}.",
            options.min_front_spacing_factor));
    }
    if (!(options.min_front_spacing_factor < 0.5 * options.max_front_spacing_factor)) {
        throw std::invalid_argument(std::format(
            "min_front_spacing_factor ({}) must be below half of max_front_spacing_factor "
            "({}); overlapping refine and coarsen bands thrash, inserting and retiring the "
            "same rung on alternate passes.",
            options.min_front_spacing_factor, options.max_front_spacing_factor));
    }
    if (!(options.initial_line_clustering >= -1.0 && options.initial_line_clustering <= 1.0)) {
        throw std::invalid_argument(std::format(
            "initial_line_clustering must lie in [-1, 1]; got {}.", options.initial_line_clustering));
    }
    if (!(options.front_spacing_growth >= 0.0 && options.front_spacing_growth <= 1.0)) {
        throw std::invalid_argument(std::format(
            "front_spacing_growth must lie in [0, 1]; got {}.", options.front_spacing_growth));
    }
    if (!(options.max_cell_aspect_ratio > 1.0)) {
        throw std::invalid_argument(std::format(
            "max_cell_aspect_ratio must be greater than 1; got {}. Healthy fronts routinely "
            "reach 2-5, so a bound at or below 1 retires the entire front.",
            options.max_cell_aspect_ratio));
    }
    if (options.num_characteristics < 3) {
        throw std::invalid_argument(std::format(
            "num_characteristics must be at least 3; got {}.", options.num_characteristics));
    }
    if (options.max_front_points != 0 &&
        options.max_front_points <= static_cast<size_t>(options.num_characteristics)) {
        throw std::invalid_argument(std::format(
            "max_front_points ({}) must exceed num_characteristics ({}), or be 0 to select "
            "the default.", options.max_front_points, options.num_characteristics));
    }
    if (!(options.solver_options.abstol > 0.0)) {
        throw std::invalid_argument(std::format(
            "solver_options.abstol must be positive; got {}.", options.solver_options.abstol));
    }
    if (options.march_scheme == MocMarchScheme::INVERSE &&
        options.mode == MocMode::DESIGN_MIN_LENGTH) {
        throw std::invalid_argument(
            "MocMarchScheme::INVERSE cannot be used with MocMode::DESIGN_MIN_LENGTH: a "
            "minimum-length contour is defined by the characteristics the DIRECT kernel "
            "absorbs at the wall, which INVERSE's prescribed fronts cannot reproduce.");
    }
    if (!(options.inverse_cfl > 0.0 && options.inverse_cfl <= 1.0)) {
        throw std::invalid_argument(std::format(
            "inverse_cfl must lie in (0, 1]; got {}.", options.inverse_cfl));
    }
    if (!(options.front_tilt_decay >= 0.0 && options.front_tilt_decay <= 1.0)) {
        throw std::invalid_argument(std::format(
            "front_tilt_decay must lie in [0, 1]; got {}.", options.front_tilt_decay));
    }
    if (!(options.max_wall_turn_per_step > 0.0)) {
        throw std::invalid_argument(std::format(
            "max_wall_turn_per_step must be positive; got {}.", options.max_wall_turn_per_step));
    }
}

ThrustCoefficient compute_thrust_coefficient(
    const MocResult& result,
    MocFlowKind flow_type,
    double ambient_pressure_ratio)
{
    const auto& ep = result.exit_plane;
    size_t n = ep.y.size();

    if (n < 2) {
        return {0.0, 0.0, 0.0, 0.0};
    }

    double momentum_integral = 0.0;
    double pressure_integral = 0.0;

    // Compute integrand components at each point
    // f_momentum = p * gamma_s * M^2 * cos^2(theta)
    // f_pressure = p
    // Total integrand: f = f_momentum + f_pressure
    // Integration measure depends on flow type:
    //   Planar:       dy
    //   Axisymmetric: 2 * y * dy  (pi cancels with pi*y_t^2 in denominator)

    for (size_t i = 0; i < n - 1; i++) {
        double f_mom_i = ep.pressure[i] * ep.gamma_s[i]
            * ep.mach[i] * ep.mach[i]
            * cos(ep.theta[i]) * cos(ep.theta[i]);
        double f_mom_next = ep.pressure[i+1] * ep.gamma_s[i+1]
            * ep.mach[i+1] * ep.mach[i+1]
            * cos(ep.theta[i+1]) * cos(ep.theta[i+1]);

        double f_pres_i = ep.pressure[i];
        double f_pres_next = ep.pressure[i+1];

        double dy = ep.y[i+1] - ep.y[i];

        if (flow_type == MocFlowKind::PLANAR) {
            // Trapezoidal rule
            momentum_integral += 0.5 * (f_mom_i + f_mom_next) * dy;
            pressure_integral += 0.5 * (f_pres_i + f_pres_next) * dy;
        } else {
            // Axisymmetric: integrand includes 2*y factor
            // Trapezoidal rule for integral of f(y)*y*dy
            momentum_integral += (f_mom_i * ep.y[i] + f_mom_next * ep.y[i+1]) * dy;
            pressure_integral += (f_pres_i * ep.y[i] + f_pres_next * ep.y[i+1]) * dy;
        }
    }

    // Normalize by throat area (p0 * A_throat)
    // Nondimensional: p0 = 1, y_throat = 1.
    // Planar: F = 2*integral (symmetry), A_throat = 2*y_t = 2.0  =>  Cf = integral
    // Axisymmetric: F = 2*pi*integral(f*y*dy), A_throat = pi*y_t^2 = pi
    //   => Cf = 2*integral(f*y*dy). The integrand already includes the 2*y factor.
    double Cf_momentum = momentum_integral;
    double Cf_pressure = pressure_integral;

    double Cf_vacuum = Cf_momentum + Cf_pressure;
    double Cf = Cf_vacuum - ambient_pressure_ratio * result.area_ratio;

    return {Cf_vacuum, Cf, Cf_momentum, Cf_pressure};
}


MocFrontShear summarize_front_shear(
    const std::vector<MocPassDiagnostics>& pass_diagnostics, size_t window)
{
    MocFrontShear shear;
    if (window == 0 || pass_diagnostics.size() < 2 * window) return shear;

    // Mean axial step per pass at each end, over the first and last `window` passes.
    // Steps are differences between consecutive passes, so a window of w passes yields
    // w-1 steps; require at least one.
    if (window < 2) return shear;

    auto mean_step = [&](size_t begin, size_t end, double MocPassDiagnostics::* member) {
        double total = 0.0;
        size_t count = 0;
        for (size_t i = begin + 1; i < end; i++) {
            const double step = pass_diagnostics[i].*member - pass_diagnostics[i - 1].*member;
            // The front's ends are rebuilt every pass and can retreat when a chain retires;
            // only forward motion is a marching step.
            if (step > 0.0) {
                total += step;
                count++;
            }
        }
        return count == 0 ? 0.0 : total / static_cast<double>(count);
    };

    const size_t n = pass_diagnostics.size();
    const double axis_early = mean_step(0, window, &MocPassDiagnostics::front_axis_x);
    const double axis_late = mean_step(n - window, n, &MocPassDiagnostics::front_axis_x);
    const double wall_early = mean_step(0, window, &MocPassDiagnostics::front_wall_x);
    const double wall_late = mean_step(n - window, n, &MocPassDiagnostics::front_wall_x);

    if (axis_early <= 0.0 || wall_early <= 0.0) return shear;

    shear.axis_growth = axis_late / axis_early;
    shear.wall_growth = wall_late / wall_early;

    const double lo = std::min(shear.axis_growth, shear.wall_growth);
    const double hi = std::max(shear.axis_growth, shear.wall_growth);
    // A stalled end (growth 0) is maximal shear, not an undefined ratio.
    shear.shear_ratio = (lo > 0.0) ? hi / lo : std::numeric_limits<double>::infinity();
    shear.valid = true;
    return shear;
}

} // namespace Goddard