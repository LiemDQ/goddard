#include <algorithm>
#include <format>
#include <limits>
#include "goddard/gas_dynamics.hpp"
#include "goddard/moc_initialization.hpp"
#include "goddard/moc_unit_processes.hpp"

namespace Goddard {

KlWallAngleFallback::KlWallAngleFallback(double R_in, double miss_in, double threshold_in)
    : std::runtime_error(std::format(
          "Kliegel-Levine series misses the wall angle by {} rad at R = {} "
          "(threshold {} rad).", miss_in, R_in, threshold_in)),
      R(R_in), miss(miss_in), threshold(threshold_in)
{}

double gquad(double gamma, double a, double b, double c, double d) {
    return (a*gamma*gamma + b*gamma + c)/d;
}

namespace {
// MocInitialization communicates failure by exception (its public methods return a plain
// std::vector<CharacteristicPoint>, with no room for an error code), so a nonzero MocThermo
// chokepoint result has to be converted back into one here. The mapping matches what
// MocNozzle::solve()'s existing exception boundary already does for every other initializer
// failure, so a chokepoint failure reached during start-line construction still surfaces as
// the same MocErrorCode it always has: TABLE_RANGE_EXCEEDED throws std::out_of_range (solve()
// maps that type to MocErrorCode::TABLE_RANGE_EXCEEDED); anything else throws
// Goddard::ConvergenceError (maps to MocErrorCode::INITIALIZATION_FAILED).
void throw_if_thermo_error(MocErrorCode code, std::string_view context) {
    if (code == MocErrorCode::NONE) return;
    if (code == MocErrorCode::TABLE_RANGE_EXCEEDED) {
        throw std::out_of_range(std::format("{}: {}", context, to_string(code)));
    }
    throw ConvergenceError(std::format("{}: {}", context, to_string(code)));
}
} // namespace

MocInitialization::MocInitialization(
    NozzleGeometry geom,
    const MocThermo& thermo,
    const MocOptions& options)
    : geometry(geom), m_thermo(&thermo), m_options(options)
{

}

void MocInitialization::set_state_from_critical_velocity_ratio(
    CharacteristicPoint& pt, double m_star, const ThroatCondition& throat)
{
    if (m_thermo->chemistry() == GasChemistry::PERFECT_GAS) {
        // The perfect-gas path keeps velocity dimensionless and equal to the Mach number,
        // so the conversion has to happen here rather than being folded into a velocity.
        throw_if_thermo_error(
            m_thermo->set_state_from_mach(pt, mach_from_critical_velocity_ratio(m_star, throat.gamma_s)),
            "Transonic start line (critical velocity ratio, perfect gas)");
        return;
    }

    // For a real gas the M*<->M algebra is not available in closed form, but its defining
    // relation is: a* is the throat speed of sound, so V follows directly and the
    // PrandtlMeyerTable can be looked up by velocity.
    throw_if_thermo_error(
        m_thermo->set_state_from_V(pt, m_star * throat.speed_of_sound),
        "Transonic start line (critical velocity ratio, real gas)");
}

std::vector<CharacteristicPoint> MocInitialization::initialize_sauer(const ThroatCondition& throat) {
    
    // predictor step: assume gamma is constant throughout
    double gamma = throat.gamma_s;
    size_t num_points = static_cast<size_t>(m_options.num_characteristics);
    double dy = 1.0/(num_points-1);
    double delt = delta();
    double alpha = sauer_alpha(gamma);

    std::vector<CharacteristicPoint> points(num_points);

    for (size_t i = 0; i < num_points; i++) {
        CharacteristicPoint& pt = points[i];
        pt.y = dy * i;
        pt.x = (gamma + 1)*alpha/(2*(3+delta()))*(1.0 - pt.y*pt.y);
        pt.theta = 0.0; // by definition
        // Sauer's series, like Hall's and Kliegel-Levine's, gives the axial velocity
        // normalized by the critical speed of sound: this is u/a* = M*, not the Mach
        // number. On the transonic line the radial component vanishes by construction,
        // so the magnitude is the axial component alone.
        double u_star = 1+alpha*pt.x + (gamma+1)*alpha*alpha*pt.y*pt.y/(2*(1+delt));
        set_state_from_critical_velocity_ratio(pt, u_star, throat);
        pt.update_Ks();
    }

    //TODO: optional corrector step with updated gamma

    return points;
}

std::vector<CharacteristicPoint> MocInitialization::initialize_kliegel_levine(const ThroatCondition& throat) {
    double gamma = throat.gamma_s;
    size_t num_points = static_cast<size_t>(m_options.num_characteristics);
    double R = KL_R();
    double alpha = sauer_alpha(gamma);

    // Lift the sonic (zero-radial-velocity) locus off itself by a rigid downstream shift so
    // it can be seeded with both characteristic families (see MocOptions::initial_line_axial_shift
    // and CharacteristicNet::add_initial_data_line's nullopt branch).
    const double x_shift = m_options.initial_line_axial_shift;

    auto station_x = [&](double y) {
        double x_guess = (gamma + 1)*alpha/(2*(3+delta()))*(1.0 - y*y); //Sauer transonic x-coordinate
        return KL_solve_transonic_x(y, gamma, R, x_guess) + x_shift;
    };

    // Where the start line meets the prescribed wall.
    //
    // The line is the sonic locus translated rigidly downstream, so its top station sits at
    // some x > 0 -- where the wall has already expanded past the throat radius. Pinning that
    // station at y = throat_radius therefore places it *inside* the flow, yet
    // CharacteristicNet::add_initial_data_line registers it as the wall point regardless, so
    // the wall march begins off the wall with an uninitialized sliver above the line. The
    // gap is geometric and fixed, so refining the grid shrinks the cells around an error
    // that does not shrink with them.
    //
    // Solve y = wall_radius(station_x(y)) instead. The fixed point converges quickly because
    // station_x varies weakly with y near the wall and the contour's slope is bounded. The
    // KL series is evaluated at r slightly greater than 1 in the process, which is correct
    // rather than an extrapolation: the throat arc's own wall lies at r > 1 for x > 0.
    const NozzleProfile& wall_profile = m_options.nozzle_profile;
    double y_wall = 1.0;
    if (wall_profile.size() >= 2) {
        for (int iteration = 0; iteration < 10; iteration++) {
            const double x_station = station_x(y_wall);
            if (!(x_station >= wall_profile.x_min() && x_station <= wall_profile.x_max())) break;
            const double y_next = wall_profile.radius_at(x_station);
            if (!(y_next > 0.0) || !std::isfinite(y_next)) break;
            const double delta_y = std::abs(y_next - y_wall);
            y_wall = y_next;
            if (delta_y < 1e-12) break;
        }
    }

    auto evaluate = [&](double y) {
        CharacteristicPoint pt{};
        pt.y = y;
        pt.x = station_x(y);
        // Both series are velocity components normalized by the *critical* speed of sound
        // a*, i.e. they are the components of M* = V/a*, not of the Mach number. See
        // Kliegel & Levine Eqs. (10) and (12), which these reduce to at the throat plane
        // and which the paper states for u/a*.
        const double u_star = KL_xMach(y, KL_z_coordinate(pt.x, gamma), gamma, R);
        // Off the sonic locus, the radial velocity component (KL_yMach) is
        // generally nonzero, so theta can no longer be pinned to 0. Both series are
        // normalized by the same a*, so their ratio gives the flow angle directly.
        const double v_star = KL_yMach(pt.x, y, gamma, R);
        pt.theta = std::atan2(v_star, u_star);
        set_state_from_critical_velocity_ratio(pt, std::hypot(u_star, v_star), throat);
        pt.update_Ks();
        return pt;
    };

    // Axial station at which the C- leaving height y reaches the axis, on a straight-ray
    // estimate, is the coordinate that would need to be uniform to equalize characteristic
    // spacing near the axis -- not y.
    //
    // Spacing points uniformly in y spaces their characteristics ~7x non-uniformly here,
    // because the near-axis region is a double zero: y -> 0 and cot(mu) -> 0 together (the
    // start line is near-sonic on the axis, so mu -> 90 deg). The C- from the lower third of
    // the line therefore all arrive within a few percent of a throat radius of each other,
    // the axis consumes them almost at once, and thereafter is resupplied only by wall
    // reflections -- so its point spacing jumps by an order of magnitude and the marching
    // front shears until pairing breaks down. The ratio is set by the flow, not the mesh, so
    // it is *independent of num_characteristics*: refining the grid halves every gap and
    // leaves the grading (measured 6.2 to 7.2 for N = 8 to 61) intact. That is why grid
    // refinement never cured the axisymmetric breakdown -- and why clustering the line
    // toward the axis to compensate was tried and abandoned (the benefit does not survive
    // refinement); interior stations are simply uniform in y.
    std::vector<CharacteristicPoint> points(num_points);
    const double last = static_cast<double>(num_points - 1);

    for (size_t i = 0; i < num_points; i++) {
        double y;
        if (i == 0) {
            y = 0.0;                       // anchor the axis point exactly
        }
        else if (i == num_points - 1) {
            y = 1.0;                       // anchor the wall point on the contour
        }
        else {
            y = static_cast<double>(i) / last;
        }
        // The interior stations are laid out on the unit interval; stretch them onto the
        // line's actual span so they stay evenly distributed when the wall end moves.
        points[i] = evaluate(y * y_wall);
    }

    // Wall-consistent correction. The wall point's theta was left as the raw series value
    // above (evaluate() has no notion of the contour), and that raw value is frequently
    // wrong by a lot: at the default r_arc = 0.382 the series recovers only ~40% of the
    // contour's own wall angle there (diagnosis.md Sec A2 -- the truncated z-dependent
    // terms of a third-order expansion do not converge for a throat this sharply curved).
    // Leaving that mismatch in place, as an earlier version of this function argued for,
    // makes K+ = theta - nu on the near-wall interior points 8-10 deg too negative, which
    // the first wall solve reads as an over-expansion and the wall Mach dips for the next
    // several wall points before recovering (Sec A1) -- refining the grid does not cure it,
    // since the mismatch is a property of the series, not of the mesh.
    //
    // Measure the raw mismatch first (before touching any theta), then either let the
    // caller fall back to the centered fan (AUTO, above threshold -- the series cannot be
    // trusted here) or spread the mismatch across the line so the wall end matches the
    // contour exactly.
    //
    // The correction is multiplicative (theta_i *= theta_wall / theta_series_wall), not a
    // uniform or linearly-graded additive offset: three variants were measured against
    // conical AR=4, r_arc=0.382, N=15/31/61, mesh control non-binding (see the package
    // report for the full table). An additive correction weighted by (y_i/y_wall) -- p=1 --
    // reduces the wall-Mach dip diagnosis.md Sec A1 documents but leaves it 7 wall points
    // deep at N=61 (peak-to-trough 0.035 in Mach) with the largest single-step regression
    // (-0.015). Squaring the weight (p=2) or scaling multiplicatively both concentrate the
    // correction much more sharply at the wall end, where the series is weakest and
    // irrotationality ties every interior theta to the same u/v polynomials the wall angle
    // was read from; either shrinks the dip to 4 points at N=61 (0.014 peak-to-trough,
    // largest single-step -0.006) and to nothing at all at N=15. Multiplicative was kept
    // over p=2 because it gave the smaller |mass_flow_error| at every N tested.
    if (num_points >= 2 && wall_profile.size() >= 2) {
        CharacteristicPoint& wall_point = points.back();
        if (wall_point.x >= wall_profile.x_min() && wall_point.x <= wall_profile.x_max()) {
            const double theta_series_wall = wall_point.theta;
            const double theta_wall = wall_profile.theta_at(wall_point.x);
            const double miss = std::abs(theta_series_wall - theta_wall);
            last_wall_bc_residual =
                (theta_wall != 0.0) ? (1.0 - theta_series_wall / theta_wall) : 0.0;

            if (miss > m_options.kl_max_wall_angle_error) {
                if (m_options.start_line == MocStartLine::AUTO) {
                    // Caller (build_start_line) catches this specific
                    // type and rebuilds the line as a centered fan instead.
                    throw KlWallAngleFallback(R, miss, m_options.kl_max_wall_angle_error);
                }
                // Forced MocStartLine::KLIEGEL_LEVINE: fail honestly rather than silently
                // stretching a correction over a mismatch this large.
                throw ConvergenceError(std::format(
                    "Kliegel-Levine series misses the wall angle by {} rad at R = {} "
                    "(threshold kl_max_wall_angle_error = {} rad); forced "
                    "MocStartLine::KLIEGEL_LEVINE does not fall back to the centered fan.",
                    miss, R, m_options.kl_max_wall_angle_error));
            }

            // theta is 0 only on the axis (v* vanishes there identically), where the
            // multiplier would be a 0/0 indeterminate; the axis point needs no correction
            // in any case (it is already exact), so it is simply skipped.
            if (std::abs(theta_series_wall) > 1e-12) {
                const double multiplier = theta_wall / theta_series_wall;
                for (CharacteristicPoint& pt : points) {
                    if (std::abs(pt.theta) <= 1e-12) continue;
                    pt.theta *= multiplier;
                    pt.update_Ks();
                }
            }
        }
    }
    // Else: no wall profile to check against (e.g. MocInitialization exercised directly,
    // without a MocNozzle solve() around it) -- nothing to correct or measure against, so
    // the series' own theta is returned unmodified, exactly as before this change.

    //TODO: corrector step with updated gamma

    return points;
}


std::vector<CharacteristicPoint> MocInitialization::initialize_centered_expansion(const ThroatCondition& throat) {
    size_t num_points = static_cast<size_t>(m_options.num_characteristics);
    std::vector<CharacteristicPoint> points;
    points.reserve(num_points);

    CharacteristicPoint sonic_point;
    sonic_point.temperature = 1.0; // should be 1.0 by definition
    sonic_point.pressure = 1.0;
    sonic_point.cantera_state = throat.state;
    sonic_point.mach = 1.0; // by definition
    sonic_point.theta = 0.0; // by definition
    sonic_point.nu = 0.0;
    sonic_point.x = 0.0;
    sonic_point.y = 1.0; // dimensionless throat radius
    sonic_point.K_minus = 0.0;
    sonic_point.K_plus = 0.0;

    std::vector<double> theta_schedule;

    if (!m_options.theta_schedule.empty()) {
        theta_schedule = m_options.theta_schedule;
    }
    else {
        double dtheta_initial = m_options.theta_max/(num_points*10);
        double dtheta = (m_options.theta_max - dtheta_initial)/(num_points - 1);
        theta_schedule.resize(num_points);
        for (size_t i = 0; i < num_points; i++) {
            theta_schedule[i] = dtheta_initial + i * dtheta;
        }
    }

    for (size_t i = 0; i < theta_schedule.size(); i++) {
        CharacteristicPoint pt;
        pt.x = sonic_point.x;
        pt.y = sonic_point.y; // dimensionless throat radius
        pt.theta = theta_schedule[i];

        throw_if_thermo_error(
            m_thermo->set_state_from_nu(pt, pt.theta, 1.0),
            "Centered-fan initialization");
        pt.update_Ks();
        points.push_back(pt);
    }
    return points;
}

double MocInitialization::KL_z_coordinate(double x, double gamma) const {
    // x is dimensionless (in throat radii), as are all net coordinates.
    double R = KL_R();
    return std::sqrt(2*R/(gamma+1))*x;
}

double MocInitialization::KL_dzdx(double, double gamma) const {
    double R = KL_R();
    return std::sqrt(2*R/(gamma+1));
}

double MocInitialization::KL_R() const {
    return geometry.downstream_wall_curvature_radius/geometry.throat_radius;
}

double MocInitialization::KL_u1(double r, double z) const {
        return 0.5*r*r - 0.25 + z;
    }

double MocInitialization::KL_u2(double r, double z, double gamma) const {
    double r2 = r*r;
    return (2*gamma+9)/24*r2*r2 - (4*gamma+15)/24*r2 + (10*gamma + 57)/288 + z*(r2 - 5.0/8.0)-(2*gamma-3)/6*z*z;
}

double MocInitialization::KL_u3(double r, double z, double gamma) const {
    double r2 = r*r;
    double r4 = r2*r2;
    
    return gquad(gamma, 556, 1737, 3069, 10368)*r4*r2
    - gquad(gamma, 388, 1161, 1881, 2304)*r4
    + gquad(gamma, 304, 831, 1242, 1728)*r2
    - gquad(gamma, 2708, 7839, 14211, 82944)
    + z*(gquad(gamma, 52, 51, 327, 384)*r4 - gquad(gamma, 52, 75, 279, 192)*r2 + gquad(gamma, 92, 180, 639, 1152))
    + z*z*(-(7*gamma-3)/8*r2 + (13*gamma-27)/48)
    + z*z*z*gquad(gamma, 4, -57, 27, 144);
}

double MocInitialization::KL_v1(double r, double z) const {
    return 0.25*r*r*r - 0.25*r + r*z;
}

double MocInitialization::KL_v2(double r, double z, double gamma) const {
    double r2 = r*r;
    
    return (gamma+3)/9*r2*r2*r 
    - (20*gamma+63)/96*r2*r 
    + (28*gamma+93)/288*r 
    + z*((2*gamma+9)/6*r2*r - (4*gamma+15)/12*r)
    + z*z*r;
}

double MocInitialization::KL_v3(double r, double z, double gamma) const {
    double r2 = r*r;
    double r3 = r*r2;
    return gquad(gamma, 6836, 23031, 30627, 82944)*r3*r2*r2
        - gquad(gamma, 3380, 11391, 15291, 13824)*r3*r2
        + gquad(gamma, 3424, 11271, 15228, 13824)*r3
        - gquad(gamma, 7100, 22311, 30249, 82944)*r
        + z*(gquad(gamma, 556, 1737, 3069, 1728)*r3*r2 - gquad(gamma, 388, 1161, 1881, 576)*r3 + gquad(gamma, 304, 831, 1242, 864)*r)
        + z*z*(gquad(gamma, 52, 51, 327, 192)*r3 - gquad(gamma, 52, 75, 279, 192)*r)
        - z*z*z*gquad(gamma, 0, 7, -3, 12)*r;
}

double MocInitialization::KL_dv1dz(double r, double) const {
    return r;
}

double MocInitialization::KL_dv2dz(double r, double z, double gamma) const {
    double r2 = r*r;
    return ((2*gamma+9)/6*r2*r - (4*gamma+15)/12*r) + 2*z*r;
}

double MocInitialization::KL_dv3dz(double r, double z, double gamma) const {
    double r2 = r*r;
    double r3 = r*r*r;

    return (gquad(gamma, 556, 1737, 3069, 1728)*r3*r2 - gquad(gamma, 388, 1161, 1881, 576)*r3 + gquad(gamma, 304, 831, 1242, 864)*r)
        + 2*z*(gquad(gamma, 52, 51, 327, 192)*r3 - gquad(gamma, 52, 75, 279, 192)*r)
        - 3*z*z*gquad(gamma, 0, 7, -3, 12)*r;
}

double MocInitialization::KL_xMach(double r, double z, double gamma, double R) const {
    double denom = 1.0/(R+1);
    double u1 = KL_u1(r, z);
    double u2 = KL_u2(r, z, gamma);
    double u3 = KL_u3(r, z, gamma);
    return 1 + denom*u1 + denom*denom * (u1 + u2) + denom*denom*denom*(u1 + 2*u2 + u3);
}

double MocInitialization::KL_axisymmetric_coeff(double gamma, double R) const {
    return std::sqrt((gamma + 1)/(2*(R+1)));
}

double MocInitialization::KL_yMach(double x, double y, double gamma, double R) const {
    double r = y; // dimensionless (throat radii)
    double z = KL_z_coordinate(x, gamma);
    double denom = 1.0/(R+1);
    double v1 = KL_v1(r, z);
    double v2 = KL_v2(r, z, gamma);
    double v3 = KL_v3(r, z, gamma);
    //TODO: this may only be valid for axisymmetric geometries. Verify.
    double factor = std::sqrt((gamma+1)/(2*(R+1)));

    // The radial-velocity series has no constant term: every v_i vanishes on the
    // axis (v_i ~ r), so v must too.
    return factor*(denom*v1 + denom*denom * (1.5*v1 + v2) + denom*denom*denom*(15.0/8*v1 + 5.0/2*v2 + v3));
}

double MocInitialization::KL_dyMachdx(double x, double y, double gamma, double R) const {
    double r = y; // dimensionless (throat radii)
    double dzdx = KL_dzdx(x, gamma);
    double z = KL_z_coordinate(x, gamma);
    double factor = std::sqrt((gamma+1)/(2*(R+1)));
    double denom = 1.0/(R+1);
    double dv1dx = KL_dv1dz(r, z)*dzdx;
    double dv2dx = KL_dv2dz(r, z, gamma)*dzdx;
    double dv3dx = KL_dv3dz(r, z, gamma)*dzdx;

    return factor*(denom*dv1dx
        + denom*denom * (1.5*dv1dx + dv2dx)
        + denom*denom*denom*(15.0/8*dv1dx + 5.0/2*dv2dx + dv3dx));
}

double MocInitialization::KL_solve_transonic_x(double y, double gamma, double R, double x_guess) const {
    // Solve for x where the radial velocity KL_yMach = 0 using Newton's method.
    // The radial velocity series is exactly proportional to r (every v_i ~ r), so
    // the residual is normalized by r to make the tolerance r-independent, and
    // on-axis stations (y = 0) are evaluated at a small finite r, where the
    // normalized residual converges to the r->0 limit of the v = 0 locus.
    double r = std::max(y, 1e-4);
    double residual = 1.0;
    double tol = 1e-6;
    double x = x_guess;
    int max_iters = 30;
    for (int i = 0; i <= max_iters; i++) {

        residual = KL_yMach(x, r, gamma, R) / r;
        if (std::abs(residual) < tol)
            return x;

        x = x - residual/(KL_dyMachdx(x, r, gamma, R) / r);
    }

    throw ConvergenceError("Newton's method for transonic line x-coordinate failed to converge.", max_iters, tol);
}

namespace {

// The centered-fan march loop, written once: MocMode::DESIGN_MIN_LENGTH and the
// ANALYSIS/DESIGN_RAO fan branch both ran this identical loop (previously duplicated).
std::vector<CharacteristicPoint> march_centered_fan(
    const MocSolveContext& ctx, const std::vector<CharacteristicPoint>& expansion_line)
{
    std::vector<CharacteristicPoint> data_line;
    data_line.reserve(expansion_line.size());
    CharacteristicPoint upstream_point;
    for (size_t i = 0; i < expansion_line.size(); i++) {
        const auto& expansion_point = expansion_line[i];
        // special case: the first point is on the centerline
        if (i == 0) {
            upstream_point = unwrap_or_throw(
                solve_fan_axis_point(ctx, expansion_point),
                "Initial data line (centered fan, axis point)");
        }
        else {
            upstream_point = unwrap_or_throw(
                solve_interior_point(ctx, expansion_point, upstream_point),
                "Initial data line (centered fan, interior point)");
        }
        data_line.push_back(upstream_point);
    }
    return data_line;
}

/**
 * Relative error of the mass flow carried across the initial data line, against the 1-D
 * critical mass flow through the throat. See MocInitDiagnostics::mass_flow_error.
 */
double start_line_mass_flow_error(
    const std::vector<CharacteristicPoint>& data_line, const MocSolveContext& ctx, double throat_radius)
{
    // The integrand needs density, which for frozen/equilibrium would take the molar mass
    // that CharacteristicPoint does not carry. Report nothing rather than something wrong.
    if (ctx.options.chemistry != GasChemistry::PERFECT_GAS) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (data_line.size() < 2) return std::numeric_limits<double>::quiet_NaN();

    const double gamma = ctx.options.gamma;
    if (!(gamma > 1.0)) return std::numeric_limits<double>::quiet_NaN();

    // rho*V normalized by its stagnation reference rho0*a0, from the isentropic relations:
    //   rho/rho0 = (1 + (g-1)/2 M^2)^(-1/(g-1)),  a/a0 = (1 + (g-1)/2 M^2)^(-1/2)
    // so rho*V/(rho0*a0) = M * (1 + (g-1)/2 M^2)^(-(g+1)/(2(g-1))).
    const double exponent = -(gamma + 1.0) / (2.0 * (gamma - 1.0));
    auto mass_flux = [&](double mach) {
        return mach * std::pow(1.0 + 0.5 * (gamma - 1.0) * mach * mach, exponent);
    };

    const bool axisymmetric = (ctx.options.flow_type == MocFlowKind::AXISYMMETRIC);
    const double r_throat = (throat_radius > 0.0) ? throat_radius : 1.0;

    double mdot = 0.0;
    for (size_t i = 1; i < data_line.size(); i++) {
        const CharacteristicPoint& below = data_line[i - 1];
        const CharacteristicPoint& above = data_line[i];
        const double dx = above.x - below.x;
        const double dy = above.y - below.y;
        const double theta_avg = 0.5 * (below.theta + above.theta);
        const double y_avg = 0.5 * (below.y + above.y);
        const double flux = 0.5 * (mass_flux(below.mach) + mass_flux(above.mach));

        // V dot n over the segment, where n is the downstream normal (dy, -dx)/|ds|; the
        // |ds| cancels against the area element, leaving the bare differences.
        const double normal_flux = std::cos(theta_avg) * dy - std::sin(theta_avg) * dx;
        const double area_weight = axisymmetric ? (2.0 * M_PI * y_avg) : 1.0;
        mdot += flux * normal_flux * area_weight;
    }

    // 1-D critical mass flow through the throat, in the same normalization. The
    // axisymmetric case is the full disc; the planar case is the half-height the solver
    // actually meshes.
    const double throat_area = axisymmetric ? (M_PI * r_throat * r_throat) : r_throat;
    const double mdot_reference = mass_flux(1.0) * throat_area;
    if (!(std::abs(mdot_reference) > 0.0)) return std::numeric_limits<double>::quiet_NaN();

    return (mdot - mdot_reference) / mdot_reference;
}

} // namespace

StartLine build_start_line(const MocSolveContext& ctx, const ThroatCondition& throat) {
    StartLine line;
    std::vector<CharacteristicPoint>& data_line = line.points;
    const size_t num_points = static_cast<size_t>(ctx.options.num_characteristics);
    data_line.reserve(num_points);
    if (!ctx.options.theta_schedule.empty()) {
        line.theta_schedule = ctx.options.theta_schedule;
    }

    // Resolve which start line this solve actually uses, once, as a single enum-valued
    // decision that both the fan/KL choice below and the theta_max derivation key off of
    // (previously a bool computed the choice while a separate ad hoc condition guarded the
    // theta_max derivation, and the two could disagree once a fallback was possible).
    // MocStartLine::AUTO keeps today's rule: a positive downstream curvature radius and
    // axisymmetric flow select the Kliegel-Levine series; anything else selects the fan.
    // The series is also not valid for planar flow when forced (validate_moc_options
    // rejects that combination outright).
    MocStartLine resolved_start_line = ctx.options.start_line;
    if (resolved_start_line == MocStartLine::AUTO) {
        resolved_start_line =
            (ctx.options.geometry.downstream_wall_curvature_radius <= 0.0 ||
             ctx.options.flow_type == MocFlowKind::PLANAR)
                ? MocStartLine::CENTERED_FAN
                : MocStartLine::KLIEGEL_LEVINE;
        if (resolved_start_line == MocStartLine::CENTERED_FAN &&
            ctx.options.flow_type == MocFlowKind::PLANAR &&
            ctx.options.geometry.downstream_wall_curvature_radius > 0.0 &&
            ctx.options.mode != MocMode::DESIGN_MIN_LENGTH)
        {
            ctx.log.info("Kliegel-Levine initialization is only valid for axisymmetric flow; "
                     "using centered-fan initialization for planar flow.");
        }
    }

    // In analysis mode the expansion is set by the wall contour, not by a design
    // theta_max (which is meaningless there and typically left unset). When the
    // centered-fan initializer will be used, derive theta_max from the profile. Shared
    // between the eager fan path and the AUTO fallback below (see the KL catch clause) so
    // the fallback gets exactly the same theta_max the eager path would have used.
    //
    // DESIGN_RAO is included alongside ANALYSIS: before this package DESIGN_RAO's default
    // positive curvature radius meant it never took the fan branch at all (use_centered_fan
    // was false whenever curvature was positive and axisymmetric), so this branch's
    // ANALYSIS-only guard was never exercised for it. The AUTO wall-angle-miss fallback
    // changes that -- DESIGN_RAO's default r_arc = 0.382 throat misses the threshold just
    // as an ANALYSIS contour at the same throat does -- and DESIGN_RAO's profile is equally
    // available here: solve() resolves the wall contour (which generates the Rao contour
    // for DESIGN_RAO) before build_start_line() runs, for every mode.
    // MocInitialization reads its wall contour from its own MocOptions copy (its constructor
    // predates MocSolveContext and takes MocOptions, not a resolved NozzleProfile), so the
    // copy's nozzle_profile has to be the *resolved* contour (ctx.wall) rather than
    // ctx.options.nozzle_profile as given: for DESIGN_RAO the latter is whatever the caller
    // passed (typically empty; the generated Rao contour lives only in ctx.wall now that
    // solve() no longer writes it back into options.nozzle_profile -- see resolve_wall_profile
    // in moc_nozzle.cpp). initialize_kliegel_levine's wall-consistent correction and its
    // y_wall fixed-point both key off this field.
    MocOptions init_options = ctx.options;
    init_options.nozzle_profile = ctx.wall;
    auto derive_fan_theta_max = [&]() {
        if (ctx.options.mode != MocMode::ANALYSIS && ctx.options.mode != MocMode::DESIGN_RAO) return;
        const NozzleProfile& wall = ctx.wall;
        if (wall.size() < 2) {
            throw std::runtime_error(
                "Centered-fan initialization requires a wall profile with at least 2 points.");
        }
        init_options.theta_max = wall.max_theta();
        if (init_options.theta_max <= 0.0) {
            throw std::runtime_error(
                "Wall profile must have a positive expansion angle for centered-fan initialization.");
        }
    };
    if (resolved_start_line == MocStartLine::CENTERED_FAN) {
        derive_fan_theta_max();
    }
    MocInitialization initializer(ctx.options.geometry, ctx.thermo, init_options);

    switch (ctx.options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            // The minimum length nozzle is a special case as the sonic line is straight, and
            // a centered expansion fan is the "exact" solution.
            // Therefore all initialization methods simplify to straight line initialization.

            // A centered expansion fan is collinear along a single C+ characteristic.
            line.family = CharacteristicFamily::PLUS;
            line.used = MocStartLine::CENTERED_FAN;

            auto expansion_line = initializer.initialize_centered_expansion(throat);
            // Mirror the fan angles (user-supplied or auto-generated inside the
            // initializer) into line.theta_schedule: the min-length wall solve reads
            // theta_wall = theta_max - line.theta_schedule[k] for each wall point.
            line.theta_schedule.resize(expansion_line.size());
            for (size_t i = 0; i < expansion_line.size(); i++) {
                line.theta_schedule[i] = expansion_line[i].theta;
            }
            //for a minimum length nozzle, the sonic line is straight.
            //the initial data line is generated through a single Prandtl-Meyer expansion at the throat.
            data_line = march_centered_fan(ctx, expansion_line);
            break;
        }
        case MocMode::DESIGN_RAO:
        case MocMode::ANALYSIS: {
            // The Kliegel-Levine transonic start line crosses many characteristics, so it is
            // seeded as a generic data line (line.family stays empty).
             // the first theta should be small to minimize approximation error.
            if (resolved_start_line == MocStartLine::KLIEGEL_LEVINE) {
                try {
                    // The Kliegel-Levine transonic line already spans axis-to-wall with full
                    // thermodynamic state and K+/K- set: it IS the initial data line. Unlike
                    // the centered fan (whose rays all emanate from the throat lip and must
                    // be marched into the flow field), re-marching these points pairwise
                    // would intersect characteristics *behind* their parents. See
                    // CharacteristicNet::add_initial_data_line's nullopt branch for how this
                    // non-collinear line is actually seeded into the net without that
                    // re-marching.
                    data_line = initializer.initialize_kliegel_levine(throat);
                    line.wall_bc_residual = initializer.last_wall_bc_residual;
                    line.used = MocStartLine::KLIEGEL_LEVINE;
                }
                catch (const KlWallAngleFallback& e) {
                    // MocOptions::start_line was AUTO (a forced KLIEGEL_LEVINE throws
                    // ConvergenceError instead, which is not this type and is left to
                    // propagate to solve()'s exception boundary as INITIALIZATION_FAILED).
                    // Rebuild as the centered fan below, with theta_max derived exactly as
                    // the eager fan path derives it.
                    ctx.log.info("Kliegel-Levine series misses the wall angle by {:.6f} rad at "
                             "R = {:.6f}; using centered-fan initialization.", e.miss, e.R);
                    resolved_start_line = MocStartLine::CENTERED_FAN;
                    derive_fan_theta_max();
                    initializer = MocInitialization(ctx.options.geometry, ctx.thermo, init_options);
                }
            }
            if (resolved_start_line == MocStartLine::CENTERED_FAN) {
                line.family = CharacteristicFamily::PLUS;
                line.used = MocStartLine::CENTERED_FAN;
                auto expansion_line = initializer.initialize_centered_expansion(throat);
                line.theta_schedule.resize(expansion_line.size());
                for (size_t i = 0; i < expansion_line.size(); i++) {
                    line.theta_schedule[i] = expansion_line[i].theta;
                }
                //for a minimum length nozzle, the sonic line is straight.
                //the initial data line is generated through a single Prandtl-Meyer expansion at the throat.
                data_line = march_centered_fan(ctx, expansion_line);
            }
            break;
        }
        case MocMode::DESIGN_CENTERLINE: {
            throw NotImplementedError("Centerline nozzle initialization not implemented.");
        }
    }

    for (size_t i = 0; i < data_line.size(); i++) {
        const auto& pt = data_line[i];
        ctx.log.debug("INIT[{}] x={:.6f} y={:.6f} theta={:.6f} nu={:.6f} mach={:.6f} mu={:.6f}",
            i, pt.x, pt.y, pt.theta, pt.nu, pt.mach, pt.mu);
        // The Kliegel-Levine and Sauer transonic-line initializers compute Mach
        // directly (not via Prandtl-Meyer inversion) and are not funneled through a
        // unit process that would otherwise validate them; check them here so a
        // breakdown of the small-perturbation series (e.g. a subsonic or non-finite
        // result near the wall) is reported as MocFailure::INITIALIZATION_FAILED
        // instead of an invalid point silently entering the net.
        MocErrorCode validity = check_point_validity(pt, ctx.options.solver_options.abstol);
        if (validity != MocErrorCode::NONE) {
            throw ConvergenceError(std::format(
                "Initial data line point {} at (x={}, y={}) failed validity check: {}.",
                i, pt.x, pt.y, to_string(validity)));
        }
    }
    return line;
}

MocInitDiagnostics measure_start_line(const StartLine& line, const MocSolveContext& ctx,
                                      double reference_spacing, double throat_radius)
{
    const std::vector<CharacteristicPoint>& data_line = line.points;
    MocInitDiagnostics diag;
    diag.points = data_line.size();
    if (data_line.size() < 2) return diag;

    // The generators return the line ordered axis to wall (index 0 on the axis for both
    // the Kliegel-Levine line and the marched centered fan).
    const CharacteristicPoint& axis_pt = data_line.front();
    const CharacteristicPoint& wall_pt = data_line.back();
    const double spacing = (reference_spacing > 0.0) ? reference_spacing : 1.0;

    diag.mach_axis = axis_pt.mach;
    diag.mach_wall = wall_pt.mach;
    diag.mach_ratio = (axis_pt.mach > 0.0) ? wall_pt.mach / axis_pt.mach : 0.0;
    diag.mu_axis = axis_pt.mu;
    diag.mu_wall = wall_pt.mu;
    {
        const double cot_axis = std::tan(axis_pt.mu) > 0.0 ? 1.0 / std::tan(axis_pt.mu) : 0.0;
        const double cot_wall = std::tan(wall_pt.mu) > 0.0 ? 1.0 / std::tan(wall_pt.mu) : 0.0;
        diag.cot_mu_ratio = (cot_wall > 0.0) ? cot_axis / cot_wall : 0.0;
    }

    diag.start_line_used = line.used;

    // K+ of the topmost *interior* point. The Kliegel-Levine line carries a separate wall
    // point (data_line.back()), so its topmost interior point is the one just below it; the
    // centered fan carries no such boundary point within data_line at all (solve() seeds
    // the throat lip anchor outside it), so its own last point already is the topmost
    // interior point.
    diag.kplus_wall_end =
        (diag.start_line_used == MocStartLine::KLIEGEL_LEVINE && data_line.size() >= 2)
            ? data_line[data_line.size() - 2].K_plus
            : wall_pt.K_plus;

    // Fit to the prescribed wall. Two things gate this. A design mode has no prescribed
    // contour to be off by. And a line that lies along one characteristic (the centered
    // fan) does not carry its own wall point at all -- solve() seeds the throat lip
    // separately as the wall anchor, and that lip is on the wall by construction -- so its
    // top point is an interior point and measuring a "gap" from it would be meaningless.
    const NozzleProfile& wall = ctx.wall;
    const bool line_owns_wall_point = !line.family.has_value();
    if (line_owns_wall_point && wall.size() >= 2 &&
        wall_pt.x >= wall.x_min() && wall_pt.x <= wall.x_max()) {
        diag.wall_gap = std::abs(wall_pt.y - wall.radius_at(wall_pt.x));
        diag.wall_gap_over_spacing = diag.wall_gap / spacing;
        diag.wall_theta_mismatch = std::abs(wall_pt.theta - wall.theta_at(wall_pt.x));

        // How far the raw (pre-correction) series missed this wall angle. The correction
        // has already overwritten wall_pt.theta (which is why wall_theta_mismatch above is
        // ~0), so the value comes from the initializer, which measured it before correcting.
        if (diag.start_line_used == MocStartLine::KLIEGEL_LEVINE) {
            diag.wall_bc_residual = line.wall_bc_residual;
        }

        // Distance to the sharpest slope break on the contour. A conical profile's arc runs
        // into a straight cone, and seeding the wall march exactly there means the first
        // wall angle query is a one-sided difference across the break.
        double worst_jump = 0.0;
        double break_x = std::numeric_limits<double>::quiet_NaN();
        for (size_t i = 2; i + 1 < wall.size(); i++) {
            const double jump = std::abs(wall.slope_at_idx(i) - wall.slope_at_idx(i - 1));
            if (jump > worst_jump) {
                worst_jump = jump;
                break_x = wall.x[i - 1];
            }
        }
        // Below this the contour is smooth to within its own discretization and there is no
        // break to be near; report a sentinel rather than a meaningless distance.
        constexpr double slope_break_threshold = 1e-6;
        diag.wall_station_to_tangency = (worst_jump > slope_break_threshold)
            ? (wall_pt.x - break_x) / spacing
            : std::numeric_limits<double>::max();
    }
    else {
        diag.wall_station_to_tangency = std::numeric_limits<double>::max();
    }

    // The shift's natural scale. The transonic region's axial extent goes as
    // sqrt(r_throat * R_curvature), so the same absolute shift is a different fraction of
    // it at every throat curvature.
    if (ctx.options.geometry.downstream_wall_curvature_radius > 0.0 &&
        ctx.options.geometry.throat_radius > 0.0)
    {
        const double R = ctx.options.geometry.downstream_wall_curvature_radius
                       / ctx.options.geometry.throat_radius;
        diag.shift_over_transonic_length = ctx.options.initial_line_axial_shift / std::sqrt(R);
    }

    // Where each point's C- reaches the axis, on a straight-ray estimate, and how unevenly
    // graded those arrivals are. A centered fan's points share one location, so its rays
    // still land at distinct stations and the measure stays meaningful.
    std::vector<double> arrivals;
    arrivals.reserve(data_line.size());
    for (const CharacteristicPoint& pt : data_line) {
        const double slope = std::tan(pt.theta - pt.mu);
        arrivals.push_back((std::abs(slope) > 1e-12) ? pt.x + pt.y / std::abs(slope) : pt.x);
    }
    double min_gap = std::numeric_limits<double>::max();
    double max_gap = 0.0;
    for (size_t i = 1; i < arrivals.size(); i++) {
        const double gap = std::abs(arrivals[i] - arrivals[i - 1]);
        min_gap = std::min(min_gap, gap);
        max_gap = std::max(max_gap, gap);
    }
    diag.axis_arrival_grading = (min_gap > 0.0 && min_gap != std::numeric_limits<double>::max())
        ? max_gap / min_gap
        : std::numeric_limits<double>::max();

    // Spacelike margin of the line's own segments: how far a segment sits from being
    // parallel to a characteristic through one of its endpoints, as a fraction of its
    // height. A centered fan is degenerate here (all points share a location) and is
    // skipped.
    double min_margin = std::numeric_limits<double>::max();
    for (size_t i = 1; i < data_line.size(); i++) {
        const CharacteristicPoint& below = data_line[i - 1];
        const CharacteristicPoint& above = data_line[i];
        const double dx = above.x - below.x;
        const double dy = above.y - below.y;
        if (dy <= 0.0) continue;
        const double margin_minus = (dy - std::tan(above.theta - above.mu) * dx) / dy;
        const double margin_plus = (dy - std::tan(below.theta + below.mu) * dx) / dy;
        min_margin = std::min({min_margin, margin_minus, margin_plus});
    }
    diag.min_spacelike_margin =
        (min_margin == std::numeric_limits<double>::max()) ? 0.0 : min_margin;

    diag.mass_flow_error = start_line_mass_flow_error(data_line, ctx, throat_radius);
    return diag;
}

} // namespace Goddard