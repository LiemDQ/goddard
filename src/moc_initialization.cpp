#include <algorithm>
#include "goddard/gas_dynamics.hpp"
#include "goddard/moc_initialization.hpp"

namespace Goddard {


double gquad(double gamma, double a, double b, double c, double d) {
    return (a*gamma*gamma + b*gamma + c)/d;
}

MocInitialization::MocInitialization(
    NozzleGeometry geom, 
    ThermodynamicContext& thermo, 
    const MocOptions& options) 
    : geometry(geom), m_thermo(thermo), m_options(options)
{

}

void MocInitialization::set_state_from_critical_velocity_ratio(
    CharacteristicPoint& pt, double m_star, const ThroatCondition& throat)
{
    GasChemistry chemistry = m_thermo.gas.has_value()
        ? m_thermo.gas->chemistry
        : GasChemistry::PERFECT_GAS;

    if (chemistry == GasChemistry::PERFECT_GAS) {
        // The perfect-gas path keeps velocity dimensionless and equal to the Mach number,
        // so the conversion has to happen here rather than being folded into a velocity.
        pt.update_thermodynamic_state_from_mach(
            m_thermo, mach_from_critical_velocity_ratio(m_star, throat.gamma_s));
        return;
    }

    // For a real gas the M*<->M algebra is not available in closed form, but its defining
    // relation is: a* is the throat speed of sound, so V follows directly and the
    // PrandtlMeyerTable can be looked up by velocity.
    pt.update_thermodynamic_state_from_V(m_thermo, m_star * throat.speed_of_sound);
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
    // estimate. This -- not y -- is the coordinate the start line must be uniform in.
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
    // refinement never cured the axisymmetric breakdown.
    auto axis_arrival = [&](double y) {
        const double x = station_x(y);
        const double u_star = KL_xMach(y, KL_z_coordinate(x, gamma), gamma, R);
        const double v_star = KL_yMach(x, y, gamma, R);
        const double mach = mach_from_critical_velocity_ratio(std::hypot(u_star, v_star), gamma);
        if (!(mach > 1.0)) return x;   // still sonic: the C- is vertical, it arrives at x
        const double theta = std::atan2(v_star, u_star);
        const double slope = std::tan(theta - mach_to_mu(mach));
        if (std::abs(slope) < 1e-12) return x;
        return x + y / std::abs(slope);
    };

    // Tabulate the arrival map once, then invert it by interpolation: far cheaper than a
    // root solve per point, and it makes the monotonicity check free.
    constexpr size_t table_size = 201;
    std::vector<double> y_table(table_size), arrival_table(table_size);
    bool monotone = true;
    for (size_t k = 0; k < table_size; k++) {
        y_table[k] = static_cast<double>(k) / static_cast<double>(table_size - 1);
        arrival_table[k] = axis_arrival(y_table[k]);
        if (k > 0 && !(arrival_table[k] > arrival_table[k - 1])) monotone = false;
    }

    std::vector<CharacteristicPoint> points(num_points);
    const double last = static_cast<double>(num_points - 1);

    for (size_t i = 0; i < num_points; i++) {
        double y;
        // Interior branches produce a normalized station in [0,1]; the two anchors are
        // written in normalized terms too so the single scaling below covers every case.
        if (i == 0) {
            y = 0.0;                       // anchor the axis point exactly
        }
        else if (i == num_points - 1) {
            y = 1.0;                       // anchor the wall point on the contour
        }
        else if (!monotone) {
            // The arrival map should be monotone for any physical throat; if the series is
            // being evaluated somewhere it is not, fall back to a uniform line rather than
            // producing a scrambled start line.
            y = static_cast<double>(i) / last;
        }
        else {
            const double target = arrival_table.front()
                + (arrival_table.back() - arrival_table.front()) * static_cast<double>(i) / last;
            const size_t k = static_cast<size_t>(
                std::lower_bound(arrival_table.begin(), arrival_table.end(), target)
                - arrival_table.begin());
            const size_t hi = std::clamp<size_t>(k, 1, table_size - 1);
            const double a0 = arrival_table[hi - 1], a1 = arrival_table[hi];
            const double w = (a1 > a0) ? (target - a0) / (a1 - a0) : 0.0;
            const double y_arrival = y_table[hi - 1] + w * (y_table[hi] - y_table[hi - 1]);
            const double y_uniform = static_cast<double>(i) / last;
            const double c = m_options.initial_line_clustering;
            // c < 0 pushes points toward the axis (y_arrival > y_uniform everywhere), c > 0
            // toward the wall. Clamped so the line stays inside the throat and ordered.
            y = std::clamp((1.0 - c) * y_uniform + c * y_arrival, 1e-6, 1.0 - 1e-6);
        }
        // The interior stations are laid out on the unit interval; stretch them onto the
        // line's actual span so they stay evenly distributed when the wall end moves.
        points[i] = evaluate(y * y_wall);
    }

    // The wall point keeps the flow angle the series gives it, and is deliberately *not*
    // overwritten with the contour's own angle. Forcing tangency there was tried and
    // measured: it leaves AR >= 4 unchanged and regresses AR = 2 from full coverage to a
    // PM inversion failure at N = 31. The reason is that a large shift on a sharply curved
    // throat puts the line's wall end downstream of the throat arc -- with the default
    // r_arc = 0.382 and a 15 degree cone the arc ends at x = 0.0989 and the default shift
    // is 0.1 -- so the contour there is the straight cone at 15 degrees while the transonic
    // solution is still turning through about 5.5. Overwriting the angle imposes that 9.5
    // degree jump on one end of the start line; keeping the series' own value leaves the
    // line self-consistent and lets the wall solve turn the flow over its first few steps.
    //
    // MocInitDiagnostics::wall_theta_mismatch reports the disagreement, which is the signal
    // that the start line has been placed outside the region the series describes.

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
        
        pt.update_thermodynamic_state_from_nu(m_thermo, pt.theta, 1.0);
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

} // namespace Goddard