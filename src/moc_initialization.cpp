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
        double machx = 1+alpha*pt.x + (gamma+1)*alpha*alpha*pt.y*pt.y/(2*(1+delt));
        pt.update_thermodynamic_state_from_mach(m_thermo, machx);
        pt.update_Ks();
    }

    //TODO: optional corrector step with updated gamma

    return points;
}

std::vector<CharacteristicPoint> MocInitialization::initialize_kliegel_levine(const ThroatCondition& throat) {
    double gamma = throat.gamma_s;
    size_t num_points = static_cast<size_t>(m_options.num_characteristics);
    double dy = 1.0/(num_points-1);
    double R = KL_R();
    double alpha = sauer_alpha(gamma);

    // Lift the sonic (zero-radial-velocity) locus off itself by a rigid downstream shift so
    // it can be seeded with both characteristic families (see MocOptions::initial_line_axial_shift
    // and CharacteristicNet::add_initial_data_line's nullopt branch).
    const double x_shift = m_options.initial_line_axial_shift;

    std::vector<CharacteristicPoint> points(num_points);
    for (size_t i = 0; i < num_points; i++) {
        CharacteristicPoint& pt = points[i];
        pt.y = dy * i;
        double x_guess = (gamma + 1)*alpha/(2*(3+delta()))*(1.0 - pt.y*pt.y); //Sauer transonic x-coordinate

        double x_sonic = KL_solve_transonic_x(pt.y, gamma, R, x_guess);
        pt.x = x_sonic + x_shift;
        double machx = KL_xMach(pt.y, KL_z_coordinate(pt.x, gamma), gamma, R);
        // Off the sonic locus, the radial velocity component (KL_yMach) is
        // generally nonzero, so theta can no longer be pinned to 0. Both machx and
        // the KL_yMach series are normalized velocity-like quantities (by a*), so
        // their ratio gives the flow angle to the same order as the series itself.
        double v = KL_yMach(pt.x, pt.y, gamma, R);
        pt.theta = std::atan2(v, machx);
        pt.update_thermodynamic_state_from_mach(m_thermo, machx);
        pt.update_Ks();
    }

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