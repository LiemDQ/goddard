#include <cmath>
#include <algorithm>
#include "goddard/moc.hpp"

namespace Goddard {


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


} // namespace Goddard