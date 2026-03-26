#include <cmath>

#include "goddard/prandtlmeyer.hpp"

namespace Goddard {

double prandtl_meyer(double mach, double gamma) {
    double gp1_gm1 = (gamma + 1.0) / (gamma - 1.0);
    double m2_minus_1 = mach * mach - 1.0;
    
    return sqrt(gp1_gm1) * atan(sqrt(m2_minus_1 / gp1_gm1)) - atan(sqrt(m2_minus_1));
}

double prandtl_meyer_derivative(double mach, double gamma) {
    double m2 = mach * mach;
    return sqrt(m2 - 1.0) / ( mach * (1.0 + (gamma - 1.0) / 2.0 * m2));
}


double mach_from_prandtl_meyer(double nu, double gamma,
                                double mach_guess,
                                double tol,
                                int max_iter) {

    double mach = (mach_guess > 1.0) ? mach_guess : 1.0 + nu;

    for (int i = 0; i < max_iter; i++) {
        double residual = prandtl_meyer(mach, gamma) - nu;
        if (std::abs(residual) < tol) return mach;

        mach -= residual / prandtl_meyer_derivative(mach, gamma);
    }
    
    // should not reach here for well-formed inputs
    return -1.0;
}

} // namespace Goddard