#pragma once

namespace Goddard {

// Prandtl-Meyer angle function
double prandtl_meyer(double mach, double gamma); 

//derivative of Prandtl-Meyer function w.r.t. mach number
double prandtl_meyer_derivative(double mach, double gamma); 

// inverse function: obtain Mach number from Prandtl-meter angle using Newton's method
double mach_from_prandtl_meyer(double nu, double gamma, 
                                double mach_guess = 0.0,
                                double tol = 1e-10,
                                int max_iter = 20);

} // namespace Goddard