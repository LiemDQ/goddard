#pragma once
#include <vector>

#include "cantera/core.h"
#include "eigen3/Eigen/Dense"

namespace Goddard {


struct EquilibriumDerivatives {
    Eigen::ArrayXd dpi_dlogT_P;
    double dlogn_dlogT_P;
    Eigen::ArrayXd dpi_dlogP_T;
    double dlogn_dlogP_T;
};

struct EquilibriumProperties {
    double dlogV_dlogT_P;
    double dlogV_dlogP_T;
    double spec_heat_p;
    double gamma_s;
};

Eigen::ArrayXXd get_stoichiometric_coeffs(Cantera::Solution& gas);
Eigen::ArrayXd get_mole_vector(Cantera::Solution& gas);
Eigen::ArrayXd get_enthalpyRT_vector(Cantera::Solution& gas);
Eigen::ArrayXd get_cpR_vector(Cantera::Solution& gas);

EquilibriumDerivatives get_thermo_equilibrium_derivatives(Cantera::Solution& gas);
EquilibriumProperties get_thermo_equilibrium_properties(Cantera::Solution& gas, const EquilibriumDerivatives& derivatives);

} //namespace Goddard