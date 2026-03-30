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

/**
 * @brief Get matrix of stoichiometric coefficients of the species contained in the `ThermoPhase` object.
 * 
 * @return 2D Eigen array of stoichiometric coefficients. Rows represent species, while columns represent elements.
 * Ordering is the same as the data file used to generated the `ThermoPhase` object. 
 */
Eigen::ArrayXXd get_stoichiometric_coeffs(const Cantera::ThermoPhase& gas);
Eigen::ArrayXd get_mole_vector(const Cantera::ThermoPhase& gas);
Eigen::ArrayXd get_enthalpyRT_vector(const Cantera::ThermoPhase& gas);
Eigen::ArrayXd get_cpR_vector(const Cantera::ThermoPhase& gas);

EquilibriumDerivatives get_thermo_equilibrium_derivatives(const Cantera::ThermoPhase& gas);
EquilibriumProperties get_thermo_equilibrium_properties(const Cantera::ThermoPhase& gas, const EquilibriumDerivatives& derivatives);
EquilibriumProperties get_thermo_equilibrium_properties(const Cantera::ThermoPhase& gas);
double get_equilibrium_gamma(const Cantera::ThermoPhase& gas);

} //namespace Goddard