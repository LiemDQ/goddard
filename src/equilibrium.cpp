#include "cantera/core.h"
#include "eigen3/Eigen/Dense"

#include "goddard/equilibrium.hpp"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;

namespace Goddard {


ArrayXXd get_stoichiometric_coeffs(const Cantera::ThermoPhase& gas){
    size_t n_elements = gas.nElements();
    size_t n_species = gas.nSpecies();

    //construct matrix of elemental stoichiometric coefficients
    ArrayXXd stoich_coeffs(n_species, n_elements);
    for(size_t i = 0; i < n_elements; i++){
        for(size_t j = 0; j < n_species; j++){
            stoich_coeffs(j,i) = gas.nAtoms(j, i);
        }
    }

    return stoich_coeffs;
}

ArrayXd get_mole_vector(const Cantera::ThermoPhase& gas){
    size_t n_species = gas.nSpecies();
    double total_moles = 1.0/gas.meanMolecularWeight();

    ArrayXd moles(n_species);
    gas.getMoleFractions(moles.data());
    moles *= total_moles;
    return moles;
}

ArrayXd get_enthalpyRT_vector(const Cantera::ThermoPhase& gas){
    size_t n_species = gas.nSpecies();

    ArrayXd std_enthalpies_RT(n_species);
    
    gas.getEnthalpy_RT(std_enthalpies_RT.data());

    return std_enthalpies_RT;
}

ArrayXd get_cpR_vector(const Cantera::ThermoPhase& gas){
    size_t n_species = gas.nSpecies();
    ArrayXd cp_R(n_species);
    
    gas.getCp_R(cp_R.data());
    return cp_R;
}

/**
 * @brief Calculate thermodynamic derivatives for a reactive gas at chemical equilibrium.
 * 
 * Cantera calculates thermodynamic derivatives assuming a fixed composition, 
 * which is not correct for reactive flows, such as ones in chemical equilibrium, as the 
 * change in molar quantity and chemical heat release/absorption must also be taken into account. 
 * See [Gordon & McBride, 1994, "Computer program for calculation of complex chemical equilibrium 
 * compositions and applications. Part 1: Analysis"](https://ntrs.nasa.gov/citations/19950013764) for more details.
 * 
 * 
*/
EquilibriumDerivatives get_thermo_equilibrium_derivatives(const Cantera::ThermoPhase& gas) {
    size_t n_elements = gas.nElements();
    size_t num_var = 2 * n_elements + 2;

    
    MatrixXd coeff_matrix(num_var, num_var);
    coeff_matrix.setZero();

    VectorXd rhs(num_var);
    rhs.setZero();

    auto std_H_RT = get_enthalpyRT_vector(gas);
    auto moles = get_mole_vector(gas);

    //indices
    size_t idx_dpi_dlogT_P = 0;
    size_t idx_dlogn_dlogT_P = idx_dpi_dlogT_P + n_elements;
    size_t idx_dpi_dlogP_T = idx_dlogn_dlogT_P + 1;
    size_t idx_dlogn_dlogP_T = idx_dpi_dlogP_T + n_elements;


    //construct matrix of elemental stoichiometric coefficients
    auto stoich_coeffs = get_stoichiometric_coeffs(gas);

    //equations for derivatives w.r.t. temperature
    //Gordon & McBride equation 2.56
    for (size_t i = 0; i < n_elements; i++) {
        for (size_t j = 0; j < n_elements; j++){
            coeff_matrix(i,j) = (stoich_coeffs.col(i)*stoich_coeffs.col(j)*moles).sum();
        }
        coeff_matrix(i, n_elements) = (stoich_coeffs.col(i) * moles).sum();
        rhs(i) = -(stoich_coeffs.col(i)*moles*std_H_RT).sum();
    }

    //TODO: temperature derivative equation 2.57 for condensed species
    
    //single term for (dpi/dlogT)_P
    for (size_t i = 0; i < n_elements; i++){
        coeff_matrix(n_elements, i) = (stoich_coeffs.col(i) * moles).sum();
    }
    //Gordon & McBride equation 2.58
    rhs(n_elements) = -(moles * std_H_RT).sum();

    // equations for derivatives w.r.t. pressure
    //Gordon & McBride equation 2.64
    for (size_t i = 0; i < n_elements; i++) {
        for (size_t j = 0; j < n_elements; j++){
            coeff_matrix(n_elements+1+i, n_elements+1+j) = (stoich_coeffs.col(i)*stoich_coeffs.col(j) * moles).sum();
        }
        coeff_matrix(n_elements+1+i, 2*n_elements+1) = (stoich_coeffs.col(i) * moles).sum();
        rhs(n_elements+1+i) = (stoich_coeffs.col(i)*moles).sum();
    }

    //TODO: pressure derivative equation 2.65 for condensed species

    for (size_t i = 0; i < n_elements; i++){
        coeff_matrix(2*n_elements+1, n_elements+1+i) = (stoich_coeffs.col(i) * moles).sum();
    }
    //Gordon & McBride equation 2.66
    rhs(2*n_elements+1) = moles.sum();

    // the derivatives are obtained from solving the system of equations. 
    VectorXd derivs = coeff_matrix.colPivHouseholderQr().solve(rhs);

    ArrayXd dpi_dlogT_P = derivs(Eigen::seq(idx_dpi_dlogT_P, idx_dpi_dlogT_P + n_elements-1));
    double dlogn_dlogT_P = derivs(idx_dlogn_dlogT_P);
    ArrayXd dpi_dlogP_T = derivs(Eigen::seq(idx_dpi_dlogP_T, idx_dpi_dlogP_T + n_elements-1));
    double dlogn_dlogP_T = derivs(idx_dlogn_dlogP_T);
    
    return {dpi_dlogT_P, dlogn_dlogT_P, dpi_dlogP_T, dlogn_dlogP_T};
}

/**
 * @brief Calculate thermodynamic properties of a reacting gas at equilibrium: volume derivatives, heat capacity and adiabatic index.
*/
EquilibriumProperties get_thermo_equilibrium_properties(const Cantera::ThermoPhase& gas, const EquilibriumDerivatives& derivs){
    auto moles = get_mole_vector(gas);
    
    auto stoich_coeffs = get_stoichiometric_coeffs(gas);
    size_t n_elements = gas.nElements();

    auto H_RT = get_enthalpyRT_vector(gas);
    auto cp_R = get_cpR_vector(gas);

    //isobaric specific heat
    //NOTE: this is on a kg basis because of the normalization of the mole vector to the gas molar mass
    double dpi_cp_contribution = 0;
    for (size_t i=0; i < n_elements; i++){
        dpi_cp_contribution += derivs.dpi_dlogT_P(i) * (moles * H_RT * stoich_coeffs.col(i)).sum();
    }
    double spec_heat_p = Cantera::GasConstant * (
        dpi_cp_contribution 
        + (moles * H_RT).sum()*derivs.dlogn_dlogT_P
        + (moles * cp_R).sum()
        + (moles * H_RT * H_RT).sum()
    );

    //volumetric derivatives
    double dlogV_dlogT_P = 1 + derivs.dlogn_dlogT_P;
    double dlogV_dlogP_T = -1 + derivs.dlogn_dlogP_T;

    //isochoric specific heat
    double mass_volume =  gas.molarVolume() / gas.meanMolecularWeight(); 
    double additional_Cpv_term = gas.pressure() * mass_volume / gas.temperature() * dlogV_dlogT_P * dlogV_dlogT_P / dlogV_dlogP_T;
    double spec_heat_v = spec_heat_p + additional_Cpv_term;
    
    //heat capacity ratios
    double gamma = spec_heat_p / spec_heat_v;
    double gamma_s = -gamma/dlogV_dlogP_T;
    
    return {dlogV_dlogT_P, dlogV_dlogP_T, spec_heat_p, gamma_s};
}

EquilibriumProperties get_thermo_equilibrium_properties(const Cantera::ThermoPhase& gas) {
    return get_thermo_equilibrium_properties(gas, get_thermo_equilibrium_derivatives(gas));
}

double get_equilibrium_gamma(const Cantera::ThermoPhase& gas) {
    return get_thermo_equilibrium_properties(gas).gamma_s;
}

} //namespace Goddard