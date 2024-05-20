#include "cantera/core.h"
#include "cantera/base/Array.h"
#include "eigen3/Eigen/Dense"

#include "goddard/thermo.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::ArrayXd;
using Eigen::ArrayXXd;

namespace Goddard {


ArrayXXd get_stoichiometric_coeffs(Cantera::Solution& gas){
    auto thermo = gas.thermo();
    size_t n_elements = thermo->nElements();
    size_t n_species = thermo->nSpecies();

    //construct matrix of elemental stoichiometric coefficients
    ArrayXXd stoich_coeffs(n_elements, n_species);
    for(size_t i = 0; i < n_elements; i++){
        for(size_t j = 0; j < n_species; j++){
            stoich_coeffs(i,j) = thermo->nAtoms(j, i);
        }
    }

    return stoich_coeffs;
}

ArrayXd get_mole_vector(Cantera::Solution& gas){
    auto thermo = gas.thermo();
    
    size_t n_species = thermo->nSpecies();
    double total_moles = 1.0/thermo->meanMolecularWeight();

    ArrayXd moles(n_species);
    thermo->getMoleFractions(moles.data());
    moles *= total_moles;
    return moles;
}

ArrayXd get_enthalpyRT_vector(Cantera::Solution& gas){
    auto thermo = gas.thermo();

    size_t n_species = thermo->nSpecies();

    ArrayXd std_enthalpies_RT(n_species);
    
    thermo->getEnthalpy_RT(std_enthalpies_RT.data());

    return std_enthalpies_RT;
}

ArrayXd get_cpR_vector(Cantera::Solution& gas){
    auto thermo = gas.thermo();

    size_t n_species = thermo->nSpecies();
    ArrayXd cp_R(n_species);
    
    thermo->getCp_R(cp_R.data());
    return cp_R;
}

ThermoDerivatives get_thermo_derivatives(Cantera::Solution& gas) {
    auto thermo = gas.thermo();

    size_t n_elements = thermo->nElements();
    size_t num_var = 2 * n_elements + 2;

    
    MatrixXd coeff_matrix(num_var, num_var);
    VectorXd rhs(num_var);


    auto std_enthalpies_RT = get_enthalpyRT_vector(gas);
    auto moles = get_mole_vector(gas);

    //indices
    size_t idx_dpi_dlogT_P = 0;
    size_t idx_dlogn_dlogT_P = idx_dpi_dlogT_P + n_elements;
    size_t idx_dpi_dlogP_T = idx_dlogn_dlogT_P + 1;
    size_t idx_dlogn_dlogP_T = idx_dpi_dlogP_T + n_elements;


    //construct matrix of elemental stoichiometric coefficients
    auto stoich_coeffs = get_stoichiometric_coeffs(gas);

    //equations for derivatives w.r.t. temperature
    //first n_elements equations
    for (size_t i = 0; i < n_elements; i++) {
        for (size_t j = 0; j < n_elements; j++){
            coeff_matrix(i,j) = (stoich_coeffs(i)*stoich_coeffs(j)*moles).sum();
        }
        coeff_matrix(i, n_elements) = (stoich_coeffs(i) * moles).sum();
        rhs(i) = -(stoich_coeffs(i)*moles*std_enthalpies_RT).sum();
    }
    
    for (size_t i = 0; i < n_elements; i++){
        coeff_matrix(n_elements, i) = (stoich_coeffs(i) * moles).sum();
    }
    rhs(n_elements) = -(moles * std_enthalpies_RT).sum();

    // equations for derivatives w.r.t. pressure
    for (size_t i = 0; i < n_elements; i++) {
        for (size_t j = 0; j < n_elements; j++){
            coeff_matrix(n_elements+1+i, n_elements+1+j) = (stoich_coeffs(i)*stoich_coeffs(j) * moles).sum();
        }
        coeff_matrix(n_elements+1+i, 2*n_elements+1) = (stoich_coeffs(i) * moles).sum();
        rhs(n_elements+1+i) = (stoich_coeffs(i)*moles).sum();
    }

    for (size_t i = 0; i < n_elements; i++){
        coeff_matrix(2*n_elements+1, n_elements+1+i) = (stoich_coeffs(i) * moles).sum();
    }
    rhs(2*n_elements+1) = moles.sum();

    VectorXd derivs = coeff_matrix.colPivHouseholderQr().solve(rhs);

    ArrayXd dpi_dlogT_P = derivs(Eigen::seq(idx_dpi_dlogT_P, idx_dpi_dlogT_P + n_elements-1));
    double dlogn_dlogT_P = derivs(idx_dlogn_dlogT_P);
    ArrayXd dpi_dlogP_T = derivs(Eigen::seq(idx_dpi_dlogP_T, idx_dpi_dlogP_T + n_elements-1));
    double dlogn_dlogP_T = derivs(idx_dlogn_dlogP_T);
    
    return {dpi_dlogT_P, dlogn_dlogT_P, dpi_dlogP_T, dlogn_dlogP_T};
}

ThermoProperties get_thermo_properties(Cantera::Solution& gas, const ThermoDerivatives& derivs){
    auto moles = get_mole_vector(gas);
    
    auto stoich_coeffs = get_stoichiometric_coeffs(gas);
    auto thermo = gas.thermo();
    
    size_t n_elements = thermo->nElements();
    size_t n_species = thermo->nSpecies();

    auto std_enthalpy_RT = get_enthalpyRT_vector(gas);
    auto cp_R = get_cpR_vector(gas);

    //isobaric specific heat
    double spec_heat_p = 0.0;

    for (size_t k=0; k < n_species; k++){
        for (size_t i=0; i < n_elements; i++){
            spec_heat_p += derivs.dpi_dlogP_T(k)*stoich_coeffs(i, k) * moles(k, 0) * std_enthalpy_RT(k,0);
        }
        spec_heat_p += moles(k,0) * std_enthalpy_RT(k,0) * derivs.dlogn_dlogT_P;
        spec_heat_p += moles(k,0) * cp_R(k,0);
        spec_heat_p += moles(k,0) * std_enthalpy_RT(k,0) * std_enthalpy_RT(k,0);
    }

    double dlogV_dlogT_P = 1 + derivs.dlogn_dlogT_P;
    double dlogV_dlogP_T = -1 + derivs.dlogn_dlogP_T;

    //isochoric specific heat
    double additional_Cpv_term = thermo->pressure() * thermo->molarVolume() / thermo->temperature() * dlogV_dlogT_P * dlogV_dlogT_P / dlogV_dlogP_T;
    double spec_heat_v = spec_heat_p + additional_Cpv_term;
    
    //heat capacity ratios
    double gamma = spec_heat_p / spec_heat_v;
    double gamma_s = -gamma/dlogV_dlogP_T;
    
    return {dlogV_dlogT_P, dlogV_dlogP_T, spec_heat_p, gamma_s};
}

} //namespace Goddard