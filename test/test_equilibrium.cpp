#include "goddard/thermo.hpp"

#include "eigen3/Eigen/Dense"
#include "cantera/core.h"
#include "gtest/gtest.h"
#include <memory>
#include <vector>
#include <cmath>
#include <algorithm>

class DerivativeTests: public ::testing::Test {
    protected:
    DerivativeTests() {
        this->sln = Cantera::newSolution("h2o2.yaml", "ohmech");
        double temp = 2400.0; //K
        double pressure = 50.0*Cantera::OneAtm;
        auto gas = sln->thermo();
        
        n_species = sln->thermo()->nSpecies();
        n_elements = sln->thermo()->nElements();

        gas->setState_TPX(temp, pressure, "H2O:1, N2:1, O2:1. AR:0.1"); //completely random composition lol
    }
    std::shared_ptr<Cantera::Solution> sln;
    size_t n_species;
    size_t n_elements;
    

};

double max_fp_error(double val, double reltol = 1e-4, double abstol = 1e-10) {
    return std::max(abs(val*reltol), abstol);
}


TEST_F(DerivativeTests, stoichiometricCoeffsAreCorrect){
    auto coeffs = Goddard::get_stoichiometric_coeffs(*sln);
    Eigen::ArrayXXd expected_coeffs{
        {0, 0, 1, 2, 1, 1, 2, 2, 0, 0},
        {2, 1, 0, 0, 1, 2, 1, 2, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 2},
    };

    long n_elements_signed = static_cast<long>(n_elements); //suppress compiler warnings
    long n_species_signed = static_cast<long>(n_species);

    EXPECT_EQ(coeffs.size(), expected_coeffs.size());

    for (int i = 0; i < n_elements_signed; i++){
        for (int j = 0; j < n_species_signed; j++){
            // EXPECT_EQ(coeffs(i,j), expected_coeffs(i,j));
        }
    }
}

TEST_F(DerivativeTests, moleVectorIsCorrect){
    auto moles = Goddard::get_mole_vector(*sln);
    
    std::vector<double> expected_moles{0., 0., 0.,  0.01219185, 0., 0.01219185, 0. , 0. , 0.00121919, 0.01219185};
    ASSERT_EQ(moles.size(), n_species);

    for (size_t i = 0; i < n_species; i++){
        double tol = max_fp_error(expected_moles[i]);
        EXPECT_NEAR(moles(i), expected_moles[i], tol);
    }
}

TEST_F(DerivativeTests, moleVectorIsCorrectAfterEquilibration){
    sln->thermo()->equilibrate("HP");
    auto moles = Goddard::get_mole_vector(*sln);
    Eigen::ArrayXd expected_moles(n_species);
    expected_moles << 9.11244396e-06, 1.13250699e-06, 2.20601379e-05, 1.21177599e-02,
        2.63718243e-04, 1.20488303e-02, 2.45871008e-06, 2.53459510e-07,
        1.21918510e-03, 1.21918510e-02;

    ASSERT_EQ(moles.size(), n_species);

    for (size_t i = 0; i < n_species; i++){
        double tol = max_fp_error(expected_moles(i));
        EXPECT_NEAR(moles(i), expected_moles(i), tol);
    }
}

TEST_F(DerivativeTests, stdEnthalpiesRTAreCorrect){
    auto enthalpies = Goddard::get_enthalpyRT_vector(*sln);
    Eigen::ArrayXd exp_enthalpies(n_species);
    exp_enthalpies << 3.35335209, 13.11402496, 14.69431338,  3.73354955,  5.37563442,
       -7.39428256,  5.8582109 ,  0.05367515,  2.18942708,  3.54038084;
    
    ASSERT_EQ(enthalpies.size(), n_species);
    for (size_t i = 0; i < n_species; i++){
        double tol = max_fp_error(exp_enthalpies[i]);
        EXPECT_NEAR(enthalpies(i), exp_enthalpies(i), tol);
    }
}

TEST_F(DerivativeTests, thermoDerivativesAreCorrect){
    EXPECT_NEAR(sln->thermo()->temperature(), 2400.0, 1e-1);
    EXPECT_NEAR(sln->thermo()->enthalpy_mass(), 23985.583414274723, 1e-1);
    
    Goddard::ThermoDerivatives derivs = Goddard::get_thermo_equilibrium_derivatives(*sln);
    ASSERT_EQ(derivs.dpi_dlogP_T.size(),n_elements);
    ASSERT_EQ(derivs.dpi_dlogT_P.size(), n_elements);

    Eigen::ArrayXd exp_dpi_dlogT_P(n_species);
    exp_dpi_dlogT_P << -1.86677477,  4.63052867, -2.18942708, -1.77019042;
    
    double exp_dlogn_dlogT_P = -1.8359417928589634e-16;

    Eigen::ArrayXd exp_dpi_dlogP_T(n_species);
    exp_dpi_dlogP_T << 0.5, 0.25, 1., 0.5;

    double exp_dlogn_dlogP_T = 9.179708964294817e-17;

    for (size_t i = 0; i < n_elements; i++){
        double tol = max_fp_error(exp_dpi_dlogT_P(i));
        EXPECT_NEAR(derivs.dpi_dlogT_P(i), exp_dpi_dlogT_P(i), tol);
    }

    EXPECT_NEAR(derivs.dlogn_dlogT_P, exp_dlogn_dlogT_P, max_fp_error(exp_dlogn_dlogT_P));

    for (size_t i = 0; i < n_elements; i++){
        double tol = max_fp_error(exp_dpi_dlogP_T(i));
        EXPECT_NEAR(derivs.dpi_dlogP_T(i), exp_dpi_dlogP_T(i), tol);
    }

    EXPECT_NEAR(derivs.dlogn_dlogP_T, exp_dlogn_dlogP_T, max_fp_error(exp_dlogn_dlogP_T));

}

TEST_F(DerivativeTests, thermoDerivativesAreCorrectAfterEquilibration){
    sln->thermo()->equilibrate("HP");

    EXPECT_NEAR(sln->thermo()->temperature(), 2367.8424567893862, 1e-1);
    EXPECT_NEAR(sln->thermo()->enthalpy_mass(), 23985.583414274723, 1e-1);
    

    Goddard::ThermoDerivatives derivs = Goddard::get_thermo_equilibrium_derivatives(*sln);
    ASSERT_EQ(derivs.dpi_dlogP_T.size(),n_elements);
    ASSERT_EQ(derivs.dpi_dlogT_P.size(), n_elements);

    double exp_dlogn_dlogT_P = 0.0198466495473248;
    double exp_dlogn_dlogP_T = -0.0006600595545953462;

    EXPECT_NEAR(derivs.dlogn_dlogT_P, exp_dlogn_dlogT_P, max_fp_error(exp_dlogn_dlogT_P));
    EXPECT_NEAR(derivs.dlogn_dlogP_T, exp_dlogn_dlogP_T, max_fp_error(exp_dlogn_dlogP_T));
}