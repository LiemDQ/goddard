#include "goddard/equilibrium.hpp"
#include "goddard/utils.hpp"
#include "goddard/numerics.hpp"

#include "eigen3/Eigen/Dense"
#include "cantera/core.h"
#include "gtest/gtest.h"
#include <memory>
#include <vector>
#include <cmath>
#include <string>
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
    std::vector<std::string> elements = {"O", "H", "Ar", "N"};
    std::vector<std::string> species = {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"};
    size_t n_species;
    size_t n_elements;
    

};

using namespace Goddard;


TEST_F(DerivativeTests, stoichiometricCoeffsAreCorrect){
    auto coeffs = Goddard::get_stoichiometric_coeffs(*sln->thermo());
    Eigen::ArrayXXd expected_coeffs{
        {0, 0, 1, 2, 1, 1, 2, 2, 0, 0},
        {2, 1, 0, 0, 1, 2, 1, 2, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 2},
    };

    long n_elements_signed = static_cast<long>(n_elements); //suppress compiler warnings
    long n_species_signed = static_cast<long>(n_species);

    EXPECT_EQ(coeffs.size(), expected_coeffs.size());

    expected_coeffs.transposeInPlace();

    for (int i = 0; i < n_elements_signed; i++){
        for (int j = 0; j < n_species_signed; j++){
            std::string stoich_id = "Species: ";
            stoich_id += species[j];
            stoich_id += ", element: ";
            stoich_id += elements[i];
            
            EXPECT_EQ(coeffs(j,i), expected_coeffs(j,i)) << stoich_id;
        }
    }
}

TEST_F(DerivativeTests, moleVectorIsCorrect){
    auto moles = Goddard::get_mole_vector(*sln->thermo());
    
    std::vector<double> expected_moles{0., 0., 0.,  0.01219185, 0., 0.01219185, 0. , 0. , 0.001219185, 0.01219185};
    ASSERT_EQ(moles.size(), n_species);

    for (size_t i = 0; i < n_species; i++){
        double tol = max_fp_error(expected_moles[i], 1e-5, 1e-7);
        EXPECT_NEAR(moles(i), expected_moles[i], tol);
    }
}

TEST_F(DerivativeTests, moleVectorIsCorrectAfterEquilibration){
    sln->thermo()->equilibrate("HP");
    auto moles = Goddard::get_mole_vector(*sln->thermo());
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
    auto enthalpies = Goddard::get_enthalpyRT_vector(*sln->thermo());
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
    
    Goddard::EquilibriumDerivatives derivs = Goddard::get_thermo_equilibrium_derivatives(*sln->thermo());
    ASSERT_EQ(derivs.dpi_dlogP_T.size(),n_elements);
    ASSERT_EQ(derivs.dpi_dlogT_P.size(), n_elements);

    Eigen::ArrayXd exp_dpi_dlogT_P(n_elements);
    exp_dpi_dlogT_P << -1.86677477,  4.63052867, -2.18942708, -1.77019042;
    
    double exp_dlogn_dlogT_P = -1.8359417928589634e-16;

    Eigen::ArrayXd exp_dpi_dlogP_T(n_elements);
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
    

    Goddard::EquilibriumDerivatives derivs = Goddard::get_thermo_equilibrium_derivatives(*sln->thermo());
    ASSERT_EQ(derivs.dpi_dlogP_T.size(),n_elements);
    ASSERT_EQ(derivs.dpi_dlogT_P.size(), n_elements);

    double exp_dlogn_dlogT_P = 0.0198466495473248;
    double exp_dlogn_dlogP_T = -0.0006600595545953462;

    EXPECT_NEAR(derivs.dlogn_dlogT_P, exp_dlogn_dlogT_P, max_fp_error(exp_dlogn_dlogT_P));
    EXPECT_NEAR(derivs.dlogn_dlogP_T, exp_dlogn_dlogP_T, max_fp_error(exp_dlogn_dlogP_T));
}

class PropertyTests: public ::testing::Test {
    protected:
    PropertyTests() {
        this->sln = Cantera::newSolution("h2o2.yaml", "ohmech");
        double temp = 2400.0; //K
        double pressure = 50.0*Cantera::OneAtm;
        auto gas = sln->thermo();
        
        n_species = sln->thermo()->nSpecies();
        n_elements = sln->thermo()->nElements();

        gas->setState_TPX(temp, pressure, "H2O:1, N2:1, O2:1. AR:0.1"); //completely random composition lol
    }
    std::shared_ptr<Cantera::Solution> sln;
    Goddard::EquilibriumProperties expected_props{0.9999999999999998,-0.9999999999999999,1604.459611106932,1.2435582661124545};
    Goddard::EquilibriumProperties expected_eq_props{1.0198466495473248,-1.0006600595545954,1796.3940426543525,1.222008621037549};
    size_t n_species;
    size_t n_elements;
};

TEST_F(PropertyTests, equilibriumPropertiesAreCorrect) {
    auto derivs = Goddard::get_thermo_equilibrium_derivatives(*sln->thermo());
    Goddard::EquilibriumProperties props = Goddard::get_thermo_equilibrium_properties(*sln->thermo(), derivs);
    
    EXPECT_NEAR(props.dlogV_dlogT_P, expected_props.dlogV_dlogT_P, max_fp_error(expected_props.dlogV_dlogT_P));
    EXPECT_NEAR(props.dlogV_dlogP_T, expected_props.dlogV_dlogP_T, max_fp_error(expected_props.dlogV_dlogP_T));
    EXPECT_NEAR(props.spec_heat_p, expected_props.spec_heat_p, max_fp_error(expected_props.spec_heat_p));
    EXPECT_NEAR(props.gamma_s, expected_props.gamma_s, max_fp_error(expected_props.gamma_s));
}

TEST_F(PropertyTests, equilibriumPropertiesAreCorrectAfterEquilibrium) {
    sln->thermo()->equilibrate("HP", "gibbs");
    auto derivs = Goddard::get_thermo_equilibrium_derivatives(*sln->thermo());
    Goddard::EquilibriumProperties props = Goddard::get_thermo_equilibrium_properties(*sln->thermo(), derivs);
    
    EXPECT_NEAR(props.dlogV_dlogT_P, expected_eq_props.dlogV_dlogT_P, max_fp_error(expected_eq_props.dlogV_dlogT_P));
    EXPECT_NEAR(props.dlogV_dlogP_T, expected_eq_props.dlogV_dlogP_T, max_fp_error(expected_eq_props.dlogV_dlogP_T));
    EXPECT_NEAR(props.spec_heat_p, expected_eq_props.spec_heat_p, max_fp_error(expected_eq_props.spec_heat_p));
    EXPECT_NEAR(props.gamma_s, expected_eq_props.gamma_s, max_fp_error(expected_eq_props.gamma_s));
}