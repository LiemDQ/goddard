#include "goddard/combustor.hpp"
#include "goddard/mixture_ratio.hpp"
#include "goddard/utils.hpp"

#include <memory>
#include <iostream>
#include "eigen3/Eigen/Dense"
#include "cantera/core.h"
#include "gtest/gtest.h"

constexpr double H2_MOLAR_MASS = 2.01588;
constexpr double O2_MOLAR_MASS = 31.99880;
constexpr size_t NUM_H2O2_SPECIES = 10;

class H2O2CombustorTests: public ::testing::Test {
    protected:
    H2O2CombustorTests() {
        fuel = Cantera::newSolution("h2o2.yaml", "ohmech");
        oxidizer = Cantera::newSolution("h2o2.yaml", "ohmech");
        products = Cantera::newSolution("h2o2.yaml", "ohmech");
        auto fuel_state = fuel->thermo();
        auto ox_state = oxidizer->thermo();

        double ox_temp = 90.15; //K
        double fuel_temp = 20.2; //K
        
        double pressure = 100.0 * Cantera::OneBar;
        fuel_state->setState_TPX(fuel_temp, pressure, "H2: 1");
        ox_state->setState_TPX(ox_temp, pressure, "O2: 1");
        
        OF_ratios = Eigen::ArrayXd(5);
        OF_ratios << 6.0, 7.0, 8.0, 9.0, 10.0;

        reference_molar_ratio = OF_ratios/O2_MOLAR_MASS * H2_MOLAR_MASS;
        reference_ox_fracs = 1- 1/(1 + reference_molar_ratio);
    }

    void SetUp() override {
        MRs = std::make_unique<Goddard::MixtureRatios>(OF_ratios, fuel, oxidizer);

        combustor = std::make_unique<Goddard::Combustor>(fuel, oxidizer, products);
    }

    std::shared_ptr<Cantera::Solution> fuel, oxidizer, products;
    Goddard::CombustionOptions options;
    Eigen::ArrayXd OF_ratios;
    Eigen::ArrayXd reference_molar_ratio;
    Eigen::ArrayXd reference_ox_fracs;

    std::unique_ptr<Goddard::Combustor> combustor;
    std::unique_ptr<Goddard::MixtureRatios> MRs;
};

using namespace Goddard; 

TEST_F(H2O2CombustorTests, mixtureRatiosAreCorrect) {



    Eigen::ArrayXd molar_ratio = MRs->molar_ratio();
    Eigen::ArrayXd fuel_mole_fracs = MRs->fuel_mole_frac();
    Eigen::ArrayXd ox_mole_fracs = MRs->oxidizer_mole_frac();

    for (int i = 0; i < fuel_mole_fracs.size(); i++) {
        EXPECT_NEAR(molar_ratio[i], reference_molar_ratio[i], max_fp_error(reference_molar_ratio[i], 1e-3));
        EXPECT_NEAR(ox_mole_fracs[i], reference_ox_fracs[i], max_fp_error(reference_ox_fracs[i], 1e-3));
        double total = fuel_mole_fracs[i] + ox_mole_fracs[i];
        EXPECT_DOUBLE_EQ(total, 1.0) 
            << "Fuel and oxidizer mole fractions must add up to 1. Fuel: " 
            << fuel_mole_fracs[i] << ", Ox: " << ox_mole_fracs[i];
    }
}

TEST_F(H2O2CombustorTests, moleFracMatrixIsCorrect) {
    
    Eigen::ArrayXXd mole_fracs = combustor->generate_mole_fraction_matrix(*MRs);
    
    //
    ASSERT_EQ(mole_fracs.cols(), NUM_H2O2_SPECIES) 
        << "Mole frac matrix cols should be equal to number of species";
    ASSERT_EQ(mole_fracs.rows(), OF_ratios.size()) 
        << "Mole frac matrix rows should be equal to number of unique compositions.";
    
    //"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"
    for (int i = 0; i < reference_ox_fracs.size(); i++) {
        //check O2
        EXPECT_NEAR(mole_fracs(i,3), reference_ox_fracs(i), max_fp_error(reference_ox_fracs(i), 1e-4))<< "i = " << i;
        //check H2
        EXPECT_NEAR(mole_fracs(i,0), 1-reference_ox_fracs(i), max_fp_error(1-reference_ox_fracs(i), 1e-4)) << "i = " << i;
    }

    for (int i = 0;  i < mole_fracs.rows(); i++){
        EXPECT_NEAR(mole_fracs.row(i).sum(), 1.0, max_fp_error(1.0, 1e-4)) << "i = " << i;
    }
}

TEST_F(H2O2CombustorTests, equilibriumIsCorrect) {
    Eigen::ArrayXd pressures = Eigen::ArrayXd(5);
    pressures << 1,2,3,4,5;
    pressures *= Cantera::OneBar * 10; //10, 20, 30, 40, 50 bar

    Eigen::ArrayXd temperatures = Eigen::ArrayXd(1);
    temperatures << 92.0;
    
    // Eigen::ArrayXXd mole_fracs = combustor->generate_mole_fraction_matrix(*MRs);

    auto results = combustor->solve(
            temperatures,
            pressures, 
            *MRs
        );

    ASSERT_EQ(results.ndim(), 3);
    ASSERT_EQ(results.shape()[0], temperatures.size());
    ASSERT_EQ(results.shape()[2], OF_ratios.size());
    ASSERT_EQ(results.size(), temperatures.size()*pressures.size()*OF_ratios.size());

    //TODO: add tests to verify combustion results

}