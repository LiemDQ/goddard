#include "goddard/combustor.hpp"
#include "goddard/mixture_ratio.hpp"
#include "goddard/numerics.hpp"
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
        auto fuel_thermo = fuel->thermo();
        auto ox_thermo = oxidizer->thermo();

        double ox_temp = 90.15; //K
        double fuel_temp = 20.2; //K
        
        double pressure = 100.0 * Cantera::OneBar;
        fuel_thermo->setState_TPX(fuel_temp, pressure, "H2: 1");
        ox_thermo->setState_TPX(ox_temp, pressure, "O2: 1");
        
        OF_ratios = Eigen::ArrayXd(5);
        OF_ratios << 6.0, 7.0, 8.0, 9.0, 10.0;

        reference_molar_ratio = OF_ratios/O2_MOLAR_MASS * H2_MOLAR_MASS;
        reference_ox_fracs = 1- 1/(1 + reference_molar_ratio);

        fuel_thermo->saveState(fuel_state);
        ox_thermo->saveState(ox_state);
    }

    void SetUp() override {
        MRs = std::make_unique<Goddard::MixtureRatios>(OF_ratios, fuel, oxidizer);

        combustor = std::make_unique<Goddard::Combustor>(products, fuel_state, ox_state);
    }

    std::shared_ptr<Cantera::Solution> fuel, oxidizer, products;
    std::vector<double> fuel_state;
    std::vector<double> ox_state;
    Goddard::CombustorOptions options{ Goddard::CombustorType::INFINITE_AREA, {70.0*Cantera::OneBar}, 0.0, 0.0};
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

    // CombustorOptions options = {CombustorType::INFINITE_AREA, {70.0*Cantera::OneBar}, 0.0, 0.0};
    // Eigen::ArrayXXd mole_fracs = combustor->generate_mole_fraction_matrix(*MRs);

    auto results = combustor->solve(
            temperatures,
            pressures,
            *MRs,
            options
        );

    ASSERT_EQ(results.ndim(), 3);
    ASSERT_EQ(results.shape()[0], temperatures.size());
    ASSERT_EQ(results.shape()[2], OF_ratios.size());
    ASSERT_EQ(results.size(), temperatures.size()*pressures.size()*OF_ratios.size());

    //TODO: add tests to verify combustion results

}

// ---- RecirculatingCombustor tests ----

class H2O2RecirculatingCombustorTests: public ::testing::Test {
    protected:
    H2O2RecirculatingCombustorTests() {
        fuel = Cantera::newSolution("h2o2.yaml", "ohmech");
        oxidizer = Cantera::newSolution("h2o2.yaml", "ohmech");
        products = Cantera::newSolution("h2o2.yaml", "ohmech");
        auto fuel_thermo = fuel->thermo();
        auto ox_thermo = oxidizer->thermo();
        auto prod_thermo = products->thermo();

        double ox_temp = 90.15; //K
        double fuel_temp = 20.2; //K
        double pressure = 100.0 * Cantera::OneBar;

        fuel_thermo->setState_TPX(fuel_temp, pressure, "H2: 1");
        ox_thermo->setState_TPX(ox_temp, pressure, "O2: 1");

        // Create a flue gas state: inert N2 diluent at ambient temperature
        prod_thermo->setState_TPX(300.0, pressure, "N2: 1");

        OF_ratios = Eigen::ArrayXd(3);
        OF_ratios << 4.0, 5.0, 6.0;

        fuel_thermo->saveState(fuel_state);
        ox_thermo->saveState(ox_state);
        prod_thermo->saveState(flue_state);
    }

    void SetUp() override {
        MRs = std::make_unique<Goddard::MixtureRatios>(OF_ratios, fuel, oxidizer);
    }

    std::shared_ptr<Cantera::Solution> fuel, oxidizer, products;
    std::vector<double> fuel_state, ox_state, flue_state;
    Goddard::CombustorOptions options{Goddard::CombustorType::INFINITE_AREA, {70.0*Cantera::OneBar}, 0.0, 0.0};
    Eigen::ArrayXd OF_ratios;
    std::unique_ptr<Goddard::MixtureRatios> MRs;
};

TEST_F(H2O2RecirculatingCombustorTests, massFracMatrixRowsSumToOne) {
    double recircRatio = 0.3;
    DilutedCombustor combustor(products, fuel_state, ox_state, flue_state);

    Eigen::ArrayXXd mass_fracs = combustor.generate_mass_fraction_matrix(*MRs, recircRatio);

    ASSERT_EQ(mass_fracs.rows(), OF_ratios.size());

    for (int i = 0; i < mass_fracs.rows(); i++) {
        EXPECT_NEAR(mass_fracs.row(i).sum(), 1.0, max_fp_error(1.0, 1e-10))
            << "Row " << i << " mass fractions must sum to 1.0";
    }
}

TEST_F(H2O2RecirculatingCombustorTests, moleFracMatrixRowsSumToOne) {
    double recircRatio = 0.3;
    DilutedCombustor combustor(products, fuel_state, ox_state, flue_state);

    Eigen::ArrayXXd mole_fracs = combustor.generate_mole_fraction_matrix(*MRs, recircRatio);

    ASSERT_EQ(mole_fracs.rows(), OF_ratios.size());

    for (int i = 0; i < mole_fracs.rows(); i++) {
        EXPECT_NEAR(mole_fracs.row(i).sum(), 1.0, max_fp_error(1.0, 1e-10))
            << "Row " << i << " mole fractions must sum to 1.0";
    }
}

TEST_F(H2O2RecirculatingCombustorTests, zeroRecirculationMatchesCombustor) {
    double recircRatio = 0.0;
    DilutedCombustor recircCombustor(products, fuel_state, ox_state, flue_state);
    Combustor baseCombustor(products, fuel_state, ox_state);

    Eigen::ArrayXXd recircMass = recircCombustor.generate_mass_fraction_matrix(*MRs, recircRatio);
    Eigen::ArrayXXd baseMass = baseCombustor.generate_mass_fraction_matrix(*MRs);

    ASSERT_EQ(recircMass.rows(), baseMass.rows());
    ASSERT_EQ(recircMass.cols(), baseMass.cols());

    for (int i = 0; i < recircMass.rows(); i++) {
        for (int j = 0; j < recircMass.cols(); j++) {
            EXPECT_NEAR(recircMass(i,j), baseMass(i,j), max_fp_error(baseMass(i,j), 1e-10))
                << "Mismatch at (" << i << ", " << j << ")";
        }
    }
}

TEST_F(H2O2RecirculatingCombustorTests, solveProducesValidEquilibriumStates) {
    // Verify that solve() with recirculating flue gas runs to completion
    // and produces equilibrium states with temperatures above the initial 300K
    // (i.e. combustion actually occurred).
    Eigen::ArrayXd pressures(2);
    pressures << 50.0 * Cantera::OneBar, 70.0 * Cantera::OneBar;

    Eigen::ArrayXd temperatures(1);
    temperatures << 300.0;

    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    DilutedCombustor combustor(sln, fuel_state, ox_state, flue_state);

    auto results = combustor.solve(temperatures, pressures, *MRs, 0.3, options);

    ASSERT_EQ(results.ndim(), 3);
    ASSERT_EQ(results.size(), temperatures.size() * pressures.size() * OF_ratios.size());

    auto thermo = sln->thermo();
    for (int i = 0; i < results.size(); i++) {
        thermo->restoreState(results.get_state(i));
        double temp = thermo->temperature();
        EXPECT_GT(temp, 300.0)
            << "Equilibrated state should be hotter than initial 300K at index " << i;
    }
}