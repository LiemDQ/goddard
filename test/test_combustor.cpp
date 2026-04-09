#include "goddard/combustor.hpp"
#include "goddard/gas.hpp"
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
    H2O2CombustorTests():
        gas(Cantera::newSolution("h2o2.yaml", "ohmech"))
    {
        fuel_comp = {{"H2", 1.0}};
        oxidizer_comp = {{"O2", 1.0}};
        fuel_temperature = 20.2;    // K
        oxidizer_temperature = 90.15; // K

        OF_ratios = Eigen::ArrayXd(5);
        OF_ratios << 6.0, 7.0, 8.0, 9.0, 10.0;

        reference_molar_ratio = OF_ratios/O2_MOLAR_MASS * H2_MOLAR_MASS;
        reference_ox_fracs = 1- 1/(1 + reference_molar_ratio);
    }

    void SetUp() override {
        combustor = std::make_unique<Goddard::Combustor>(gas, fuel_comp, oxidizer_comp);
    }

    Goddard::Gas gas;
    Goddard::Composition fuel_comp;
    Goddard::Composition oxidizer_comp;
    double fuel_temperature;
    double oxidizer_temperature;

    Goddard::CombustorOptions options{ Goddard::CombustorType::INFINITE_AREA,
        Goddard::MixtureRatioType::OF_RATIO, {70.0*Cantera::OneBar}, 0.0, 0.0};
    Eigen::ArrayXd OF_ratios;
    Eigen::ArrayXd reference_molar_ratio;
    Eigen::ArrayXd reference_ox_fracs;

    std::unique_ptr<Goddard::Combustor> combustor;
};

using namespace Goddard;

TEST_F(H2O2CombustorTests, moleFracMatrixIsCorrect) {
    Eigen::ArrayXXd mole_fracs = combustor->generate_mole_fraction_matrix(
        OF_ratios, MixtureRatioType::OF_RATIO);

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

    auto results = combustor->solve(
            temperatures,
            pressures,
            OF_ratios,
            options
        );

    ASSERT_EQ(results.ndim(), 3);
    ASSERT_EQ(results.shape()[0], temperatures.size());
    ASSERT_EQ(results.shape()[2], OF_ratios.size());
    ASSERT_EQ(results.size(), temperatures.size()*pressures.size()*OF_ratios.size());
}

// ---- RecirculatingCombustor tests ----

class H2O2RecirculatingCombustorTests: public ::testing::Test {
    protected:
    H2O2RecirculatingCombustorTests():
        gas(Cantera::newSolution("h2o2.yaml", "ohmech"))
    {
        fuel_comp = {{"H2", 1.0}};
        oxidizer_comp = {{"O2", 1.0}};
        flue_comp = {{"N2", 1.0}};
        fuel_temperature = 20.2;
        oxidizer_temperature = 90.15;
        flue_temperature = 300.0;

        OF_ratios = Eigen::ArrayXd(3);
        OF_ratios << 4.0, 5.0, 6.0;
    }

    Goddard::Gas gas;
    Goddard::Composition fuel_comp, oxidizer_comp, flue_comp;
    double fuel_temperature, oxidizer_temperature, flue_temperature;

    Goddard::CombustorOptions options{Goddard::CombustorType::INFINITE_AREA,
        Goddard::MixtureRatioType::OF_RATIO, {70.0*Cantera::OneBar}, 0.0, 0.0};
    Eigen::ArrayXd OF_ratios;
};

TEST_F(H2O2RecirculatingCombustorTests, massFracMatrixRowsSumToOne) {
    double recircRatio = 0.3;
    DilutedCombustor combustor(gas, fuel_comp, oxidizer_comp, flue_comp);

    Eigen::ArrayXXd mass_fracs = combustor.generate_mass_fraction_matrix(
        OF_ratios, MixtureRatioType::OF_RATIO, recircRatio);

    ASSERT_EQ(mass_fracs.rows(), OF_ratios.size());

    for (int i = 0; i < mass_fracs.rows(); i++) {
        EXPECT_NEAR(mass_fracs.row(i).sum(), 1.0, max_fp_error(1.0, 1e-10))
            << "Row " << i << " mass fractions must sum to 1.0";
    }
}

TEST_F(H2O2RecirculatingCombustorTests, moleFracMatrixRowsSumToOne) {
    double recircRatio = 0.3;
    DilutedCombustor combustor(gas, fuel_comp, oxidizer_comp, flue_comp);

    Eigen::ArrayXXd mole_fracs = combustor.generate_mole_fraction_matrix(
        OF_ratios, MixtureRatioType::OF_RATIO, recircRatio);

    ASSERT_EQ(mole_fracs.rows(), OF_ratios.size());

    for (int i = 0; i < mole_fracs.rows(); i++) {
        EXPECT_NEAR(mole_fracs.row(i).sum(), 1.0, max_fp_error(1.0, 1e-10))
            << "Row " << i << " mole fractions must sum to 1.0";
    }
}

TEST_F(H2O2RecirculatingCombustorTests, zeroRecirculationMatchesCombustor) {
    double recircRatio = 0.0;
    DilutedCombustor recircCombustor(gas, fuel_comp, oxidizer_comp, flue_comp);
    Combustor baseCombustor(gas, fuel_comp, oxidizer_comp);

    Eigen::ArrayXXd recircMass = recircCombustor.generate_mass_fraction_matrix(
        OF_ratios, MixtureRatioType::OF_RATIO, recircRatio);
    Eigen::ArrayXXd baseMass = baseCombustor.generate_mass_fraction_matrix(
        OF_ratios, MixtureRatioType::OF_RATIO);

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
    Eigen::ArrayXd pressures(2);
    pressures << 50.0 * Cantera::OneBar, 70.0 * Cantera::OneBar;

    Eigen::ArrayXd temperatures(1);
    temperatures << 300.0;

    Gas solve_gas(Cantera::newSolution("h2o2.yaml", "ohmech"));
    DilutedCombustor combustor(solve_gas, fuel_comp, oxidizer_comp, flue_comp);

    auto results = combustor.solve(temperatures, pressures, OF_ratios, 0.3, options);

    ASSERT_EQ(results.ndim(), 3);
    ASSERT_EQ(results.size(), temperatures.size() * pressures.size() * OF_ratios.size());

    auto thermo = solve_gas.thermo();
    for (int i = 0; i < results.size(); i++) {
        thermo->restoreState(results.get_state(i));
        double temp = thermo->temperature();
        EXPECT_GT(temp, 300.0)
            << "Equilibrated state should be hotter than initial 300K at index " << i;
    }
}
