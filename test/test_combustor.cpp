#include "goddard/combustor.hpp"
#include "goddard/gas.hpp"
#include "goddard/numerics.hpp"
#include "goddard/utils.hpp"
#include "goddard/rocket.hpp"
#include "goddard/error.hpp"
#include "goddard/config.h"
#include "goddard/format.hpp"
#include "status_test_helpers.hpp"
#include <limits>
#include <algorithm>
#include <memory>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>
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
        Goddard::MixtureRatioType::OF_RATIO, 0.0, 0.0};
    Eigen::ArrayXd OF_ratios;
    Eigen::ArrayXd reference_molar_ratio;
    Eigen::ArrayXd reference_ox_fracs;

    std::unique_ptr<Goddard::Combustor> combustor;
};

using namespace Goddard;
using Goddard::test_helpers::failed_with;

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
        Goddard::MixtureRatioType::OF_RATIO, 0.0, 0.0};
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

TEST_F(H2O2RecirculatingCombustorTests, negativeFractionsAreRejectedInReleaseBuilds) {
    // A negative recirculation ratio subtracts flue gas, and the N2 of the flue then has a negative
    // mole and mass fraction in the H2/O2 feed. That is an invalid input, not a solver failure.
    DilutedCombustor combustor(gas, fuel_comp, oxidizer_comp, flue_comp);
    const double negative_recirculation = -0.5;
    EXPECT_THROW(combustor.generate_mole_fraction_matrix(
                     OF_ratios, MixtureRatioType::OF_RATIO, negative_recirculation),
                 std::invalid_argument);
    EXPECT_THROW(combustor.generate_mass_fraction_matrix(
                     OF_ratios, MixtureRatioType::OF_RATIO, negative_recirculation),
                 std::invalid_argument);

    Eigen::ArrayXd negative_ratio(1);
    negative_ratio << -2.0;
    EXPECT_THROW(combustor.generate_mass_fraction_matrix(
                     negative_ratio, MixtureRatioType::OF_RATIO, 0.3),
                 std::invalid_argument);
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

TEST_F(H2O2CombustorTests, stringCompositionConstructorParsesSpecies) {
    // The Gas argument is moved into the combustor before the compositions are parsed.
    Gas string_gas(Cantera::newSolution("h2o2.yaml", "ohmech"));
    Combustor string_combustor(string_gas, "H2:1.0", "O2:1.0");

    Eigen::ArrayXXd from_strings = string_combustor.generate_mole_fraction_matrix(
        OF_ratios, MixtureRatioType::OF_RATIO);
    Eigen::ArrayXXd from_maps = combustor->generate_mole_fraction_matrix(
        OF_ratios, MixtureRatioType::OF_RATIO);

    ASSERT_EQ(from_strings.rows(), from_maps.rows());
    ASSERT_EQ(from_strings.cols(), from_maps.cols());
    for (long i = 0; i < from_maps.rows(); i++) {
        for (long j = 0; j < from_maps.cols(); j++) {
            EXPECT_NEAR(from_strings(i, j), from_maps(i, j), 1e-12);
        }
    }
}

TEST_F(H2O2CombustorTests, malformedOrUnknownCompositionStringsAreInvalidArguments) {
    // An unknown species and a malformed entry are bad arguments, not Cantera errors.
    EXPECT_THROW(Combustor(gas, "H2:1", "O2:1, XX:1"), std::invalid_argument);
    EXPECT_THROW(Combustor(gas, "H2:abc", "O2:1"), std::invalid_argument);
    EXPECT_THROW(DilutedCombustor(gas, "H2:1", "O2:1", "XX:1"), std::invalid_argument);
    EXPECT_THROW(PhaseSpecification(300.0, 1.0e5, "H2:abc"), std::invalid_argument);
    EXPECT_NO_THROW(PhaseSpecification(300.0, 1.0e5, "H2:1, O2:2"));
}

TEST_F(H2O2CombustorTests, streamCombustorValidatesPressuresAndEveryMixtureRatio) {
    Eigen::ArrayXd good_pressures(1);
    good_pressures << 10.0 * Cantera::OneBar;
    Eigen::ArrayXd ratios(1);
    ratios << 6.0;

    Gas fuel(Cantera::newSolution("h2o2.yaml", "ohmech"));
    fuel.set_state_TPX(300.0, good_pressures(0), "H2:1");
    Gas oxidizer(Cantera::newSolution("h2o2.yaml", "ohmech"));
    oxidizer.set_state_TPX(300.0, good_pressures(0), "O2:1");
    Gas products(Cantera::newSolution("h2o2.yaml", "ohmech"));
    Combustor stream_combustor(products, fuel, oxidizer);

    for (double bad : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN()}) {
        Eigen::ArrayXd bad_pressures(2);
        bad_pressures << good_pressures(0), bad;
        EXPECT_THROW(stream_combustor.solve(bad_pressures, ratios, options), std::invalid_argument)
            << "pressure " << bad;
    }

    // The bad ratio comes last: it is rejected before the good one is solved.
    Eigen::ArrayXd bad_ratios(2);
    bad_ratios << 6.0, -1.0;
    EXPECT_THROW(stream_combustor.solve(good_pressures, bad_ratios, options), std::invalid_argument);
    // An equivalence ratio is not implemented on this path, also when it is not the first value
    CombustorOptions phi = options;
    phi.mixture_type = MixtureRatioType::PHI_RATIO;
    EXPECT_THROW(stream_combustor.solve(good_pressures, ratios, phi), NotImplementedError);
}

TEST_F(H2O2CombustorTests, streamCombustorReportsFailedPointsInTheStatuses) {
    Eigen::ArrayXd pressures(2);
    pressures << 10.0 * Cantera::OneBar, 20.0 * Cantera::OneBar;
    Eigen::ArrayXd ratios(2);
    ratios << 4.0, 6.0;

    Gas oxidizer(Cantera::newSolution("h2o2.yaml", "ohmech"));
    oxidizer.set_state_TPX(300.0, pressures(0), "O2:1");
    Gas products(Cantera::newSolution("h2o2.yaml", "ohmech"));

    // A fuel stream at 1e5 K has an enthalpy no combustion temperature can match, at every point
    Gas hot_fuel(Cantera::newSolution("h2o2.yaml", "ohmech"));
    hot_fuel.set_state_TPX(1.0e5, pressures(0), "H2:1");
    Combustor hot_combustor(products, hot_fuel, oxidizer);
    ThermoArray hot_states = hot_combustor.solve(pressures, ratios, options);
    EXPECT_EQ(hot_states.num_failed(), 4);
    EXPECT_FALSE(hot_states.ok().any());
    EXPECT_TRUE(failed_with(hot_states.status(), FailureReason::THERMO_FAILURE));
    EXPECT_NE(hot_states.status_at(0).stage.find("combustion at O/F = 4"), std::string::npos)
        << hot_states.status_at(0).describe();
    EXPECT_TRUE(hot_states.temperature().isNaN().all());

    // The same combustor with a cold fuel stream solves every point
    Gas cold_fuel(Cantera::newSolution("h2o2.yaml", "ohmech"));
    cold_fuel.set_state_TPX(300.0, pressures(0), "H2:1");
    Combustor cold_combustor(products, cold_fuel, oxidizer);
    ThermoArray cold_states = cold_combustor.solve(pressures, ratios, options);
    EXPECT_EQ(cold_states.num_failed(), 0);
    EXPECT_TRUE(cold_states.temperature().isFinite().all());
}

TEST_F(H2O2CombustorTests, emptyInputGridsAreRejectedBeforeIndexing) {
    // An empty pressure or mixture-ratio grid has no combustion state to compute. Every solve
    // entry point must say so instead of reading pressures[0] out of bounds.
    const Eigen::ArrayXd empty(0);
    Eigen::ArrayXd temperatures(1);
    temperatures << 300.0;
    Eigen::ArrayXd pressures(1);
    pressures << 10.0 * Cantera::OneBar;
    Eigen::ArrayXd mixture_ratios(1);
    mixture_ratios << 6.0;

    EXPECT_THROW(combustor->solve(300.0, 300.0, empty, mixture_ratios, options),
                 std::invalid_argument);
    EXPECT_THROW(combustor->solve(300.0, 300.0, pressures, empty, options),
                 std::invalid_argument);
    EXPECT_THROW(combustor->solve(temperatures, empty, mixture_ratios, options),
                 std::invalid_argument);
    EXPECT_THROW(combustor->solve(temperatures, pressures, empty, options),
                 std::invalid_argument);
    EXPECT_THROW(combustor->solve(empty, pressures, mixture_ratios, options),
                 std::invalid_argument);

    Gas fuel(Cantera::newSolution("h2o2.yaml", "ohmech"));
    fuel.set_state_TPX(300.0, pressures(0), "H2:1");
    Gas oxidizer(Cantera::newSolution("h2o2.yaml", "ohmech"));
    oxidizer.set_state_TPX(300.0, pressures(0), "O2:1");
    Gas products(Cantera::newSolution("h2o2.yaml", "ohmech"));
    Combustor stream_combustor(products, fuel, oxidizer);
    EXPECT_THROW(stream_combustor.solve(empty, mixture_ratios, options), std::invalid_argument);
    EXPECT_THROW(stream_combustor.solve(pressures, empty, options), std::invalid_argument);

    DilutedCombustor diluted(gas, fuel_comp, oxidizer_comp, Composition{{"N2", 1.0}});
    EXPECT_THROW(diluted.solve(300.0, 300.0, 300.0, empty, mixture_ratios, 0.1, options),
                 std::invalid_argument);
    EXPECT_THROW(diluted.solve(300.0, 300.0, 300.0, pressures, empty, 0.1, options),
                 std::invalid_argument);
    EXPECT_THROW(diluted.solve(temperatures, empty, mixture_ratios, 0.1, options),
                 std::invalid_argument);
    EXPECT_THROW(diluted.solve(temperatures, pressures, empty, 0.1, options),
                 std::invalid_argument);

}

TEST_F(H2O2CombustorTests, outOfRangeMixtureRatiosAreRejectedInReleaseBuilds) {
    // A negative O/F ratio maps to a fuel mass fraction outside [0, 1]: there is no such mixture.
    // The check must be an exception, not an assert that vanishes under NDEBUG.
    Eigen::ArrayXd negative_ratio(2);
    negative_ratio << 6.0, -2.0;
    EXPECT_THROW(combustor->generate_mole_fraction_matrix(negative_ratio, MixtureRatioType::OF_RATIO),
                 std::invalid_argument);
    EXPECT_THROW(combustor->generate_mass_fraction_matrix(negative_ratio, MixtureRatioType::OF_RATIO),
                 std::invalid_argument);

    Eigen::ArrayXd pressures(1);
    pressures << 10.0 * Cantera::OneBar;
    EXPECT_THROW(combustor->solve(300.0, 300.0, pressures, negative_ratio, options),
                 std::invalid_argument);

    Eigen::ArrayXd fuel_fraction(1);
    fuel_fraction << 1.5;
    EXPECT_THROW(combustor->generate_mass_fraction_matrix(fuel_fraction, MixtureRatioType::FUEL_FRAC),
                 std::invalid_argument);

    // The reactant-stream path checks the same domain.
    Gas fuel(Cantera::newSolution("h2o2.yaml", "ohmech"));
    fuel.set_state_TPX(300.0, pressures(0), "H2:1");
    Gas oxidizer(Cantera::newSolution("h2o2.yaml", "ohmech"));
    oxidizer.set_state_TPX(300.0, pressures(0), "O2:1");
    Gas products(Cantera::newSolution("h2o2.yaml", "ohmech"));
    Combustor stream_combustor(products, fuel, oxidizer);
    EXPECT_THROW(stream_combustor.solve(pressures, negative_ratio, options), std::invalid_argument);
}

TEST_F(H2O2CombustorTests, nonpositiveTemperaturesAndPressuresAreRejected) {
    Eigen::ArrayXd temperatures(1);
    temperatures << 300.0;
    Eigen::ArrayXd pressures(1);
    pressures << 10.0 * Cantera::OneBar;
    Eigen::ArrayXd bad_pressures(2);
    bad_pressures << 10.0 * Cantera::OneBar, 0.0;
    Eigen::ArrayXd bad_temperatures(2);
    bad_temperatures << 300.0, -10.0;
    Eigen::ArrayXd ratios(1);
    ratios << 6.0;

    EXPECT_THROW(combustor->solve(0.0, 300.0, pressures, ratios, options), std::invalid_argument);
    EXPECT_THROW(combustor->solve(300.0, -1.0, pressures, ratios, options), std::invalid_argument);
    EXPECT_THROW(combustor->solve(300.0, 300.0, bad_pressures, ratios, options), std::invalid_argument);
    EXPECT_THROW(combustor->solve(bad_temperatures, pressures, ratios, options), std::invalid_argument);
    EXPECT_THROW(combustor->solve(temperatures, bad_pressures, ratios, options), std::invalid_argument);
    EXPECT_THROW(combustor->solve(300.0, std::numeric_limits<double>::quiet_NaN(), pressures, ratios, options),
                 std::invalid_argument);

    DilutedCombustor diluted(gas, fuel_comp, oxidizer_comp, Composition{{"N2", 1.0}});
    EXPECT_THROW(diluted.solve(300.0, 300.0, 0.0, pressures, ratios, 0.1, options), std::invalid_argument);
    EXPECT_THROW(diluted.solve(temperatures, bad_pressures, ratios, 0.1, options), std::invalid_argument);
}

TEST_F(H2O2CombustorTests, unsolvableOperatingPointsAreReportedInTheArrayStatuses) {
    // A 1e5 K fuel stream has an enthalpy that no combustion temperature can match. Nothing is
    // thrown: every point records the failure and its properties are NaN.
    Eigen::ArrayXd pressures(2);
    pressures << 10.0 * Cantera::OneBar, 20.0 * Cantera::OneBar;
    Eigen::ArrayXd ratios(2);
    ratios << 4.0, 6.0;

    ThermoArray states = combustor->solve(1.0e5, 300.0, pressures, ratios, options);
    EXPECT_EQ(states.num_failed(), 4);
    EXPECT_FALSE(states.ok().any());
    EXPECT_TRUE(failed_with(states.status(), FailureReason::THERMO_FAILURE));
    EXPECT_NE(states.status_at(0).stage.find("combustion at O/F = 4"), std::string::npos)
        << states.status_at(0).describe();
    EXPECT_TRUE(states.temperature().isNaN().all());

    // A good fuel temperature on the same combustor solves every point
    ThermoArray good = combustor->solve(300.0, 300.0, pressures, ratios, options);
    EXPECT_EQ(good.num_failed(), 0);
    EXPECT_TRUE(good.temperature().isFinite().all());
}

TEST_F(H2O2CombustorTests, compositionOverloadsRejectACombustorBuiltFromStreams) {
    Eigen::ArrayXd temperatures(1);
    temperatures << 300.0;
    Eigen::ArrayXd pressures(1);
    pressures << 10.0 * Cantera::OneBar;
    Eigen::ArrayXd ratios(1);
    ratios << 6.0;

    Gas fuel(Cantera::newSolution("h2o2.yaml", "ohmech"));
    fuel.set_state_TPX(300.0, pressures(0), "H2:1");
    Gas oxidizer(Cantera::newSolution("h2o2.yaml", "ohmech"));
    oxidizer.set_state_TPX(300.0, pressures(0), "O2:1");
    Gas products(Cantera::newSolution("h2o2.yaml", "ohmech"));
    Combustor stream_combustor(products, fuel, oxidizer);

    EXPECT_THROW(stream_combustor.solve(temperatures, pressures, ratios, options), std::invalid_argument);
    EXPECT_THROW(stream_combustor.solve(300.0, 300.0, pressures, ratios, options), std::invalid_argument);
}

// ---- Isochoric combustion ----

TEST_F(H2O2CombustorTests, isochoricConservesInternalEnergyAndVolume) {
    Eigen::ArrayXd temperatures(1);
    temperatures << 300.0;
    Eigen::ArrayXd pressures(1);
    pressures << 1.0 * Cantera::OneBar;
    Eigen::ArrayXd mixture_ratios(1);
    mixture_ratios << 8.0;

    CombustorOptions isochoric = options;
    isochoric.process = CombustionProcess::ISOCHORIC;
    ThermoArray burnt = combustor->solve(temperatures, pressures, mixture_ratios, isochoric);

    auto reactants = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    Eigen::ArrayXXd mole_fracs = combustor->generate_mole_fraction_matrix(
        mixture_ratios, MixtureRatioType::OF_RATIO);
    Eigen::ArrayXd reactant_mole_fracs = mole_fracs.row(0);
    reactants->setMoleFractions(reactant_mole_fracs.data());
    reactants->setState_TP(temperatures(0), pressures(0));
    const double u_reactants = reactants->intEnergy_mass();
    const double v_reactants = 1.0 / reactants->density();

    auto products = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    for (int i = 0; i < burnt.size(); i++) {
        products->restoreState(burnt.get_state(i));
        EXPECT_NEAR(products->intEnergy_mass(), u_reactants, 1e-6 * std::abs(u_reactants)) << "i = " << i;
        EXPECT_NEAR(1.0 / products->density(), v_reactants, 1e-6 * v_reactants) << "i = " << i;
        EXPECT_GT(products->temperature(), 3000.0) << "i = " << i;
    }
}

TEST_F(H2O2CombustorTests, isochoricIsHotterThanIsobaric) {
    Eigen::ArrayXd temperatures(1);
    temperatures << 300.0;
    Eigen::ArrayXd pressures(1);
    pressures << 1.0 * Cantera::OneBar;
    Eigen::ArrayXd mixture_ratios(1);
    mixture_ratios << 8.0;

    CombustorOptions isochoric = options;
    isochoric.process = CombustionProcess::ISOCHORIC;
    ThermoArray burnt_uv = combustor->solve(temperatures, pressures, mixture_ratios, isochoric);
    ThermoArray burnt_hp = combustor->solve(temperatures, pressures, mixture_ratios, options);

    auto thermo = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    thermo->restoreState(burnt_uv.get_state(0));
    const double T_uv = thermo->temperature();
    const double P_uv = thermo->pressure();
    thermo->restoreState(burnt_hp.get_state(0));
    const double T_hp = thermo->temperature();

    // With no expansion work done, constant-volume combustion reaches a higher temperature.
    EXPECT_GT(T_uv, T_hp);
    EXPECT_GT(P_uv, pressures(0));
}


// ---- Storage order of combustor results ----

TEST_F(H2O2CombustorTests, adiabaticSolveStoresEntriesByIndex) {
    const double reactant_temperature = 300.0;
    Eigen::ArrayXd pressures(3);
    pressures << 10.0 * Cantera::OneBar, 20.0 * Cantera::OneBar, 30.0 * Cantera::OneBar;
    // Kept below O/F 8: Cantera's "gibbs" HP solver fails to converge for hotter H2/O2 flames here.
    Eigen::ArrayXd mixture_ratios(3);
    mixture_ratios << 4.0, 5.0, 6.0;

    ThermoArray states = combustor->solve(reactant_temperature, reactant_temperature,
        pressures, mixture_ratios, options);
    ASSERT_EQ(states.shape(), (std::vector<long>{mixture_ratios.size(), pressures.size()}));

    Eigen::ArrayXXd T = states.temperature();
    Eigen::ArrayXXd P = states.pressure();
    auto thermo = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    for (long i = 0; i < mixture_ratios.size(); i++) {
        for (long j = 0; j < pressures.size(); j++) {
            Eigen::ArrayXd single_pressure(1);
            single_pressure << pressures(j);
            Eigen::ArrayXd single_ratio(1);
            single_ratio << mixture_ratios(i);
            ThermoArray single = combustor->solve(reactant_temperature, reactant_temperature,
                single_pressure, single_ratio, options);
            thermo->restoreState(single.get_state(0));

            EXPECT_NEAR(P(i, j), pressures(j), 1e-9 * pressures(j)) << "i = " << i << ", j = " << j;
            EXPECT_NEAR(T(i, j), thermo->temperature(), 1e-9 * thermo->temperature())
                << "i = " << i << ", j = " << j;
            EXPECT_EQ(states.get_state(states.flat_index(i, j)), single.get_state(0))
                << "i = " << i << ", j = " << j;
        }
    }
}

// ---- Reactant streams given as Gas objects ----

namespace {

std::string reactant_file() { return std::string(DATA_DIR) + "/nasa9_reactants.yaml"; }
std::string nasa9_gas_file() { return std::string(DATA_DIR) + "/nasa9_gas.yaml"; }

} // namespace

TEST_F(H2O2CombustorTests, reactantGasPathMatchesLegacyIsobaricPath) {
    const double reactant_temperature = 300.0;
    Eigen::ArrayXd pressures(1);
    pressures << 70.0 * Cantera::OneBar;
    // Both paths call Cantera's "gibbs" HP solver, which fails to converge on `h2o2.yaml` above
    // O/F 7 at this pressure, on the legacy path as well as on this one.
    Eigen::ArrayXd of_ratios(4);
    of_ratios << 4.0, 5.0, 6.0, 7.0;

    Gas fuel(Cantera::newSolution("h2o2.yaml", "ohmech"));
    fuel.set_state_TPX(reactant_temperature, pressures(0), "H2:1");
    Gas oxidizer(Cantera::newSolution("h2o2.yaml", "ohmech"));
    oxidizer.set_state_TPX(reactant_temperature, pressures(0), "O2:1");

    Gas products(Cantera::newSolution("h2o2.yaml", "ohmech"));
    Combustor stream_combustor(products, fuel, oxidizer);
    ThermoArray from_streams = stream_combustor.solve(pressures, of_ratios, options);
    ThermoArray from_legacy = combustor->solve(reactant_temperature, reactant_temperature,
        pressures, of_ratios, options);

    ASSERT_EQ(from_streams.shape(), from_legacy.shape());
    auto stream_thermo = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    auto legacy_thermo = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    for (int i = 0; i < from_legacy.size(); i++) {
        stream_thermo->restoreState(from_streams.get_state(i));
        legacy_thermo->restoreState(from_legacy.get_state(i));

        EXPECT_NEAR(stream_thermo->temperature(), legacy_thermo->temperature(),
            1e-6 * legacy_thermo->temperature()) << "i = " << i;
        EXPECT_NEAR(stream_thermo->pressure(), legacy_thermo->pressure(),
            1e-6 * legacy_thermo->pressure()) << "i = " << i;

        std::vector<double> Y_streams(stream_thermo->nSpecies());
        std::vector<double> Y_legacy(legacy_thermo->nSpecies());
        stream_thermo->getMassFractions(Y_streams.data());
        legacy_thermo->getMassFractions(Y_legacy.data());
        for (size_t k = 0; k < Y_legacy.size(); k++) {
            EXPECT_NEAR(Y_streams[k], Y_legacy[k], 1e-6 * std::max(Y_legacy[k], 1e-6))
                << "i = " << i << ", species " << legacy_thermo->speciesName(k);
        }
    }
}

TEST(MultiSpeciesStreamTests, compositionPathBlendsEnthalpyByMassLikeStreamPath) {
    // H2 burnt in air. The composition path blends the stream enthalpies with the fuel mass
    // fraction from `ThermoPhase::mixtureFraction`, whose basis argument describes the stream
    // compositions, which are mole fractions. Read as mass fractions, air (O2:0.21, N2:0.79)
    // carries 21 % O by mass instead of 23.3 %, and at O/F 40 the fuel fraction comes out as
    // about 0.0216 instead of 1/41 = 0.0244. The fuel is preheated to 1000 K
    // (h ~ 10 MJ/kg against ~0 for the air), so that error shifts the mixture enthalpy by
    // about 30 kJ/kg and the flame temperature by roughly 20 K. The reactant-stream path blends by
    // mass fraction 1/(1+O/F) directly and is the reference.
    const double fuel_temperature = 1000.0;     // K
    const double oxidizer_temperature = 300.0;  // K
    const double pressure = 10.0 * Cantera::OneBar;
    const std::string fuel_composition = "H2:1";
    const std::string oxidizer_composition = "O2:0.21, N2:0.79";

    Eigen::ArrayXd pressures(1);
    pressures << pressure;
    Eigen::ArrayXd of_ratios(1);
    of_ratios << 40.0;  // lean H2/air, stoichiometric O/F is 34.3
    const double fuel_mass_fraction = 1.0 / (1.0 + of_ratios(0));
    const CombustorOptions options{CombustorType::INFINITE_AREA, MixtureRatioType::OF_RATIO, 0.0, 0.0};

    Gas fuel(Cantera::newSolution("h2o2.yaml", "ohmech"));
    fuel.set_state_TPX(fuel_temperature, pressure, fuel_composition);
    Gas oxidizer(Cantera::newSolution("h2o2.yaml", "ohmech"));
    oxidizer.set_state_TPX(oxidizer_temperature, pressure, oxidizer_composition);
    const double mixture_enthalpy = fuel_mass_fraction * fuel.enthalpy_mass()
        + (1.0 - fuel_mass_fraction) * oxidizer.enthalpy_mass();

    Gas stream_products(Cantera::newSolution("h2o2.yaml", "ohmech"));
    ThermoArray from_streams = Combustor(stream_products, fuel, oxidizer)
        .solve(pressures, of_ratios, options);

    Gas composition_products(Cantera::newSolution("h2o2.yaml", "ohmech"));
    ThermoArray from_compositions = Combustor(composition_products, fuel_composition,
        oxidizer_composition).solve(fuel_temperature, oxidizer_temperature, pressures, of_ratios,
        options);

    // With zero recirculation the diluted combustor burns the same fresh mixture.
    Gas diluted_products(Cantera::newSolution("h2o2.yaml", "ohmech"));
    ThermoArray from_diluted = DilutedCombustor(diluted_products, fuel_composition,
        oxidizer_composition, "N2:1").solve(fuel_temperature, oxidizer_temperature,
        oxidizer_temperature, pressures, of_ratios, 0.0, options);

    auto reference = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    reference->restoreState(from_streams.get_state(0));
    const double T_reference = reference->temperature();
    EXPECT_NEAR(reference->enthalpy_mass(), mixture_enthalpy,
        max_fp_error(mixture_enthalpy, 1e-6, 1.0));

    auto composition = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    composition->restoreState(from_compositions.get_state(0));
    EXPECT_NEAR(composition->enthalpy_mass(), mixture_enthalpy,
        max_fp_error(mixture_enthalpy, 1e-6, 1.0)) << "first law: chamber h = mass-blended h";
    EXPECT_NEAR(composition->temperature(), T_reference, 1e-6 * T_reference);

    auto diluted = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    diluted->restoreState(from_diluted.get_state(0));
    EXPECT_NEAR(diluted->enthalpy_mass(), mixture_enthalpy,
        max_fp_error(mixture_enthalpy, 1e-6, 1.0)) << "first law: chamber h = mass-blended h";
    EXPECT_NEAR(diluted->temperature(), T_reference, 1e-6 * T_reference);
}

TEST_F(H2O2CombustorTests, reactantGasPathMatchesLegacyIsochoricPath) {
    const double reactant_temperature = 300.0;
    Eigen::ArrayXd pressures(1);
    pressures << 1.0 * Cantera::OneBar;
    Eigen::ArrayXd mixture_ratios(1);
    mixture_ratios << 8.0;

    CombustorOptions isochoric = options;
    isochoric.process = CombustionProcess::ISOCHORIC;

    Gas fuel(Cantera::newSolution("h2o2.yaml", "ohmech"));
    fuel.set_state_TPX(reactant_temperature, pressures(0), "H2:1");
    Gas oxidizer(Cantera::newSolution("h2o2.yaml", "ohmech"));
    oxidizer.set_state_TPX(reactant_temperature, pressures(0), "O2:1");

    Gas products(Cantera::newSolution("h2o2.yaml", "ohmech"));
    Combustor stream_combustor(products, fuel, oxidizer);
    ThermoArray from_streams = stream_combustor.solve(pressures, mixture_ratios, isochoric);
    ThermoArray from_legacy = combustor->solve(reactant_temperature, reactant_temperature,
        pressures, mixture_ratios, isochoric);

    auto stream_thermo = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    auto legacy_thermo = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    stream_thermo->restoreState(from_streams.get_state(0));
    legacy_thermo->restoreState(from_legacy.get_state(0));

    EXPECT_NEAR(stream_thermo->temperature(), legacy_thermo->temperature(),
        1e-6 * legacy_thermo->temperature());
    EXPECT_NEAR(stream_thermo->pressure(), legacy_thermo->pressure(),
        1e-6 * legacy_thermo->pressure());
    EXPECT_GT(stream_thermo->pressure(), pressures(0));
}

TEST_F(H2O2CombustorTests, reactantGasPathRejectsUnsupportedInputs) {
    Eigen::ArrayXd pressures(1);
    pressures << 70.0 * Cantera::OneBar;

    Gas fuel(Cantera::newSolution("h2o2.yaml", "ohmech"));
    fuel.set_state_TPX(300.0, pressures(0), "H2:1");
    Gas oxidizer(Cantera::newSolution("h2o2.yaml", "ohmech"));
    oxidizer.set_state_TPX(300.0, pressures(0), "O2:1");
    Gas products(Cantera::newSolution("h2o2.yaml", "ohmech"));
    Combustor stream_combustor(products, fuel, oxidizer);

    CombustorOptions phi = options;
    phi.mixture_type = MixtureRatioType::PHI_RATIO;
    EXPECT_THROW(stream_combustor.solve(pressures, OF_ratios, phi), NotImplementedError);

    // A finite-area combustor with a valid mass flux is accepted: the combustor itself just
    // produces the ordinary injector-face HP state (the finite-area chamber solve happens on
    // the nozzle side, not here). A single, convergence-safe O/F ratio is used here: `OF_ratios`
    // spans 6-10, and Cantera's "gibbs" HP solver fails to converge on `h2o2.yaml` above O/F 7 at
    // this pressure (see `reactantGasPathMatchesLegacyIsobaricPath`).
    Eigen::ArrayXd single_of_ratio(1);
    single_of_ratio << 6.0;
    CombustorOptions finite_area = options;
    finite_area.type = CombustorType::FINITE_MASS_FLUX;
    finite_area.mass_flux = 1000.0;
    EXPECT_NO_THROW(stream_combustor.solve(pressures, single_of_ratio, finite_area));

    // mass_flux <= 0 is rejected.
    CombustorOptions zero_mass_flux = finite_area;
    zero_mass_flux.mass_flux = 0.0;
    EXPECT_THROW(stream_combustor.solve(pressures, OF_ratios, zero_mass_flux), std::invalid_argument);

    // contraction_ratio <= 1 is rejected.
    CombustorOptions invalid_contraction = options;
    invalid_contraction.type = CombustorType::FINITE_CONTRACTION_RATIO;
    invalid_contraction.contraction_ratio = 1.0;
    EXPECT_THROW(stream_combustor.solve(pressures, OF_ratios, invalid_contraction), std::invalid_argument);

    // ISOCHORIC combined with a finite-area type is rejected.
    CombustorOptions isochoric_finite = finite_area;
    isochoric_finite.process = CombustionProcess::ISOCHORIC;
    EXPECT_THROW(stream_combustor.solve(pressures, OF_ratios, isochoric_finite), std::invalid_argument);

    // The product species carry no carbon, so a hydrocarbon fuel cannot be burnt in this phase.
    Gas methane = Gas::create_from_species(reactant_file(), "reactants", {"CH4"});
    methane.set_state_TPX(300.0, pressures(0), "CH4:1");
    Combustor carbon_combustor(products, methane, oxidizer);
    EXPECT_THROW(carbon_combustor.solve(pressures, OF_ratios, options), std::invalid_argument);

    // The product-species constructors do not accept the reactant-stream solve overload.
    EXPECT_THROW(combustor->solve(pressures, OF_ratios, options), NotImplementedError);
}

// RP-1311 example 8: H2(L)/O2(L) at O/F 5.55157 and 53.3172 bar, on the NASA9 data CEA itself uses.
class CryogenicRocketTests : public ::testing::Test {
protected:
    static constexpr double H2_BOILING_POINT = 20.27;      // K
    static constexpr double O2_BOILING_POINT = 90.17;      // K
    static constexpr double OF_RATIO = 5.55157;
    static constexpr double CHAMBER_PRESSURE = 53.3172e5;   // Pa
    static constexpr double CEA_CHAMBER_TEMPERATURE = 3383.84;  // K
    static constexpr double CEA_CHAMBER_MOLECULAR_WEIGHT = 12.7157;  // kg/kmol

    static Gas make_products() {
        return Gas::create_from_elements(nasa9_gas_file(), "gas", {"H", "O"});
    }

    static Gas make_fuel() {
        Gas fuel = Gas::create_from_species(reactant_file(), "reactants", {"H2(L)"});
        fuel.set_state_TPX(H2_BOILING_POINT, CHAMBER_PRESSURE, "H2(L):1");
        return fuel;
    }

    static Gas make_oxidizer() {
        Gas oxidizer = Gas::create_from_species(reactant_file(), "reactants", {"O2(L)"});
        oxidizer.set_state_TPX(O2_BOILING_POINT, CHAMBER_PRESSURE, "O2(L):1");
        return oxidizer;
    }

    CombustorOptions options{CombustorType::INFINITE_AREA, MixtureRatioType::OF_RATIO, 0.0, 0.0};
};

TEST_F(CryogenicRocketTests, ChamberMatchesCEAExampleEight) {
    Gas products = make_products();
    Combustor combustor(products, make_fuel(), make_oxidizer());

    Eigen::ArrayXd pressures(1);
    pressures << CHAMBER_PRESSURE;
    Eigen::ArrayXd mixture_ratios(1);
    mixture_ratios << OF_RATIO;

    ThermoArray states = combustor.solve(pressures, mixture_ratios, options);
    ASSERT_EQ(states.size(), 1);

    products.restore_state(states.get_state(0));
    EXPECT_NEAR(products.temperature(), CEA_CHAMBER_TEMPERATURE, 0.25);  // observed 0.09 K
    EXPECT_NEAR(products.molecular_weight(), CEA_CHAMBER_MOLECULAR_WEIGHT, 0.002);  // observed 5.6e-4
    EXPECT_NEAR(products.pressure(), CHAMBER_PRESSURE, 1e-6 * CHAMBER_PRESSURE);
}

TEST_F(CryogenicRocketTests, RocketSolverReproducesTheChamber) {
    RocketPoint point;
    point.fuel = PhaseSpecification(H2_BOILING_POINT, CHAMBER_PRESSURE, Composition{{"H2(L)", 1.0}});
    point.oxidizer = PhaseSpecification(O2_BOILING_POINT, CHAMBER_PRESSURE, Composition{{"O2(L)", 1.0}});
    point.mixture_ratio = OF_RATIO;
    point.chamber_pressure = CHAMBER_PRESSURE;
    point.combustor = options;
    point.nozzle.chemistry = GasChemistry::EQUILIBRIUM;
    point.expansion_ratios = {5.0};

    RocketSolver solver(make_products(), reactant_file());
    const RocketResult result = solver.solve(point);

    ASSERT_TRUE(result.status.ok()) << result.status.describe();
    EXPECT_NEAR(result.chamber.thermo.temperature, CEA_CHAMBER_TEMPERATURE, 0.25);
    EXPECT_NEAR(result.chamber.thermo.molecular_weight, CEA_CHAMBER_MOLECULAR_WEIGHT, 0.002);
    EXPECT_NEAR(result.chamber.thermo.pressure, CHAMBER_PRESSURE, 1e-6 * CHAMBER_PRESSURE);
}
