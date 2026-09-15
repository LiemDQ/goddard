#include "goddard/combustor.hpp"
#include "goddard/gas.hpp"
#include "goddard/numerics.hpp"
#include "goddard/utils.hpp"
#include "goddard/problem.hpp"
#include "goddard/config.h"
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

TEST(RocketProblemIsochoric, ChamberIsConstantVolumeState) {
    const double fuel_temperature = 300.0;
    const double oxidizer_temperature = 300.0;
    const double initial_pressure = 10.0 * Cantera::OneBar;
    const double of_ratio = 6.0;

    ChemicalParameters chem_params;
    chem_params.thermo_file = std::string(DATA_DIR) + "/h2o2.yaml";
    chem_params.species = {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"};
    chem_params.cantera_fuel_state = PhaseSpecification(fuel_temperature, initial_pressure, "H2:1");
    chem_params.cantera_oxidizer_state = PhaseSpecification(oxidizer_temperature, initial_pressure, "O2:1");
    chem_params.mixture_type = MixtureRatioType::OF_RATIO;
    chem_params.OF_ratios = {of_ratio};

    RocketCaseParameters case_params;
    case_params.name = "isochoric";
    case_params.problem_type = "rocket";
    case_params.combustor_options.type = CombustorType::INFINITE_AREA;
    case_params.combustor_options.process = CombustionProcess::ISOCHORIC;
    case_params.combustor_options.pressures = {initial_pressure};
    case_params.nozzle_options.chemistry = GasChemistry::EQUILIBRIUM;
    case_params.nozzle_options.expansion_type = ExpansionType::SUPERSONIC_AREA_RATIO;
    case_params.nozzle_options.expansion_ratios = {10.0};

    RocketProblem problem(chem_params, {case_params}, "ohmech");
    RocketProblemResults results = problem.solve();

    // Reference: the same constant-volume combustion computed directly with the combustor.
    Gas gas(Cantera::newSolution("h2o2.yaml", "ohmech"));
    Combustor combustor(gas, Composition{{"H2", 1.0}}, Composition{{"O2", 1.0}});
    Eigen::ArrayXd pressures(1);
    pressures << initial_pressure;
    Eigen::ArrayXd mixture_ratios(1);
    mixture_ratios << of_ratio;
    ThermoArray reference = combustor.solve(fuel_temperature, oxidizer_temperature,
        pressures, mixture_ratios, case_params.combustor_options);
    auto thermo = gas.thermo();
    thermo->restoreState(reference.get_state(0));

    const RocketStation& chamber = results.chamber(0, "isochoric");
    EXPECT_NEAR(chamber.thermo.pressure, thermo->pressure(), 1e-6 * thermo->pressure());
    EXPECT_NEAR(chamber.thermo.temperature, thermo->temperature(), 1e-6 * thermo->temperature());
    EXPECT_GT(chamber.thermo.pressure, 5.0 * initial_pressure);

    const RocketStation& throat = results.throat(0, "isochoric");
    EXPECT_TRUE(throat.converged);
    EXPECT_LT(throat.thermo.pressure, chamber.thermo.pressure);

    std::vector<RocketStation> exits = results.exits(0, "isochoric");
    ASSERT_EQ(exits.size(), 1u);
    EXPECT_LT(exits[0].thermo.pressure, throat.thermo.pressure);

    std::string report = results.report("isochoric");
    EXPECT_NE(report.find("CONSTANT-VOLUME COMBUSTOR"), std::string::npos);
    EXPECT_NE(report.find("Pinitial"), std::string::npos);
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

TEST(RocketProblemIndexing, StationsMatchMixtureRatioAndPressure) {
    const double reactant_temperature = 300.0;
    const std::vector<double> pressures = {20.0 * Cantera::OneBar, 50.0 * Cantera::OneBar};
    const std::vector<double> of_ratios = {4.0, 8.0};

    ChemicalParameters chem_params;
    chem_params.thermo_file = std::string(DATA_DIR) + "/h2o2.yaml";
    chem_params.species = {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"};
    chem_params.cantera_fuel_state = PhaseSpecification(reactant_temperature, pressures[0], "H2:1");
    chem_params.cantera_oxidizer_state = PhaseSpecification(reactant_temperature, pressures[0], "O2:1");
    chem_params.mixture_type = MixtureRatioType::OF_RATIO;
    chem_params.OF_ratios = of_ratios;

    RocketCaseParameters case_params;
    case_params.name = "sweep";
    case_params.problem_type = "rocket";
    case_params.combustor_options.pressures = pressures;
    case_params.nozzle_options.chemistry = GasChemistry::EQUILIBRIUM;
    case_params.nozzle_options.expansion_type = ExpansionType::SUPERSONIC_AREA_RATIO;
    case_params.nozzle_options.expansion_ratios = {5.0};

    RocketProblem problem(chem_params, {case_params}, "ohmech");
    RocketProblemResults results = problem.solve();

    std::vector<RocketStation> chambers = results.stations_of_type(StationType::CHAMBER, "sweep");
    ASSERT_EQ(chambers.size(), of_ratios.size() * pressures.size());

    Gas gas(Cantera::newSolution("h2o2.yaml", "ohmech"));
    Combustor combustor(gas, Composition{{"H2", 1.0}}, Composition{{"O2", 1.0}});
    auto thermo = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    for (const RocketStation& chamber : chambers) {
        Eigen::ArrayXd single_pressure(1);
        single_pressure << pressures[chamber.pressure_index];
        Eigen::ArrayXd single_ratio(1);
        single_ratio << of_ratios[chamber.of_index];
        ThermoArray reference = combustor.solve(reactant_temperature, reactant_temperature,
            single_pressure, single_ratio, case_params.combustor_options);
        thermo->restoreState(reference.get_state(0));

        EXPECT_NEAR(chamber.thermo.pressure, thermo->pressure(), 1e-9 * thermo->pressure())
            << "of_index = " << chamber.of_index << ", pressure_index = " << chamber.pressure_index;
        EXPECT_NEAR(chamber.thermo.temperature, thermo->temperature(), 1e-9 * thermo->temperature())
            << "of_index = " << chamber.of_index << ", pressure_index = " << chamber.pressure_index;
    }
}
