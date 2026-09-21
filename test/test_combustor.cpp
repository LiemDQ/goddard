#include "goddard/combustor.hpp"
#include "goddard/gas.hpp"
#include "goddard/numerics.hpp"
#include "goddard/utils.hpp"
#include "goddard/problem.hpp"
#include "goddard/error.hpp"
#include "goddard/config.h"
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

// ---- Reactant streams given as Gas objects ----

namespace {

/** Species names of a phase, for reuse as a `ChemicalParameters::species` set. */
std::unordered_set<std::string> species_set(const Gas& gas) {
    std::vector<std::string> names = gas.species_names();
    return {names.begin(), names.end()};
}

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

    // CombustorType::NONE is not implemented.
    CombustorOptions none_type = options;
    none_type.type = CombustorType::NONE;
    EXPECT_THROW(stream_combustor.solve(pressures, OF_ratios, none_type), NotImplementedError);

    // The product species carry no carbon, so a hydrocarbon fuel cannot be burnt in this phase.
    Gas methane = Gas::create_from_species(reactant_file(), "reactants", {"CH4"});
    methane.set_state_TPX(300.0, pressures(0), "CH4:1");
    Combustor carbon_combustor(products, methane, oxidizer);
    EXPECT_THROW(carbon_combustor.solve(pressures, OF_ratios, options), FmtError);

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

    CombustorOptions options{CombustorType::INFINITE_AREA, MixtureRatioType::OF_RATIO,
        {CHAMBER_PRESSURE}, 0.0, 0.0};
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
    EXPECT_NEAR(products.temperature(), CEA_CHAMBER_TEMPERATURE, 3.0);
    EXPECT_NEAR(products.molecular_weight(), CEA_CHAMBER_MOLECULAR_WEIGHT, 0.01);
    EXPECT_NEAR(products.pressure(), CHAMBER_PRESSURE, 1e-6 * CHAMBER_PRESSURE);
}

TEST_F(CryogenicRocketTests, RocketProblemReproducesTheChamber) {
    ChemicalParameters chem_params;
    chem_params.thermo_file = nasa9_gas_file();
    chem_params.species = species_set(make_products());
    chem_params.reactant_file = reactant_file();
    chem_params.cantera_fuel_state =
        PhaseSpecification(H2_BOILING_POINT, CHAMBER_PRESSURE, Composition{{"H2(L)", 1.0}});
    chem_params.cantera_oxidizer_state =
        PhaseSpecification(O2_BOILING_POINT, CHAMBER_PRESSURE, Composition{{"O2(L)", 1.0}});
    chem_params.mixture_type = MixtureRatioType::OF_RATIO;
    chem_params.OF_ratios = {OF_RATIO};

    RocketCaseParameters case_params;
    case_params.name = "ex8";
    case_params.problem_type = "rocket";
    case_params.combustor_options = options;
    case_params.nozzle_options.chemistry = GasChemistry::EQUILIBRIUM;
    case_params.nozzle_options.expansion_type = ExpansionType::SUPERSONIC_AREA_RATIO;
    case_params.nozzle_options.expansion_ratios = {5.0};

    RocketProblem problem(chem_params, {case_params}, "gas");
    RocketProblemResults results = problem.solve();

    const RocketStation& chamber = results.chamber(0, "ex8");
    EXPECT_NEAR(chamber.thermo.temperature, CEA_CHAMBER_TEMPERATURE, 3.0);
    EXPECT_NEAR(chamber.thermo.molecular_weight, CEA_CHAMBER_MOLECULAR_WEIGHT, 0.01);
    EXPECT_NEAR(chamber.thermo.pressure, CHAMBER_PRESSURE, 1e-6 * CHAMBER_PRESSURE);
}

// ---- RocketProblemResults for a finite-area combustor ----
//
// Nozzle::solve_finite_area_chamber() and Nozzle::solve_stations() are throwing stubs in this
// worktree (implemented on the nozzle-numerics side), so these tests build a
// RocketProblemCaseResult by hand: real infinite-area combustion solves stand in for the
// injector face and the "inf" stagnation state, and a real Nozzle solve from the stagnation
// state provides the combustion-end, throat and exit stations.

TEST(FiniteAreaCombustorResults, StationsAccessorsAndReportUseStagnationState) {
    const double injector_pressure = 50.0 * Cantera::OneBar;
    // P_inf < P_inj, as the finite-area chamber momentum balance requires.
    const double stagnation_pressure = 47.8 * Cantera::OneBar;
    const double of_ratio = 6.0;
    const double reactant_temperature = 300.0;

    Eigen::ArrayXd of_ratios(1);
    of_ratios << of_ratio;

    Gas combustor_gas(Cantera::newSolution("h2o2.yaml", "ohmech"));
    Combustor infinite_combustor(combustor_gas, Composition{{"H2", 1.0}}, Composition{{"O2", 1.0}});
    CombustorOptions infinite_area_options{CombustorType::INFINITE_AREA, MixtureRatioType::OF_RATIO,
        {}, 0.0, 0.0};

    Eigen::ArrayXd injector_pressures(1);
    injector_pressures << injector_pressure;
    ThermoArray injector_states = infinite_combustor.solve(
        reactant_temperature, reactant_temperature, injector_pressures, of_ratios, infinite_area_options);

    Eigen::ArrayXd stagnation_pressures(1);
    stagnation_pressures << stagnation_pressure;
    ThermoArray stagnation_states = infinite_combustor.solve(
        reactant_temperature, reactant_temperature, stagnation_pressures, of_ratios, infinite_area_options);

    std::vector<double> stagnation_state = stagnation_states.get_state(0);

    Gas nozzle_gas(Cantera::newSolution("h2o2.yaml", "ohmech"));
    NozzleOptions nozzle_opts;
    nozzle_opts.chemistry = GasChemistry::EQUILIBRIUM;
    nozzle_opts.expansion_type = ExpansionType::SUPERSONIC_AREA_RATIO;
    nozzle_opts.expansion_ratios = {5.0};

    Nozzle nozzle(nozzle_gas, stagnation_state, nozzle_opts);
    NozzleResults exits = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, nozzle_opts.expansion_ratios);
    NozzleResults comb_end_result = nozzle.solve(ExpansionType::SUBSONIC_AREA_RATIO, 2.0);

    Gas station_gas(Cantera::newSolution("h2o2.yaml", "ohmech"), GasChemistry::EQUILIBRIUM);
    auto equilibrium_station = [&station_gas](const std::vector<double>& state) {
        station_gas.restore_state(state);
        const ExpansionProperties props = station_gas.expansion_properties();
        return NozzleStation{true, props.gamma_s, props.dlogV_dlogP_T, props.dlogV_dlogT_P,
            state, props.pinned_transition};
    };

    FiniteAreaChamber fac{
        equilibrium_station(injector_states.get_state(0)),
        equilibrium_station(stagnation_state),
        comb_end_result.expansions[0],
        exits.throat,
        injector_pressure,
        stagnation_pressure,
        /*contraction_ratio=*/2.0,
        /*mass_flux=*/500.0,
        /*iterations=*/1
    };

    RocketProblemCaseResult case_result{
        "rocket",
        injector_states,
        GasChemistry::EQUILIBRIUM,
        {NozzleResults{fac.injector, fac.throat, exits.expansions}},
        {of_ratio},
        {injector_pressure},
        {5.0},
        ExpansionType::SUPERSONIC_AREA_RATIO,
        CombustionProcess::ISOBARIC,
        CombustorType::FINITE_CONTRACTION_RATIO,
        {fac}
    };

    std::unordered_map<std::string, RocketProblemCaseResult> case_results;
    case_results.emplace("fac_case", std::move(case_result));

    Gas results_gas(Cantera::newSolution("h2o2.yaml", "ohmech"));
    RocketProblemResults results(std::move(case_results), results_gas);

    // Station order and types for the single operating point: chamber (injector), stagnation,
    // combustion end, throat, one exit.
    const std::vector<RocketStation>& stations = results.stations();
    ASSERT_EQ(stations.size(), 5u);
    EXPECT_EQ(stations[0].type, StationType::CHAMBER);
    EXPECT_EQ(stations[1].type, StationType::STAGNATION);
    EXPECT_EQ(stations[2].type, StationType::COMBUSTION_END);
    EXPECT_EQ(stations[3].type, StationType::THROAT);
    EXPECT_EQ(stations[4].type, StationType::EXIT);

    // Accessors.
    const RocketStation& chamber_station = results.chamber(0, "fac_case");
    EXPECT_EQ(chamber_station.type, StationType::CHAMBER);
    EXPECT_NEAR(chamber_station.thermo.pressure, injector_pressure, 1e-6 * injector_pressure);

    const RocketStation& stagnation_station = results.stagnation(0, "fac_case");
    EXPECT_EQ(stagnation_station.type, StationType::STAGNATION);
    EXPECT_NEAR(stagnation_station.thermo.pressure, stagnation_pressure, 1e-6 * stagnation_pressure);

    const RocketStation& comb_end_station = results.combustion_end(0, "fac_case");
    EXPECT_EQ(comb_end_station.type, StationType::COMBUSTION_END);
    EXPECT_DOUBLE_EQ(comb_end_station.area_ratio, 2.0);

    // performance() must use the stagnation station, not the chamber/injector one: the two
    // differ here, since P_inf != P_inj.
    const RocketStation& throat_station = results.throat(0, "fac_case");
    std::vector<RocketStation> exit_stations = results.exits(0, "fac_case");
    ASSERT_EQ(exit_stations.size(), 1u);

    RocketPerformance expected = RocketProblemResults::calculate_performance(
        stagnation_station.thermo, throat_station.thermo, exit_stations[0].thermo);
    RocketPerformance actual = results.performance(0, 0, "fac_case");
    EXPECT_DOUBLE_EQ(actual.cstar, expected.cstar);
    EXPECT_DOUBLE_EQ(actual.CF, expected.CF);
    EXPECT_DOUBLE_EQ(actual.isp, expected.isp);

    RocketPerformance chamber_based = RocketProblemResults::calculate_performance(
        chamber_station.thermo, throat_station.thermo, exit_stations[0].thermo);
    EXPECT_NE(actual.cstar, chamber_based.cstar);

    // Report layout: CEA's finite-area combustor page layout (see data/cea_results/h2gas_fac.output).
    std::string report = results.report("fac_case");
    EXPECT_NE(report.find("FINITE AREA COMBUSTOR"), std::string::npos);
    EXPECT_NE(report.find("INJECTOR"), std::string::npos);
    EXPECT_NE(report.find("COMB END"), std::string::npos);
    EXPECT_NE(report.find("Pinj/P"), std::string::npos);
    EXPECT_NE(report.find("Pinf/P"), std::string::npos);
}

TEST(FiniteAreaCombustorResults, InfiniteAreaHasNoCombustionEndStation) {
    const double reactant_temperature = 300.0;
    const double pressure = 50.0 * Cantera::OneBar;
    const double of_ratio = 6.0;

    ChemicalParameters chem_params;
    chem_params.thermo_file = std::string(DATA_DIR) + "/h2o2.yaml";
    chem_params.species = {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"};
    chem_params.cantera_fuel_state = PhaseSpecification(reactant_temperature, pressure, "H2:1");
    chem_params.cantera_oxidizer_state = PhaseSpecification(reactant_temperature, pressure, "O2:1");
    chem_params.mixture_type = MixtureRatioType::OF_RATIO;
    chem_params.OF_ratios = {of_ratio};

    RocketCaseParameters case_params;
    case_params.name = "infinite";
    case_params.problem_type = "rocket";
    case_params.combustor_options.type = CombustorType::INFINITE_AREA;
    case_params.combustor_options.pressures = {pressure};
    case_params.nozzle_options.chemistry = GasChemistry::EQUILIBRIUM;
    case_params.nozzle_options.expansion_type = ExpansionType::SUPERSONIC_AREA_RATIO;
    case_params.nozzle_options.expansion_ratios = {5.0};

    RocketProblem problem(chem_params, {case_params}, "ohmech");
    RocketProblemResults results = problem.solve();

    EXPECT_THROW(results.combustion_end(0, "infinite"), std::runtime_error);

    const RocketStation& stag = results.stagnation(0, "infinite");
    const RocketStation& chamber = results.chamber(0, "infinite");
    EXPECT_EQ(stag.type, StationType::CHAMBER);
    EXPECT_NEAR(stag.thermo.pressure, chamber.thermo.pressure, 1e-9 * chamber.thermo.pressure);
}

