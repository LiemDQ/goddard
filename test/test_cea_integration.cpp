#include "gtest/gtest.h"
#include "cea_loader.hpp"
#include "cea_comparer.hpp"

#include "goddard/combustor.hpp"
#include "goddard/mixture_ratio.hpp"
#include "goddard/problem.hpp"
#include "goddard/utils.hpp"
#include "goddard/thermo.hpp"
#include "goddard/config.h"

#include <algorithm>
#include <memory>
#include <iostream>
#include <string>
#include <unordered_set>

class CEAIntegrationTests : public ::testing::Test {
protected:
    CEAIntegrationTests() {
        // Initialize with H2O2 mechanism (use full path to project's data file)
        std::string h2o2_path = std::string(DATA_DIR) + "/h2o2.yaml";
        fuel = Cantera::newSolution(h2o2_path, "ohmech");
        oxidizer = Cantera::newSolution(h2o2_path, "ohmech");
        products = Cantera::newSolution(h2o2_path, "ohmech");

        // Initialize CEA data loader
        cea_loader = std::make_unique<CEADataLoader>();

        // Construct path to CEA data directory
        cea_data_dir = std::string(DATA_DIR) + "/cea_results";
    }

    // Helper to create ChemicalParameters from CEA conditions
    Goddard::ChemicalParameters createChemParamsFromCEA(const CEAConditions& conditions);

    // Helper to create RocketCaseParameters
    Goddard::RocketCaseParameters createCaseParams(
        const std::string& name,
        const CEAConditions& conditions,
        Goddard::NozzleChemistryType chemistry);

    std::shared_ptr<Cantera::Solution> fuel, oxidizer, products;
    std::unique_ptr<CEADataLoader> cea_loader;
    std::string cea_data_dir;
};

Goddard::ChemicalParameters CEAIntegrationTests::createChemParamsFromCEA(
    const CEAConditions& conditions) {

    Goddard::ChemicalParameters params;
    params.thermo_file = std::string(DATA_DIR) + "/h2o2.yaml";

    // Species set for H2/O2 combustion (all species from h2o2.yaml)
    // Include AR and N2 to match the phase definition in h2o2.yaml
    params.species = {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"};

    // Set fuel state (H2)
    // Note: CEA conditions may specify cryogenic temps (liquid H2/O2) which
    // aren't supported by ideal gas model. Use standard temp (298.15K) instead.
    // The CEA fuel_energy/oxidizer_energy fields contain the actual enthalpy.
    auto fuel_thermo = fuel->thermo();
    double fuel_temp = conditions.fuel_temp;
    params.cantera_fuel_state = Goddard::ThermodynamicState(fuel_temp, 101325.0, "H2:1");
    // double pressure_Pa = conditions.pressure_psia * 6894.76;
    
    // double fuel_MW = fuel_thermo->meanMolecularWeight();
    // double fuel_energy = conditions.fuel_energy * 1000.0 / fuel_MW;
    // fuel_thermo->setState_HP(fuel_energy, pressure_Pa);

    // Set oxidizer state (O2)
    auto ox_thermo = oxidizer->thermo();
    double ox_temp = conditions.oxidizer_temp;
    params.cantera_oxidizer_state = Goddard::ThermodynamicState(ox_temp, 101325.0, "O2:1");
    // double ox_MW = ox_thermo->meanMolecularWeight();
    // double ox_energy = conditions.oxidizer_energy * 1000.0 / ox_MW;
    // ox_thermo->setState_HP(ox_energy, pressure_Pa);

    // O/F ratio from CEA conditions
    params.OF_ratios = {conditions.of_ratio};

    return params;
}

Goddard::RocketCaseParameters CEAIntegrationTests::createCaseParams(
    const std::string& name,
    const CEAConditions& conditions,
    Goddard::NozzleChemistryType chemistry) {

    Goddard::RocketCaseParameters case_params;
    case_params.name = name;
    case_params.problem_type = "rocket";

    // Combustor: infinite area at CEA pressure (convert psia to Pa)
    case_params.combustor_options.type = Goddard::CombustorType::INFINITE_AREA;
    case_params.combustor_options.pressures = {conditions.pressure_psia * 6894.76};

    // Nozzle configuration
    case_params.nozzle_options.chemistry = chemistry;
    case_params.nozzle_options.expansion_type = Goddard::ExpansionType::SUPERSONIC_AREA_RATIO;

    // Use area ratios from CEA if available, otherwise use default
    if (!conditions.area_ratios.empty()) {
        case_params.nozzle_options.expansion_ratios = conditions.area_ratios;
    } else {
        // Default expansion ratios matching h2.json (15.0, 35.0)
        case_params.nozzle_options.expansion_ratios = {15.0, 35.0};
    }

    return case_params;
}

using namespace Goddard;

TEST_F(CEAIntegrationTests, LoadCEADataFiles) {
    // Test loading all CEA data files
    auto all_cea_data = cea_loader->load_from_directory(cea_data_dir);
    
    EXPECT_GT(all_cea_data.size(), 0) << "No CEA data files found in " << cea_data_dir;
    
    for (const auto& [filename, cea_result] : all_cea_data) {
        EXPECT_NE(cea_result, nullptr) << "Failed to load " << filename;
        
        if (cea_result) {
            std::cout << "Loaded CEA case: " << filename 
                     << " (" << cea_result->conditions.fuel 
                     << " + " << cea_result->conditions.oxidizer
                     << ", O/F=" << cea_result->conditions.of_ratio << ")" << std::endl;
            
            // Verify basic structure
            EXPECT_FALSE(cea_result->conditions.fuel.empty());
            EXPECT_FALSE(cea_result->conditions.oxidizer.empty());
            EXPECT_GT(cea_result->conditions.of_ratio, 0);
            EXPECT_GT(cea_result->equilibrium_states.size(), 0);
        }
    }
}


// Simple demonstration test showing the workflow
TEST_F(CEAIntegrationTests, DemonstrateWorkflow) {
    std::cout << "\n=== CEA Integration Test Workflow Demo ===" << std::endl;

    // Step 1: Load CEA data
    std::cout << "1. Loading CEA data files..." << std::endl;
    auto all_cea_data = cea_loader->load_from_directory(cea_data_dir);
    std::cout << "   Loaded " << all_cea_data.size() << " CEA cases" << std::endl;

    // Step 2: Pick a test case
    std::cout << "2. Processing available test cases..." << std::endl;
    for (const auto& [filename, cea_result] : all_cea_data) {
        std::cout << "   - " << filename << ": "
                 << cea_result->conditions.fuel << " + " << cea_result->conditions.oxidizer
                 << " (O/F=" << cea_result->conditions.of_ratio << ")" << std::endl;

        // Step 3: Extract test conditions
        const auto* chamber = cea_result->find_chamber_state();
        const auto* throat = cea_result->find_throat_state();
        auto exits = cea_result->find_exit_states();

        std::cout << "     Chamber: " << (chamber ? "Y" : "N")
                 << ", Throat: " << (throat ? "Y" : "N")
                 << ", Exits: " << exits.size() << std::endl;

        // Step 4: This is where you'd run Goddard simulation and compare
        // (Implementation depends on your specific Goddard API)
    }

    std::cout << "3. Workflow complete - ready for detailed comparisons!" << std::endl;

    // This test always passes - it's just demonstrating the workflow
    SUCCEED();
}


// ============================================================================
// RocketProblem Integration Tests
// ============================================================================
//
// KNOWN ISSUES identified during test development:
// 1. problem.cpp line 68-69: M_fuel incorrectly uses cantera_oxidizer_state
//    (should be cantera_fuel_state)
// 2. select_species in speciate.cpp filters species from root_node but doesn't
//    update phase species lists, causing "species not found" errors
// 3. State vectors saved from one Solution may be incompatible when restored
//    to a Solution created via select_species due to different state sizes/ordering
//
// These tests are currently skipped pending fixes to the core simulation code.
// ============================================================================

TEST_F(CEAIntegrationTests, RocketProblemChamberMatchesCEA) {
    // SKIP: RocketProblem has known issues - see comments above
    // TODO: Re-enable once problem.cpp and speciate.cpp issues are fixed
    // GTEST_SKIP() << "RocketProblem has known issues with state handling and species selection. "
    //              << "See comments in test file for details.";

    // Load H2 CEA reference data
    std::string h2_json_path = cea_data_dir + "/h2gas.json";
    auto cea_result = cea_loader->load_from_file(h2_json_path);
    ASSERT_NE(cea_result, nullptr) << "Failed to load CEA data from " << h2_json_path;

    // Get CEA chamber state
    const auto* cea_chamber = cea_result->find_chamber_state();
    ASSERT_NE(cea_chamber, nullptr) << "No chamber state found in CEA data";

    // Skip if CEA data incomplete
    if (cea_chamber->temperature_k == 0 || cea_chamber->pressure_bar == 0) {
        GTEST_SKIP() << "CEA chamber data incomplete";
    }

    std::cout << "\n=== RocketProblem vs CEA Chamber Comparison ===" << std::endl;
    std::cout << "CEA Chamber: T=" << cea_chamber->temperature_k << "K, P="
              << cea_chamber->pressure_bar << "bar, MW=" << cea_chamber->molecular_weight
              << "g/mol" << std::endl;

    // Create RocketProblem from CEA conditions
    auto chem_params = createChemParamsFromCEA(cea_result->conditions);
    auto case_params = createCaseParams("H2_O2_equilibrium",
        cea_result->conditions, NozzleChemistryType::EQUILIBRIUM);
    
    // Create and solve RocketProblem
    RocketProblem problem(chem_params, {case_params}, "ohmech");
    auto results = problem.solve();

    // Extract chamber (inlet) state - index 0 is the first state
    auto thermo_states = results.extract_thermo_info("H2_O2_equilibrium", 0);
    ASSERT_GE(thermo_states.size(), 1) << "No thermo states extracted";

    const auto& goddard_chamber = thermo_states[0]; // inlet = chamber
    // const auto& goddard_chamber = results.get_chamber_state("H2_O2_equilibrium", 0);

    std::cout << "Goddard Chamber: T=" << goddard_chamber.temperature << "K, P="
              << goddard_chamber.pressure / 1e5 << "bar, MW=" << goddard_chamber.molecular_weight
              << "g/mol" << std::endl;

    // Compare using CEATestUtils
    auto comparisons = CEATestUtils::compare_chamber_states(goddard_chamber, *cea_chamber);
    CEATestUtils::print_comparison_summary(comparisons, "H2+O2 Chamber");

    // Check all comparisons pass
    for (const auto& result : comparisons) {
        EXPECT_TRUE(result.passed) << result.message;
    }
}

TEST_F(CEAIntegrationTests, RocketProblemThroatMatchesCEA) {
    // Previously skipped due to state handling issues - now fixed

    // Load H2 CEA reference data
    std::string h2_json_path = cea_data_dir + "/h2gas.json";
    auto cea_result = cea_loader->load_from_file(h2_json_path);
    ASSERT_NE(cea_result, nullptr) << "Failed to load CEA data from " << h2_json_path;

    // Get CEA throat state
    const auto* cea_throat = cea_result->find_throat_state();
    ASSERT_NE(cea_throat, nullptr) << "No throat state found in CEA data";

    // Skip if CEA data incomplete
    if (cea_throat->temperature_k == 0) {
        GTEST_SKIP() << "CEA throat data incomplete";
    }

    std::cout << "\n=== RocketProblem vs CEA Throat Comparison ===" << std::endl;
    std::cout << "CEA Throat: T=" << cea_throat->temperature_k << "K, P="
              << cea_throat->pressure_bar << "bar, Mach=" << cea_throat->mach_number
              << std::endl;

    // Create RocketProblem from CEA conditions
    auto chem_params = createChemParamsFromCEA(cea_result->conditions);
    auto case_params = createCaseParams("H2_O2_equilibrium",
        cea_result->conditions, NozzleChemistryType::EQUILIBRIUM);

    // Create and solve RocketProblem
    RocketProblem problem(chem_params, {case_params}, "ohmech");
    auto results = problem.solve();

    // Extract throat state - index 1 is the throat
    auto thermo_states = results.extract_thermo_info("H2_O2_equilibrium", 0);
    ASSERT_GE(thermo_states.size(), 2) << "Not enough thermo states (need throat)";

    const auto& goddard_throat = thermo_states[1]; // throat

    std::cout << "Goddard Throat: T=" << goddard_throat.temperature << "K, P="
              << goddard_throat.pressure / 1e5 << "bar" << std::endl;

    // Compare using CEATestUtils
    auto comparisons = CEATestUtils::compare_chamber_states(goddard_throat, *cea_throat);
    CEATestUtils::print_comparison_summary(comparisons, "H2+O2 Throat");

    // Check all comparisons pass
    for (const auto& result : comparisons) {
        EXPECT_TRUE(result.passed) << result.message;
    }
}

TEST_F(CEAIntegrationTests, RocketProblemExitMatchesCEA) {
    // Previously skipped due to state handling issues - now fixed

    // Load H2 CEA reference data
    std::string h2_json_path = cea_data_dir + "/h2gas.json";
    auto cea_result = cea_loader->load_from_file(h2_json_path);
    ASSERT_NE(cea_result, nullptr) << "Failed to load CEA data from " << h2_json_path;

    // Get CEA exit states
    auto cea_exits = cea_result->find_exit_states();
    ASSERT_GT(cea_exits.size(), 0) << "No exit states found in CEA data";

    std::cout << "\n=== RocketProblem vs CEA Exit Comparison ===" << std::endl;
    std::cout << "CEA has " << cea_exits.size() << " exit state(s)" << std::endl;

    // Create RocketProblem from CEA conditions
    auto chem_params = createChemParamsFromCEA(cea_result->conditions);
    auto case_params = createCaseParams("H2_O2_equilibrium",
        cea_result->conditions, NozzleChemistryType::EQUILIBRIUM);

    // Create and solve RocketProblem
    RocketProblem problem(chem_params, {case_params}, "ohmech");
    auto results = problem.solve();

    // Extract all states
    auto thermo_states = results.extract_thermo_info("H2_O2_equilibrium", 0);
    // States: [0]=inlet/chamber, [1]=throat, [2+]=exits
    size_t num_exits = thermo_states.size() > 2 ? thermo_states.size() - 2 : 0;
    std::cout << "Goddard has " << num_exits << " exit state(s)" << std::endl;

    // Compare each exit state
    for (size_t i = 0; i < std::min(num_exits, cea_exits.size()); i++) {
        const auto& goddard_exit = thermo_states[i + 2];
        const auto* cea_exit = cea_exits[i];

        if (cea_exit->temperature_k == 0) continue;

        std::cout << "\nExit " << i << ":" << std::endl;
        std::cout << "  CEA: T=" << cea_exit->temperature_k << "K, P="
                  << cea_exit->pressure_bar << "bar" << std::endl;
        std::cout << "  Goddard: T=" << goddard_exit.temperature << "K, P="
                  << goddard_exit.pressure / 1e5 << "bar" << std::endl;

        auto comparisons = CEATestUtils::compare_chamber_states(goddard_exit, *cea_exit);
        CEATestUtils::print_comparison_summary(comparisons, "Exit " + std::to_string(i));

        for (const auto& result : comparisons) {
            EXPECT_TRUE(result.passed) << "Exit " << i << ": " << result.message;
        }
    }
}