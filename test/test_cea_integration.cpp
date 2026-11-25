#include "gtest/gtest.h"
#include "cea_loader.hpp"
#include "cea_comparer.hpp"

#include "goddard/combustor.hpp"
#include "goddard/rocket.hpp"
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

    // Solve the rocket point of the CEA conditions with the given nozzle chemistry.
    Goddard::RocketResult solve_point(const CEAConditions& conditions,
                                      Goddard::GasChemistry chemistry);

    std::shared_ptr<Cantera::Solution> fuel, oxidizer, products;
    std::unique_ptr<CEADataLoader> cea_loader;
    std::string cea_data_dir;
};

Goddard::RocketResult CEAIntegrationTests::solve_point(
    const CEAConditions& conditions, Goddard::GasChemistry chemistry) {

    // Every species of h2o2.yaml's ohmech phase, including AR and N2.
    const Goddard::Gas product_gas(std::string(DATA_DIR) + "/h2o2.yaml", "ohmech");

    // CEA conditions may specify cryogenic temperatures, which the ideal-gas reactants take as
    // given; the CEA fuel_energy/oxidizer_energy fields hold the actual enthalpies.
    Goddard::RocketPoint point;
    point.fuel = Goddard::PhaseSpecification(conditions.fuel_temp, 101325.0, "H2:1");
    point.oxidizer = Goddard::PhaseSpecification(conditions.oxidizer_temp, 101325.0, "O2:1");
    point.mixture_ratio = conditions.of_ratio;
    // Infinite-area combustor at the CEA pressure (psia to Pa)
    point.chamber_pressure = conditions.pressure_psia * 6894.76;
    point.nozzle.chemistry = chemistry;
    point.expansion_type = Goddard::ExpansionType::SUPERSONIC_AREA_RATIO;
    // The CEA area ratios if given, otherwise those of h2.json
    point.expansion_ratios = conditions.area_ratios.empty()
        ? std::vector<double>{15.0, 35.0} : conditions.area_ratios;

    Goddard::RocketSolver solver(product_gas);
    return solver.solve(point);
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
// RocketSolver Integration Tests
// ============================================================================

TEST_F(CEAIntegrationTests, RocketSolverChamberMatchesCEA) {
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

    std::cout << "\n=== RocketSolver vs CEA Chamber Comparison ===" << std::endl;
    std::cout << "CEA Chamber: T=" << cea_chamber->temperature_k << "K, P="
              << cea_chamber->pressure_bar << "bar, MW=" << cea_chamber->molecular_weight
              << "g/mol" << std::endl;

    // Solve the rocket point of the CEA conditions
    const RocketResult rocket = solve_point(cea_result->conditions, GasChemistry::EQUILIBRIUM);

    // Extract chamber state
    const auto& goddard_chamber = rocket.chamber.thermo;

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

TEST_F(CEAIntegrationTests, RocketSolverThroatMatchesCEA) {
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

    std::cout << "\n=== RocketSolver vs CEA Throat Comparison ===" << std::endl;
    std::cout << "CEA Throat: T=" << cea_throat->temperature_k << "K, P="
              << cea_throat->pressure_bar << "bar, Mach=" << cea_throat->mach_number
              << std::endl;

    // Solve the rocket point of the CEA conditions
    const RocketResult rocket = solve_point(cea_result->conditions, GasChemistry::EQUILIBRIUM);

    // Extract throat state
    const auto& goddard_throat = rocket.throat.thermo;

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

TEST_F(CEAIntegrationTests, RocketSolverExitMatchesCEA) {
    // Previously skipped due to state handling issues - now fixed

    // Load H2 CEA reference data
    std::string h2_json_path = cea_data_dir + "/h2gas.json";
    auto cea_result = cea_loader->load_from_file(h2_json_path);
    ASSERT_NE(cea_result, nullptr) << "Failed to load CEA data from " << h2_json_path;

    // Get CEA exit states
    auto cea_exits = cea_result->find_exit_states();
    ASSERT_GT(cea_exits.size(), 0) << "No exit states found in CEA data";

    std::cout << "\n=== RocketSolver vs CEA Exit Comparison ===" << std::endl;
    std::cout << "CEA has " << cea_exits.size() << " exit state(s)" << std::endl;

    // Solve the rocket point of the CEA conditions
    const RocketResult rocket = solve_point(cea_result->conditions, GasChemistry::EQUILIBRIUM);

    // Extract exit states
    const auto& goddard_exits = rocket.exits;
    std::cout << "Goddard has " << goddard_exits.size() << " exit state(s)" << std::endl;

    // Compare each exit state
    for (size_t i = 0; i < std::min(goddard_exits.size(), cea_exits.size()); i++) {
        const auto& goddard_exit = goddard_exits[i].thermo;
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


TEST_F(CEAIntegrationTests, FrozenRocketSolverChamberMatchesCEA) {

    // Load H2 CEA reference data
    std::string h2_json_path = cea_data_dir + "/h2gas.json";
    auto cea_result = cea_loader->load_from_file(h2_json_path);
    ASSERT_NE(cea_result, nullptr) << "Failed to load CEA data from " << h2_json_path;

    // Get CEA chamber state
    const auto* cea_chamber = cea_result->find_frozen_chamber_state();
    ASSERT_NE(cea_chamber, nullptr) << "No chamber state found in CEA data";

    // Skip if CEA data incomplete
    if (cea_chamber->temperature_k == 0 || cea_chamber->pressure_bar == 0) {
        GTEST_SKIP() << "CEA chamber data incomplete";
    }

    std::cout << "\n=== RocketSolver vs CEA Chamber Comparison ===" << std::endl;
    std::cout << "CEA Chamber: T=" << cea_chamber->temperature_k << "K, P="
              << cea_chamber->pressure_bar << "bar, MW=" << cea_chamber->molecular_weight
              << "g/mol" << std::endl;

    // Solve the rocket point of the CEA conditions
    const RocketResult rocket = solve_point(cea_result->conditions, GasChemistry::FROZEN);

    // Extract chamber state
    const auto& goddard_chamber = rocket.chamber.thermo;

    std::cout << "Goddard Chamber: T=" << goddard_chamber.temperature << "K, P="
              << goddard_chamber.pressure / 1e5 << "bar, MW=" << goddard_chamber.molecular_weight
              << "g/mol" << std::endl;

    // Compare using CEATestUtils
    auto comparisons = CEATestUtils::compare_chamber_states(goddard_chamber, *cea_chamber);
    CEATestUtils::print_comparison_summary(comparisons, "H2+O2 Frozen Chamber");

    // Check all comparisons pass
    for (const auto& result : comparisons) {
        EXPECT_TRUE(result.passed) << result.message;
    }
}


TEST_F(CEAIntegrationTests, FrozenRocketSolverThroatMatchesCEA) {
    // Previously skipped due to state handling issues - now fixed

    // Load H2 CEA reference data
    std::string h2_json_path = cea_data_dir + "/h2gas.json";
    auto cea_result = cea_loader->load_from_file(h2_json_path);
    ASSERT_NE(cea_result, nullptr) << "Failed to load CEA data from " << h2_json_path;

    // Get CEA throat state
    const auto* cea_throat = cea_result->find_frozen_throat_state();
    ASSERT_NE(cea_throat, nullptr) << "No throat state found in CEA data";

    // Skip if CEA data incomplete
    if (cea_throat->temperature_k == 0) {
        GTEST_SKIP() << "CEA throat data incomplete";
    }

    std::cout << "\n=== RocketSolver vs CEA Throat Comparison ===" << std::endl;
    std::cout << "CEA Throat: T=" << cea_throat->temperature_k << "K, P="
              << cea_throat->pressure_bar << "bar, Mach=" << cea_throat->mach_number
              << std::endl;

    // Solve the rocket point of the CEA conditions
    const RocketResult rocket = solve_point(cea_result->conditions, GasChemistry::FROZEN);

    // Extract throat state
    const auto& goddard_throat = rocket.throat.thermo;

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

TEST_F(CEAIntegrationTests, FrozenRocketSolverExitMatchesCEA) {
    // Load H2 CEA reference data
    std::string h2_json_path = cea_data_dir + "/h2gas.json";
    auto cea_result = cea_loader->load_from_file(h2_json_path);
    ASSERT_NE(cea_result, nullptr) << "Failed to load CEA data from " << h2_json_path;

    // Get CEA exit states
    auto cea_exits = cea_result->find_frozen_exit_states();
    ASSERT_GT(cea_exits.size(), 0) << "No exit states found in CEA data";

    std::cout << "\n=== RocketSolver vs CEA Exit Comparison ===" << std::endl;
    std::cout << "CEA has " << cea_exits.size() << " exit state(s)" << std::endl;

    // Solve the rocket point of the CEA conditions
    const RocketResult rocket = solve_point(cea_result->conditions, GasChemistry::FROZEN);

    // Extract exit states
    const auto& goddard_exits = rocket.exits;
    std::cout << "Goddard has " << goddard_exits.size() << " exit state(s)" << std::endl;

    // Compare each exit state
    for (size_t i = 0; i < std::min(goddard_exits.size(), cea_exits.size()); i++) {
        const auto& goddard_exit = goddard_exits[i].thermo;
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
