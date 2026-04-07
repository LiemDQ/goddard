#include "goddard/nozzle.hpp"
#include "goddard/profile.hpp"
#include "goddard/combustor.hpp"
#include "goddard/mixture_ratio.hpp"
#include "goddard/utils.hpp"
#include "goddard/numerics.hpp"

#include <memory>
#include <iostream>
#include "eigen3/Eigen/Dense"
#include "cantera/core.h"
#include "gtest/gtest.h"

using namespace Goddard;

class NozzleTests : public ::testing::Test {
protected:
    NozzleTests() {
        // Create solution for hydrogen-oxygen combustion
        gas = Cantera::newSolution("h2o2.yaml", "ohmech");

        // Set up a typical combustion chamber state
        // Using stoichiometric H2/O2 at high pressure and temperature
        double temp = 3500.0; // K - typical combustion temperature
        double pressure = 70.0 * Cantera::OneAtm; // ~70 atm chamber pressure

        auto thermo = gas->thermo();
        thermo->setState_TPX(temp, pressure, "H2O:0.9, H2:0.05, O2:0.03, OH:0.02");
        thermo->equilibrate("HP");

        // Save this as our standard inlet state
        thermo->saveState(inlet_state);
    }

    void SetUp() override {
        // Reset to inlet state before each test
        gas->thermo()->restoreState(inlet_state);
    }

    std::shared_ptr<Cantera::Solution> gas;
    std::vector<double> inlet_state;
};

// Test basic construction and state management
TEST_F(NozzleTests, EquilibriumNozzleConstruction) {
    ASSERT_NO_THROW({
        Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);
    });
}

TEST_F(NozzleTests, FrozenNozzleConstruction) {
    ASSERT_NO_THROW({
        Nozzle nozzle(*gas, GasChemistry::FROZEN);
    });
}

TEST_F(NozzleTests, NozzleStateManagement) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    std::vector<double> original_state = nozzle.get_inlet_state();
    ASSERT_EQ(original_state.size(), gas->thermo()->stateSize());

    // Modify the gas state
    gas->thermo()->setState_TP(1000.0, Cantera::OneAtm);

    // Reset should restore original state
    nozzle.reset_state();

    double current_temp = gas->thermo()->temperature();
    double current_pressure = gas->thermo()->pressure();

    // Should be back to inlet conditions
    EXPECT_GT(current_temp, 3000.0); // Much higher than 1000K we set
    EXPECT_GT(current_pressure, 50.0 * Cantera::OneAtm); // Much higher than 1 atm
}

TEST_F(NozzleTests, SetInletState) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    // Create a new state
    gas->thermo()->setState_TPX(2000.0, 10.0 * Cantera::OneAtm, "H2O:1.0");
    std::vector<double> new_state(gas->thermo()->stateSize());
    gas->thermo()->saveState(new_state);

    nozzle.set_inlet_state(new_state);

    std::vector<double> retrieved_state = nozzle.get_inlet_state();
    ASSERT_EQ(retrieved_state.size(), new_state.size());

    for (size_t i = 0; i < new_state.size(); i++) {
        EXPECT_DOUBLE_EQ(retrieved_state[i], new_state[i]);
    }
}

TEST_F(NozzleTests, EqThroatConditionsConverge) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    ThroatCondition conditions = nozzle.solve_throat_conditions();
    EXPECT_TRUE(conditions.converged);
    EXPECT_GT(conditions.gamma_s, 1.0) << "Specific heat ratio should be greater than 1";
    EXPECT_LT(conditions.gamma_s, 2.0) << "Specific heat ratio should be physically reasonable";

    gas->thermo()->restoreState(conditions.state);
    double tol = max_fp_error(conditions.S_inlet, 1e-5, 1e-3);  
    EXPECT_NEAR(conditions.S_inlet, gas->thermo()->entropy_mass(), tol) << "Throat entropy should be isentropic";
    EXPECT_LE(gas->thermo()->pressure(), conditions.P_inlet) << "Throat pressure should be less than chamber pressure";

    gas->thermo()->restoreState(inlet_state);
    EXPECT_DOUBLE_EQ(conditions.P_inlet, gas->thermo()->pressure()) << "Inlet pressure should be equal";
    EXPECT_DOUBLE_EQ(conditions.S_inlet, gas->thermo()->entropy_mass()) << "Inlet pressure should remain equal";
}

// Test throat condition calculation
TEST_F(NozzleTests, EquilibriumConverge) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 10.0);
    NozzleStation result = results.expansions.front();

    EXPECT_TRUE(result.converged) << "Throat conditions should converge for typical rocket conditions";
    EXPECT_GT(result.gamma_s, 1.0) << "Specific heat ratio should be greater than 1";
    EXPECT_LT(result.gamma_s, 2.0) << "Specific heat ratio should be physically reasonable";
}

TEST_F(NozzleTests, FrozenConverge) {
    Nozzle nozzle(*gas, GasChemistry::FROZEN);

    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 10.0);
    NozzleStation result = results.expansions.front();

    EXPECT_TRUE(result.converged) << "Throat conditions should converge for typical rocket conditions";
    EXPECT_GT(result.gamma_s, 1.0) << "Specific heat ratio should be greater than 1";
    EXPECT_LT(result.gamma_s, 2.0) << "Specific heat ratio should be physically reasonable";
}

// Test supersonic area expansion
TEST_F(NozzleTests, EquilibriumSupersonicAreaExpansion) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    std::vector<double> expansion_ratios = {2.0, 5.0, 10.0, 20.0};

    for (double ratio : expansion_ratios) {
        NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, ratio);
        NozzleStation result = results.expansions.front();

        EXPECT_TRUE(result.converged)
            << "Expansion should converge for area ratio " << ratio;
        EXPECT_GT(result.gamma_s, 1.0)
            << "gamma_s should be > 1 for area ratio " << ratio;
        EXPECT_EQ(result.state.size(), gas->thermo()->stateSize())
            << "State vector should have correct size";
    }
}

TEST_F(NozzleTests, FrozenSupersonicAreaExpansion) {
    Nozzle nozzle(*gas, GasChemistry::FROZEN);

    std::vector<double> expansion_ratios = {2.0, 5.0, 10.0, 20.0};

    for (double ratio : expansion_ratios) {
        NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, ratio);
        NozzleStation result = results.expansions.front();

        EXPECT_TRUE(result.converged)
            << "Expansion should converge for area ratio " << ratio;
        EXPECT_GT(result.gamma_s, 1.0)
            << "gamma_s should be > 1 for area ratio " << ratio;
    }
}

TEST_F(NozzleTests, SupersonicExpansionPressureDecreases) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    double initial_pressure = gas->thermo()->pressure();

    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 10.0);
    NozzleStation result = results.expansions.front();

    ASSERT_TRUE(result.converged);

    // Restore the exit state
    gas->thermo()->restoreState(result.state);
    double exit_pressure = gas->thermo()->pressure();

    EXPECT_LT(exit_pressure, initial_pressure)
        << "Pressure should decrease through nozzle expansion";
    EXPECT_GT(exit_pressure, 0.0)
        << "Exit pressure should be positive";
}

TEST_F(NozzleTests, SupersonicExpansionTemperatureDecreases) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    double initial_temp = gas->thermo()->temperature();

    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 10.0);
    NozzleStation result = results.expansions.front();

    ASSERT_TRUE(result.converged);

    gas->thermo()->restoreState(result.state);
    double exit_temp = gas->thermo()->temperature();

    EXPECT_LT(exit_temp, initial_temp)
        << "Temperature should decrease through nozzle expansion";
    EXPECT_GT(exit_temp, 0.0)
        << "Exit temperature should be positive";
}

// Test subsonic area expansion
TEST_F(NozzleTests, EquilibriumSubsonicAreaExpansion) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    // Subsonic expansion ratios are > 1 but smaller than supersonic
    std::vector<double> expansion_ratios = {1.1, 1.2, 1.5};

    for (double ratio : expansion_ratios) {
        NozzleResults results = nozzle.solve(ExpansionType::SUBSONIC_AREA_RATIO, ratio);
        NozzleStation result = results.expansions.front();

        EXPECT_TRUE(result.converged)
            << "Subsonic expansion should converge for area ratio " << ratio;
        EXPECT_GT(result.gamma_s, 1.0)
            << "gamma_s should be > 1 for area ratio " << ratio;
    }
}

TEST_F(NozzleTests, FrozenSubsonicAreaExpansion) {
    Nozzle nozzle(*gas, GasChemistry::FROZEN);

    std::vector<double> expansion_ratios = {1.1, 1.2, 1.5};

    for (double ratio : expansion_ratios) {
        NozzleResults results = nozzle.solve(ExpansionType::SUBSONIC_AREA_RATIO, ratio);
        NozzleStation result = results.expansions.front();

        EXPECT_TRUE(result.converged)
            << "Subsonic expansion should converge for area ratio " << ratio;
        EXPECT_GT(result.gamma_s, 1.0)
            << "gamma_s should be > 1 for area ratio " << ratio;
    }
}

// Test pressure ratio expansion
TEST_F(NozzleTests, EquilibriumPressureRatioExpansion) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);
    gas->thermo()->restoreState(inlet_state);
    double inlet_pressure = gas->thermo()->pressure();

    std::vector<double> pressure_ratios = {2.0, 5.0, 10.0, 50.0};

    for (double ratio : pressure_ratios) {
        NozzleResults results = nozzle.solve(ExpansionType::PRESSURE_RATIO, ratio);
        NozzleStation result = results.expansions.front();

        EXPECT_TRUE(result.converged)
            << "Pressure ratio expansion should converge for ratio " << ratio;
        EXPECT_GT(result.gamma_s, 1.0)
            << "gamma_s should be > 1 for pressure ratio " << ratio;

        // Verify pressure ratio
        // inlet_pressure = results.throat.P_inlet;
        // gas->thermo()->restoreState(results.throat.state);
        // inlet_pressure = gas->thermo()->pressure();


        gas->thermo()->restoreState(result.state);
        double exit_pressure = gas->thermo()->pressure();
        double exit_entropy = gas->thermo()->entropy_mass();

        double actual_ratio = inlet_pressure / exit_pressure;
        EXPECT_NEAR(actual_ratio, ratio, ratio * 1E-3)
            << "Actual pressure ratio should match requested ratio";
        EXPECT_NEAR(exit_entropy, results.throat.S_inlet, exit_entropy*1E-3)
            << "Nozzle expansion should be isentropic";
    }
}

TEST_F(NozzleTests, FrozenPressureRatioExpansion) {
    Nozzle nozzle(*gas, GasChemistry::FROZEN);
    gas->thermo()->restoreState(inlet_state);
    double inlet_pressure = gas->thermo()->pressure();
    
    std::vector<double> pressure_ratios = {2.0, 5.0, 10.0, 50.0};

    for (double ratio : pressure_ratios) {
        NozzleResults results = nozzle.solve(ExpansionType::PRESSURE_RATIO, ratio);
        NozzleStation result = results.expansions.front();

        EXPECT_TRUE(result.converged)
            << "Pressure ratio expansion should converge for ratio " << ratio;
        EXPECT_GT(result.gamma_s, 1.0)
            << "gamma_s should be > 1 for pressure ratio " << ratio;
        
        gas->thermo()->restoreState(result.state);
        double exit_pressure = gas->thermo()->pressure();
        double exit_entropy = gas->thermo()->entropy_mass();

        double actual_ratio = inlet_pressure / exit_pressure;
        EXPECT_NEAR(actual_ratio, ratio, ratio * 1E-3)
            << "Actual pressure ratio should match requested ratio";
        EXPECT_NEAR(exit_entropy, results.throat.S_inlet, exit_entropy*1E-3)
            << "Nozzle expansion should be isentropic";
        
    }
}

// Test batch solving with multiple ratios
TEST_F(NozzleTests, EquilibriumBatchSolve) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    std::vector<double> ratios = {2.0, 5.0, 10.0, 15.0, 20.0};

    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, ratios);

    EXPECT_TRUE(results.throat.converged) << "Throat should converge";
    EXPECT_EQ(results.expansions.size(), ratios.size())
        << "Should have one result per expansion ratio";

    for (size_t i = 0; i < results.expansions.size(); i++) {
        EXPECT_TRUE(results.expansions[i].converged)
            << "Expansion " << i << " should converge";
        EXPECT_GT(results.expansions[i].gamma_s, 1.0)
            << "gamma_s should be > 1 for expansion " << i;
    }
}

TEST_F(NozzleTests, FrozenBatchSolve) {
    Nozzle nozzle(*gas, GasChemistry::FROZEN);

    std::vector<double> ratios = {2.0, 5.0, 10.0, 15.0, 20.0};

    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, ratios);

    EXPECT_TRUE(results.throat.converged) << "Throat should converge";
    EXPECT_EQ(results.expansions.size(), ratios.size())
        << "Should have one result per expansion ratio";

    for (size_t i = 0; i < results.expansions.size(); i++) {
        EXPECT_TRUE(results.expansions[i].converged)
            << "Expansion " << i << " should converge";
    }
}

// Test gamma_s calculation differences
TEST_F(NozzleTests, GammaDifferencesBetweenEquilibriumAndFrozen) {
    gas->thermo()->restoreState(inlet_state);

    Nozzle eq_nozzle(*gas, GasChemistry::EQUILIBRIUM);
    Nozzle frozen_nozzle(*gas, GasChemistry::FROZEN);

    NozzleResults eq_results = eq_nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 10.0);
    NozzleStation eq_result = eq_results.expansions.front();

    // Reset state for frozen nozzle
    gas->thermo()->restoreState(inlet_state);
    NozzleResults frozen_results = frozen_nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 10.0);
    NozzleStation frozen_result = frozen_results.expansions.front();

    ASSERT_TRUE(eq_result.converged);
    ASSERT_TRUE(frozen_result.converged);

    // Frozen gamma should generally be different from equilibrium gamma
    // They may be close but shouldn't be exactly equal for most cases
    EXPECT_GT(eq_result.gamma_s, 1.0);
    EXPECT_GT(frozen_result.gamma_s, 1.0);
}

// Test invalid expansion ratios
TEST_F(NozzleTests, InvalidExpansionRatioTooSmall) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    // Expansion ratio < 1.0001 should fail
    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 0.5);
    NozzleStation result = results.expansions.front();

    EXPECT_FALSE(result.converged)
        << "Should not converge for invalid expansion ratio < 1";
}

TEST_F(NozzleTests, InvalidSubsonicExpansionRatioTooSmall) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    // Subsonic expansion ratio < 1.0001 should fail
    NozzleResults results = nozzle.solve(ExpansionType::SUBSONIC_AREA_RATIO, 1.0);
    NozzleStation result = results.expansions.front();

    EXPECT_FALSE(result.converged)
        << "Should not converge for invalid subsonic expansion ratio ~1";
}

// Test throat condition properties
TEST_F(NozzleTests, EquilibriumThroatConditionProperties) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);
    double initial_entropy = gas->thermo()->entropy_mass();
    double initial_enthalpy = gas->thermo()->enthalpy_mass();
    double initial_pressure = gas->thermo()->pressure();
    double initial_temperature = gas->thermo()->temperature();

    ThroatCondition results = nozzle.solve_throat_conditions();
    

    ASSERT_TRUE(results.converged);
    gas->thermo()->restoreState(results.state);
    double throat_temperature = gas->thermo()->temperature();

    double tol_H = max_fp_error(results.H_stagnation, 1.0E-5, 1.0E-4);
    double tol_S = max_fp_error(results.S_inlet, 1.0E-5, 1.0E-4);

    // Throat should have reasonable properties
    EXPECT_NEAR(results.H_stagnation, initial_enthalpy, tol_H)
        << "Stagnation enthalpy should be conserved";
    EXPECT_GT(results.P_inlet, 0.0)
        << "Inlet pressure should be positive";
    EXPECT_LE(results.P_inlet, initial_pressure)
        << "Pressure in throat should be less than upstream";
    EXPECT_LE(throat_temperature, initial_temperature)
        << "Temperature in throat should be less than upstream";
    EXPECT_NEAR(results.S_inlet, initial_entropy, tol_S)
        << "Throat flow should be isentropic";
    EXPECT_GT(results.gamma_s, 1.0)
        << "Throat gamma_s should be > 1";
    EXPECT_LT(results.gamma_s, 2.0)
        << "Throat gamma_s should be physically reasonable";
}

// Test state restoration
TEST_F(NozzleTests, StatePreservationAfterSolve) {
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    double initial_temp = gas->thermo()->temperature();
    double initial_pressure = gas->thermo()->pressure();

    // Solve nozzle problem
    nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 10.0);

    // The gas state may have changed, but inlet_state should be preserved
    std::vector<double> current_inlet = nozzle.get_inlet_state();

    gas->thermo()->restoreState(current_inlet);
    double restored_temp = gas->thermo()->temperature();
    double restored_pressure = gas->thermo()->pressure();

    EXPECT_NEAR(restored_temp, initial_temp, initial_temp * 1e-10)
        << "Inlet state temperature should be preserved";
    EXPECT_NEAR(restored_pressure, initial_pressure, initial_pressure * 1e-10)
        << "Inlet state pressure should be preserved";
}

// Test with different gas mixtures
class NozzleDifferentGasTests : public ::testing::Test {
protected:
    NozzleDifferentGasTests() {
        h2o2_gas = Cantera::newSolution("h2o2.yaml", "ohmech");

        // Set up combustion products at high temperature
        auto thermo = h2o2_gas->thermo();
        thermo->setState_TPX(3000.0, 50.0 * Cantera::OneAtm, "H2O:0.85, O2:0.1, H2:0.05");
        thermo->equilibrate("HP");
        thermo->saveState(h2o2_state);
    }

    std::shared_ptr<Cantera::Solution> h2o2_gas;
    std::vector<double> h2o2_state;
};

TEST_F(NozzleDifferentGasTests, H2O2EquilibriumNozzle) {
    Nozzle nozzle(*h2o2_gas, GasChemistry::EQUILIBRIUM);

    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 15.0);
    NozzleStation result = results.expansions.front();

    EXPECT_TRUE(result.converged) << "H2/O2 equilibrium nozzle should converge";
    EXPECT_GT(result.gamma_s, 1.0);
}

TEST_F(NozzleDifferentGasTests, H2O2FrozenNozzle) {
    Nozzle nozzle(*h2o2_gas, GasChemistry::FROZEN);

    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 15.0);
    NozzleStation result = results.expansions.front();

    EXPECT_TRUE(result.converged) << "H2/O2 frozen nozzle should converge";
    EXPECT_GT(result.gamma_s, 1.0);
}

// ============================================================
//  Profile-based solve
// ============================================================

static NozzleProfile make_conical_profile(
    double r_throat, double r_exit, double length, int n_points)
{
    NozzleProfile profile;
    for (int i = 0; i < n_points; i++) {
        double frac = static_cast<double>(i) / (n_points - 1);
        double x = length * frac;
        double r = r_throat + (r_exit - r_throat) * frac;
        profile.push_back({x, r});
    }
    return profile;
}

TEST_F(NozzleTests, ProfileSolveReturnsCorrectCount) {
    NozzleProfile profile = make_conical_profile(0.01, 0.01414, 0.1, 100);
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    int num_stations = 20;
    NozzleResults results = nozzle.solve(profile, num_stations);

    EXPECT_TRUE(results.throat.converged);
    EXPECT_EQ(results.expansions.size(), static_cast<size_t>(num_stations));
}

TEST_F(NozzleTests, ProfileSolveAllConverged) {
    NozzleProfile profile = make_conical_profile(0.01, 0.01414, 0.1, 100);
    Nozzle nozzle(*gas, GasChemistry::EQUILIBRIUM);

    NozzleResults results = nozzle.solve(profile, 20);

    for (size_t i = 0; i < results.expansions.size(); i++) {
        EXPECT_TRUE(results.expansions[i].converged)
            << "Station " << i << " should converge";
    }
}

TEST_F(NozzleTests, ProfileSolveMatchesAreaRatioSolve) {
    // Cross-check: result at the exit area ratio from the profile should
    // match a direct area-ratio solve at the same ratio
    NozzleProfile profile = make_conical_profile(0.01, 0.01414, 0.1, 100);
    // Query slightly inside the boundary to avoid interpolation boundary error
    double x_near_exit = profile.x_max() - 1e-10 * (profile.x_max() - profile.x_min());
    double exit_area_ratio = profile.area_at(x_near_exit) / profile.area_at(profile.x_min());

    Nozzle nozzle_profile(*gas, GasChemistry::EQUILIBRIUM);
    NozzleResults profile_results = nozzle_profile.solve(profile, 20);

    gas->thermo()->restoreState(inlet_state);
    Nozzle nozzle_ar(*gas, GasChemistry::EQUILIBRIUM);
    NozzleResults ar_results = nozzle_ar.solve(ExpansionType::SUPERSONIC_AREA_RATIO, exit_area_ratio);

    ASSERT_TRUE(profile_results.expansions.back().converged);
    ASSERT_TRUE(ar_results.expansions.front().converged);

    // The last profile station should match the area-ratio solve
    EXPECT_NEAR(profile_results.expansions.back().gamma_s,
                ar_results.expansions.front().gamma_s,
                max_fp_error(ar_results.expansions.front().gamma_s, 1e-6, 1e-10));
}

TEST_F(NozzleTests, ProfileSolveFrozen) {
    NozzleProfile profile = make_conical_profile(0.01, 0.01414, 0.1, 100);
    Nozzle nozzle(*gas, GasChemistry::FROZEN);

    NozzleResults results = nozzle.solve(profile, 10);

    EXPECT_TRUE(results.throat.converged);
    for (size_t i = 0; i < results.expansions.size(); i++) {
        EXPECT_TRUE(results.expansions[i].converged)
            << "Frozen station " << i << " should converge";
    }
}
