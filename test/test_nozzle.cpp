#include "goddard/nozzle.hpp"
#include "goddard/profile.hpp"
#include "goddard/combustor.hpp"
#include "goddard/utils.hpp"
#include "goddard/numerics.hpp"
#include "goddard/error.hpp"
#include "goddard/gas.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/config.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <memory>
#include <iostream>
#include "eigen3/Eigen/Dense"
#include "cantera/core.h"
#include "cantera/base/logger.h"
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
        options.chemistry = GasChemistry::EQUILIBRIUM;

    }

    void SetUp() override {
        // Reset to inlet state before each test
        gas->thermo()->restoreState(inlet_state);
    }

    std::shared_ptr<Cantera::Solution> gas;
    std::vector<double> inlet_state;
    NozzleOptions options;
};

// Test basic construction and state management
TEST_F(NozzleTests, EquilibriumNozzleConstruction) {
    ASSERT_NO_THROW({
        Nozzle nozzle(*gas, options);
    });
}

TEST_F(NozzleTests, FrozenNozzleConstruction) {
    options.chemistry = GasChemistry::FROZEN;
    ASSERT_NO_THROW({
        Nozzle nozzle(*gas, options);
    });
}

TEST_F(NozzleTests, NozzleStateManagement) {
    Nozzle nozzle(*gas, options);

    std::vector<double> original_state = nozzle.get_inlet_state();
    ASSERT_EQ(original_state.size(), gas->thermo()->stateSize());

    // Solve throat to modify nozzle's internal state
    nozzle.solve_throat_conditions();

    // Reset should restore original inlet state
    nozzle.reset_state();

    std::vector<double> restored_state = nozzle.get_inlet_state();
    ASSERT_EQ(restored_state.size(), original_state.size());
    for (size_t i = 0; i < original_state.size(); i++) {
        EXPECT_DOUBLE_EQ(restored_state[i], original_state[i]);
    }
}

TEST_F(NozzleTests, SetInletState) {
    Nozzle nozzle(*gas, options);

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
    Nozzle nozzle(*gas, options);

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

TEST_F(NozzleTests, ThroatNonConvergenceThrows) {
    Nozzle nozzle(*gas, options);
    EXPECT_THROW(nozzle.solve_throat_conditions(1e-14), ConvergenceError);
}

// Test throat condition calculation
TEST_F(NozzleTests, EquilibriumConverge) {
    Nozzle nozzle(*gas, options);

    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 10.0);
    NozzleStation result = results.expansions.front();

    EXPECT_TRUE(result.converged) << "Throat conditions should converge for typical rocket conditions";
    EXPECT_GT(result.gamma_s, 1.0) << "Specific heat ratio should be greater than 1";
    EXPECT_LT(result.gamma_s, 2.0) << "Specific heat ratio should be physically reasonable";
}

TEST_F(NozzleTests, FrozenConverge) {
    options.chemistry = GasChemistry::FROZEN;
    Nozzle nozzle(*gas, options);

    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 10.0);
    NozzleStation result = results.expansions.front();

    EXPECT_TRUE(result.converged) << "Throat conditions should converge for typical rocket conditions";
    EXPECT_GT(result.gamma_s, 1.0) << "Specific heat ratio should be greater than 1";
    EXPECT_LT(result.gamma_s, 2.0) << "Specific heat ratio should be physically reasonable";
}

// Test supersonic area expansion
TEST_F(NozzleTests, EquilibriumSupersonicAreaExpansion) {
    Nozzle nozzle(*gas, options);

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
    options.chemistry = GasChemistry::FROZEN;
    Nozzle nozzle(*gas, options);

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
    Nozzle nozzle(*gas, options);

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
    Nozzle nozzle(*gas, options);

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
    Nozzle nozzle(*gas, options);

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
    options.chemistry = GasChemistry::FROZEN;
    Nozzle nozzle(*gas, options);

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
    Nozzle nozzle(*gas, options);
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
    options.chemistry = GasChemistry::FROZEN;
    Nozzle nozzle(*gas, options);
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
    Nozzle nozzle(*gas, options);

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
    options.chemistry = GasChemistry::FROZEN;
    Nozzle nozzle(*gas, options);

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

    Nozzle eq_nozzle(*gas, options);
    NozzleOptions frozen_options = options;
    frozen_options.chemistry = GasChemistry::FROZEN;
    Nozzle frozen_nozzle(*gas, frozen_options);

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
    Nozzle nozzle(*gas, options);

    // Expansion ratio < 1.0001 should fail
    EXPECT_THROW(nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 0.5), std::invalid_argument);
}

TEST_F(NozzleTests, InvalidSubsonicExpansionRatioTooSmall) {
    Nozzle nozzle(*gas, options);

    // Subsonic expansion ratio < 1.0001 should fail
    EXPECT_THROW(nozzle.solve(ExpansionType::SUBSONIC_AREA_RATIO, 1.0), std::invalid_argument);
}

// Test throat condition properties
TEST_F(NozzleTests, EquilibriumThroatConditionProperties) {
    Nozzle nozzle(*gas, options);
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
    Nozzle nozzle(*gas, options);

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
    NozzleOptions options;
};

TEST_F(NozzleDifferentGasTests, H2O2EquilibriumNozzle) {
    
    Nozzle nozzle(*h2o2_gas, options);

    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 15.0);
    NozzleStation result = results.expansions.front();

    EXPECT_TRUE(result.converged) << "H2/O2 equilibrium nozzle should converge";
    EXPECT_GT(result.gamma_s, 1.0);
}

TEST_F(NozzleDifferentGasTests, H2O2FrozenNozzle) {
    options.chemistry = GasChemistry::FROZEN;
    Nozzle nozzle(*h2o2_gas, options);

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
    Nozzle nozzle(*gas, options);

    int num_stations = 20;
    NozzleResults results = nozzle.solve(profile, num_stations);

    EXPECT_TRUE(results.throat.converged);
    EXPECT_EQ(results.expansions.size(), static_cast<size_t>(num_stations));
}

TEST_F(NozzleTests, ProfileSolveAllConverged) {
    NozzleProfile profile = make_conical_profile(0.01, 0.01414, 0.1, 100);
    Nozzle nozzle(*gas, options);

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

    Nozzle nozzle_profile(*gas, options);
    NozzleResults profile_results = nozzle_profile.solve(profile, 20);

    gas->thermo()->restoreState(inlet_state);
    Nozzle nozzle_ar(*gas, options);
    NozzleResults ar_results = nozzle_ar.solve(ExpansionType::SUPERSONIC_AREA_RATIO, exit_area_ratio);

    ASSERT_TRUE(profile_results.expansions.back().converged);
    ASSERT_TRUE(ar_results.expansions.front().converged);

    // The last profile station should match the area-ratio solve
    EXPECT_NEAR(profile_results.expansions.back().gamma_s,
                ar_results.expansions.front().gamma_s,
                max_fp_error(ar_results.expansions.front().gamma_s, 1e-6, 1e-10));
}

TEST_F(NozzleTests, ProfileSolveFrozen) {
    options.chemistry = GasChemistry::FROZEN;
    NozzleProfile profile = make_conical_profile(0.01, 0.01414, 0.1, 100);
    Nozzle nozzle(*gas, options);

    NozzleResults results = nozzle.solve(profile, 10);

    EXPECT_TRUE(results.throat.converged);
    for (size_t i = 0; i < results.expansions.size(); i++) {
        EXPECT_TRUE(results.expansions[i].converged)
            << "Frozen station " << i << " should converge";
    }
}


// ---------------------------------------------------------------------------------------------
// Perfect-gas finite-area combustor relations
// ---------------------------------------------------------------------------------------------

TEST(FiniteAreaGasDynamicsTests, MachFromAreaRatioInvertsAreaMachRelation) {
    const double gamma = 1.4;
    for (double mach : {0.02, 0.3, 0.7, 0.95}) {
        const double area_ratio = area_mach_relation(mach, gamma);
        EXPECT_NEAR(mach_from_area_ratio(area_ratio, gamma, false), mach, 1e-9 * mach)
            << "subsonic root at M = " << mach;
    }
    for (double mach : {1.05, 2.0, 4.0, 12.0}) {
        const double area_ratio = area_mach_relation(mach, gamma);
        EXPECT_NEAR(mach_from_area_ratio(area_ratio, gamma, true), mach, 1e-9 * mach)
            << "supersonic root at M = " << mach;
    }
    EXPECT_DOUBLE_EQ(mach_from_area_ratio(1.0, gamma, false), 1.0);
    EXPECT_THROW(mach_from_area_ratio(0.99, gamma, true), std::invalid_argument);
}

TEST(FiniteAreaGasDynamicsTests, PressureLossLimits) {
    const double gamma = 1.4;
    // phi(0) = 1: a chamber with stagnant gas loses no stagnation pressure.
    EXPECT_DOUBLE_EQ(finite_area_pressure_loss(0.0, gamma), 1.0);
    // phi(1) = (gamma+1) / ((gamma+1)/2)^(gamma/(gamma-1)) = 1.2679 at gamma = 1.4.
    EXPECT_NEAR(finite_area_pressure_loss(1.0, gamma), 1.2679, 1e-4);
    // Small-M expansion phi = 1 + gamma/2 M^2 + O(M^4).
    const double mach = 0.01;
    EXPECT_NEAR(finite_area_pressure_loss(mach, gamma), 1.0 + 0.5 * gamma * mach * mach,
        std::pow(mach, 4));
}

// ---------------------------------------------------------------------------------------------
// Finite-area combustor chamber
// ---------------------------------------------------------------------------------------------

namespace {
std::string fac_reactant_file() { return std::string(DATA_DIR) + "/nasa9_reactants.yaml"; }
std::string fac_products_file() { return std::string(DATA_DIR) + "/nasa9_gas.yaml"; }
} // namespace

/**
 * Gaseous H2/O2 at 298.15 K, O/F 5.55157, P_inj = 53.3172 bar. Reference values are from the
 * `cea` Python package (`cea.RocketSolver(..., iac=False)`) with all NASA9 H/O gas products.
 */
class FiniteAreaCombustorTests : public ::testing::Test {
protected:
    static constexpr double OF_RATIO = 5.55157;
    static constexpr double INJECTOR_PRESSURE = 53.3172e5;  // Pa
    static constexpr double REACTANT_TEMPERATURE = 298.15;  // K
    // Gordon-McBride and Cantera NASA9 data agree to about this relative level.
    static constexpr double CEA_RELTOL = 2e-4;

    FiniteAreaCombustorTests()
        : products(Gas::create_from_elements(fac_products_file(), "gas", {"H", "O"})) {
        Gas fuel = Gas::create_from_species(fac_reactant_file(), "reactants", {"H2"});
        fuel.set_state_TPX(REACTANT_TEMPERATURE, INJECTOR_PRESSURE, "H2:1");
        Gas oxidizer = Gas::create_from_species(fac_reactant_file(), "reactants", {"O2"});
        oxidizer.set_state_TPX(REACTANT_TEMPERATURE, INJECTOR_PRESSURE, "O2:1");

        // Injector state: equilibrium at the propellant enthalpy and P_inj. The propellant
        // enthalpy is almost exactly zero (elements in their reference states at 298.15 K), where
        // Cantera's "gibbs" HP solver cannot converge because its enthalpy criterion is relative.
        // `Combustor::solve` uses that solver, so this solves with "vcs" instead.
        const double fuel_fraction = 1.0 / (1.0 + OF_RATIO);
        const double enthalpy = fuel_fraction * fuel.enthalpy_mass()
            + (1.0 - fuel_fraction) * oxidizer.enthalpy_mass();
        const std::vector<std::string> elements = products.element_names();
        Eigen::ArrayXd element_moles(static_cast<long>(elements.size()));
        for (size_t m = 0; m < elements.size(); m++) {
            element_moles(static_cast<long>(m)) =
                fuel_fraction * stream_element_moles(fuel, elements[m])
                + (1.0 - fuel_fraction) * stream_element_moles(oxidizer, elements[m]);
        }
        products.set_element_moles(element_moles, 3500.0, INJECTOR_PRESSURE);
        products.equilibrate_TP(3500.0, INJECTOR_PRESSURE);
        products.set_state_HP(enthalpy, INJECTOR_PRESSURE);
        products.equilibrate("HP", "vcs");
        injector_state = products.save_state();
        products.restore_state(injector_state);
        products.chemistry = GasChemistry::EQUILIBRIUM;
        injector_enthalpy = products.enthalpy_mass();
    }

    /** Element amount [kmol/kg] of `element` in a reactant stream. */
    static double stream_element_moles(const Gas& stream, const std::string& element) {
        const std::vector<std::string> names = stream.element_names();
        const Eigen::ArrayXd moles = stream.element_moles();
        for (size_t m = 0; m < names.size(); m++) {
            if (names[m] == element) {
                return moles(static_cast<long>(m));
            }
        }
        return 0.0;
    }

    double pressure_of(const std::vector<double>& state) {
        products.restore_state(state);
        return products.pressure();
    }

    double temperature_of(const std::vector<double>& state) {
        products.restore_state(state);
        return products.temperature();
    }

    /** Equilibrium Mach number of a station on the isentrope with stagnation enthalpy h_inj. */
    double mach_of(const std::vector<double>& state) {
        products.restore_state(state);
        return products.isenthalpic_velocity(injector_enthalpy) / products.speed_of_sound();
    }

    Gas products;
    std::vector<double> injector_state;
    double injector_enthalpy = 0.0;
    NozzleOptions options;
};

TEST_F(FiniteAreaCombustorTests, ContractionRatioMatchesCEA) {
    Nozzle nozzle(products, injector_state, options);
    const FiniteAreaChamber chamber = nozzle.solve_finite_area_chamber(
        injector_state, CombustorType::FINITE_CONTRACTION_RATIO, 1.58);

    // CEA, ac_at = 1.58.
    EXPECT_NEAR(chamber.injector_pressure, INJECTOR_PRESSURE, 1e-9 * INJECTOR_PRESSURE);
    EXPECT_NEAR(chamber.stagnation_pressure, 49.16230e5, CEA_RELTOL * 49.16230e5);
    EXPECT_NEAR(pressure_of(chamber.stagnation_state), 49.16230e5, CEA_RELTOL * 49.16230e5);
    EXPECT_NEAR(pressure_of(chamber.combustion_end.state), 44.62730e5, CEA_RELTOL * 44.62730e5);
    EXPECT_NEAR(pressure_of(chamber.throat.state), 28.32907e5, CEA_RELTOL * 28.32907e5);
    // The injector temperature itself differs from CEA (3498.17 K) by 1.5 K from the species
    // data, so temperatures are compared as drops from the injector.
    const double T_injector = temperature_of(injector_state);
    EXPECT_NEAR(T_injector, 3498.17, 5e-4 * 3498.17);
    EXPECT_NEAR(T_injector - temperature_of(chamber.stagnation_state), 3498.17 - 3488.69, 0.1);
    EXPECT_NEAR(T_injector - temperature_of(chamber.combustion_end.state), 3498.17 - 3454.61, 0.1);
    EXPECT_NEAR(T_injector - temperature_of(chamber.throat.state), 3498.17 - 3297.73, 0.5);
    EXPECT_NEAR(mach_of(chamber.combustion_end.state), 0.41317, 1e-3 * 0.41317);
    EXPECT_DOUBLE_EQ(chamber.contraction_ratio, 1.58);
    EXPECT_NEAR(chamber.throat.P_inlet, chamber.stagnation_pressure,
        1e-12 * chamber.stagnation_pressure);

    // c* = P_inf / (rho_t a_t) = P_inf A_t / mdot; CEA reports 2389.852 m/s.
    products.restore_state(chamber.throat.state);
    const double cstar_value = chamber.stagnation_pressure
        / (products.density() * chamber.throat.speed_of_sound);
    EXPECT_NEAR(cstar_value, 2389.852, CEA_RELTOL * 2389.852);

    // Momentum balance at the combustion end: P_inj = P_c + rho_c u_c^2.
    products.restore_state(chamber.combustion_end.state);
    const double velocity = products.isenthalpic_velocity(injector_enthalpy);
    EXPECT_NEAR(products.pressure() + products.density() * velocity * velocity,
        INJECTOR_PRESSURE, 1e-5 * INJECTOR_PRESSURE);

    // In contraction mode P_inj,calc is nearly proportional to P_inf, so the secant iteration
    // needs very few momentum-balance evaluations.
    EXPECT_LE(chamber.iterations, 4);
}

TEST_F(FiniteAreaCombustorTests, MassFluxMatchesCEA) {
    Nozzle nozzle(products, injector_state, options);
    const FiniteAreaChamber chamber = nozzle.solve_finite_area_chamber(
        injector_state, CombustorType::FINITE_MASS_FLUX, 1333.9);

    // CEA, mdot = 1333.9 kg/(m^2 s).
    EXPECT_NEAR(chamber.contraction_ratio, 1.534894, CEA_RELTOL * 1.534894);
    EXPECT_NEAR(chamber.stagnation_pressure, 48.92327e5, CEA_RELTOL * 48.92327e5);
    EXPECT_NEAR(pressure_of(chamber.combustion_end.state), 44.09505e5, CEA_RELTOL * 44.09505e5);
    EXPECT_NEAR(pressure_of(chamber.throat.state), 28.19171e5, CEA_RELTOL * 28.19171e5);
    EXPECT_NEAR(mach_of(chamber.combustion_end.state), 0.428209, 1e-3 * 0.428209);
    EXPECT_DOUBLE_EQ(chamber.mass_flux, 1333.9);

    // Continuity: mdot / A_c = rho_t a_t / (A_c / A_t).
    products.restore_state(chamber.throat.state);
    EXPECT_NEAR(products.density() * chamber.throat.speed_of_sound / chamber.contraction_ratio,
        1333.9, 1e-9 * 1333.9);
}

TEST_F(FiniteAreaCombustorTests, MassFluxInvertsContractionRatio) {
    // 1.01 puts the combustion end near M = 0.9, where the subsonic station is hardest to solve.
    for (double contraction_ratio : {1.01, 2.0}) {
        Nozzle contraction_nozzle(products, injector_state, options);
        const FiniteAreaChamber by_contraction = contraction_nozzle.solve_finite_area_chamber(
            injector_state, CombustorType::FINITE_CONTRACTION_RATIO, contraction_ratio, 1e-8);

        Nozzle mass_flux_nozzle(products, injector_state, options);
        const FiniteAreaChamber by_mass_flux = mass_flux_nozzle.solve_finite_area_chamber(
            injector_state, CombustorType::FINITE_MASS_FLUX, by_contraction.mass_flux, 1e-8);

        EXPECT_NEAR(by_mass_flux.contraction_ratio, contraction_ratio, 1e-5 * contraction_ratio);
        EXPECT_NEAR(by_mass_flux.stagnation_pressure, by_contraction.stagnation_pressure,
            1e-6 * by_contraction.stagnation_pressure);
    }
}

TEST_F(FiniteAreaCombustorTests, LargeContractionRatioRecoversInfiniteArea) {
    Nozzle infinite_nozzle(products, injector_state, options);
    const NozzleResults infinite = infinite_nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 10.0);

    Nozzle fac_nozzle(products, injector_state, options);
    const FiniteAreaChamber chamber = fac_nozzle.solve_finite_area_chamber(
        injector_state, CombustorType::FINITE_CONTRACTION_RATIO, 1000.0);
    const std::vector<NozzleStation> exits = fac_nozzle.solve_stations(
        chamber.throat, ExpansionType::SUPERSONIC_AREA_RATIO, {10.0});
    ASSERT_EQ(exits.size(), 1u);

    // At M_c ~ 6e-4 the stagnation-pressure loss is gamma/2 M_c^2 ~ 2e-7.
    EXPECT_NEAR(chamber.stagnation_pressure, INJECTOR_PRESSURE, 1e-5 * INJECTOR_PRESSURE);
    EXPECT_NEAR(pressure_of(chamber.throat.state), pressure_of(infinite.throat.state),
        1e-5 * pressure_of(infinite.throat.state));
    EXPECT_NEAR(temperature_of(chamber.throat.state), temperature_of(infinite.throat.state),
        1e-5 * temperature_of(infinite.throat.state));
    EXPECT_NEAR(pressure_of(exits[0].state), pressure_of(infinite.expansions[0].state),
        1e-5 * pressure_of(infinite.expansions[0].state));
    EXPECT_NEAR(temperature_of(exits[0].state), temperature_of(infinite.expansions[0].state),
        1e-5 * temperature_of(infinite.expansions[0].state));
}

namespace {

/** Cantera logger that records warning messages instead of printing them. */
class WarningRecorder : public Cantera::Logger {
public:
    explicit WarningRecorder(std::vector<std::string>* messages) : m_messages(messages) {}
    void warn(const std::string& /*warning*/, const std::string& msg) override {
        m_messages->push_back(msg);
    }

private:
    std::vector<std::string>* m_messages;
};

/** Warnings emitted by Cantera::warn_user while `action` runs. */
std::vector<std::string> recorded_warnings(const std::function<void()>& action) {
    std::vector<std::string> messages;
    Cantera::setLogger(std::make_unique<WarningRecorder>(&messages));
    try {
        action();
    } catch (...) {
        Cantera::setLogger(std::make_unique<Cantera::Logger>());
        throw;
    }
    Cantera::setLogger(std::make_unique<Cantera::Logger>());
    return messages;
}

bool any_contains(const std::vector<std::string>& messages, const std::string& text) {
    return std::any_of(messages.begin(), messages.end(),
        [&](const std::string& message) { return message.find(text) != std::string::npos; });
}

} // namespace

TEST_F(FiniteAreaCombustorTests, UnresolvedCombustionEndVelocityWarns) {
    const std::string warning = "below the solver's resolution";

    Nozzle realistic_nozzle(products, injector_state, options);
    const std::vector<std::string> realistic = recorded_warnings([&] {
        realistic_nozzle.solve_finite_area_chamber(
            injector_state, CombustorType::FINITE_CONTRACTION_RATIO, 1.58);
    });
    EXPECT_FALSE(any_contains(realistic, warning));

    // M_c ~ 6e-4: the kinetic energy is comparable to the station solve's enthalpy error.
    Nozzle huge_nozzle(products, injector_state, options);
    const std::vector<std::string> huge = recorded_warnings([&] {
        huge_nozzle.solve_finite_area_chamber(
            injector_state, CombustorType::FINITE_CONTRACTION_RATIO, 1000.0);
    });
    EXPECT_TRUE(any_contains(huge, warning));
}

TEST_F(FiniteAreaCombustorTests, FrozenExpansionContinuesFromEquilibriumChamber) {
    options.chemistry = GasChemistry::FROZEN;
    options.frozen_NFZ = 1;
    Nozzle nozzle(products, injector_state, options);
    const FiniteAreaChamber chamber = nozzle.solve_finite_area_chamber(
        injector_state, CombustorType::FINITE_CONTRACTION_RATIO, 1.58);

    // The chamber is in equilibrium regardless of the nozzle chemistry, so it matches CEA.
    EXPECT_NEAR(chamber.stagnation_pressure, 49.16230e5, CEA_RELTOL * 49.16230e5);

    const std::vector<NozzleStation> exits = nozzle.solve_stations(
        chamber.throat, ExpansionType::SUPERSONIC_AREA_RATIO, {10.0});
    ASSERT_EQ(exits.size(), 1u);

    // Frozen at the throat: the exit keeps the throat composition and the stagnation entropy.
    products.restore_state(chamber.throat.state);
    const std::vector<double> throat_composition = products.mole_fractions();
    products.restore_state(chamber.stagnation_state);
    const double stagnation_entropy = products.entropy_mass();
    const double throat_pressure = pressure_of(chamber.throat.state);
    products.restore_state(exits[0].state);
    const std::vector<double> exit_composition = products.mole_fractions();
    for (size_t k = 0; k < throat_composition.size(); k++) {
        EXPECT_NEAR(exit_composition[k], throat_composition[k], 1e-12);
    }
    EXPECT_NEAR(products.entropy_mass(), stagnation_entropy, 1e-5 * std::abs(stagnation_entropy));
    EXPECT_LT(products.pressure(), throat_pressure);
}

TEST_F(FiniteAreaCombustorTests, InvalidInputsThrow) {
    Nozzle nozzle(products, injector_state, options);
    EXPECT_THROW(nozzle.solve_finite_area_chamber(injector_state, CombustorType::INFINITE_AREA, 2.0),
        std::invalid_argument);
    EXPECT_THROW(nozzle.solve_finite_area_chamber(injector_state, CombustorType::NONE, 2.0),
        std::invalid_argument);
    EXPECT_THROW(nozzle.solve_finite_area_chamber(
        injector_state, CombustorType::FINITE_CONTRACTION_RATIO, 1.0), std::invalid_argument);
    EXPECT_THROW(nozzle.solve_finite_area_chamber(
        injector_state, CombustorType::FINITE_MASS_FLUX, 0.0), std::invalid_argument);

    options.chemistry = GasChemistry::FROZEN;
    options.frozen_NFZ = 0;
    Nozzle frozen_nozzle(products, injector_state, options);
    EXPECT_THROW(frozen_nozzle.solve_finite_area_chamber(
        injector_state, CombustorType::FINITE_CONTRACTION_RATIO, 2.0), NotImplementedError);
}

TEST_F(FiniteAreaCombustorTests, MassFluxAboveThermalChokingThrows) {
    Nozzle nozzle(products, injector_state, options);
    // The throat mass flux of the infinite-area chamber is about P_inj / c* = 2231 kg/(m^2 s);
    // a constant-area chamber cannot pass more than roughly that divided by phi(1) ~ 1.25.
    try {
        nozzle.solve_finite_area_chamber(injector_state, CombustorType::FINITE_MASS_FLUX, 3000.0);
        FAIL() << "Expected std::invalid_argument for a thermally choked chamber";
    } catch (const std::invalid_argument& error) {
        EXPECT_NE(std::string(error.what()).find("Maximum mass flux"), std::string::npos)
            << error.what();
    }
}
