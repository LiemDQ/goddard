#include "goddard/gas.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/numerics.hpp"

#include <memory>
#include "cantera/core.h"
#include "gtest/gtest.h"

using namespace Goddard;

class GasTests : public ::testing::Test {
protected:
    GasTests() {
        sol = Cantera::newSolution("h2o2.yaml", "ohmech");

        double temp = 3500.0;
        double pressure = 70.0 * Cantera::OneAtm;

        auto thermo = sol->thermo();
        thermo->setState_TPX(temp, pressure, "H2O:0.9, H2:0.05, O2:0.03, OH:0.02");
        thermo->equilibrate("HP");

        thermo->saveState(inlet_state);
    }

    void SetUp() override {
        sol->thermo()->restoreState(inlet_state);
    }

    std::shared_ptr<Cantera::Solution> sol;
    std::vector<double> inlet_state;
};

// Construction

TEST_F(GasTests, ConstructFromSharedPtr) {
    ASSERT_NO_THROW(Gas(sol, GasChemistry::FROZEN));
}

TEST_F(GasTests, ConstructFromReference) {
    ASSERT_NO_THROW(Gas(*sol, GasChemistry::EQUILIBRIUM));
}

// Reference state defaults

TEST_F(GasTests, DefaultStagnationEnthalpy) {
    Gas gas(sol, GasChemistry::FROZEN);
    EXPECT_DOUBLE_EQ(gas.get_stagnation_enthalpy(), sol->thermo()->enthalpy_mass());
}

TEST_F(GasTests, DefaultReferenceEntropy) {
    Gas gas(sol, GasChemistry::FROZEN);
    EXPECT_DOUBLE_EQ(gas.get_reference_entropy(), sol->thermo()->entropy_mass());
}

// Basic properties match Cantera directly

TEST_F(GasTests, BasicPropertiesMatchCantera) {
    Gas gas(sol, GasChemistry::FROZEN);
    auto thermo = sol->thermo();

    EXPECT_DOUBLE_EQ(gas.temperature(), thermo->temperature());
    EXPECT_DOUBLE_EQ(gas.pressure(), thermo->pressure());
    EXPECT_DOUBLE_EQ(gas.density(), thermo->density());
    EXPECT_DOUBLE_EQ(gas.enthalpy_mass(), thermo->enthalpy_mass());
    EXPECT_DOUBLE_EQ(gas.entropy_mass(), thermo->entropy_mass());
    EXPECT_DOUBLE_EQ(gas.cp_mass(), thermo->cp_mass());
    EXPECT_DOUBLE_EQ(gas.cv_mass(), thermo->cv_mass());
    EXPECT_DOUBLE_EQ(gas.molecular_weight(), thermo->meanMolecularWeight());
}

// Gamma_s

TEST_F(GasTests, FrozenGammaIsCpOverCv) {
    Gas gas(sol, GasChemistry::FROZEN);
    double expected = sol->thermo()->cp_mass() / sol->thermo()->cv_mass();
    EXPECT_DOUBLE_EQ(gas.gamma_s(), expected);
}

TEST_F(GasTests, EquilibriumGammaMatchesEquilibriumFunction) {
    Gas gas(sol, GasChemistry::EQUILIBRIUM);
    double expected = get_thermo_equilibrium_properties(*sol->thermo()).gamma_s;
    EXPECT_DOUBLE_EQ(gas.gamma_s(), expected);
}

TEST_F(GasTests, FrozenAndEquilibriumGammaDiffer) {
    Gas frozen(sol, GasChemistry::FROZEN);
    Gas eq(sol, GasChemistry::EQUILIBRIUM);

    // For a reactive mixture these should generally differ
    EXPECT_NE(frozen.gamma_s(), eq.gamma_s());
}

// Speed of sound

TEST_F(GasTests, SpeedOfSoundMatchesFreeFunction) {
    Gas gas(sol, GasChemistry::FROZEN);
    double expected = gas_sonic_velocity(*sol->thermo(), gas.gamma_s());
    EXPECT_DOUBLE_EQ(gas.speed_of_sound(), expected);
}

// Stagnation properties

TEST_F(GasTests, StagnationEnthalpyMatchesFreeFunction) {
    Gas gas(sol, GasChemistry::FROZEN);
    double velocity = 500.0;
    double expected = gas_stagnation_enthalpy(*sol->thermo(), velocity);
    EXPECT_DOUBLE_EQ(gas.stagnation_enthalpy(velocity), expected);
}

TEST_F(GasTests, StagnationPressureMatchesFreeFunction) {
    Gas gas(sol, GasChemistry::FROZEN);
    double velocity = 500.0;
    double expected = gas_stagnation_pressure(*sol->thermo(), velocity);
    EXPECT_NEAR(gas.stagnation_pressure(velocity), expected,
                max_fp_error(expected, 1e-10, 1e-6));
}

TEST_F(GasTests, IsenthalpicVelocityMatchesFreeFunction) {
    Gas gas(sol, GasChemistry::FROZEN);
    double H_stag = gas.enthalpy_mass() + 100000.0;
    double expected = gas_isenthalpic_velocity(*sol->thermo(), H_stag);
    EXPECT_DOUBLE_EQ(gas.isenthalpic_velocity(H_stag), expected);
}

TEST_F(GasTests, IsenthalpicVelocityNoArgUsesStoredEnthalpy) {
    Gas gas(sol, GasChemistry::FROZEN);
    // At construction, H_stag = enthalpy_mass, so velocity should be ~0
    // (within floating point precision of sqrt(2*(H - H)) = sqrt(0))
    EXPECT_NEAR(gas.isenthalpic_velocity(), 0.0, 1e-6);

    // Set a higher stagnation enthalpy and verify it uses the stored value
    double H_stag = gas.enthalpy_mass() + 100000.0;
    gas.set_stagnation_enthalpy(H_stag);
    double expected = gas_isenthalpic_velocity(*sol->thermo(), H_stag);
    EXPECT_DOUBLE_EQ(gas.isenthalpic_velocity(), expected);
}

// Expansion properties

TEST_F(GasTests, EquilibriumExpansionPropertiesMatchFunction) {
    Gas gas(sol, GasChemistry::EQUILIBRIUM);
    auto expected = get_thermo_equilibrium_properties(*sol->thermo());
    auto result = gas.expansion_properties();

    EXPECT_DOUBLE_EQ(result.dlogV_dlogT_P, expected.dlogV_dlogT_P);
    EXPECT_DOUBLE_EQ(result.dlogV_dlogP_T, expected.dlogV_dlogP_T);
    EXPECT_DOUBLE_EQ(result.spec_heat_p, expected.spec_heat_p);
    EXPECT_DOUBLE_EQ(result.gamma_s, expected.gamma_s);
}

TEST_F(GasTests, FrozenExpansionPropertiesAreIdealGas) {
    Gas gas(sol, GasChemistry::FROZEN);
    auto result = gas.expansion_properties();

    EXPECT_DOUBLE_EQ(result.dlogV_dlogT_P, 1.0);
    EXPECT_DOUBLE_EQ(result.dlogV_dlogP_T, -1.0);
    EXPECT_DOUBLE_EQ(result.spec_heat_p, sol->thermo()->cp_mass());
    EXPECT_DOUBLE_EQ(result.gamma_s, sol->thermo()->cp_mass() / sol->thermo()->cv_mass());
}

// Snapshot

TEST_F(GasTests, SnapshotPopulatesAllFields) {
    Gas gas(sol, GasChemistry::FROZEN);
    auto info = gas.snapshot();

    EXPECT_DOUBLE_EQ(info.temperature, gas.temperature());
    EXPECT_DOUBLE_EQ(info.pressure, gas.pressure());
    EXPECT_DOUBLE_EQ(info.density, gas.density());
    EXPECT_DOUBLE_EQ(info.enthalpy, gas.enthalpy_mass());
    EXPECT_DOUBLE_EQ(info.entropy, gas.entropy_mass());
    EXPECT_DOUBLE_EQ(info.molecular_weight, gas.molecular_weight());
    EXPECT_DOUBLE_EQ(info.cp, gas.cp_mass());
    EXPECT_DOUBLE_EQ(info.gamma_s, gas.gamma_s());
    EXPECT_DOUBLE_EQ(info.speed_of_sound, gas.speed_of_sound());
    EXPECT_DOUBLE_EQ(info.stagnation_enthalpy, gas.get_stagnation_enthalpy());
    EXPECT_FALSE(info.composition.empty());
}

// Save/restore roundtrip

TEST_F(GasTests, SaveRestoreRoundtrip) {
    Gas gas(sol, GasChemistry::FROZEN);

    std::vector<double> state;
    gas.copy_state(state);
    double T_original = gas.temperature();

    // Modify state
    gas.set_state_TP(2000.0, 1e5);
    EXPECT_NE(gas.temperature(), T_original);

    // Restore
    gas.restore_state(state);
    EXPECT_DOUBLE_EQ(gas.temperature(), T_original);
}
