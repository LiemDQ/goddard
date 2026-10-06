#include "goddard/gas.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/numerics.hpp"
#include "goddard/utils.hpp"
#include "goddard/config.h"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
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

TEST_F(GasTests, EquilibriumStagnationPressureConservesEntropyAndTotalEnthalpy) {
    // The fixture state is at equilibrium. Decelerating it isentropically with shifting
    // equilibrium must reach the stagnation enthalpy h + v^2/2 at the returned pressure.
    Gas gas(sol, GasChemistry::EQUILIBRIUM);
    const double velocity = 1500.0;
    const double P_static = gas.pressure();
    const double P_stagnation = gas.stagnation_pressure(velocity);

    EXPECT_DOUBLE_EQ(gas.pressure(), P_static) << "stagnation_pressure must not modify the state";

    Gas stagnation = gas.clone();
    stagnation.equilibrate_SP(gas.entropy_mass(), P_stagnation);
    // The solver stops on a 1e-6 relative pressure step, i.e. an enthalpy error of up to
    // 1e-6 * P/rho ~ 2 J/kg here; 10 J/kg is ~1e-5 of the kinetic energy.
    EXPECT_NEAR(stagnation.enthalpy_mass(), gas.stagnation_enthalpy(velocity), 10.0);

    // Recombination during compression releases heat, so the equilibrium and frozen
    // stagnation pressures of this dissociated mixture differ well beyond solver tolerance.
    Gas frozen(sol, GasChemistry::FROZEN);
    const double P_stagnation_frozen = frozen.stagnation_pressure(velocity);
    EXPECT_GT(std::abs(P_stagnation - P_stagnation_frozen) / P_stagnation_frozen, 1e-4);
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

// Copy semantics

TEST_F(GasTests, PerfectGasCopyPreservesGamma) {
    Gas gas(1.4);
    Gas copied(gas); // NOLINT(performance-unnecessary-copy-initialization): the copy is under test
    EXPECT_DOUBLE_EQ(copied.gamma_s(), 1.4);

    Gas assigned(1.2);
    assigned = gas;
    EXPECT_DOUBLE_EQ(assigned.gamma_s(), 1.4);
    EXPECT_EQ(assigned.chemistry, GasChemistry::PERFECT_GAS);
}

TEST_F(GasTests, CopySharesSolution) {
    Gas gas(sol, GasChemistry::FROZEN);
    Gas copied = gas;
    EXPECT_EQ(copied.solution(), gas.solution());

    copied.set_state_TP(2000.0, 1e5);
    EXPECT_DOUBLE_EQ(gas.temperature(), 2000.0);
}

TEST_F(GasTests, CloneIsIndependentDeepCopy) {
    Gas gas(sol, GasChemistry::EQUILIBRIUM);
    gas.set_stagnation_enthalpy(1.234e6);
    gas.set_reference_entropy(5.678e3);

    Gas cloned = gas.clone();
    EXPECT_NE(cloned.solution(), gas.solution());
    EXPECT_EQ(cloned.chemistry, GasChemistry::EQUILIBRIUM);
    EXPECT_DOUBLE_EQ(cloned.temperature(), gas.temperature());
    EXPECT_DOUBLE_EQ(cloned.pressure(), gas.pressure());
    std::vector<double> X_original = gas.mole_fractions();
    std::vector<double> X_cloned = cloned.mole_fractions();
    ASSERT_EQ(X_cloned.size(), X_original.size());
    for (size_t i = 0; i < X_original.size(); i++) {
        EXPECT_DOUBLE_EQ(X_cloned[i], X_original[i]);
    }
    EXPECT_DOUBLE_EQ(cloned.get_stagnation_enthalpy(), 1.234e6);
    EXPECT_DOUBLE_EQ(cloned.get_reference_entropy(), 5.678e3);

    double T_original = gas.temperature();
    cloned.set_state_TP(1000.0, 1e5);
    EXPECT_DOUBLE_EQ(gas.temperature(), T_original);
}

TEST_F(GasTests, ClonePerfectGas) {
    Gas gas(1.3);
    Gas cloned = gas.clone();
    EXPECT_FALSE(cloned.has_cantera_sln());
    EXPECT_DOUBLE_EQ(cloned.gamma_s(), 1.3);
}

// Element bookkeeping and reactant streams

TEST(GasElementTests, ElementMolesMatchElementalMassFractions) {
    Gas gas(Cantera::newSolution("h2o2.yaml", "ohmech"));
    gas.set_state_TPX(1200.0, 3.0 * Cantera::OneAtm, "H2:2, O2:1");

    auto thermo = gas.thermo();
    Eigen::ArrayXd moles = gas.element_moles();
    std::vector<std::string> names = gas.element_names();
    ASSERT_EQ(names.size(), thermo->nElements());
    ASSERT_EQ(static_cast<size_t>(moles.size()), thermo->nElements());

    for (size_t m = 0; m < thermo->nElements(); m++) {
        // kmol of element m per kg of mixture.
        const double expected =
            thermo->elementalMassFraction(m) / thermo->atomicWeight(m);
        EXPECT_NEAR(moles(static_cast<long>(m)), expected, 1e-12 * std::max(expected, 1e-3))
            << "element " << names[m];
    }
}

// The NASA9 reactant database holds condensed reactants such as H2(L) as ordinary species. Their
// polynomials are evaluated in an ideal-gas phase, which gives the right enthalpy and element
// amounts; density and entropy of such a stream are meaningless.
class ReactantGasTests : public ::testing::Test {
protected:
    static double species_enthalpy_mass(const Gas& gas, const std::string& name, double T) {
        auto thermo = gas.thermo();
        thermo->setState_TP(T, Cantera::OneAtm);
        std::vector<double> h_RT(thermo->nSpecies());
        thermo->getEnthalpy_RT_ref(h_RT.data());
        const size_t k = thermo->speciesIndex(name);
        return h_RT[k] * Cantera::GasConstant * T / thermo->molecularWeight(k);
    }

    std::string reactant_file = std::string(DATA_DIR) + "/nasa9_reactants.yaml";
};

TEST_F(ReactantGasTests, CryogenicStreamsCarryEnthalpyAndElements) {
    const double H2_boiling_point = 20.27;
    const double O2_boiling_point = 90.17;

    Gas fuel = Gas::create_from_species(reactant_file, "reactants", {"H2(L)"});
    fuel.set_state_TPX(H2_boiling_point, Cantera::OneAtm, "H2(L):1");
    Gas oxidizer = Gas::create_from_species(reactant_file, "reactants", {"O2(L)"});
    oxidizer.set_state_TPX(O2_boiling_point, Cantera::OneAtm, "O2(L):1");

    EXPECT_NEAR(fuel.enthalpy_mass(),
        species_enthalpy_mass(fuel, "H2(L)", H2_boiling_point), 1e-9 * 1e6);
    EXPECT_NEAR(oxidizer.enthalpy_mass(),
        species_enthalpy_mass(oxidizer, "O2(L)", O2_boiling_point), 1e-9 * 1e6);
    // Liquid hydrogen and oxygen have negative formation enthalpies.
    EXPECT_LT(fuel.enthalpy_mass(), 0.0);
    EXPECT_LT(oxidizer.enthalpy_mass(), 0.0);

    // 2 kmol of H per kmol of H2(L), i.e. 2/M kmol per kg.
    fuel.set_state_TPX(H2_boiling_point, Cantera::OneAtm, "H2(L):1");
    ASSERT_EQ(fuel.element_names(), (std::vector<std::string>{"H"}));
    EXPECT_NEAR(fuel.element_moles()(0), 2.0 / fuel.molecular_weight(), 1e-12);

    oxidizer.set_state_TPX(O2_boiling_point, Cantera::OneAtm, "O2(L):1");
    ASSERT_EQ(oxidizer.element_names(), (std::vector<std::string>{"O"}));
    EXPECT_NEAR(oxidizer.element_moles()(0), 2.0 / oxidizer.molecular_weight(), 1e-12);
}

TEST_F(ReactantGasTests, AssignedEnthalpyReactantIsTemperatureIndependent) {
    // CEA gives reactants such as RP-1 an assigned enthalpy: constant-cp species with cp0 = 0,
    // so h(T) is the formation enthalpy at every temperature.
    Gas fuel = Gas::create_from_species(reactant_file, "reactants", {"RP-1"});
    fuel.set_state_TPX(298.15, Cantera::OneAtm, "RP-1:1");
    const double h_reference = fuel.enthalpy_mass();
    fuel.set_state_TPX(500.0, Cantera::OneAtm, "RP-1:1");
    EXPECT_NEAR(fuel.enthalpy_mass(), h_reference, 1e-9 * std::abs(h_reference));
}

// Species sets arrive unordered (Python passes a set), so the phase must not take their iteration
// order: species indices would then differ between runs and processes.
TEST(GasSpeciesOrderTests, CreateFromSpeciesUsesFileOrder) {
    const std::string infile = "h2o2.yaml";
    const std::vector<std::string> forward{"H2", "O2", "H2O", "OH", "H"};
    const std::vector<std::string> reverse(forward.rbegin(), forward.rend());

    // Expected: the requested species in the order of the file's `species` section.
    const std::unordered_set<std::string> requested(forward.begin(), forward.end());
    const Cantera::AnyMap root_node = load_root_node(infile);
    std::vector<std::string> file_order;
    for (const Cantera::AnyMap& node : root_node.at("species").asVector<Cantera::AnyMap>()) {
        if (requested.count(node.at("name").asString())) {
            file_order.push_back(node.at("name").asString());
        }
    }
    ASSERT_EQ(file_order.size(), forward.size());
    ASSERT_NE(file_order, forward) << "the request must not already be in file order";

    // Sets built by inserting in different orders, as a binding would receive them.
    const std::unordered_set<std::string> set_forward(forward.begin(), forward.end());
    const std::unordered_set<std::string> set_reverse(reverse.begin(), reverse.end());

    const std::vector<std::string> from_set_forward(set_forward.begin(), set_forward.end());
    const std::vector<std::string> from_set_reverse(set_reverse.begin(), set_reverse.end());

    const std::vector<std::vector<std::string>> requests{
        forward, reverse, from_set_forward, from_set_reverse};
    for (const std::vector<std::string>& request : requests) {
        Gas gas = Gas::create_from_species(infile, "gas", request);
        EXPECT_EQ(gas.species_names(), file_order);
    }
}

TEST(GasSpeciesOrderTests, CreateFromSpeciesRejectsUnknownSpecies) {
    EXPECT_THROW(Gas::create_from_species("h2o2.yaml", "gas", {"H2", "NOT_A_SPECIES"}),
                 std::invalid_argument);
}
