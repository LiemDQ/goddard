#include "goddard/combustor.hpp"
#include "goddard/condensed.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/error.hpp"
#include "goddard/gas.hpp"
#include "goddard/global.hpp"
#include "goddard/thermoarray.hpp"
#include "goddard/utils.hpp"

#include "cantera/core.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <string>
#include <vector>

using namespace Goddard;
using ::testing::Contains;
using ::testing::Not;
using ::testing::UnorderedElementsAreArray;

namespace {

/** Cantera's shipped NASA (TM-4513) data files, used until the CEA-derived data lands (WP-E). */
constexpr const char* GAS_FILE = "nasa_gas.yaml";
constexpr const char* CONDENSED_FILE = "nasa_condensed.yaml";

/**
 * Build an ideal-gas `Gas` from named species of `nasa_gas.yaml`.
 *
 * `Gas::create_from_species` cannot be used here: `create_phase_node` declares `kinetics: gas`
 * with `reactions: declared-species`, and `nasa_gas.yaml` has no `reactions` section, so Cantera
 * raises an InputFileError. The phase node below leaves kinetics out.
 */
Gas make_gas(const std::string& name, const std::vector<std::string>& species) {
    setup_defaults();
    Cantera::AnyMap root = load_root_node(GAS_FILE);

    Cantera::AnyMap phase;
    phase["name"] = name;
    phase["thermo"] = "ideal-gas";
    phase["species"] = species;
    phase["state"]["T"] = 298.15;
    phase["state"]["P"] = Cantera::OneAtm;

    return Gas(Cantera::newSolution(phase, root), GasChemistry::EQUILIBRIUM);
}

/** Gas phase spanning the elements of every condensed species used in these tests. */
Gas make_product_gas() {
    return make_gas("products", {"H2", "O2", "H2O", "AL", "ALO", "Be", "BeO", "C", "CO", "CO2"});
}

double element_amount(const Gas& gas, const std::string& element) {
    const std::vector<std::string> names = gas.element_names();
    auto found = std::find(names.begin(), names.end(), element);
    EXPECT_NE(found, names.end()) << "element " << element << " missing from phase";
    return gas.element_moles()(found - names.begin());
}

} // namespace

class CondensedSetTests : public ::testing::Test {
protected:
    CondensedSetTests() : gas(make_product_gas()), root(load_root_node(CONDENSED_FILE)) {}

    Gas gas;
    Cantera::AnyMap root;
};

// ---- CondensedPhaseSet construction ----

TEST_F(CondensedSetTests, ConstructFromNamedSpecies) {
    CondensedPhaseSet set(root, {"C(gr)", "H2O(L)", "AL2O3(a)"}, *gas.thermo());

    EXPECT_EQ(set.size(), 3u);
    EXPECT_THAT(set.names(), UnorderedElementsAreArray({"C(gr)", "H2O(L)", "AL2O3(a)"}));
    EXPECT_EQ(set.moles, std::vector<double>({0.0, 0.0, 0.0}));
    EXPECT_EQ(set.pinned_group, -1);

    const size_t graphite = set.species_index("C(gr)");
    EXPECT_NEAR(set.species[graphite].molar_mass, 12.011, 1e-3);
    EXPECT_DOUBLE_EQ(set.species[graphite].T_min, 200.0);
    EXPECT_DOUBLE_EQ(set.species[graphite].T_max, 5000.0);
}

TEST_F(CondensedSetTests, UnknownSpeciesThrows) {
    EXPECT_THROW(CondensedPhaseSet(root, {"NotASpecies(s)"}, *gas.thermo()), FmtError);
}

TEST_F(CondensedSetTests, SpeciesWithElementMissingFromGasThrows) {
    Gas ho_gas = make_gas("ho", {"H2", "O2", "H2O"});
    EXPECT_THROW(CondensedPhaseSet(root, {"AL2O3(a)"}, *ho_gas.thermo()), FmtError);
}

TEST_F(CondensedSetTests, ElementAtomsUseGasElementOrder) {
    CondensedPhaseSet set(root, {"AL2O3(a)"}, *gas.thermo());

    const std::vector<std::string> elements = gas.element_names();
    const Eigen::ArrayXd& atoms = set.species[0].element_atoms;
    ASSERT_EQ(atoms.size(), static_cast<long>(elements.size()));

    for (size_t m = 0; m < elements.size(); m++) {
        double expected = 0.0;
        if (elements[m] == "Al") expected = 2.0;
        if (elements[m] == "O") expected = 3.0;
        EXPECT_DOUBLE_EQ(atoms(static_cast<long>(m)), expected) << "element " << elements[m];
    }
}

// ---- Polymorph groups ----

TEST_F(CondensedSetTests, WaterGroupHasFreezingTransition) {
    CondensedPhaseSet set(root, {"H2O(L)", "H2O(s)"}, *gas.thermo());

    ASSERT_EQ(set.groups.size(), 1u);
    EXPECT_EQ(set.groups[0].members.size(), 2u);
    // Members are sorted by increasing T_min, so ice comes first even though it was added second.
    EXPECT_EQ(set.species[set.groups[0].members[0]].name, "H2O(s)");
    EXPECT_EQ(set.species[set.groups[0].members[1]].name, "H2O(L)");
    ASSERT_EQ(set.groups[0].transition_temperatures.size(), 1u);
    EXPECT_DOUBLE_EQ(set.groups[0].transition_temperatures[0], 273.15);
}

TEST_F(CondensedSetTests, AluminaGroupHasMeltingTransition) {
    CondensedPhaseSet set(root, {"AL2O3(a)", "AL2O3(L)"}, *gas.thermo());

    ASSERT_EQ(set.groups.size(), 1u);
    ASSERT_EQ(set.groups[0].transition_temperatures.size(), 1u);
    EXPECT_DOUBLE_EQ(set.groups[0].transition_temperatures[0], 2327.0);
}

TEST_F(CondensedSetTests, BerylliaGroupHasTwoTransitions) {
    CondensedPhaseSet set(root, {"BeO(L)", "BeO(a)", "BeO(b)"}, *gas.thermo());

    ASSERT_EQ(set.groups.size(), 1u);
    ASSERT_EQ(set.groups[0].members.size(), 3u);
    EXPECT_EQ(set.species[set.groups[0].members[0]].name, "BeO(a)");
    EXPECT_EQ(set.species[set.groups[0].members[1]].name, "BeO(b)");
    EXPECT_EQ(set.species[set.groups[0].members[2]].name, "BeO(L)");

    ASSERT_EQ(set.groups[0].transition_temperatures.size(), 2u);
    EXPECT_DOUBLE_EQ(set.groups[0].transition_temperatures[0], 2373.001);
    EXPECT_DOUBLE_EQ(set.groups[0].transition_temperatures[1], 2821.22);
}

TEST_F(CondensedSetTests, DistinctCompositionsFormSeparateGroups) {
    CondensedPhaseSet set(root, {"H2O(s)", "H2O(L)", "C(gr)"}, *gas.thermo());

    ASSERT_EQ(set.groups.size(), 2u);
    const size_t graphite = set.species_index("C(gr)");
    const PolymorphGroup& carbon = set.groups[static_cast<size_t>(set.species[graphite].group)];
    EXPECT_EQ(carbon.members, std::vector<size_t>({graphite}));
    EXPECT_TRUE(carbon.transition_temperatures.empty());
}

// ---- Compatible species ----

TEST_F(CondensedSetTests, CompatibleSpeciesRespectsGasElements) {
    const std::vector<std::string> compatible =
        CondensedPhaseSet::compatible_species(root, *gas.thermo());

    EXPECT_THAT(compatible, Contains("C(gr)"));
    EXPECT_THAT(compatible, Contains("H2O(L)"));
    EXPECT_THAT(compatible, Contains("AL2O3(a)"));
    EXPECT_THAT(compatible, Contains("BeO(b)"));
    // The gas phase has no Si, so an aluminosilicate must be excluded.
    EXPECT_THAT(compatible, Not(Contains("AL2SiO5(an)")));

    // Every compatible species can actually be built.
    EXPECT_NO_THROW(CondensedPhaseSet(root, compatible, *gas.thermo()));
}

TEST_F(CondensedSetTests, CompatibleSpeciesOfSmallerGasIsSubset) {
    Gas ho_gas = make_gas("ho", {"H2", "O2", "H2O"});
    const std::vector<std::string> compatible =
        CondensedPhaseSet::compatible_species(root, *ho_gas.thermo());

    EXPECT_THAT(compatible, Contains("H2O(s)"));
    EXPECT_THAT(compatible, Not(Contains("C(gr)")));
    EXPECT_THAT(compatible, Not(Contains("AL2O3(a)")));
}

// ---- Temperature ranges ----

TEST_F(CondensedSetTests, InRangeIsAClosedInterval) {
    CondensedPhaseSet set(root, {"H2O(s)", "H2O(L)"}, *gas.thermo());
    const size_t ice = set.species_index("H2O(s)");
    const size_t liquid = set.species_index("H2O(L)");

    EXPECT_TRUE(set.in_range(ice, 200.0));
    EXPECT_TRUE(set.in_range(ice, 273.15));
    EXPECT_FALSE(set.in_range(ice, 273.16));
    EXPECT_TRUE(set.in_range(liquid, 273.15));
    EXPECT_FALSE(set.in_range(liquid, 273.0));
}

TEST_F(CondensedSetTests, OfferedAtSelectsTheInRangePolymorph) {
    CondensedPhaseSet set(root, {"H2O(s)", "H2O(L)", "C(gr)"}, *gas.thermo());
    const size_t ice = set.species_index("H2O(s)");
    const size_t liquid = set.species_index("H2O(L)");
    const size_t graphite = set.species_index("C(gr)");

    const std::vector<size_t> cold = set.offered_at(250.0);
    EXPECT_THAT(cold, UnorderedElementsAreArray({ice, graphite}));

    const std::vector<size_t> warm = set.offered_at(300.0);
    EXPECT_THAT(warm, UnorderedElementsAreArray({liquid, graphite}));
}

TEST_F(CondensedSetTests, OfferedAtATransitionKeepsOnlyTheLowerPolymorph) {
    CondensedPhaseSet set(root, {"H2O(s)", "H2O(L)"}, *gas.thermo());
    const size_t ice = set.species_index("H2O(s)");

    // Both polymorphs are in range at the transition, but offering both makes MultiPhaseEquil
    // oscillate, so only the lower one is offered.
    EXPECT_TRUE(set.in_range(set.species_index("H2O(L)"), 273.15));
    EXPECT_EQ(set.offered_at(273.15), std::vector<size_t>({ice}));
}

TEST_F(CondensedSetTests, OfferedAtSkipsOutOfRangeGroups) {
    CondensedPhaseSet set(root, {"AL2O3(a)", "AL2O3(L)"}, *gas.thermo());
    // The alumina data starts at 300 K.
    EXPECT_TRUE(set.offered_at(250.0).empty());
}

// ---- Reference-state data ----

TEST_F(CondensedSetTests, ReferenceStateDataMatchesNasaPolynomial) {
    CondensedPhaseSet set(root, {"AL2O3(L)"}, *gas.thermo());
    const double T = 3000.0;

    // AL2O3(L) is a one-interval NASA7 fit with only a1, a6 and a7 nonzero, so the reference
    // state follows directly from the data file.
    const double a1 = 23.148241;
    const double a6 = -2.114052e+05;
    const double a7 = -138.60205;

    EXPECT_NEAR(set.cp_R(0, T), a1, 1e-10);
    EXPECT_NEAR(set.enthalpy_RT(0, T), a1 + a6 / T, 1e-10);
    EXPECT_NEAR(set.entropy_R(0, T), a1 * std::log(T) + a7, 1e-10);
}

TEST_F(CondensedSetTests, LatentHeatSeparatesPolymorphsAtTheTransition) {
    CondensedPhaseSet set(root, {"H2O(s)", "H2O(L)"}, *gas.thermo());
    const size_t ice = set.species_index("H2O(s)");
    const size_t liquid = set.species_index("H2O(L)");
    const double T = 273.15;

    const double latent_heat =
        (set.enthalpy_RT(liquid, T) - set.enthalpy_RT(ice, T)) * Cantera::GasConstant * T;
    // Heat of fusion of water is ~6.0 MJ/kmol.
    EXPECT_NEAR(latent_heat, 6.0e6, 3e5);
}

// ---- clone ----

TEST_F(CondensedSetTests, CloneIsIndependent) {
    CondensedPhaseSet set(root, {"C(gr)", "H2O(L)"}, *gas.thermo());
    set.moles = {0.01, 0.02};
    set.pinned_group = 1;

    auto copy = set.clone();
    EXPECT_EQ(copy->names(), set.names());
    EXPECT_EQ(copy->moles, set.moles);
    EXPECT_EQ(copy->pinned_group, 1);
    EXPECT_NE(copy->species[0].phase.get(), set.species[0].phase.get());

    copy->moles[0] = 0.5;
    EXPECT_DOUBLE_EQ(set.moles[0], 0.01);
}

// ---- add_species ----

TEST_F(CondensedSetTests, AddSpeciesAppendsAndRegroups) {
    CondensedPhaseSet set(root, {"H2O(L)"}, *gas.thermo());
    set.moles[0] = 0.03;

    set.add_species(root, {"H2O(s)", "C(gr)"}, *gas.thermo());

    EXPECT_EQ(set.size(), 3u);
    EXPECT_DOUBLE_EQ(set.moles[0], 0.03); // existing amounts survive
    EXPECT_DOUBLE_EQ(set.moles[1], 0.0);
    EXPECT_EQ(set.groups.size(), 2u);
    ASSERT_EQ(set.groups[set.species[0].group].transition_temperatures.size(), 1u);
    EXPECT_DOUBLE_EQ(set.groups[set.species[0].group].transition_temperatures[0], 273.15);
}

// ---- Gas integration ----

class CondensedGasTests : public ::testing::Test {
protected:
    CondensedGasTests() : gas(make_product_gas()) {
        gas.set_state_TPX(400.0, Cantera::OneAtm, "H2O:1.0");
    }

    Gas gas;
};

TEST_F(CondensedGasTests, GasOnlyByDefault) {
    EXPECT_FALSE(gas.has_condensed_candidates());
    EXPECT_FALSE(gas.has_condensed_phases());
    EXPECT_TRUE(gas.condensed_species_names().empty());
    EXPECT_TRUE(gas.condensed_moles().empty());
    EXPECT_DOUBLE_EQ(gas.gas_mass_fraction(), 1.0);
    EXPECT_DOUBLE_EQ(gas.mixture_molecular_weight(), gas.molecular_weight());
    EXPECT_FALSE(gas.at_phase_transition());
    EXPECT_EQ(gas.pinned_polymorphs(), std::make_pair(-1L, -1L));
}

TEST_F(CondensedGasTests, AddCondensedSpecies) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)", "H2O(s)"});

    EXPECT_TRUE(gas.has_condensed_candidates());
    EXPECT_FALSE(gas.has_condensed_phases()); // no moles yet
    EXPECT_EQ(gas.condensed_species_names(), std::vector<std::string>({"H2O(L)", "H2O(s)"}));
    EXPECT_EQ(gas.condensed_moles(), std::vector<double>({0.0, 0.0}));
}

TEST_F(CondensedGasTests, AddCondensedSpeciesTwiceAppends) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)"});
    gas.add_condensed_species(CONDENSED_FILE, {"C(gr)"});

    EXPECT_EQ(gas.condensed_species_names(), std::vector<std::string>({"H2O(L)", "C(gr)"}));
}

TEST_F(CondensedGasTests, AddAllCondensedSpecies) {
    gas.add_all_condensed_species(CONDENSED_FILE);

    const std::vector<std::string> names = gas.condensed_species_names();
    EXPECT_THAT(names, Contains("C(gr)"));
    EXPECT_THAT(names, Contains("AL2O3(L)"));
    EXPECT_THAT(names, Not(Contains("AL2SiO5(an)")));
    EXPECT_EQ(gas.condensed_moles().size(), names.size());
}

TEST_F(CondensedGasTests, SetCondensedMolesChecksSize) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)", "H2O(s)"});

    EXPECT_THROW(gas.set_condensed_moles({0.01}), std::invalid_argument);
    EXPECT_NO_THROW(gas.set_condensed_moles({0.01, 0.0}));
    EXPECT_TRUE(gas.has_condensed_phases());
    EXPECT_EQ(gas.condensed_moles(), std::vector<double>({0.01, 0.0}));
}

TEST_F(CondensedGasTests, GasMassFractionAndMixtureMolecularWeight) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)"});
    const double n_liquid = 0.01; // kmol per kg of mixture
    gas.set_condensed_moles({n_liquid});

    const double gas_molar_mass = gas.thermo()->meanMolecularWeight(); // pure H2O gas
    const double w_gas = 1.0 - n_liquid * gas_molar_mass;
    EXPECT_DOUBLE_EQ(gas.gas_mass_fraction(), w_gas);

    const double gas_moles = w_gas / gas_molar_mass;
    EXPECT_DOUBLE_EQ(gas.mixture_molecular_weight(), 1.0 / (gas_moles + n_liquid));
    // CEA's "M" counts the gas moles against the whole kg of mixture.
    EXPECT_DOUBLE_EQ(gas.molecular_weight(), 1.0 / gas_moles);
    EXPECT_GT(gas.molecular_weight(), gas_molar_mass);
}

TEST_F(CondensedGasTests, MixtureMassFractionsSumToOne) {
    gas.set_state_TPX(400.0, Cantera::OneAtm, "H2O:0.8, H2:0.1, O2:0.1");
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)", "C(gr)"});
    gas.set_condensed_moles({0.005, 0.002});

    const std::vector<double> fractions = gas.mixture_mass_fractions();
    ASSERT_EQ(fractions.size(), gas.num_species() + 2);

    double total = 0.0;
    for (double y : fractions) total += y;
    EXPECT_NEAR(total, 1.0, 1e-12);

    const std::vector<std::string> names = gas.condensed_species_names();
    const size_t n_gas = gas.num_species();
    EXPECT_NEAR(fractions[n_gas], 0.005 * 18.015, 1e-3);   // H2O(L)
    EXPECT_NEAR(fractions[n_gas + 1], 0.002 * 12.011, 1e-3); // C(gr)
}

TEST_F(CondensedGasTests, ElementMolesConservedWhenMassMovesIntoACondensedPhase) {
    gas.set_state_TPX(400.0, Cantera::OneAtm, "H2O:1.0");
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)"});

    const Eigen::ArrayXd gas_only = gas.element_moles();
    // Pure water vapour: b_H = 2/M, b_O = 1/M.
    EXPECT_NEAR(element_amount(gas, "H"), 2.0 / gas.molecular_weight(), 1e-12);
    EXPECT_NEAR(element_amount(gas, "O"), 1.0 / gas.molecular_weight(), 1e-12);

    // Condensing part of the water changes neither the elements nor their amounts.
    gas.set_condensed_moles({0.02});
    const Eigen::ArrayXd with_condensed = gas.element_moles();

    ASSERT_EQ(with_condensed.size(), gas_only.size());
    for (long m = 0; m < gas_only.size(); m++) {
        EXPECT_NEAR(with_condensed(m), gas_only(m), 1e-12) << "element index " << m;
    }
}

TEST_F(CondensedGasTests, CondensedPerSpeciesDataCoversPresentSpeciesOnly) {
    gas.set_state_TP(3000.0, Cantera::OneAtm);
    gas.add_condensed_species(CONDENSED_FILE, {"AL2O3(L)", "C(gr)"});
    gas.set_condensed_moles({0.0, 0.003}); // only graphite present

    EXPECT_EQ(gas.condensed_enthalpy_RT().size(), 1);
    EXPECT_EQ(gas.condensed_cp_R().size(), 1);
    ASSERT_EQ(gas.condensed_molar_masses().size(), 1);
    EXPECT_NEAR(gas.condensed_molar_masses()(0), 12.011, 1e-3);

    const Eigen::ArrayXXd coeffs = gas.condensed_stoich_coeffs();
    ASSERT_EQ(coeffs.rows(), 1);
    ASSERT_EQ(coeffs.cols(), static_cast<long>(gas.element_names().size()));

    const std::vector<std::string> elements = gas.element_names();
    for (long m = 0; m < coeffs.cols(); m++) {
        EXPECT_DOUBLE_EQ(coeffs(0, m), elements[static_cast<size_t>(m)] == "C" ? 1.0 : 0.0);
    }

    gas.set_condensed_moles({0.001, 0.003});
    EXPECT_EQ(gas.condensed_molar_masses().size(), 2);
    EXPECT_EQ(gas.condensed_stoich_coeffs().rows(), 2);
}

// ---- State round trips ----

TEST_F(CondensedGasTests, SaveAndRestoreExtendedState) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)", "C(gr)"});
    gas.set_condensed_moles({0.01, 0.002});

    const std::vector<double> state = gas.save_state();
    EXPECT_EQ(state.size(), gas.thermo()->stateSize() + 2);

    gas.set_state_TPX(800.0, 5.0 * Cantera::OneAtm, "H2:1.0");
    gas.set_condensed_moles({0.0, 0.0});

    gas.restore_state(state);
    EXPECT_NEAR(gas.temperature(), 400.0, 1e-9);
    EXPECT_EQ(gas.condensed_moles(), std::vector<double>({0.01, 0.002}));
}

TEST_F(CondensedGasTests, RestoreGasOnlyStateZerosCondensedMoles) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)"});

    std::vector<double> cantera_state(gas.thermo()->stateSize());
    gas.thermo()->saveState(cantera_state);

    gas.set_condensed_moles({0.01});
    gas.restore_state(cantera_state);

    EXPECT_EQ(gas.condensed_moles(), std::vector<double>({0.0}));
    EXPECT_FALSE(gas.has_condensed_phases());
}

TEST_F(CondensedGasTests, RestoreStateOfWrongLengthThrows) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)", "C(gr)"});
    std::vector<double> bad(gas.thermo()->stateSize() + 1);
    EXPECT_THROW(gas.restore_state(bad), std::invalid_argument);
}

TEST_F(CondensedGasTests, GasOnlyStateVectorLengthUnchanged) {
    EXPECT_EQ(gas.save_state().size(), gas.thermo()->stateSize());
}

TEST_F(CondensedGasTests, CopiesShareTheCondensedSetAndCloneDoesNot) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)"});
    gas.set_condensed_moles({0.01});

    Gas shallow = gas;          // shares the set, like the Solution
    Gas deep = gas.clone();     // independent copy

    EXPECT_EQ(deep.condensed_moles(), std::vector<double>({0.01}));

    gas.set_condensed_moles({0.05});
    EXPECT_EQ(shallow.condensed_moles(), std::vector<double>({0.05}));
    EXPECT_EQ(deep.condensed_moles(), std::vector<double>({0.01}));

    deep.set_condensed_moles({0.2});
    EXPECT_EQ(gas.condensed_moles(), std::vector<double>({0.05}));
}

// ---- Gas-only equilibrium entry points ----

class GasOnlyEquilibriumTests : public ::testing::Test {
protected:
    GasOnlyEquilibriumTests() : gas("h2o2.yaml", "ohmech", GasChemistry::EQUILIBRIUM) {
        setup_defaults();
        gas.set_state_TPX(3500.0, 70.0 * Cantera::OneAtm, "H2O:0.9, H2:0.05, O2:0.03, OH:0.02");
    }

    Gas gas;
};

TEST_F(GasOnlyEquilibriumTests, EquilibrateHPMatchesCantera) {
    Gas reference = gas.clone();

    const double H = gas.enthalpy_mass();
    const double P = gas.pressure();

    gas.equilibrate_HP(H, P);

    reference.thermo()->setState_HP(H, P);
    reference.thermo()->equilibrate("HP", "gibbs");

    EXPECT_DOUBLE_EQ(gas.temperature(), reference.temperature());
    EXPECT_DOUBLE_EQ(gas.pressure(), reference.pressure());
    EXPECT_DOUBLE_EQ(gas.molecular_weight(), reference.molecular_weight());
}

TEST_F(GasOnlyEquilibriumTests, EquilibrateTPAndSPMatchCantera) {
    Gas reference = gas.clone();

    gas.equilibrate_TP(3000.0, 10.0 * Cantera::OneAtm);
    reference.thermo()->setState_TP(3000.0, 10.0 * Cantera::OneAtm);
    reference.thermo()->equilibrate("TP", "gibbs");
    EXPECT_DOUBLE_EQ(gas.enthalpy_mass(), reference.enthalpy_mass());

    const double S = gas.entropy_mass();
    gas.equilibrate_SP(S, Cantera::OneAtm);
    reference.thermo()->setState_SP(S, Cantera::OneAtm);
    reference.thermo()->equilibrate("SP", "gibbs");
    EXPECT_DOUBLE_EQ(gas.temperature(), reference.temperature());
}

TEST_F(GasOnlyEquilibriumTests, SetElementMolesReproducesTheElementAmounts) {
    const Eigen::ArrayXd target = gas.element_moles();

    gas.set_element_moles(target, 1000.0, Cantera::OneAtm);

    EXPECT_DOUBLE_EQ(gas.temperature(), 1000.0);
    const Eigen::ArrayXd result = gas.element_moles();
    ASSERT_EQ(result.size(), target.size());
    // Only the element ratios are preserved: the composition is renormalized to one kg of gas.
    const double scale = result(0) / target(0);
    for (long m = 0; m < target.size(); m++) {
        EXPECT_NEAR(result(m), scale * target(m), 1e-9 * std::abs(result(m)) + 1e-12);
    }
}

TEST_F(GasOnlyEquilibriumTests, SetElementMolesPicksTheDiatomicBasisSpecies) {
    Eigen::ArrayXd amounts = Eigen::ArrayXd::Zero(static_cast<long>(gas.thermo()->nElements()));
    const std::vector<std::string> elements = gas.element_names();
    for (size_t m = 0; m < elements.size(); m++) {
        if (elements[m] == "H") amounts(static_cast<long>(m)) = 2.0;
        if (elements[m] == "O") amounts(static_cast<long>(m)) = 1.0;
    }

    gas.set_element_moles(amounts, 300.0, Cantera::OneAtm);

    const std::vector<double> moles = gas.mole_fractions();
    const std::vector<std::string> names = gas.species_names();
    for (size_t j = 0; j < names.size(); j++) {
        if (names[j] == "H2") EXPECT_NEAR(moles[j], 2.0 / 3.0, 1e-12);
        else if (names[j] == "O2") EXPECT_NEAR(moles[j], 1.0 / 3.0, 1e-12);
        else EXPECT_NEAR(moles[j], 0.0, 1e-12) << names[j];
    }
}

TEST_F(GasOnlyEquilibriumTests, SetElementMolesRejectsWrongSize) {
    Eigen::ArrayXd amounts = Eigen::ArrayXd::Zero(1);
    EXPECT_THROW(gas.set_element_moles(amounts, 300.0, Cantera::OneAtm), std::invalid_argument);
}

TEST_F(GasOnlyEquilibriumTests, MixtureOverloadsMatchTheThermoPhaseOverloads) {
    gas.equilibrate_HP(gas.enthalpy_mass(), gas.pressure());

    const ExpansionProperties from_gas = get_thermo_equilibrium_properties(gas);
    const ExpansionProperties from_thermo = get_thermo_equilibrium_properties(*gas.thermo());

    EXPECT_DOUBLE_EQ(from_gas.gamma_s, from_thermo.gamma_s);
    EXPECT_DOUBLE_EQ(from_gas.spec_heat_p, from_thermo.spec_heat_p);
    EXPECT_DOUBLE_EQ(from_gas.dlogV_dlogT_P, from_thermo.dlogV_dlogT_P);
    EXPECT_DOUBLE_EQ(from_gas.dlogV_dlogP_T, from_thermo.dlogV_dlogP_T);
    EXPECT_DOUBLE_EQ(get_equilibrium_gamma(gas), get_equilibrium_gamma(*gas.thermo()));

    // The mixture fields are filled for a gas-only state.
    EXPECT_NEAR(from_gas.gas_moles, 1.0 / gas.molecular_weight(), 1e-15);
    EXPECT_DOUBLE_EQ(from_gas.total_moles, from_gas.gas_moles);
    EXPECT_NEAR(from_gas.density, gas.density(), 1e-9 * gas.density());
    EXPECT_NEAR(from_gas.speed_of_sound, gas.speed_of_sound(), 1e-9 * gas.speed_of_sound());
    EXPECT_DOUBLE_EQ(from_gas.frozen_spec_heat_p, gas.cp_mass());
    EXPECT_DOUBLE_EQ(from_gas.frozen_gamma, gas.cp_mass() / gas.cv_mass());
    EXPECT_FALSE(from_gas.pinned_transition);
    EXPECT_GT(from_gas.spec_heat_v, 0.0);
    EXPECT_LT(from_gas.spec_heat_v, from_gas.spec_heat_p);
}

TEST_F(GasOnlyEquilibriumTests, FrozenPropertiesUseFixedComposition) {
    const ExpansionProperties props = get_frozen_properties(gas);

    EXPECT_DOUBLE_EQ(props.dlogV_dlogT_P, 1.0);
    EXPECT_DOUBLE_EQ(props.dlogV_dlogP_T, -1.0);
    EXPECT_DOUBLE_EQ(props.spec_heat_p, gas.cp_mass());
    EXPECT_DOUBLE_EQ(props.spec_heat_v, gas.cv_mass());
    EXPECT_DOUBLE_EQ(props.gamma_s, gas.cp_mass() / gas.cv_mass());
    EXPECT_DOUBLE_EQ(props.frozen_gamma, props.gamma_s);
    EXPECT_NEAR(props.density, gas.density(), 1e-9 * gas.density());
}

TEST_F(GasOnlyEquilibriumTests, DerivativesMatchTheThermoPhaseOverload) {
    const EquilibriumDerivatives from_gas = get_thermo_equilibrium_derivatives(gas);
    const EquilibriumDerivatives from_thermo = get_thermo_equilibrium_derivatives(*gas.thermo());

    EXPECT_DOUBLE_EQ(from_gas.dlogn_dlogT_P, from_thermo.dlogn_dlogT_P);
    EXPECT_DOUBLE_EQ(from_gas.dlogn_dlogP_T, from_thermo.dlogn_dlogP_T);
    EXPECT_EQ(from_gas.dn_condensed_dlogT_P.size(), 0);
    EXPECT_EQ(from_gas.dn_condensed_dlogP_T.size(), 0);
    EXPECT_FALSE(from_gas.pinned_transition);
}

// ---- Paths that wait on later work packages ----

TEST_F(CondensedGasTests, UnsupportedEquilibriumRequestsAreRejected) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)"});

    EXPECT_THROW(gas.equilibrate("TP", "vcs"), std::invalid_argument);
    EXPECT_THROW(gas.equilibrate("UV"), NotImplementedError);
    EXPECT_THROW(gas.equilibrate("TV"), std::invalid_argument);
    EXPECT_THROW(gas.set_state_UV(1.0e6, 1.0), NotImplementedError);
}

TEST_F(CondensedGasTests, SetElementMolesClearsTheCondensedState) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)"});
    gas.set_condensed_moles({0.01});

    const Eigen::ArrayXd amounts = gas.element_moles();
    gas.set_element_moles(amounts, 400.0, Cantera::OneAtm);

    EXPECT_EQ(gas.condensed_moles(), std::vector<double>({0.0}));
    EXPECT_FALSE(gas.at_phase_transition());
    EXPECT_DOUBLE_EQ(gas.temperature(), 400.0);
}

TEST_F(CondensedGasTests, EquilibriumPropertyOverloadsThrowWithCondensedPhasesPresent) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)"});

    // Candidates alone are harmless: the state is still a gas.
    EXPECT_NO_THROW(get_thermo_equilibrium_derivatives(gas));

    gas.set_condensed_moles({0.01});
    EXPECT_THROW(get_thermo_equilibrium_derivatives(gas), NotImplementedError);
    EXPECT_THROW(get_thermo_equilibrium_properties(gas), NotImplementedError);
    EXPECT_THROW(get_frozen_properties(gas), NotImplementedError);
    EXPECT_THROW(get_equilibrium_gamma(gas), NotImplementedError);
}

TEST_F(CondensedGasTests, ThermoArrayStoresCondensedAmountsPerLocation) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)", "C(gr)"});
    gas.set_condensed_moles({0.004, 0.001});

    ThermoArray states(gas, {2, 3});
    EXPECT_EQ(states.size(), 6);
    EXPECT_EQ(states.num_condensed(), 2u);
    EXPECT_EQ(states.condensed_species_names(), std::vector<std::string>({"H2O(L)", "C(gr)"}));

    // Every entry starts from the amounts of the Gas it was built from.
    const size_t cantera_size = gas.thermo()->stateSize();
    for (int loc = 0; loc < states.size(); loc++) {
        EXPECT_EQ(states.get_condensed_moles(loc), std::vector<double>({0.004, 0.001}));
        EXPECT_EQ(states.get_state(loc).size(), cantera_size + 2);
    }

    std::vector<double> state = states.get_state(3);
    state[cantera_size] = 0.02;
    state[cantera_size + 1] = 0.0;
    states.set_state(3, state);
    EXPECT_EQ(states.get_condensed_moles(3), std::vector<double>({0.02, 0.0}));
    EXPECT_EQ(states.get_condensed_moles(2), std::vector<double>({0.004, 0.001}));
    EXPECT_EQ(states.get_state(3), state);

    // A bare Cantera state vector empties that entry.
    states.set_state(3, std::vector<double>(state.begin(), state.begin() + static_cast<long>(cantera_size)));
    EXPECT_EQ(states.get_condensed_moles(3), std::vector<double>({0.0, 0.0}));

    EXPECT_THROW(states.set_state(0, std::vector<double>(cantera_size + 1)), std::invalid_argument);
    EXPECT_THROW(states.set_condensed_moles(0, {0.1}), std::invalid_argument);
}

TEST_F(CondensedGasTests, ThermoArrayFromGasWithoutCandidatesBehavesAsBefore) {
    ThermoArray states(gas, {2, 3});
    EXPECT_EQ(states.num_condensed(), 0u);
    EXPECT_TRUE(states.condensed_species_names().empty());
    EXPECT_NO_THROW(states.get_state(0));
}

TEST_F(CondensedGasTests, CombustorFromReactantGasesThrows) {
    Gas fuel = make_gas("fuel", {"H2"});
    Gas oxidizer = make_gas("oxidizer", {"O2"});
    Combustor combustor(gas, fuel, oxidizer);

    Eigen::ArrayXd pressures(1);
    pressures << Cantera::OneAtm;
    Eigen::ArrayXd ratios(1);
    ratios << 8.0;

    EXPECT_THROW(combustor.solve(pressures, ratios), NotImplementedError);
}

// ---------------------------------------------------------------------------------------------
// Multiphase equilibrium against CEA (RP-1311) on the NASA9 data Goddard and CEA share.
// ---------------------------------------------------------------------------------------------

namespace {

/** CEA-derived NASA9 data converted to Cantera YAML (work package E). */
constexpr const char* NASA9_GAS = "nasa9_gas.yaml";
constexpr const char* NASA9_CONDENSED = "nasa9_condensed.yaml";
constexpr const char* NASA9_REACTANTS = "nasa9_reactants.yaml";

constexpr double BAR = 1.0e5;

/** Product gas spanning `elements`, frozen so no accidental call needs the WP-A derivatives. */
Gas make_products(const std::vector<std::string>& elements) {
    setup_defaults();
    return Gas::create_from_elements(NASA9_GAS, "products", elements, GasChemistry::FROZEN);
}

/** Reactant stream: a phase of the named species from the reactant database at 298.15 K. */
Gas make_reactants(const std::vector<std::string>& species, const std::string& mass_fractions) {
    setup_defaults();
    Gas reactants = Gas::create_from_species(NASA9_REACTANTS, "reactants", species,
                                             GasChemistry::FROZEN);
    reactants.set_state_TPY(298.15, Cantera::OneAtm, mass_fractions);
    return reactants;
}

/** Element amounts of `from` [kmol/kg] re-indexed onto the element order of `to`. */
Eigen::ArrayXd map_elements(const Gas& from, const Gas& to) {
    const Eigen::ArrayXd source = from.element_moles();
    const std::vector<std::string> source_names = from.element_names();
    const std::vector<std::string> target_names = to.element_names();

    Eigen::ArrayXd target = Eigen::ArrayXd::Zero(static_cast<long>(target_names.size()));
    for (size_t m = 0; m < source_names.size(); m++) {
        if (source(static_cast<long>(m)) == 0.0) continue;
        auto found = std::find(target_names.begin(), target_names.end(), source_names[m]);
        EXPECT_NE(found, target_names.end()) << "element " << source_names[m] << " missing";
        if (found == target_names.end()) continue;
        target(found - target_names.begin()) += source(static_cast<long>(m));
    }
    return target;
}

/** Mole fraction of a condensed species among *all* species, CEA's convention. */
double condensed_mole_fraction(const Gas& gas, const std::string& name) {
    const std::vector<std::string> names = gas.condensed_species_names();
    auto found = std::find(names.begin(), names.end(), name);
    if (found == names.end()) return 0.0;
    return gas.condensed_moles()[static_cast<size_t>(found - names.begin())]
        * gas.mixture_molecular_weight();
}

/** Mass and element balances of a converged multiphase state. */
void expect_balances(const Gas& gas, const Eigen::ArrayXd& elements_before) {
    double mixture_mass = 0.0;
    for (double y : gas.mixture_mass_fractions()) mixture_mass += y;
    EXPECT_NEAR(mixture_mass, 1.0, 1e-10);
    EXPECT_GT(gas.gas_mass_fraction(), 0.0);
    EXPECT_LE(gas.gas_mass_fraction(), 1.0);

    const Eigen::ArrayXd elements_after = gas.element_moles();
    ASSERT_EQ(elements_after.size(), elements_before.size());
    for (long m = 0; m < elements_before.size(); m++) {
        EXPECT_NEAR(elements_after(m), elements_before(m),
                    1e-7 * std::abs(elements_before(m)) + 1e-14)
            << "element index " << m;
    }
}

/** CH4/O2 products with every compatible condensed candidate attached. */
Gas make_methane_oxygen_products() {
    Gas products = make_products({"C", "H", "O"});
    products.add_all_condensed_species(NASA9_CONDENSED);
    return products;
}

/** Element amounts per kg for a CH4/O2 mixture at the given oxidizer-to-fuel mass ratio. */
Eigen::ArrayXd methane_oxygen_elements(const Gas& products, double OF_ratio) {
    const double fuel_fraction = 1.0 / (1.0 + OF_ratio);
    Gas reactants = make_reactants({"CH4", "O2"},
        "CH4:" + std::to_string(fuel_fraction) + ", O2:" + std::to_string(1.0 - fuel_fraction));
    return map_elements(reactants, products);
}

} // namespace

class MultiphaseTPTests : public ::testing::Test {};

TEST_F(MultiphaseTPTests, MethaneOxygenGraphiteAt1000K) {
    Gas products = make_methane_oxygen_products();
    const Eigen::ArrayXd elements = methane_oxygen_elements(products, 0.5);

    products.set_element_moles(elements, 1000.0, BAR);
    const Eigen::ArrayXd before = products.element_moles();
    products.equilibrate_TP(1000.0, BAR);

    EXPECT_NEAR(condensed_mole_fraction(products, "C(gr)"), 0.18459, 5e-3);
    EXPECT_NEAR(products.molecular_weight(), 10.676, 0.05);
    EXPECT_DOUBLE_EQ(products.temperature(), 1000.0);
    EXPECT_NEAR(products.pressure(), BAR, 1e-6 * BAR);
    expect_balances(products, before);
}

TEST_F(MultiphaseTPTests, MethaneOxygenGraphiteAt1500K) {
    Gas products = make_methane_oxygen_products();
    const Eigen::ArrayXd elements = methane_oxygen_elements(products, 0.5);

    products.set_element_moles(elements, 1500.0, BAR);
    const Eigen::ArrayXd before = products.element_moles();
    products.equilibrate_TP(1500.0, BAR);

    EXPECT_NEAR(condensed_mole_fraction(products, "C(gr)"), 0.16554, 5e-3);
    EXPECT_NEAR(products.molecular_weight(), 9.639, 0.05);
    expect_balances(products, before);
}

TEST_F(MultiphaseTPTests, MethaneOxygenGraphiteAtFiftyBar) {
    Gas products = make_methane_oxygen_products();
    const Eigen::ArrayXd elements = methane_oxygen_elements(products, 1.0);

    products.set_element_moles(elements, 1000.0, 50.0 * BAR);
    const Eigen::ArrayXd before = products.element_moles();
    products.equilibrate_TP(1000.0, 50.0 * BAR);

    EXPECT_NEAR(condensed_mole_fraction(products, "C(gr)"), 0.09619, 5e-3);
    EXPECT_NEAR(products.molecular_weight(), 18.087, 0.05);
    expect_balances(products, before);
}

TEST_F(MultiphaseTPTests, MethaneOxygenHasNoGraphiteWhenOxygenRich) {
    Gas products = make_methane_oxygen_products();
    const Eigen::ArrayXd elements = methane_oxygen_elements(products, 2.0);

    products.set_element_moles(elements, 1500.0, BAR);
    const Eigen::ArrayXd before = products.element_moles();
    products.equilibrate_TP(1500.0, BAR);

    EXPECT_DOUBLE_EQ(condensed_mole_fraction(products, "C(gr)"), 0.0);
    EXPECT_FALSE(products.has_condensed_phases());
    EXPECT_DOUBLE_EQ(products.gas_mass_fraction(), 1.0);
    expect_balances(products, before);
}

// ---- HP: constant enthalpy and pressure ----

class MultiphaseHPTests : public ::testing::Test {};

TEST_F(MultiphaseHPTests, MethaneOxygenAtUnitMixtureRatioDepositsGraphite) {
    Gas products = make_methane_oxygen_products();
    Gas reactants = make_reactants({"CH4", "O2"}, "CH4:0.5, O2:0.5");
    const Eigen::ArrayXd elements = map_elements(reactants, products);
    const double enthalpy = reactants.enthalpy_mass();

    products.set_element_moles(elements, 2000.0, BAR);
    const Eigen::ArrayXd before = products.element_moles();
    products.equilibrate_HP(enthalpy, BAR);

    EXPECT_NEAR(products.temperature(), 1039.24, 2.0);
    EXPECT_NEAR(condensed_mole_fraction(products, "C(gr)"), 0.03137, 3e-3);
    EXPECT_NEAR(products.enthalpy_mass(), enthalpy, 1e-6 * std::abs(enthalpy));
    EXPECT_FALSE(products.at_phase_transition());
    EXPECT_GT(products.last_equilibrium_solve_count(), 0);
    expect_balances(products, before);
}

TEST_F(MultiphaseHPTests, MethaneOxygenAtMixtureRatioThreeHasNoGraphite) {
    Gas products = make_methane_oxygen_products();
    Gas reactants = make_reactants({"CH4", "O2"}, "CH4:0.25, O2:0.75");
    const Eigen::ArrayXd elements = map_elements(reactants, products);
    const double enthalpy = reactants.enthalpy_mass();

    products.set_element_moles(elements, 2000.0, BAR);
    const Eigen::ArrayXd before = products.element_moles();
    products.equilibrate_HP(enthalpy, BAR);

    EXPECT_NEAR(products.temperature(), 3025.45, 3.0);
    EXPECT_FALSE(products.has_condensed_phases());
    expect_balances(products, before);
}

TEST_F(MultiphaseHPTests, MethaneOxygenAtTwentyBar) {
    Gas products = make_methane_oxygen_products();
    Gas reactants = make_reactants({"CH4", "O2"}, "CH4:0.25, O2:0.75");
    const Eigen::ArrayXd elements = map_elements(reactants, products);

    products.set_element_moles(elements, 2000.0, 20.0 * BAR);
    products.equilibrate_HP(reactants.enthalpy_mass(), 20.0 * BAR);

    EXPECT_NEAR(products.temperature(), 3398.53, 3.0);
}

// ---- The water dew point (RP-1311 example 14 style) ----

namespace {

/** H2/O2 products at 100 : 60 moles with the water candidates attached. */
Gas make_hydrogen_oxygen_products() {
    Gas products = make_products({"H", "O"});
    products.add_all_condensed_species(NASA9_CONDENSED);
    return products;
}

Eigen::ArrayXd hydrogen_oxygen_elements(const Gas& products) {
    setup_defaults();
    Gas reactants = Gas::create_from_species(NASA9_REACTANTS, "reactants", {"H2", "O2"},
                                             GasChemistry::FROZEN);
    reactants.set_state_TPX(298.15, Cantera::OneAtm, "H2:100, O2:60");
    return map_elements(reactants, products);
}

} // namespace

class DewPointTests : public ::testing::Test {
protected:
    DewPointTests() : products(make_hydrogen_oxygen_products()),
                      elements(hydrogen_oxygen_elements(products)),
                      pressure(0.05 * Cantera::OneAtm) {}

    /** Solve at `T` from a cold start and return the mole fraction of `name` among all species. */
    double solve(double T, const std::string& name) {
        products.set_element_moles(elements, T, pressure);
        const Eigen::ArrayXd before = products.element_moles();
        products.equilibrate_TP(T, pressure);
        expect_balances(products, before);
        return condensed_mole_fraction(products, name);
    }

    Gas products;
    Eigen::ArrayXd elements;
    double pressure;
};

TEST_F(DewPointTests, NothingCondensesAboveTheDewPoint) {
    products.set_element_moles(elements, 305.0, pressure);
    products.equilibrate_TP(305.0, pressure);
    EXPECT_FALSE(products.has_condensed_phases());
}

TEST_F(DewPointTests, LiquidWaterJustBelowTheDewPoint) {
    // Within about a degree of the dew point the liquid fraction climbs by roughly 0.25 per kelvin,
    // so this station is far more sensitive to the thermodynamic data than any other; CEA reports
    // 0.2488 here. What it really exercises is convergence, see `NearTheDewPointNeedsManySteps`.
    const double liquid = solve(304.0, "H2O(L)");
    EXPECT_GT(liquid, 0.1);
    EXPECT_LT(liquid, 0.5);
}

TEST_F(DewPointTests, NearTheDewPointTheStepLimitDoesNotChangeTheAnswer) {
    // A barely-present condensed phase converges slowly, and Cantera's default step limit of 1000
    // is not enough from every starting guess. The restart ladder covers that: the answer is the
    // same with either limit, it just takes another rung.
    const double reference = solve(304.0, "H2O(L)");

    products.equilibrium_options.max_steps = 1000;
    const double restricted = solve(304.0, "H2O(L)");
    EXPECT_NEAR(restricted, reference, 1e-6);
}

TEST_F(DewPointTests, LiquidWaterAt300K) {
    EXPECT_NEAR(solve(300.0, "H2O(L)"), 0.6995, 0.01);
}

TEST_F(DewPointTests, IceBelowFreezing) {
    EXPECT_NEAR(solve(270.0, "H2O(cr)"), 0.8998, 5e-3);
    EXPECT_NEAR(solve(250.0, "H2O(cr)"), 0.9077, 5e-3);
}

TEST_F(DewPointTests, ExactlyAtTheFreezingPointEitherPolymorphIsAccepted) {
    products.set_element_moles(elements, 273.15, pressure);
    ASSERT_NO_THROW(products.equilibrate_TP(273.15, pressure));

    const double water = condensed_mole_fraction(products, "H2O(cr)")
        + condensed_mole_fraction(products, "H2O(L)");
    EXPECT_GT(water, 0.85);
}

// ---- Aluminized ammonium perchlorate: the Al2O3 melting transition ----

namespace {

/** Products of an ammonium-perchlorate / aluminium propellant, with the alumina candidates. */
Gas make_propellant_products() {
    Gas products = make_products({"N", "H", "Cl", "O", "Al"});
    products.add_condensed_species(NASA9_CONDENSED, {"AL2O3(a)", "AL2O3(L)"});
    return products;
}

/**
 * Solve the propellant at constant enthalpy and pressure for a given aluminium mass fraction.
 * @return the reactant enthalpy [J/kg] the products were equilibrated to.
 */
double solve_propellant(Gas& products, double aluminium_fraction, double P, double T_guess) {
    Gas reactants = make_reactants({"NH4CLO4(I)", "AL(cr)"},
        std::format("NH4CLO4(I):{:.12g}, AL(cr):{:.12g}", 1.0 - aluminium_fraction,
                    aluminium_fraction));

    const double enthalpy = reactants.enthalpy_mass();
    products.set_element_moles(map_elements(reactants, products), T_guess, P);
    products.equilibrate_HP(enthalpy, P);
    return enthalpy;
}

} // namespace

class PropellantTests : public ::testing::Test {};

TEST_F(PropellantTests, MoltenAluminaAtThreePressures) {
    Gas products = make_propellant_products();

    double previous_temperature = 1e9;
    for (double P : {34.47 * BAR, 3.447 * BAR, 0.3447 * BAR}) {
        solve_propellant(products, 0.2, P, 3000.0);

        EXPECT_GT(condensed_mole_fraction(products, "AL2O3(L)"), 0.0) << "P = " << P;
        EXPECT_DOUBLE_EQ(condensed_mole_fraction(products, "AL2O3(a)"), 0.0) << "P = " << P;
        EXPECT_GT(products.temperature(), 2327.0) << "P = " << P;
        // Dropping the pressure shifts the equilibrium towards dissociation, cooling the flame.
        EXPECT_LT(products.temperature(), previous_temperature);
        previous_temperature = products.temperature();
    }
}

TEST_F(PropellantTests, ColdAndWarmStartsAgree) {
    Gas warm = make_propellant_products();
    solve_propellant(warm, 0.2, 34.47 * BAR, 3000.0);
    const double from_warm = warm.temperature();

    Gas cold = make_propellant_products();
    solve_propellant(cold, 0.2, 34.47 * BAR, 1200.0);

    EXPECT_NEAR(cold.temperature(), from_warm, 1.0);
    EXPECT_NEAR(cold.gas_mass_fraction(), warm.gas_mass_fraction(), 1e-4);
}

TEST_F(PropellantTests, LittleAluminiumLeavesSolidAlumina) {
    Gas products = make_propellant_products();
    solve_propellant(products, 0.05, 34.47 * BAR, 3000.0);

    EXPECT_LT(products.temperature(), 2327.0);
    EXPECT_GT(condensed_mole_fraction(products, "AL2O3(a)"), 0.0);
    EXPECT_DOUBLE_EQ(condensed_mole_fraction(products, "AL2O3(L)"), 0.0);
    EXPECT_FALSE(products.at_phase_transition());
}

TEST_F(PropellantTests, PinnedAtTheAluminaMeltingPoint) {
    // Between the two cases above lies a band of aluminium fractions whose flame enthalpy falls
    // inside the heat of fusion of alumina. Bisect on the aluminium fraction to land in it.
    Gas products = make_propellant_products();
    const double P = 34.47 * BAR;

    double solid_side = 0.05;   // ends below the melting point
    double liquid_side = 0.2;   // ends above it
    double target_enthalpy = 0.0;
    for (int iteration = 0; iteration < 40 && !products.at_phase_transition(); iteration++) {
        const double fraction = 0.5 * (solid_side + liquid_side);
        target_enthalpy = solve_propellant(products, fraction, P, 2500.0);
        if (products.at_phase_transition()) break;
        if (products.temperature() < 2327.0) {
            solid_side = fraction;
        } else {
            liquid_side = fraction;
        }
    }

    ASSERT_TRUE(products.at_phase_transition())
        << "no aluminium fraction in [" << solid_side << ", " << liquid_side << "] pinned the melt";
    EXPECT_DOUBLE_EQ(products.temperature(), 2327.0);

    const std::vector<std::string> names = products.condensed_species_names();
    const std::vector<double> moles = products.condensed_moles();
    const size_t solid = static_cast<size_t>(
        std::find(names.begin(), names.end(), "AL2O3(a)") - names.begin());
    const size_t liquid = static_cast<size_t>(
        std::find(names.begin(), names.end(), "AL2O3(L)") - names.begin());
    EXPECT_GT(moles[solid], 0.0);
    EXPECT_GT(moles[liquid], 0.0);

    // Both polymorphs are reported as the pinned pair, lower-temperature one first.
    const std::pair<long, long> pinned = products.pinned_polymorphs();
    EXPECT_GE(pinned.first, 0);
    EXPECT_GE(pinned.second, 0);

    // The split is what makes the mixture enthalpy match the target exactly.
    EXPECT_NEAR(products.enthalpy_mass(), target_enthalpy, 1e-9 * std::abs(target_enthalpy));
}

// ---- RP-1311 example 13: N2H4/Be with H2O2, whose beryllia passes two transitions ----

namespace {

/**
 * Products of the RP-1311 example 13 propellant.
 *
 * Every N/H/Be/O species of the database except gaseous Be(OH)2, which CEA's product list for this
 * example does not contain. Including it moves 3.8 % of the beryllium out of the condensate and
 * costs 16 K of flame temperature, which would swamp the solver differences this case measures.
 */
Gas make_beryllium_products() {
    const Gas all = make_products({"N", "H", "Be", "O"});
    std::vector<std::string> species;
    for (const std::string& name : all.species_names()) {
        if (name != "Be(OH)2") species.push_back(name);
    }

    setup_defaults();
    Gas products = Gas::create_from_species(NASA9_GAS, "products", species, GasChemistry::FROZEN);
    products.add_all_condensed_species(NASA9_CONDENSED);
    return products;
}

/** 67 % fuel (80 % N2H4(L), 20 % Be(a)) and 33 % H2O2(L) by mass, all at 298.15 K. */
Gas make_beryllium_reactants() {
    return make_reactants({"N2H4(L)", "Be(a)", "H2O2(L)"},
                          "N2H4(L):0.536, Be(a):0.134, H2O2(L):0.33");
}

} // namespace

class BerylliumRocketTests : public ::testing::Test {
protected:
    BerylliumRocketTests()
        : products(make_beryllium_products()), reactants(make_beryllium_reactants()) {}

    /** Equilibrate the chamber at constant enthalpy and return its entropy [J/(kg.K)]. */
    double solve_chamber() {
        const double P = 206.8419 * BAR;
        products.set_element_moles(map_elements(reactants, products), 3000.0, P);
        products.equilibrate_HP(reactants.enthalpy_mass(), P);
        return products.entropy_mass();
    }

    Gas products;
    Gas reactants;
};

TEST_F(BerylliumRocketTests, ChamberMatchesCea) {
    const Eigen::ArrayXd before = [&] {
        const double P = 206.8419 * BAR;
        products.set_element_moles(map_elements(reactants, products), 3000.0, P);
        return products.element_moles();
    }();

    products.equilibrate_HP(reactants.enthalpy_mass(), 206.8419 * BAR);

    EXPECT_NEAR(products.temperature(), 3018.89, 2.0);
    EXPECT_NEAR(products.molecular_weight(), 16.6222, 5e-3);
    EXPECT_NEAR(condensed_mole_fraction(products, "BeO(L)"), 0.19796, 2e-3);
    EXPECT_FALSE(products.at_phase_transition());
    expect_balances(products, before);
}

TEST_F(BerylliumRocketTests, ExpansionFollowsTheBerylliaPolymorphs) {
    const double entropy = solve_chamber();

    struct Station {
        double pressure_bar;
        double temperature;
        const char* species;
        double mole_fraction;
        bool pinned;
    };
    // CEA's stations for this expansion. The first two sit exactly on the BeO(b)/BeO(L) melting
    // point, where the entropy falls inside the latent heat and the temperature stops moving.
    const Station stations[] = {
        {127.2265, 2851.0,  "BeO(L)", 0.18656, true},
        {68.9473,  2851.0,  "BeO(L)", 0.04510, true},
        {20.6842,  2453.58, "BeO(b)", 0.19862, false},
        {6.8947,   2066.65, "BeO(a)", 0.19886, false},
    };

    for (const Station& station : stations) {
        products.equilibrate_SP(entropy, station.pressure_bar * BAR);

        EXPECT_NEAR(products.temperature(), station.temperature, 2.0)
            << "at " << station.pressure_bar << " bar";
        EXPECT_NEAR(condensed_mole_fraction(products, station.species), station.mole_fraction, 2e-3)
            << "at " << station.pressure_bar << " bar";
        EXPECT_EQ(products.at_phase_transition(), station.pinned)
            << "at " << station.pressure_bar << " bar";
        EXPECT_NEAR(products.entropy_mass(), entropy, 1e-8 * std::abs(entropy));
    }
}

TEST_F(BerylliumRocketTests, PinnedStationsSplitTheBerylliaGroup) {
    const double entropy = solve_chamber();
    products.equilibrate_SP(entropy, 127.2265 * BAR);

    ASSERT_TRUE(products.at_phase_transition());
    EXPECT_NEAR(products.temperature(), 2851.0, 1e-6);
    EXPECT_NEAR(condensed_mole_fraction(products, "BeO(L)"), 0.18656, 2e-3);
    EXPECT_NEAR(condensed_mole_fraction(products, "BeO(b)"), 0.01170, 2e-3);
    EXPECT_NEAR(products.molecular_weight(), 16.6427, 5e-3);

    products.equilibrate_SP(entropy, 68.9473 * BAR);
    EXPECT_NEAR(products.temperature(), 2851.0, 1e-6);
    EXPECT_NEAR(condensed_mole_fraction(products, "BeO(L)"), 0.04510, 2e-3);
    EXPECT_NEAR(condensed_mole_fraction(products, "BeO(b)"), 0.15288, 2e-3);
}

TEST_F(BerylliumRocketTests, DeepExpansionStations) {
    const double entropy = solve_chamber();

    products.equilibrate_SP(entropy, 0.6895 * BAR);
    EXPECT_NEAR(products.temperature(), 1395.67, 2.0);

    products.equilibrate_SP(entropy, 0.2068 * BAR);
    EXPECT_NEAR(products.temperature(), 1119.74, 2.0);
}

TEST_F(BerylliumRocketTests, WarmStartsCostFewerSolvesThanColdStarts) {
    const double entropy = solve_chamber();
    const std::vector<double> pressures = {127.2265, 68.9473, 20.6842, 6.8947, 0.6895, 0.2068};

    int warm_total = 0;
    for (double pressure_bar : pressures) {
        products.equilibrate_SP(entropy, pressure_bar * BAR);   // continues from the last station
        warm_total += products.last_equilibrium_solve_count();
    }

    int cold_total = 0;
    for (double pressure_bar : pressures) {
        // Reset the temperature guess to the default the solver falls back on.
        products.set_state_TP(products.equilibrium_options.T_default, pressure_bar * BAR);
        products.equilibrate_SP(entropy, pressure_bar * BAR);
        cold_total += products.last_equilibrium_solve_count();
    }

    EXPECT_LT(warm_total, cold_total);
}

// ---- State round trips, frozen expansion and batch operations ----

class CondensedStateTests : public ::testing::Test {};

TEST_F(CondensedStateTests, PinnedStateSurvivesSaveRestoreAndClone) {
    Gas products = make_propellant_products();
    const double P = 34.47 * BAR;

    double solid_side = 0.05;
    double liquid_side = 0.2;
    for (int iteration = 0; iteration < 40 && !products.at_phase_transition(); iteration++) {
        const double fraction = 0.5 * (solid_side + liquid_side);
        solve_propellant(products, fraction, P, 2500.0);
        if (products.at_phase_transition()) break;
        (products.temperature() < 2327.0 ? solid_side : liquid_side) = fraction;
    }
    ASSERT_TRUE(products.at_phase_transition());

    const std::vector<double> state = products.save_state();
    const std::vector<double> moles = products.condensed_moles();
    const std::pair<long, long> pinned = products.pinned_polymorphs();
    EXPECT_EQ(state.size(), products.thermo()->stateSize() + moles.size());

    // The pinned flag is not in the state vector; it is re-derived from the two coexisting
    // polymorphs, so a round trip restores it.
    products.set_state_TP(3000.0, P);
    products.set_condensed_moles(std::vector<double>(moles.size(), 0.0));
    products.clear_phase_transition();
    ASSERT_FALSE(products.at_phase_transition());

    products.restore_state(state);
    EXPECT_EQ(products.condensed_moles(), moles);
    EXPECT_TRUE(products.at_phase_transition());
    EXPECT_EQ(products.pinned_polymorphs(), pinned);

    Gas copy = products.clone();
    EXPECT_EQ(copy.condensed_moles(), moles);
    EXPECT_TRUE(copy.at_phase_transition());
    EXPECT_EQ(copy.pinned_polymorphs(), pinned);
    EXPECT_NEAR(copy.enthalpy_mass(), products.enthalpy_mass(),
                1e-12 * std::abs(products.enthalpy_mass()));

    copy.set_condensed_moles(std::vector<double>(moles.size(), 0.0));
    EXPECT_EQ(products.condensed_moles(), moles);
}

TEST_F(CondensedStateTests, PhaseTransitionCanBeSetByName) {
    Gas products = make_propellant_products();
    products.set_state_TP(2327.0, 34.47 * BAR);
    products.set_condensed_moles({0.004, 0.002});

    EXPECT_FALSE(products.at_phase_transition());
    products.set_phase_transition("AL2O3(a)", "AL2O3(L)");
    EXPECT_TRUE(products.at_phase_transition());
    EXPECT_EQ(products.pinned_polymorphs(), std::make_pair(0L, 1L));

    products.clear_phase_transition();
    EXPECT_FALSE(products.at_phase_transition());
    products.set_phase_transition(0L, 1L);
    EXPECT_TRUE(products.at_phase_transition());

    EXPECT_THROW(products.set_phase_transition("AL2O3(a)", "C(gr)"), std::invalid_argument);
    EXPECT_THROW(products.set_phase_transition(0L, 0L), std::invalid_argument);
}

TEST_F(CondensedStateTests, FrozenExpansionKeepsTheMixtureEnthalpy) {
    Gas products = make_propellant_products();
    solve_propellant(products, 0.2, 34.47 * BAR, 3000.0);
    ASSERT_TRUE(products.has_condensed_phases());

    const std::vector<double> moles = products.condensed_moles();
    const double entropy = products.entropy_mass();
    const double chamber_temperature = products.temperature();

    products.set_state_SP(entropy, 20.0 * BAR);
    EXPECT_NEAR(products.entropy_mass(), entropy, 1e-9 * std::abs(entropy));
    EXPECT_LT(products.temperature(), chamber_temperature);
    // Frozen: neither the gas composition nor the condensed amounts move.
    EXPECT_EQ(products.condensed_moles(), moles);

    const double enthalpy = products.enthalpy_mass();
    products.set_state_TP(2500.0, 20.0 * BAR);
    products.set_state_HP(enthalpy, 20.0 * BAR);
    EXPECT_NEAR(products.enthalpy_mass(), enthalpy, 1e-9 * std::abs(enthalpy));
}

TEST_F(CondensedStateTests, FrozenExpansionStopsAtACondensedTemperatureRange) {
    Gas products = make_propellant_products();
    solve_propellant(products, 0.2, 34.47 * BAR, 3000.0);
    ASSERT_GT(products.condensed_moles()[1], 0.0); // AL2O3(L), whose data starts at 2327 K

    // Expanding far enough drives the frozen temperature below the melting point of alumina, which
    // a frozen expansion cannot represent because the liquid is not allowed to solidify.
    EXPECT_THROW(products.set_state_SP(products.entropy_mass(), 0.001 * BAR), FmtError);
}

TEST_F(CondensedStateTests, ThermoArrayEquilibratesEachLocationLikeGas) {
    Gas gas = make_methane_oxygen_products();
    Gas reactants = make_reactants({"CH4", "O2"}, "CH4:0.5, O2:0.5");
    gas.set_element_moles(map_elements(reactants, gas), 2000.0, BAR);

    ThermoArray states(gas, {2});
    ASSERT_EQ(states.num_condensed(), 3u);
    states.equilibrate("HP", "gibbs", 1e-9, 20000);

    Gas reference = gas.clone();
    reference.equilibrate_HP(reference.enthalpy_mass(), reference.pressure());

    for (int loc = 0; loc < states.size(); loc++) {
        EXPECT_NEAR(states.temperature()(loc, 0), reference.temperature(), 1e-9)
            << "location " << loc;
        const std::vector<double> moles = states.get_condensed_moles(loc);
        ASSERT_EQ(moles.size(), reference.condensed_moles().size());
        for (size_t k = 0; k < moles.size(); k++) {
            EXPECT_NEAR(moles[k], reference.condensed_moles()[k], 1e-12) << "location " << loc;
        }
    }
}

TEST_F(CondensedStateTests, ThermoArrayRejectsTheVcsSolverWithCondensedSpecies) {
    Gas gas = make_methane_oxygen_products();
    Gas reactants = make_reactants({"CH4", "O2"}, "CH4:0.5, O2:0.5");
    gas.set_element_moles(map_elements(reactants, gas), 2000.0, BAR);

    ThermoArray states(gas, {1});
    EXPECT_THROW(states.equilibrate("TP", "vcs"), std::invalid_argument);
}
