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

    const double molar_mass = gas.molecular_weight(); // pure H2O gas
    const double w_gas = 1.0 - n_liquid * molar_mass;
    EXPECT_DOUBLE_EQ(gas.gas_mass_fraction(), w_gas);

    const double total_moles = w_gas / molar_mass + n_liquid;
    EXPECT_DOUBLE_EQ(gas.mixture_molecular_weight(), 1.0 / total_moles);
    // The gas-phase molecular weight is untouched.
    EXPECT_DOUBLE_EQ(gas.molecular_weight(), molar_mass);
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

TEST_F(CondensedGasTests, EquilibriumEntryPointsThrowWithCandidates) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)"});

    EXPECT_THROW(gas.equilibrate_TP(400.0, Cantera::OneAtm), NotImplementedError);
    EXPECT_THROW(gas.equilibrate_HP(gas.enthalpy_mass(), Cantera::OneAtm), NotImplementedError);
    EXPECT_THROW(gas.equilibrate_SP(gas.entropy_mass(), Cantera::OneAtm), NotImplementedError);

    Eigen::ArrayXd amounts = gas.element_moles();
    EXPECT_THROW(gas.set_element_moles(amounts, 400.0, Cantera::OneAtm), NotImplementedError);
}

TEST_F(CondensedGasTests, EquilibriumPropertyOverloadsAcceptCondensedPhases) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)"});

    // Candidates alone are harmless: the state is still a gas.
    EXPECT_NO_THROW(get_thermo_equilibrium_derivatives(gas));

    // With condensed phases present the system grows one unknown and one row per species (WP-A).
    gas.set_condensed_moles({0.01});
    const EquilibriumDerivatives derivs = get_thermo_equilibrium_derivatives(gas);
    EXPECT_EQ(derivs.dn_condensed_dlogT_P.size(), 1);
    EXPECT_EQ(derivs.dn_condensed_dlogP_T.size(), 1);
    EXPECT_FALSE(derivs.pinned_transition);

    EXPECT_NO_THROW(get_thermo_equilibrium_properties(gas));
    EXPECT_NO_THROW(get_frozen_properties(gas));
    EXPECT_NO_THROW(get_equilibrium_gamma(gas));
}

TEST_F(CondensedGasTests, ThermoArrayCarriesTheCandidatesButCannotStoreThem) {
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)", "C(gr)"});

    ThermoArray states(gas, {2, 3});
    EXPECT_EQ(states.size(), 6);
    EXPECT_EQ(states.num_condensed(), 2u);
    EXPECT_EQ(states.condensed_species_names(), std::vector<std::string>({"H2O(L)", "C(gr)"}));

    EXPECT_THROW(states.get_state(0), NotImplementedError);
    EXPECT_THROW(states.set_state(0, std::vector<double>{}), NotImplementedError);
}

TEST_F(CondensedGasTests, ThermoArrayFromGasWithoutCandidatesBehavesAsBefore) {
    ThermoArray states(gas, {2, 3});
    EXPECT_EQ(states.num_condensed(), 0u);
    EXPECT_TRUE(states.condensed_species_names().empty());
    EXPECT_NO_THROW(states.get_state(0));
}

TEST_F(CondensedGasTests, CombustorFromReactantGasesWaitsOnTheMultiphaseSolver) {
    // The reactant-stream path itself is implemented (work package C), but it equilibrates through
    // `Gas`, which cannot yet handle candidate condensed species.
    gas.add_condensed_species(CONDENSED_FILE, {"H2O(L)"});

    Gas fuel = make_gas("fuel", {"H2"});
    fuel.set_state_TPX(300.0, Cantera::OneAtm, "H2:1");
    Gas oxidizer = make_gas("oxidizer", {"O2"});
    oxidizer.set_state_TPX(300.0, Cantera::OneAtm, "O2:1");
    Combustor combustor(gas, fuel, oxidizer);

    Eigen::ArrayXd pressures(1);
    pressures << Cantera::OneAtm;
    Eigen::ArrayXd ratios(1);
    ratios << 8.0;

    EXPECT_THROW(combustor.solve(pressures, ratios), NotImplementedError);

    CombustorOptions isochoric;
    isochoric.process = CombustionProcess::ISOCHORIC;
    EXPECT_THROW(combustor.solve(pressures, ratios, isochoric), NotImplementedError);
}
