// Nozzle expansion and rocket results with condensed products (work package D).
//
// The fixtures mirror those of `test_condensed.cpp`: an aluminized ammonium-perchlorate
// propellant whose alumina condenses, and the RP-1311 example 13 beryllium rocket whose throat
// sits exactly on the BeO(b)/BeO(L) melting point.

#include "goddard/combustor.hpp"
#include "goddard/error.hpp"
#include "goddard/gas.hpp"
#include "goddard/global.hpp"
#include "goddard/kinetic_nozzle.hpp"
#include "goddard/moc_thermo.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/rocket.hpp"
#include "goddard/profile.hpp"
#include "goddard/shocks.hpp"
#include "goddard/thermoarray.hpp"
#include "status_test_helpers.hpp"

#include "cantera/core.h"
#include "gtest/gtest.h"

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace Goddard;
using Goddard::test_helpers::failed_with;
using Goddard::test_helpers::is_ok;

namespace {

/** CEA-derived NASA9 data converted to Cantera YAML (work package E). */
constexpr const char* NASA9_GAS = "nasa9_gas.yaml";
constexpr const char* NASA9_CONDENSED = "nasa9_condensed.yaml";
constexpr const char* NASA9_REACTANTS = "nasa9_reactants.yaml";

constexpr double BAR = 1.0e5;

/** Product gas spanning `elements`. */
Gas make_products(const std::vector<std::string>& elements) {
    setup_defaults();
    return Gas::create_from_elements(NASA9_GAS, "products", elements);
}

/** Reactant stream: a phase of the named species from the reactant database at 298.15 K. */
Gas make_reactants(const std::vector<std::string>& species, const std::string& mass_fractions) {
    setup_defaults();
    Gas reactants = Gas::create_from_species(NASA9_REACTANTS, "reactants", species);
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

/** Amount of one candidate condensed species [kmol per kg of mixture]; 0 if it is not a candidate. */
double condensed_moles_of(const Gas& gas, const std::string& name) {
    const std::vector<std::string> names = gas.condensed_species_names();
    auto found = std::find(names.begin(), names.end(), name);
    if (found == names.end()) return 0.0;
    return gas.condensed_moles()[static_cast<size_t>(found - names.begin())];
}

/** Products of an ammonium-perchlorate / aluminium propellant, with the alumina candidates. */
Gas make_propellant_products() {
    Gas products = make_products({"N", "H", "Cl", "O", "Al"});
    products.add_condensed_species(NASA9_CONDENSED, {"AL2O3(a)", "AL2O3(L)"});
    return products;
}

/** 80 % ammonium perchlorate, 20 % aluminium by mass, at 298.15 K. */
Gas make_propellant_reactants() {
    return make_reactants({"NH4CLO4(I)", "AL(cr)"}, "NH4CLO4(I):0.8, AL(cr):0.2");
}

/**
 * Burn the propellant in an infinite-area combustor.
 *
 * The propellant is a single monopropellant stream, so it is handed to the combustor as both the
 * fuel and the oxidizer and burnt at a fuel fraction of one.
 *
 * @param products Product mixture, whose condensed candidate set the combustor shares.
 * @param pressure Chamber pressure [Pa].
 * @return Chamber state vector, gas state followed by the condensed amounts.
 */
std::vector<double> burn_propellant(const Gas& products, double pressure) {
    Gas fuel = make_propellant_reactants();
    Gas combustor_gas = products;

    CombustorOptions options;
    options.mixture_type = MixtureRatioType::FUEL_FRAC;

    Eigen::ArrayXd pressures(1);
    pressures << pressure;
    Eigen::ArrayXd fuel_fractions(1);
    fuel_fractions << 1.0;

    Combustor combustor(combustor_gas, fuel, fuel);
    ThermoArray states = combustor.solve(pressures, fuel_fractions, options);
    return states.get_state(0);
}

/**
 * Products of the RP-1311 example 13 propellant.
 *
 * Every N/H/Be/O species of the database except gaseous Be(OH)2, which CEA's product list for this
 * example does not contain.
 */
Gas make_beryllium_products() {
    const Gas all = make_products({"N", "H", "Be", "O"});
    std::vector<std::string> species;
    for (const std::string& name : all.species_names()) {
        if (name != "Be(OH)2") species.push_back(name);
    }

    setup_defaults();
    Gas products = Gas::create_from_species(NASA9_GAS, "products", species);
    products.add_all_condensed_species(NASA9_CONDENSED);
    return products;
}

/** 67 % fuel (80 % N2H4(L), 20 % Be(a)) and 33 % H2O2(L) by mass, all at 298.15 K. */
Gas make_beryllium_reactants() {
    return make_reactants({"N2H4(L)", "Be(a)", "H2O2(L)"},
                          "N2H4(L):0.536, Be(a):0.134, H2O2(L):0.33");
}

} // namespace


// ---- Equilibrium expansion of an aluminized propellant ----

class PropellantNozzleTests : public ::testing::Test {
protected:
    PropellantNozzleTests()
        : products(make_propellant_products()),
          chamber_state(burn_propellant(products, 70.0 * BAR)) {}

    /** A mixture restored to `state`, sharing the candidate set of `products`. */
    Gas at_state(const std::vector<double>& state) {
        Gas gas = products;
        gas.restore_state(state);
        return gas;
    }

    Gas products;
    std::vector<double> chamber_state;
};

TEST_F(PropellantNozzleTests, EquilibriumExpansionCarriesTheAlumina) {
    double chamber_enthalpy = 0.0;
    double chamber_entropy = 0.0;
    double chamber_pressure = 0.0;
    {
        // `at_state` hands out mixtures that share one Cantera `Solution`, so the chamber values
        // are read out before any other station is restored.
        Gas chamber = at_state(chamber_state);
        chamber_enthalpy = chamber.enthalpy_mass();
        chamber_entropy = chamber.entropy_mass();
        chamber_pressure = chamber.pressure();
        ASSERT_LT(chamber.gas_mass_fraction(), 1.0) << "no alumina condensed in the chamber";
    }

    Gas nozzle_gas = products;
    NozzleOptions options;
    options.chemistry = GasChemistry::EQUILIBRIUM;

    Nozzle nozzle(nozzle_gas, chamber_state, options);
    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO,
                                         std::vector<double>{5.0, 20.0});
    ASSERT_TRUE(is_ok(results.status)) << results.first_failure().message;

    EXPECT_NEAR(results.throat.H_stagnation, chamber_enthalpy, 1e-10 * std::abs(chamber_enthalpy));
    EXPECT_NEAR(results.throat.S_inlet, chamber_entropy, 1e-10 * std::abs(chamber_entropy));

    // The throat is sonic: the reported Mach number is the one recomputed from its state.
    Gas throat = at_state(results.throat.state);
    const double throat_mach = throat.isenthalpic_velocity(results.throat.H_stagnation)
        / throat.speed_of_sound(CompositionModel::EQUILIBRIUM);
    EXPECT_NEAR(throat_mach, 1.0, 1e-4);
    EXPECT_NEAR(throat_mach, results.throat.mach, 1e-12);

    ASSERT_EQ(results.expansions.size(), 2u);

    double previous_pressure = chamber_pressure;
    for (const std::vector<double>* state :
         {&results.throat.state, &results.expansions[0].state, &results.expansions[1].state}) {
        Gas station = at_state(*state);

        EXPECT_LT(station.gas_mass_fraction(), 1.0) << "no condensate at P = " << station.pressure();
        EXPECT_GT(condensed_moles_of(station, "AL2O3(a)") + condensed_moles_of(station, "AL2O3(L)"),
                  0.0);

        const double velocity = station.isenthalpic_velocity(results.throat.H_stagnation);
        EXPECT_NEAR(station.enthalpy_mass() + 0.5 * velocity * velocity, chamber_enthalpy,
                    1e-8 * std::abs(chamber_enthalpy));
        EXPECT_NEAR(station.entropy_mass(), chamber_entropy, 1e-8 * std::abs(chamber_entropy));

        // Expanding to a larger area drops the pressure.
        EXPECT_LT(station.pressure(), previous_pressure);
        previous_pressure = station.pressure();
    }
}

TEST_F(PropellantNozzleTests, FrozenExpansionHoldsTheCondensedAmounts) {
    Gas nozzle_gas = products;
    NozzleOptions options;
    options.chemistry = GasChemistry::FROZEN;
    options.frozen_NFZ = 1;   // equilibrium up to the throat, frozen beyond it

    Nozzle nozzle(nozzle_gas, chamber_state, options);
    // Area ratio 20 would drive the frozen temperature below the melting point of alumina, which
    // `FrozenExpansionStopsWhereTheAluminaWouldFreeze` covers.
    NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO,
                                         std::vector<double>{2.0, 5.0});

    Gas throat = at_state(results.throat.state);
    const std::vector<double> throat_moles = throat.condensed_moles();
    ASSERT_GT(condensed_moles_of(throat, "AL2O3(L)"), 0.0);

    for (const NozzleStation& station : results.expansions) {
        Gas exit = at_state(station.state);
        EXPECT_EQ(exit.condensed_moles(), throat_moles);
        EXPECT_FALSE(station.thermo.pinned_transition);
    }
}

TEST_F(PropellantNozzleTests, FrozenExpansionStopsWhereTheAluminaWouldFreeze) {
    Gas nozzle_gas = products;
    NozzleOptions options;
    options.chemistry = GasChemistry::FROZEN;
    options.frozen_NFZ = 1;

    Nozzle nozzle(nozzle_gas, chamber_state, options);
    // Far enough to drive the frozen temperature below the 2327 K melting point of alumina,
    // which a frozen expansion cannot represent: CEA stops there too.
    NozzleResults results = nozzle.solve(ExpansionType::PRESSURE_RATIO, 1000.0);
    ASSERT_EQ(results.expansions.size(), 1u);
    const SolveStatus& status = results.expansions.front().status;
    EXPECT_TRUE(failed_with(status, FailureReason::OUT_OF_DATA_RANGE, "Frozen expansion"));
    EXPECT_NE(status.message.find("AL2O3(L)"), std::string::npos) << status.message;
    EXPECT_FALSE(results.status.ok());
}


TEST_F(PropellantNozzleTests, FailedThroatMakesEveryStationNotComputed) {
    // Liquid alumina just above its 2327 K melting point, frozen from the chamber on: the frozen
    // throat is colder than the lower bound of the AL2O3(L) data, which a frozen expansion
    // cannot represent.
    std::vector<double> cool_chamber_state;
    {
        Gas chamber = at_state(chamber_state);
        chamber.equilibrate_TP(2350.0, chamber.pressure());
        ASSERT_GT(condensed_moles_of(chamber, "AL2O3(L)"), 0.0);
        cool_chamber_state = chamber.save_state();
    }

    NozzleOptions options;
    options.chemistry = GasChemistry::FROZEN;
    options.frozen_NFZ = 0;
    Nozzle nozzle(products, cool_chamber_state, options);
    const NozzleResults results = nozzle.solve(ExpansionType::SUPERSONIC_AREA_RATIO, 10.0);

    EXPECT_TRUE(failed_with(results.throat.status, FailureReason::OUT_OF_DATA_RANGE, "AL2O3(L)"));
    EXPECT_TRUE(std::isnan(results.throat.velocity));
    ASSERT_EQ(results.expansions.size(), 1u);
    const NozzleStation& exit = results.expansions.front();
    EXPECT_EQ(exit.status.reason, FailureReason::NOT_COMPUTED);
    EXPECT_NE(exit.status.message.find(results.throat.status.message), std::string::npos)
        << "NOT_COMPUTED must repeat the upstream message: " << exit.status.message;
    EXPECT_TRUE(std::isnan(exit.mach));
    EXPECT_TRUE(std::isnan(exit.area_ratio));
    EXPECT_TRUE(exit.state.empty());

    EXPECT_FALSE(results.status.ok());
    EXPECT_EQ(results.first_failure().reason, FailureReason::OUT_OF_DATA_RANGE);

    // Stations solved from the failed throat are not computed either
    const std::vector<NozzleStation> downstream = nozzle.solve_stations(
        results.throat, ExpansionType::SUPERSONIC_AREA_RATIO, {3.0, 6.0});
    ASSERT_EQ(downstream.size(), 2u);
    for (const NozzleStation& station : downstream) {
        EXPECT_EQ(station.status.reason, FailureReason::NOT_COMPUTED);
    }
}

TEST_F(PropellantNozzleTests, OneFailedStationLeavesTheOtherStationsValid) {
    Gas nozzle_gas = products;
    NozzleOptions options;
    options.chemistry = GasChemistry::FROZEN;
    options.frozen_NFZ = 1;

    // The pressure ratio of 1000 drives the frozen temperature below the melting point of
    // alumina (see the previous test); the stations around it do not.
    Nozzle nozzle(nozzle_gas, chamber_state, options);
    const NozzleResults results = nozzle.solve(ExpansionType::PRESSURE_RATIO, std::vector<double>{2.0, 1000.0, 5.0});
    ASSERT_EQ(results.expansions.size(), 3u);

    EXPECT_TRUE(is_ok(results.expansions[0].status));
    EXPECT_TRUE(is_ok(results.expansions[2].status));
    EXPECT_TRUE(failed_with(results.expansions[1].status, FailureReason::OUT_OF_DATA_RANGE, "Frozen expansion"));
    EXPECT_GT(results.expansions[0].mach, 1.0);
    EXPECT_GT(results.expansions[2].mach, results.expansions[0].mach);
    EXPECT_TRUE(std::isnan(results.expansions[1].mach));
    EXPECT_TRUE(std::isnan(results.expansions[1].thermo.temperature));
    EXPECT_TRUE(results.expansions[1].state.empty());

    EXPECT_TRUE(is_ok(results.throat.status));
    EXPECT_FALSE(results.status.ok());
    EXPECT_EQ(results.first_failure().reason, FailureReason::OUT_OF_DATA_RANGE);
    EXPECT_EQ(results.first_failure().message, results.expansions[1].status.message);

    // The station after the failed one is the station a solve without the failing ratio gives:
    // nothing the failed solve left behind leaks into it.
    Nozzle reference_nozzle(nozzle_gas, chamber_state, options);
    const NozzleResults reference =
        reference_nozzle.solve(ExpansionType::PRESSURE_RATIO, std::vector<double>{2.0, 5.0});
    ASSERT_TRUE(is_ok(reference.status));
    EXPECT_NEAR(results.expansions[2].thermo.pressure, reference.expansions[1].thermo.pressure,
                1e-9 * reference.expansions[1].thermo.pressure);
    EXPECT_NEAR(results.expansions[2].thermo.temperature, reference.expansions[1].thermo.temperature,
                1e-9 * reference.expansions[1].thermo.temperature);
    EXPECT_NEAR(results.expansions[2].mach, reference.expansions[1].mach,
                1e-9 * reference.expansions[1].mach);
    EXPECT_NEAR(results.expansions[0].mach, reference.expansions[0].mach,
                1e-9 * reference.expansions[0].mach);
}


// ---- RP-1311 example 13: a throat pinned at the beryllia melting point ----

class BerylliumNozzleTests : public ::testing::Test {
protected:
    BerylliumNozzleTests() : products(make_beryllium_products()) {
        Gas reactants = make_beryllium_reactants();
        const double P = 206.8419 * BAR;
        products.set_element_moles(map_elements(reactants, products), 3000.0, P);
        products.equilibrate_HP(reactants.enthalpy_mass(), P);
        chamber_state = products.save_state();
    }

    Gas at_state(const std::vector<double>& state) {
        Gas gas = products;
        gas.restore_state(state);
        return gas;
    }

    Gas products;
    std::vector<double> chamber_state;
};

TEST_F(BerylliumNozzleTests, ThroatIsPinnedAtTheBerylliaMeltingPoint) {
    Gas chamber = at_state(chamber_state);
    EXPECT_NEAR(chamber.temperature(), 3018.89, 0.05);  // CEA; observed 0.02 K

    Gas nozzle_gas = products;
    NozzleOptions options;
    options.chemistry = GasChemistry::EQUILIBRIUM;

    Nozzle nozzle(nozzle_gas, chamber_state, options);
    const ThroatCondition throat = nozzle.solve_throat_conditions();

    EXPECT_TRUE(throat.pinned_transition);

    Gas throat_gas = at_state(throat.state);
    EXPECT_NEAR(throat_gas.temperature(), 2851.0, 1e-6);
    EXPECT_TRUE(throat_gas.at_phase_transition());
    // CEA's throat pressure for this expansion; observed error 1.6e-6.
    EXPECT_NEAR(throat_gas.pressure(), 127.2265 * BAR, 1e-5 * 127.2265 * BAR);
    // An isentropic change at a pinned temperature is isothermal, so gamma_s reduces to the
    // isothermal compressibility.
    EXPECT_NEAR(throat.gamma_s, -1.0 / throat.dlV_dlP_T, 1e-10);
}

TEST_F(BerylliumNozzleTests, PressureRatioStationFollowsTheSolidPolymorph) {
    Gas nozzle_gas = products;
    NozzleOptions options;
    options.chemistry = GasChemistry::EQUILIBRIUM;

    Nozzle nozzle(nozzle_gas, chamber_state, options);
    NozzleResults results = nozzle.solve(ExpansionType::PRESSURE_RATIO, 10.0);

    ASSERT_EQ(results.expansions.size(), 1u);
    Gas station = at_state(results.expansions[0].state);

    EXPECT_NEAR(station.pressure(), 20.68419 * BAR, 1e-6 * 20.68419 * BAR);
    EXPECT_NEAR(station.temperature(), 2453.58, 0.1);  // CEA; observed 0.05 K
    EXPECT_GT(condensed_moles_of(station, "BeO(b)"), 0.0);
    EXPECT_FALSE(results.expansions[0].thermo.pinned_transition);
}


// ---- Rocket solver results ----

TEST(RocketSolverCondensedTests, ReportsCondensedProductsAndTheMixtureMolecularWeight) {
    setup_defaults();

    Gas products = make_products({"N", "H", "Cl", "O", "Al"});
    products.add_all_condensed_species(NASA9_CONDENSED);

    // A reactant composition is a mole-fraction map; the propellant is specified by mass, so
    // convert it through the reactant phase itself.
    Gas propellant = make_propellant_reactants();
    Composition composition;
    {
        const std::vector<double> fractions = propellant.mole_fractions();
        const std::vector<std::string> names = propellant.species_names();
        for (size_t i = 0; i < names.size(); i++) {
            composition[names[i]] = fractions[i];
        }
    }

    RocketPoint point;
    point.label = "ap_al";
    point.fuel = PhaseSpecification(298.15, Cantera::OneAtm, composition);
    point.oxidizer = point.fuel;
    point.mixture_ratio = 1.0;
    point.combustor.mixture_type = MixtureRatioType::FUEL_FRAC;
    point.chamber_pressure = 70.0 * BAR;
    point.nozzle.chemistry = GasChemistry::EQUILIBRIUM;
    point.expansion_ratios = {5.0};

    RocketSolver solver(products, NASA9_REACTANTS);
    const RocketResult result = solver.solve(point);

    ASSERT_TRUE(result.status.ok()) << result.status.describe();
    EXPECT_LT(result.chamber.thermo.gas_mass_fraction, 1.0);
    EXPECT_GT(result.chamber.thermo.gas_mass_fraction, 0.0);
    EXPECT_LT(result.chamber.thermo.mixture_molecular_weight, result.chamber.thermo.molecular_weight);

    const std::string report = rocket_report(result);
    EXPECT_NE(report.find("MW, MOL WT"), std::string::npos) << report;
    EXPECT_NE(report.find("AL2O3(L)"), std::string::npos) << report;
}


// ---- Solvers that do not support condensed phases ----

class CondensedGuardTests : public ::testing::Test {
protected:
    CondensedGuardTests() : products(make_propellant_products()) {
        products.restore_state(burn_propellant(products, 70.0 * BAR));
    }

    Gas products;
};

TEST_F(CondensedGuardTests, KineticNozzleRejectsCondensedPhases) {
    ASSERT_TRUE(products.has_condensed_phases());

    NozzleProfile profile;
    for (int i = 0; i < 20; i++) {
        const double fraction = static_cast<double>(i) / 19.0;
        profile.push_back({10.0 * fraction, 1.0 + 0.4 * fraction});
    }

    EXPECT_THROW(KineticNozzle(products, profile, 0.01, 1.0), NotImplementedError);
}

TEST_F(CondensedGuardTests, MethodOfCharacteristicsRejectsCondensedPhases) {
    ASSERT_TRUE(products.has_condensed_phases());

    ThroatCondition throat;
    throat.speed_of_sound = products.speed_of_sound(CompositionModel::FROZEN);
    throat.H_stagnation = products.enthalpy_mass();
    throat.P_inlet = products.pressure();
    throat.S_inlet = products.entropy_mass();
    throat.gamma_s = products.gamma_s(CompositionModel::FROZEN);
    throat.dlV_dlP_T = -1.0;
    throat.dlV_dlT_P = 1.0;
    throat.state = products.save_state();

    Gas gas = products;
    EXPECT_THROW(MocThermo::tabulated(gas, GasChemistry::FROZEN, throat), NotImplementedError);
}

TEST_F(CondensedGuardTests, ShockSolverRejectsCondensedPhases) {
    ASSERT_TRUE(products.has_condensed_phases());

    ShockSolver solver(products, GasChemistry::FROZEN);
    EXPECT_THROW(solver.normal_shock(2.0), std::invalid_argument);
    EXPECT_THROW(solver.normal_shock_from_velocity(3000.0), std::invalid_argument);
    EXPECT_THROW(solver.reflected_shock(2.0), std::invalid_argument);
    EXPECT_THROW(solver.reflected_shock(2.0, GasChemistry::FROZEN, GasChemistry::EQUILIBRIUM),
                 std::invalid_argument);
    EXPECT_THROW(solver.oblique_shock_from_wave_angle(2.0, 0.6), std::invalid_argument);
    EXPECT_THROW(solver.oblique_shock_from_deflection(2.0, 0.2), std::invalid_argument);
    EXPECT_THROW(solver.max_deflection(2.0), std::invalid_argument);
}
