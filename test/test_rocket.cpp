#include <cmath>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "cantera/core.h"
#include "eigen3/Eigen/Dense"
#include "gtest/gtest.h"

#include "goddard/combustor.hpp"
#include "goddard/config.h"
#include "goddard/error.hpp"
#include "goddard/format.hpp"
#include "goddard/gas.hpp"
#include "goddard/numerics.hpp"
#include "goddard/rocket.hpp"
#include "goddard/utils.hpp"
#include "status_test_helpers.hpp"

using namespace Goddard;
using Goddard::test_helpers::failed_with;

namespace {

/** Every species of h2o2.yaml's `ohmech` phase. */
Gas h2o2_products() {
    return Gas(std::string(DATA_DIR) + "/h2o2.yaml", "ohmech");
}

/** H2/O2 from 300 K gaseous reactants (product-species path), equilibrium expansion. */
RocketPoint h2o2_point(double of_ratio, double chamber_pressure, std::vector<double> area_ratios) {
    RocketPoint point;
    point.fuel = PhaseSpecification(300.0, chamber_pressure, "H2:1");
    point.oxidizer = PhaseSpecification(300.0, chamber_pressure, "O2:1");
    point.mixture_ratio = of_ratio;
    point.chamber_pressure = chamber_pressure;
    point.nozzle.chemistry = GasChemistry::EQUILIBRIUM;
    point.expansion_type = ExpansionType::SUPERSONIC_AREA_RATIO;
    point.expansion_ratios = std::move(area_ratios);
    return point;
}

/** Chamber state of a point computed directly with a `Combustor`. */
std::vector<double> direct_chamber_state(const RocketPoint& point) {
    Combustor combustor(h2o2_products(), point.fuel.composition, point.oxidizer.composition);
    const ThermoArray states = combustor.solve(point.fuel.T, point.oxidizer.T,
        Eigen::ArrayXd::Constant(1, point.chamber_pressure),
        Eigen::ArrayXd::Constant(1, point.mixture_ratio), point.combustor);
    return states.get_state(0);
}

} // namespace

// ---- Sweeps: one result per point, in order ----

TEST(RocketSolverSweep, ResultsFollowThePointsAndMatchDirectCombustion) {
    // A sweep in a deliberately unsorted order: the results must follow the points, not a grid.
    const std::vector<std::pair<double, double>> operating_points = {
        {8.0, 20.0e5}, {4.0, 50.0e5}, {8.0, 50.0e5}, {4.0, 20.0e5}};
    std::vector<RocketPoint> points;
    points.reserve(operating_points.size());
    for (const auto& [of_ratio, pressure] : operating_points) {
        points.push_back(h2o2_point(of_ratio, pressure, {5.0}));
    }

    RocketSolver solver(h2o2_products());
    const std::vector<RocketResult> results = solver.solve(points);
    ASSERT_EQ(results.size(), points.size());

    Gas reference = h2o2_products();
    for (std::size_t i = 0; i < points.size(); i++) {
        ASSERT_TRUE(results[i].status.ok()) << results[i].status.describe();
        EXPECT_EQ(results[i].point.mixture_ratio, points[i].mixture_ratio);
        EXPECT_EQ(results[i].point.chamber_pressure, points[i].chamber_pressure);

        reference.restore_state(direct_chamber_state(points[i]));
        EXPECT_NEAR(results[i].chamber.thermo.pressure, reference.pressure(),
                    1e-9 * reference.pressure()) << "point " << i;
        EXPECT_NEAR(results[i].chamber.thermo.temperature, reference.temperature(),
                    1e-9 * reference.temperature()) << "point " << i;
    }
}

TEST(RocketSolverSweep, ReusedSolversDoNotCarryStateBetweenPoints) {
    // The solver keeps the combustor and the nozzle of the previous point. A point solved after
    // a different one must give what it gives alone.
    RocketPoint equilibrium = h2o2_point(6.0, 70.0e5, {10.0, 40.0});
    RocketPoint frozen = h2o2_point(4.0, 20.0e5, {10.0, 40.0});
    frozen.nozzle.chemistry = GasChemistry::FROZEN;
    frozen.nozzle.frozen_NFZ = 1;
    frozen.fuel.T = 400.0;

    RocketSolver solver(h2o2_products());
    const RocketResult alone = solver.solve(equilibrium);
    const std::vector<RocketResult> sweep = solver.solve({equilibrium, frozen, equilibrium});
    ASSERT_EQ(sweep.size(), 3u);

    for (const RocketResult* result : {&sweep[0], &sweep[2]}) {
        ASSERT_TRUE(result->status.ok()) << result->status.describe();
        EXPECT_NEAR(result->chamber.thermo.temperature, alone.chamber.thermo.temperature,
                    1e-9 * alone.chamber.thermo.temperature);
        ASSERT_EQ(result->performance.size(), 2u);
        for (std::size_t e = 0; e < 2; e++) {
            EXPECT_NEAR(result->exits[e].thermo.temperature, alone.exits[e].thermo.temperature,
                        1e-7 * alone.exits[e].thermo.temperature) << "exit " << e;
            EXPECT_NEAR(result->performance[e].ivac, alone.performance[e].ivac,
                        1e-7 * alone.performance[e].ivac) << "exit " << e;
        }
    }
    // The point in between used its own nozzle chemistry.
    EXPECT_EQ(sweep[1].exits[1].thermo.composition_model, CompositionModel::FROZEN);
}

TEST(RocketSolverSweep, PointsWithTheSameLabelAreKept) {
    RocketPoint first = h2o2_point(4.0, 20.0e5, {5.0});
    RocketPoint second = h2o2_point(8.0, 20.0e5, {5.0});
    first.label = second.label = "same";

    RocketSolver solver(h2o2_products());
    const std::vector<RocketResult> results = solver.solve({first, second});
    ASSERT_EQ(results.size(), 2u);
    EXPECT_NE(results[0].chamber.thermo.temperature, results[1].chamber.thermo.temperature);
    EXPECT_EQ(results[0].point.label, "same");
    EXPECT_EQ(results[1].point.label, "same");
}

// ---- Stations and performance ----

TEST(RocketSolverStations, PerformanceRecoversTheRequestedAreaRatioAndPressureRatio) {
    // Area ratio from mass flux continuity, rho_t a_t / (rho_e u_e): the nozzle solved for the
    // requested area ratio, so this recovers it to the station solver tolerance.
    const std::vector<double> area_ratios = {5.0, 20.0};
    RocketSolver solver(h2o2_products());
    for (double pressure : {20.0e5, 50.0e5}) {
        const RocketResult result = solver.solve(h2o2_point(6.0, pressure, area_ratios));
        ASSERT_TRUE(result.status.ok()) << result.status.describe();
        EXPECT_NEAR(result.chamber.thermo.pressure, pressure, 1e-9 * pressure);
        EXPECT_LT(result.throat.thermo.pressure, result.chamber.thermo.pressure);
        EXPECT_GT(result.throat.thermo.pressure, 0.5 * result.chamber.thermo.pressure);
        ASSERT_EQ(result.performance.size(), area_ratios.size());
        for (std::size_t e = 0; e < area_ratios.size(); e++) {
            const RocketPerformance& performance = result.performance[e];
            EXPECT_NEAR(performance.area_ratio, area_ratios[e], 1e-3 * area_ratios[e]);
            EXPECT_NEAR(performance.pressure_ratio, pressure / result.exits[e].thermo.pressure,
                        1e-9 * performance.pressure_ratio);
            // Exit velocity from the energy balance, reported as the specific impulse
            EXPECT_NEAR(performance.isp, result.exits[e].velocity, 1e-6 * performance.isp);
        }
    }
}

TEST(RocketSolverStations, InfiniteAreaChamberIsTheStagnationStateAndTheCombustionEnd) {
    RocketSolver solver(h2o2_products());
    const RocketResult result = solver.solve(h2o2_point(6.0, 50.0e5, {5.0}));
    ASSERT_TRUE(result.status.ok()) << result.status.describe();

    for (const NozzleStation* station : {&result.stagnation, &result.combustion_end}) {
        EXPECT_EQ(station->thermo.pressure, result.chamber.thermo.pressure);
        EXPECT_EQ(station->thermo.temperature, result.chamber.thermo.temperature);
        EXPECT_EQ(station->velocity, 0.0);
    }
    EXPECT_TRUE(std::isinf(result.contraction_ratio));
    EXPECT_EQ(result.mass_flux, 0.0);
}

TEST(RocketSolverStations, StationsCarryTheStagnationEnthalpyOfTheirPoint) {
    // The expansion is adiabatic, so every station of a point shares the stagnation enthalpy of
    // its chamber; the points differ in mixture ratio, so their chamber enthalpies differ.
    RocketSolver solver(h2o2_products());
    const std::vector<RocketResult> results = solver.solve(
        {h2o2_point(4.0, 20.0e5, {5.0, 20.0}), h2o2_point(6.0, 50.0e5, {5.0, 20.0})});
    ASSERT_NE(results[0].chamber.thermo.enthalpy, results[1].chamber.thermo.enthalpy);

    for (const RocketResult& result : results) {
        const double h0 = result.chamber.thermo.enthalpy;
        const double tolerance = max_fp_error(h0, 1e-9, 1e-6);
        EXPECT_NEAR(result.chamber.thermo.stagnation_enthalpy, h0, tolerance);
        EXPECT_NEAR(result.throat.thermo.stagnation_enthalpy, h0, tolerance);
        for (const NozzleStation& exit : result.exits) {
            EXPECT_NEAR(exit.thermo.stagnation_enthalpy, h0, tolerance);
            // h + u^2/2 = h0 at every exit
            EXPECT_NEAR(exit.thermo.enthalpy + 0.5 * exit.velocity * exit.velocity, h0,
                        1e-6 * std::abs(h0) + 1.0);
        }
    }
}

TEST(RocketSolverStations, IsochoricChamberIsTheConstantVolumeState) {
    const double initial_pressure = 10.0e5;
    RocketPoint point = h2o2_point(6.0, initial_pressure, {10.0});
    point.combustor.process = CombustionProcess::ISOCHORIC;

    RocketSolver solver(h2o2_products());
    const RocketResult result = solver.solve(point);
    ASSERT_TRUE(result.status.ok()) << result.status.describe();

    Gas reference = h2o2_products();
    reference.restore_state(direct_chamber_state(point));
    EXPECT_NEAR(result.chamber.thermo.pressure, reference.pressure(), 1e-6 * reference.pressure());
    EXPECT_NEAR(result.chamber.thermo.temperature, reference.temperature(),
                1e-6 * reference.temperature());
    // Constant-volume combustion raises the pressure well above the initial one
    EXPECT_GT(result.chamber.thermo.pressure, 5.0 * initial_pressure);
    EXPECT_LT(result.throat.thermo.pressure, result.chamber.thermo.pressure);
    ASSERT_EQ(result.exits.size(), 1u);
    EXPECT_LT(result.exits[0].thermo.pressure, result.throat.thermo.pressure);

    const std::string report = rocket_report(result);
    EXPECT_NE(report.find("CONSTANT-VOLUME COMBUSTOR"), std::string::npos) << report;
    EXPECT_NE(report.find("Pinitial"), std::string::npos) << report;
}

TEST(RocketSolverStations, FiniteAreaChamberSatisfiesTheMomentumBalance) {
    // Finite-area chamber at contraction ratio 2: the gas accelerates along the chamber, so
    // P_inj = P_c + rho_c u_c^2 at the combustion end, and the nozzle expands from the
    // stagnation state "inf" at P_inf < P_inj.
    const double injector_pressure = 50.0e5;
    const double contraction_ratio = 2.0;
    RocketPoint point = h2o2_point(6.0, injector_pressure, {5.0});
    point.combustor.type = CombustorType::FINITE_CONTRACTION_RATIO;
    point.combustor.contraction_ratio = contraction_ratio;

    RocketSolver solver(h2o2_products());
    const RocketResult result = solver.solve(point);
    ASSERT_TRUE(result.status.ok()) << result.status.describe();

    EXPECT_NEAR(result.chamber.thermo.pressure, injector_pressure, 1e-9 * injector_pressure);
    EXPECT_LT(result.stagnation.thermo.pressure, injector_pressure);
    EXPECT_DOUBLE_EQ(result.contraction_ratio, contraction_ratio);

    const NozzleStation& end = result.combustion_end;
    const double rho_u2 = end.thermo.density * end.velocity * end.velocity;
    EXPECT_NEAR(end.thermo.pressure + rho_u2, injector_pressure, 1e-5 * injector_pressure);
    EXPECT_NEAR(result.mass_flux, end.thermo.density * end.velocity, 1e-5 * result.mass_flux);

    // Performance is measured from the stagnation state, not from the injector face.
    ASSERT_EQ(result.performance.size(), 1u);
    const RocketPerformance expected = rocket_performance(
        result.stagnation.thermo, result.throat.thermo, result.exits[0].thermo);
    EXPECT_DOUBLE_EQ(result.performance[0].cstar, expected.cstar);
    EXPECT_DOUBLE_EQ(result.performance[0].isp, expected.isp);
    const RocketPerformance injector_based = rocket_performance(
        result.chamber.thermo, result.throat.thermo, result.exits[0].thermo);
    EXPECT_NE(result.performance[0].cstar, injector_based.cstar);

    // CEA's finite-area page layout (see data/cea_results/h2gas_fac.output)
    const std::string report = rocket_report(result);
    for (const char* label : {"FINITE AREA COMBUSTOR", "INJECTOR", "COMB END", "Pinj/P", "Pinf/P"}) {
        EXPECT_NE(report.find(label), std::string::npos) << label << "\n" << report;
    }
}

// ---- Report ----

TEST(RocketSolverReport, FrozenReportPrintsTheCompositionOfTheFrozenFlow) {
    // Frozen at the throat (frozen_NFZ = 1): the exits keep the throat composition, which differs
    // from the chamber composition because the flow is in equilibrium up to the throat. The
    // report's single composition block must show the frozen (exit) composition.
    RocketPoint point = h2o2_point(6.0, 50.0e5, {5.0, 20.0});
    point.nozzle.chemistry = GasChemistry::FROZEN;
    point.nozzle.frozen_NFZ = 1;
    point.label = "frozen";

    RocketSolver solver(h2o2_products());
    const RocketResult result = solver.solve(point);
    ASSERT_TRUE(result.status.ok()) << result.status.describe();
    const std::string report = rocket_report(result);

    EXPECT_NE(report.find("CASE = frozen"), std::string::npos) << report;
    EXPECT_NE(report.find("FROZEN AT THE THROAT"), std::string::npos) << report;

    // The report pads each species name to 16 characters before its mass fraction.
    auto composition_entry = [](const std::string& species, double mass_fraction) {
        return species + std::string(16 - species.size(), ' ')
            + format_mass_fraction(mass_fraction, 10);
    };
    const std::string species = "OH";
    const double frozen_fraction = result.exits.back().thermo.composition.at(species);
    const double chamber_fraction = result.chamber.thermo.composition.at(species);
    EXPECT_NEAR(frozen_fraction, result.throat.thermo.composition.at(species), 1e-12)
        << "the exits keep the throat composition";
    ASSERT_NE(composition_entry(species, frozen_fraction),
              composition_entry(species, chamber_fraction))
        << "the chamber and frozen OH fractions must differ for this test to discriminate";
    EXPECT_NE(report.find(composition_entry(species, frozen_fraction)), std::string::npos) << report;
    EXPECT_EQ(report.find(composition_entry(species, chamber_fraction)), std::string::npos) << report;
}

// ---- A failed point does not abort the sweep ----

class RocketSolverPartialFailureTests : public ::testing::Test {
protected:
    // At an area ratio of 1e5 the equilibrium exit falls outside what the backend can solve; the
    // exit at an area ratio of 10 solves. Both mixture ratios behave alike. The sweep takes a
    // few seconds, so the tests share one solve.
    static void SetUpTestSuite() {
        RocketSolver solver(h2o2_products());
        shared_results = std::make_unique<std::vector<RocketResult>>(solver.solve(
            {h2o2_point(4.0, 70.0e5, {10.0, 1.0e5}), h2o2_point(6.0, 70.0e5, {10.0, 1.0e5})}));
    }
    static void TearDownTestSuite() { shared_results.reset(); }

    const std::vector<RocketResult>& results() const { return *shared_results; }

    static std::unique_ptr<std::vector<RocketResult>> shared_results;
};

std::unique_ptr<std::vector<RocketResult>> RocketSolverPartialFailureTests::shared_results;

TEST_F(RocketSolverPartialFailureTests, onlyTheUnsolvableExitFailsAndTheStatusNamesIt) {
    ASSERT_EQ(results().size(), 2u);
    for (const RocketResult& result : results()) {
        EXPECT_TRUE(result.chamber.status.ok());
        EXPECT_TRUE(result.throat.status.ok());
        ASSERT_EQ(result.exits.size(), 2u);
        EXPECT_TRUE(result.exits[0].status.ok());
        EXPECT_GT(result.exits[0].thermo.temperature, 500.0);

        const NozzleStation& failed = result.exits[1];
        EXPECT_EQ(failed.status.reason, FailureReason::THERMO_FAILURE) << failed.status.describe();
        EXPECT_TRUE(std::isnan(failed.thermo.temperature));
        // A failed exit does not report the requested area ratio as if it had been solved
        EXPECT_TRUE(std::isnan(failed.area_ratio));

        EXPECT_TRUE(failed_with(result.status, FailureReason::THERMO_FAILURE));
        EXPECT_EQ(result.status.message, failed.status.message);
    }
}

TEST_F(RocketSolverPartialFailureTests, performanceOfTheFailedExitIsNaNAndTheGoodExitIsNot) {
    for (const RocketResult& result : results()) {
        const RocketPerformance& good = result.performance[0];
        EXPECT_TRUE(good.status.ok()) << good.status.describe();
        EXPECT_GT(good.isp, 2000.0);
        EXPECT_GT(good.cstar, 1500.0);

        const RocketPerformance& bad = result.performance[1];
        EXPECT_EQ(bad.status.reason, FailureReason::THERMO_FAILURE) << bad.status.describe();
        for (double value : {bad.isp, bad.ivac, bad.CF, bad.cstar, bad.mach_number, bad.area_ratio}) {
            EXPECT_TRUE(std::isnan(value));
        }
    }
}

TEST_F(RocketSolverPartialFailureTests, reportMarksTheFailedStationsAndStillPrintsTheTable) {
    const std::string report = rocket_report(results());
    EXPECT_NE(report.find("STATIONS NOT SOLVED (NaN ABOVE)"), std::string::npos) << report;
    EXPECT_EQ(report.find("OPERATING POINT NOT SOLVED"), std::string::npos) << report;

    // The failed exit is the last column of every row: NaN, never a formatted garbage value or 0
    for (const std::string label : {"RHO, KG/CU M", "T, K", "H2O "}) {
        const std::size_t start = report.find("\n" + label);
        ASSERT_NE(start, std::string::npos) << label << "\n" << report;
        const std::string row = report.substr(start + 1, report.find('\n', start + 1) - start - 1);
        EXPECT_EQ(row.substr(row.size() - 4), " nan") << row;
    }
}

TEST(RocketSolverFailures, failedCombustionFailsThePointAndTheSweepGoesOn) {
    // A 1e5 K fuel stream has an enthalpy no equilibrium state of the products can reach.
    RocketPoint hot = h2o2_point(6.0, 70.0e5, {10.0});
    hot.fuel.T = 1.0e5;
    RocketSolver solver(h2o2_products());
    const std::vector<RocketResult> results = solver.solve({hot, h2o2_point(6.0, 70.0e5, {10.0})});

    const RocketResult& failed = results[0];
    ASSERT_FALSE(failed.status.ok());
    EXPECT_FALSE(failed.chamber.status.ok());
    for (const SolveStatus* status : {&failed.stagnation.status, &failed.throat.status,
                                      &failed.exits[0].status}) {
        EXPECT_EQ(status->reason, FailureReason::NOT_COMPUTED) << status->describe();
    }
    EXPECT_TRUE(std::isnan(failed.performance[0].isp));
    EXPECT_NE(rocket_report(failed).find("OPERATING POINT NOT SOLVED"), std::string::npos);

    EXPECT_TRUE(results[1].status.ok()) << results[1].status.describe();
}

TEST(RocketSolverFailures, finiteAreaChamberChokedByTheMassFluxFailsEveryDownstreamStation) {
    RocketPoint choked = h2o2_point(6.0, 70.0e5, {5.0, 10.0});
    choked.combustor.type = CombustorType::FINITE_MASS_FLUX;
    // The throat mass flux at 70 bar is about 3e3 kg/(m^2 s)
    choked.combustor.mass_flux = 1.0e5;
    RocketSolver solver(h2o2_products());
    const RocketResult result = solver.solve(choked);

    EXPECT_TRUE(failed_with(result.status, FailureReason::NO_SOLUTION, "thermally chokes"));
    // The injector face does not depend on the chamber solve; everything computed from it does
    EXPECT_TRUE(result.chamber.status.ok());
    EXPECT_GT(result.chamber.thermo.temperature, 2000.0);
    for (const SolveStatus* status : {&result.stagnation.status, &result.combustion_end.status,
                                      &result.throat.status}) {
        EXPECT_TRUE(failed_with(*status, FailureReason::NO_SOLUTION, "thermally chokes"))
            << status->describe();
    }
    ASSERT_EQ(result.exits.size(), 2u);
    for (const NozzleStation& exit : result.exits) {
        EXPECT_TRUE(failed_with(exit.status, FailureReason::NOT_COMPUTED, "not computed because"));
        // The reason is given once, not as a chain of "not computed because"
        const std::string needle = "not computed because";
        EXPECT_EQ(exit.status.message.find(needle, exit.status.message.find(needle) + 1),
                  std::string::npos) << exit.status.message;
        EXPECT_TRUE(std::isnan(exit.area_ratio));
    }
    EXPECT_EQ(result.performance[0].status.reason, FailureReason::NO_SOLUTION);
    EXPECT_TRUE(std::isnan(result.performance[0].isp));
    EXPECT_TRUE(std::isnan(result.contraction_ratio));
}

TEST(RocketSolverFailures, finiteAreaPressureRatioAboveTheStagnationPressureFailsOnlyThatExit) {
    RocketPoint finite = h2o2_point(6.0, 70.0e5, {1.0, 20.0});
    finite.combustor.type = CombustorType::FINITE_CONTRACTION_RATIO;
    finite.combustor.contraction_ratio = 2.0;
    finite.expansion_type = ExpansionType::PRESSURE_RATIO;
    RocketSolver solver(h2o2_products());
    const RocketResult result = solver.solve(finite);

    // P_inj/P_inf is about 1.1 here, so the station at P_inj/P = 1 would be above the
    // stagnation pressure; it fails alone
    ASSERT_EQ(result.exits.size(), 2u);
    EXPECT_TRUE(failed_with(result.exits[0].status, FailureReason::NO_SOLUTION, "stagnation pressure"))
        << result.exits[0].status.describe();
    EXPECT_TRUE(std::isnan(result.exits[0].area_ratio));
    EXPECT_TRUE(result.exits[1].status.ok());
    EXPECT_GT(result.exits[1].area_ratio, 1.0);
    EXPECT_TRUE(result.stagnation.status.ok());
    EXPECT_TRUE(result.throat.status.ok());
}

// ---- Invalid points are rejected before anything is solved ----

TEST(RocketSolverValidation, invalidPointsThrowInvalidArgument) {
    RocketSolver solver(h2o2_products());
    const RocketPoint good = h2o2_point(6.0, 70.0e5, {10.0});

    RocketPoint unset_fuel = good;
    unset_fuel.fuel = PhaseSpecification();
    EXPECT_THROW(solver.solve(unset_fuel), std::invalid_argument);

    RocketPoint cold_oxidizer = good;
    cold_oxidizer.oxidizer.T = 0.0;
    EXPECT_THROW(solver.solve(cold_oxidizer), std::invalid_argument);

    RocketPoint no_reactant_pressure = good;
    no_reactant_pressure.fuel.P = 0.0;
    EXPECT_THROW(solver.solve(no_reactant_pressure), std::invalid_argument);

    RocketPoint no_chamber_pressure = good;
    no_chamber_pressure.chamber_pressure = NaN;
    EXPECT_THROW(solver.solve(no_chamber_pressure), std::invalid_argument);

    RocketPoint negative_ratio = good;
    negative_ratio.mixture_ratio = -1.0;
    EXPECT_THROW(solver.solve(negative_ratio), std::invalid_argument);

    RocketPoint fuel_fraction_above_one = good;
    fuel_fraction_above_one.combustor.mixture_type = MixtureRatioType::FUEL_FRAC;
    fuel_fraction_above_one.mixture_ratio = 1.5;
    EXPECT_THROW(solver.solve(fuel_fraction_above_one), std::invalid_argument);

    // Without a reactant file the reactants must be product species
    RocketPoint not_a_product = good;
    not_a_product.fuel = PhaseSpecification(300.0, 70.0e5, Composition{{"CH4", 1.0}});
    EXPECT_THROW(solver.solve(not_a_product), std::invalid_argument);

    RocketPoint contraction_ratio_of_one = good;
    contraction_ratio_of_one.combustor.type = CombustorType::FINITE_CONTRACTION_RATIO;
    contraction_ratio_of_one.combustor.contraction_ratio = 1.0;
    EXPECT_THROW(solver.solve(contraction_ratio_of_one), std::invalid_argument);

    // An area ratio of 1 is not a supersonic station
    RocketPoint sonic_exit = good;
    sonic_exit.expansion_ratios = {1.0};
    EXPECT_THROW(solver.solve(sonic_exit), std::invalid_argument);

    // A pressure ratio below 1 is above the chamber pressure
    RocketPoint below_one = good;
    below_one.expansion_type = ExpansionType::PRESSURE_RATIO;
    below_one.expansion_ratios = {20.0, 0.5};
    EXPECT_THROW(solver.solve(below_one), std::invalid_argument);

    // A freezing station beyond the last station (the chamber, the throat and one exit are
    // stations 0 to 2) would never freeze the flow
    RocketPoint beyond_last = good;
    beyond_last.nozzle.chemistry = GasChemistry::FROZEN;
    beyond_last.nozzle.frozen_NFZ = 3;
    EXPECT_THROW(solver.solve(beyond_last), std::invalid_argument);

    RocketPoint kinetic = good;
    kinetic.nozzle.chemistry = GasChemistry::KINETIC;
    EXPECT_THROW(solver.solve(kinetic), std::invalid_argument);

    RocketPoint frozen_finite_area = good;
    frozen_finite_area.combustor.type = CombustorType::FINITE_CONTRACTION_RATIO;
    frozen_finite_area.combustor.contraction_ratio = 2.0;
    frozen_finite_area.nozzle.chemistry = GasChemistry::FROZEN;
    frozen_finite_area.nozzle.frozen_NFZ = 0;
    EXPECT_THROW(solver.solve(frozen_finite_area), NotImplementedError);

    solver.include_transport = true;
    EXPECT_THROW(solver.solve(good), NotImplementedError);
}

TEST(RocketSolverValidation, everyPointIsCheckedBeforeAnyIsSolvedAndTheErrorNamesIt) {
    RocketPoint bad = h2o2_point(6.0, 70.0e5, {10.0});
    bad.expansion_ratios = {10.0, 0.5};
    RocketSolver solver(h2o2_products());
    try {
        solver.solve({h2o2_point(6.0, 70.0e5, {10.0}), h2o2_point(4.0, 70.0e5, {10.0}), bad});
        FAIL() << "the bad last point was not rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_EQ(std::string(error.what()).rfind("point 2: ", 0), 0u) << error.what();
    }
}

TEST(RocketSolverValidation, solverWorksOnItsOwnCopyOfTheProducts) {
    Gas products = h2o2_products();
    products.set_state_TPX(500.0, 2.0e5, "N2:1");
    RocketSolver solver(products);
    ASSERT_TRUE(solver.solve(h2o2_point(6.0, 70.0e5, {10.0})).status.ok());
    EXPECT_DOUBLE_EQ(products.temperature(), 500.0);
    EXPECT_DOUBLE_EQ(products.pressure(), 2.0e5);
}
