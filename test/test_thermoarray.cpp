#include "goddard/thermoarray.hpp"
#include "goddard/thermo.hpp"
#include "goddard/utils.hpp"
#include "goddard/numerics.hpp"
#include "goddard/error.hpp"
#include "goddard/global.hpp"
#include "status_test_helpers.hpp"

#include "cantera/core.h"
#include "gtest/gtest.h"
#include "eigen3/Eigen/Dense"
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>
#include <iostream>


constexpr size_t NUM_H2O2_SPECIES = 10;

TEST(ThermoArrayConstruction, createThermoArrayWithoutShape) {
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    auto array = Goddard::ThermoArray(sln, 5);

    ASSERT_FALSE(array.is_shape_set());

}

TEST(ThermoArrayConstruction, createThermoArrayWithShape) {

    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    std::vector<long> shape = {5,5,5};
    auto array = std::make_unique<Goddard::ThermoArray>(sln, shape);

    ASSERT_TRUE(array->is_shape_set());
    ASSERT_EQ(array->shape(), shape);
    EXPECT_EQ(array->size(), 125);

}


class ThermoArray2DTests: public ::testing::Test {
    protected:
    ThermoArray2DTests():
        sln(Cantera::newSolution("h2o2.yaml", "ohmech"))
    {
        sln->thermo()->setMoleFractionsByName("H2:1, N2:1, O2:1, AR:0.1");
        sln->thermo()->setTemperature(100.0);
        sln->thermo()->setPressure(1.0*Cantera::OneBar);
        pressures = Eigen::ArrayXd(9);
        pressures << 1,2,3,4,5,6,7,8,9;
        pressures *= Cantera::OneBar;

        temperatures = Eigen::ArrayXd(6);
        temperatures << 1,2,3,4,5,6;
        temperatures *=  100.0;

    

        //See https://cantera.org/dev/cxx/d7/dfa/classCantera_1_1IdealGasPhase.html#a59983094a3e5c4391305b2dd93e83ea0
        densities = pressures*sln->thermo()->meanMolecularWeight()/(Cantera::GasConstant*100.0);
        shape = {temperatures.size(), pressures.size()};
    }

    void SetUp() override {
        array = std::make_unique<Goddard::ThermoArray>(sln, shape);
    }
    std::shared_ptr<Cantera::Solution> sln;
    std::vector<long> shape;
    Eigen::ArrayXd pressures;
    Eigen::ArrayXd temperatures;
    Eigen::ArrayXd densities;
    // Eigen::ArrayXd enthalpies; 
    // Eigen::ArrayXXd compositions;
    std::unique_ptr<Goddard::ThermoArray> array;
};

TEST_F(ThermoArray2DTests, constructorWorks){
    ASSERT_EQ(array->shape(), shape);
    EXPECT_EQ(array->size(), temperatures.size() * pressures.size());
    EXPECT_EQ(array->get_state(0).size(), sln->thermo()->stateSize());
    
}

TEST_F(ThermoArray2DTests, broadcastTD) {
    ASSERT_NO_THROW(array->TD(temperatures, densities));
    
    int loc = 0;
    
    ASSERT_EQ(array->size(), densities.size()*temperatures.size()) << "SolutionArray must have same size as thermodynamic states.";
    ASSERT_GE(array->get_state(0).size(), 2) << "State vector should have at least 2 entries";

    for (int j = 0; j < densities.size(); j++) {
        for (int i = 0; i < temperatures.size(); i++){
            std::vector<double> state = array->get_state(loc);
            EXPECT_DOUBLE_EQ(state[0], temperatures(i)) << "i = " << i << ", j = " << j;
            EXPECT_DOUBLE_EQ(state[1], densities(j)) << "i = " << i << ", j = " << j; 
            loc++;
        }
    }
}

TEST_F(ThermoArray2DTests, broadcastTP) {
    
    
    ASSERT_NO_THROW(array->TP(temperatures, pressures));
    
    int loc = 0;
    
    ASSERT_EQ(array->size(), pressures.size()*temperatures.size()) << "SolutionArray must have same size as thermodynamic states.";
    ASSERT_GE(array->get_state(0).size(), 2) << "State vector should have at least 2 entries";

    auto MW = sln->thermo()->meanMolecularWeight();

    for (int j = 0; j < pressures.size(); j++) {
        for (int i = 0; i < temperatures.size(); i++){
            std::vector<double> state = array->get_state(loc);
            EXPECT_DOUBLE_EQ(state[0], temperatures(i)) << "i = " << i << ", j = " << j;
            
            //the state vector stores pressure implicitly as density, so needs to be calculated.
            double density = state[1];
            double calculated_pressure = Goddard::ideal_gas_D_to_P(density, temperatures(i), MW);
            EXPECT_DOUBLE_EQ(calculated_pressure, pressures(j)) << "i = " << i << ", j = " << j; 
            loc++;
        }
    }
}


TEST_F(ThermoArray2DTests, equilibrationHP) {
    ASSERT_NO_THROW(array->TP(temperatures, pressures));
    
    array->equilibrate("HP", "gibbs", Goddard::DEFAULT_RELTOL);

    /**
     * TODO: test results of equilibration.
     */

    
    int loc = 0;
    
    ASSERT_EQ(array->size(), pressures.size()*temperatures.size()) << "SolutionArray must have same size as thermodynamic states.";
    ASSERT_GE(array->get_state(0).size(), 2) << "State vector should have at least 2 entries";

    std::shared_ptr<Cantera::Solution> ref_sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    std::shared_ptr<Cantera::ThermoPhase> ref_thermo = ref_sln->thermo();
    ref_thermo->setMoleFractionsByName("H2:1, N2:1, O2:1, AR:0.1");

    // auto MW = sln->thermo()->meanMolecularWeight();
    std::vector<double> prev_state(ref_thermo->stateSize());
    ref_thermo->saveState(prev_state);

    std::vector<double> ref_state(ref_thermo->stateSize());


    for (int j = 0; j < pressures.size(); j++) {
        for (int i = 0; i < temperatures.size(); i++){
            std::vector<double> state = array->get_state(loc);

            ref_thermo->setState_TP(temperatures(i), pressures(j));
            ref_thermo->equilibrate("HP", "gibbs", Goddard::DEFAULT_RELTOL);
            ref_thermo->saveState(ref_state);
            
            for (size_t k = 0; k < state.size(); k++){
                EXPECT_NEAR(state[k], ref_state[k], Goddard::max_fp_error(ref_state[k])) <<
                    "State vector should be identical. Indices (i,j,k): (" << i << "," << j << "," << k <<")" ;
            }

            ref_thermo->restoreState(prev_state);
            loc++;
        }
    }
}


class ThermoArray3DTests: public ::testing::Test {
    protected:
    
    ThermoArray3DTests():
        sln(Cantera::newSolution("h2o2.yaml", "ohmech"))
    {
        sln->thermo()->setTemperature(100.0);
        sln->thermo()->setPressure(1.0*Cantera::OneBar);
        pressures = Eigen::ArrayXd(NUM_PRESSURES);
        pressures << 1,2,3,4,5;
        pressures *= Cantera::OneBar;

        temperatures = Eigen::ArrayXd(NUM_TEMPERATURES);
        temperatures << 1,2,3,4,5;
        temperatures *=  100.0;

        compositions = Eigen::ArrayXXd(NUM_COMPOSITIONS, NUM_H2O2_SPECIES);
        for (size_t i = 0; i < NUM_COMPOSITIONS; i++) {
            double h2 = static_cast<double>(i);
            double o2 = static_cast<double>(NUM_COMPOSITIONS - i);
            //"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"
            compositions.row(i) << h2, 0.0e-10, 0.0e-10, o2, 0.0e-10, 0.1, 0.0e-10, 0.5, 1.0, 0.0e-10;
        }
    }

    void SetUp() override {
        shape = {temperatures.size(), pressures.size(), compositions.rows()};
        array = std::make_unique<Goddard::ThermoArray>(sln, shape);
    }
    std::shared_ptr<Cantera::Solution> sln;
    std::vector<long> shape;
    Eigen::ArrayXd pressures;
    Eigen::ArrayXd temperatures;
    Eigen::ArrayXXd compositions;
    std::unique_ptr<Goddard::ThermoArray> array;
    const size_t NUM_PRESSURES = 5;
    const size_t NUM_TEMPERATURES = 5;
    const size_t NUM_COMPOSITIONS = 5;
    
};

TEST_F(ThermoArray3DTests, broadcastTPX){
    ASSERT_NO_THROW(array->TPX(temperatures, pressures, compositions));

    std::shared_ptr<Cantera::Solution> ref_sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    std::shared_ptr<Cantera::ThermoPhase> ref_thermo = ref_sln->thermo();

    
    ASSERT_EQ(array->size(), pressures.size()*temperatures.size()*compositions.rows()) 
        << "SolutionArray must have same size as thermodynamic states.";
    ASSERT_GE(array->get_state(0).size(), 2) 
        << "State vector should have at least 2 entries";
    
    
    int loc = 0;
    std::vector<double> ref_state(ref_thermo->stateSize());

    for (int k = 0; k < compositions.rows(); k++) {

        Eigen::ArrayXd composition = compositions.row(k);
        ref_thermo->setMoleFractions(composition.data());

        for (int j = 0; j < pressures.size(); j++) {
            for (int i = 0; i < temperatures.size(); i++){
                std::vector<double> state = array->get_state(loc);

                ref_thermo->setTemperature(temperatures(i));
                ref_thermo->setPressure(pressures(j));
                ref_thermo->saveState(ref_state);

                for (size_t l = 0; l < state.size(); l++) {
                    EXPECT_DOUBLE_EQ(state[l], ref_state[l])
                        << "State vectors should be identical. "
                        << "Indices (i,j,k,l): (" << i << "," << j << "," << k << "," << l << ")";
                }
                loc++;
            }
        }
    }
}


TEST_F(ThermoArray3DTests, equilibrateTPX){
    ASSERT_NO_THROW(array->TPX(temperatures, pressures, compositions));
    ASSERT_NO_THROW(array->equilibrate("HP", "gibbs", Goddard::DEFAULT_RELTOL)) << "Composition:\n" << compositions;


    std::shared_ptr<Cantera::Solution> ref_sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    std::shared_ptr<Cantera::ThermoPhase> ref_thermo = ref_sln->thermo();

    
    ASSERT_EQ(array->size(), pressures.size()*temperatures.size()*compositions.rows()) 
        << "SolutionArray must have same size as thermodynamic states.";
    ASSERT_GE(array->get_state(0).size(), 2) 
        << "State vector should have at least 2 entries";
    
    
    int loc = 0;
    std::vector<double> ref_state(ref_thermo->stateSize());

    for (int k = 0; k < compositions.rows(); k++) {

        Eigen::ArrayXd composition = compositions.row(k);

        for (int j = 0; j < pressures.size(); j++) {
            for (int i = 0; i < temperatures.size(); i++){

                std::vector<double> state = array->get_state(loc);

                EXPECT_NO_THROW(
                    ref_thermo->setMoleFractions(composition.data());
                    ref_thermo->setTemperature(temperatures(i));
                    ref_thermo->setPressure(pressures(j));
                    ref_thermo->equilibrate("HP", "gibbs", Goddard::DEFAULT_RELTOL);
                ) << "Composition: " << composition.data();
                ref_thermo->saveState(ref_state);

                for (size_t l = 0; l < state.size(); l++) {
                    EXPECT_NEAR(state[l], ref_state[l], Goddard::max_fp_error(ref_state[l], 1.0e-6, 1.0E-7)) 
                        << "State vectors should be identical. " 
                        << "Indices (i,j,k,l): (" << i << "," << j << "," << k << "," << l << ")";
                }
                loc++;
            }
        }
    }
}

// ---- Storage order and property getters ----

TEST(ThermoArrayIndexing, flatIndexFirstDimensionFastest) {
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    Goddard::ThermoArray array(sln, std::vector<long>{2, 3, 4});

    EXPECT_EQ(array.flat_index(0, 0, 0), 0);
    EXPECT_EQ(array.flat_index(1, 0, 0), 1);
    EXPECT_EQ(array.flat_index(0, 1, 0), 2);
    EXPECT_EQ(array.flat_index(0, 0, 1), 6);
    EXPECT_EQ(array.flat_index(1, 2, 3), 1 + 2*2 + 3*6);
    EXPECT_THROW(array.flat_index(2, 0, 0), std::out_of_range);
    EXPECT_THROW(array.flat_index(0, 0, 4), std::out_of_range);

    Goddard::ThermoArray array_2d(sln, std::vector<long>{2, 3});
    EXPECT_THROW(array_2d.flat_index(0, 0, 1), std::out_of_range);
}

TEST(ThermoArrayIndexing, emptyShapeLeavesShapeUnset) {
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    Goddard::ThermoArray array(sln, std::vector<long>{});
    EXPECT_FALSE(array.is_shape_set());

    Eigen::ArrayXd temperatures(2);
    temperatures << 500.0, 600.0;
    Eigen::ArrayXd pressures(3);
    pressures << 1e5, 2e5, 3e5;
    ASSERT_NO_THROW(array.TP(temperatures, pressures));
    EXPECT_TRUE(array.is_shape_set());
    EXPECT_EQ(array.shape(), (std::vector<long>{2, 3}));
}

TEST(ThermoArrayIndexing, setterDimensionMismatchThrows) {
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    Goddard::ThermoArray array(sln, std::vector<long>{2, 3, 1});
    Eigen::ArrayXd temperatures(2);
    temperatures << 500.0, 600.0;
    Eigen::ArrayXd pressures(3);
    pressures << 1e5, 2e5, 3e5;
    EXPECT_THROW(array.TP(temperatures, pressures), std::length_error);
}

TEST(ThermoArrayGetters, oneDimensionalArray) {
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    sln->thermo()->setState_TPX(300.0, 1e5, "H2:2, O2:1");
    Goddard::ThermoArray array(sln, std::vector<long>{3});

    auto scratch = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    std::vector<double> state(scratch->stateSize());
    for (int loc = 0; loc < 3; loc++) {
        scratch->setState_TPX(1000.0 + 100.0 * loc, 1e5, "H2:2, O2:1");
        scratch->saveState(state);
        array.set_state(loc, state);
    }

    Eigen::ArrayXXd T = array.temperature();
    ASSERT_EQ(T.rows(), 3);
    ASSERT_EQ(T.cols(), 1);
    for (int loc = 0; loc < 3; loc++) {
        EXPECT_DOUBLE_EQ(T(loc, 0), 1000.0 + 100.0 * loc);
    }
}

TEST_F(ThermoArray2DTests, gettersReturnEntriesByIndex) {
    ASSERT_NO_THROW(array->TP(temperatures, pressures));

    Eigen::ArrayXXd T = array->temperature();
    Eigen::ArrayXXd P = array->pressure();
    Eigen::ArrayXXd H = array->enthalpy_mass();
    ASSERT_EQ(T.rows(), temperatures.size());
    ASSERT_EQ(T.cols(), pressures.size());

    auto ref_thermo = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    ref_thermo->setMoleFractionsByName("H2:1, N2:1, O2:1, AR:0.1");
    for (int j = 0; j < pressures.size(); j++) {
        for (int i = 0; i < temperatures.size(); i++) {
            ref_thermo->setState_TP(temperatures(i), pressures(j));
            EXPECT_DOUBLE_EQ(T(i, j), temperatures(i)) << "i = " << i << ", j = " << j;
            EXPECT_NEAR(P(i, j), pressures(j), 1e-9 * pressures(j)) << "i = " << i << ", j = " << j;
            EXPECT_NEAR(H(i, j), ref_thermo->enthalpy_mass(), Goddard::max_fp_error(ref_thermo->enthalpy_mass()))
                << "i = " << i << ", j = " << j;
        }
    }

    EXPECT_THROW(array->temperature(1), std::out_of_range);
}

TEST_F(ThermoArray3DTests, gettersReturnEntriesBySlice) {
    ASSERT_NO_THROW(array->TPX(temperatures, pressures, compositions));

    auto ref_thermo = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    for (int k = 0; k < compositions.rows(); k++) {
        Eigen::ArrayXd composition = compositions.row(k);
        ref_thermo->setMoleFractions(composition.data());
        const double ref_molecular_weight = ref_thermo->meanMolecularWeight();

        Eigen::ArrayXXd T = array->temperature(k);
        Eigen::ArrayXXd P = array->pressure(k);
        Eigen::ArrayXXd MW = array->mean_molecular_weight(k);
        ASSERT_EQ(T.rows(), temperatures.size());
        ASSERT_EQ(T.cols(), pressures.size());
        for (int j = 0; j < pressures.size(); j++) {
            for (int i = 0; i < temperatures.size(); i++) {
                EXPECT_DOUBLE_EQ(T(i, j), temperatures(i)) << "i,j,k = " << i << "," << j << "," << k;
                EXPECT_NEAR(P(i, j), pressures(j), 1e-9 * pressures(j)) << "i,j,k = " << i << "," << j << "," << k;
                EXPECT_NEAR(MW(i, j), ref_molecular_weight, 1e-12 * ref_molecular_weight)
                    << "i,j,k = " << i << "," << j << "," << k;
            }
        }
    }

    EXPECT_THROW(array->temperature(static_cast<int>(compositions.rows())), std::out_of_range);
}


// ---- State isolation (Cantera 3.2 SolutionArray::getState buffering, Cantera issue #2067) ----

TEST_F(ThermoArray2DTests, lastEntryCorrectImmediatelyAfterEquilibrate) {
    ASSERT_NO_THROW(array->TP(temperatures, pressures));
    array->equilibrate("HP", "gibbs", Goddard::DEFAULT_RELTOL);

    const long i = temperatures.size() - 1;
    const long j = pressures.size() - 1;
    std::vector<double> state = array->get_state(array->flat_index(i, j));

    auto ref_thermo = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    ref_thermo->setMoleFractionsByName("H2:1, N2:1, O2:1, AR:0.1");
    ref_thermo->setState_TP(temperatures(i), pressures(j));
    ref_thermo->equilibrate("HP", "gibbs", Goddard::DEFAULT_RELTOL);
    std::vector<double> ref_state(ref_thermo->stateSize());
    ref_thermo->saveState(ref_state);

    ASSERT_EQ(state.size(), ref_state.size());
    for (size_t l = 0; l < state.size(); l++) {
        EXPECT_NEAR(state[l], ref_state[l], Goddard::max_fp_error(ref_state[l])) << "l = " << l;
    }
}

TEST(ThermoArrayIsolation, sizeOneArrayAfterEquilibrate) {
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    sln->thermo()->setMoleFractionsByName("H2:2, O2:1");
    Goddard::ThermoArray array(sln, std::vector<long>{1, 1});

    Eigen::ArrayXd temperatures(1);
    temperatures << 1000.0;
    Eigen::ArrayXd pressures(1);
    pressures << 1e5;
    array.TP(temperatures, pressures);
    array.equilibrate("HP", "gibbs", Goddard::DEFAULT_RELTOL);
    const std::vector<double> burnt = array.get_state(0);
    EXPECT_GT(burnt[0], 2500.0);

    sln->thermo()->setState_TP(300.0, 1e5);
    EXPECT_EQ(array.get_state(0), burnt);
    EXPECT_DOUBLE_EQ(array.temperature()(0, 0), burnt[0]);
}

TEST_F(ThermoArray2DTests, externalSolutionChangesDoNotLeak) {
    ASSERT_NO_THROW(array->TP(temperatures, pressures));
    const int last = array->size() - 1;
    const std::vector<double> last_state = array->get_state(last);

    sln->thermo()->setState_TP(555.0, 2e5);
    EXPECT_EQ(array->get_state(last), last_state);
    EXPECT_DOUBLE_EQ(array->get_state(0)[0], temperatures(0));
}

TEST_F(ThermoArray2DTests, operationsDoNotChangeConstructorSolution) {
    std::vector<double> original(sln->thermo()->stateSize());
    sln->thermo()->saveState(original);

    ASSERT_NO_THROW(array->TP(temperatures, pressures));
    array->equilibrate("HP", "gibbs", Goddard::DEFAULT_RELTOL);
    array->temperature();
    array->get_state(array->size() - 1);

    std::vector<double> after(sln->thermo()->stateSize());
    sln->thermo()->saveState(after);
    EXPECT_EQ(after, original);
}

// ============================================================
//  Failure reporting
// ============================================================

namespace {

using Goddard::FailureReason;
using Goddard::SolveError;
using Goddard::SolveStatus;
using Goddard::test_helpers::failed_with;
using Goddard::test_helpers::is_ok;
using Goddard::test_helpers::throws_solve_error;

/**
 * 2 x 2 array of H2/O2 at T = {300 K, 1e5 K} by P = {1 bar, 2 bar}. The hot row (i = 1) is far
 * outside the temperature range of the thermodynamic data, so a constant-enthalpy equilibrium
 * from it fails in the backend; the cold row (i = 0) equilibrates.
 */
class ThermoArrayFailureTests : public ::testing::Test {
protected:
    ThermoArrayFailureTests() : sln(Cantera::newSolution("h2o2.yaml", "ohmech")) {
        Goddard::setup_defaults();
        sln->thermo()->setMoleFractionsByName("H2:2, O2:1");
        array = std::make_unique<Goddard::ThermoArray>(sln, std::vector<long>{2, 2});
        temperatures = Eigen::ArrayXd(2);
        temperatures << 300.0, 1.0e5;
        pressures = Eigen::ArrayXd(2);
        pressures << 1.0e5, 2.0e5;
    }

    std::shared_ptr<Cantera::Solution> sln;
    std::unique_ptr<Goddard::ThermoArray> array;
    Eigen::ArrayXd temperatures;
    Eigen::ArrayXd pressures;
};

}

TEST_F(ThermoArrayFailureTests, freshArrayHasNoFailures) {
    EXPECT_EQ(array->num_failed(), 0);
    EXPECT_TRUE(is_ok(array->status()));
    EXPECT_TRUE(is_ok(array->status_at(0)));
    EXPECT_TRUE(array->ok().all());
    array->TP(temperatures, pressures);
    EXPECT_EQ(array->num_failed(), 0);
    EXPECT_TRUE(array->temperature().isFinite().all());
}

TEST_F(ThermoArrayFailureTests, failedEquilibriumIsRecordedPerEntryAndOthersSurvive) {
    array->TP(temperatures, pressures);
    ASSERT_NO_THROW(array->equilibrate("HP"));

    EXPECT_EQ(array->num_failed(), 2);
    const auto ok = array->ok();
    ASSERT_EQ(ok.rows(), 2);
    ASSERT_EQ(ok.cols(), 2);
    EXPECT_TRUE(ok(0, 0));
    EXPECT_TRUE(ok(0, 1));
    EXPECT_FALSE(ok(1, 0));
    EXPECT_FALSE(ok(1, 1));

    const int failed_loc = array->flat_index(1, 0);
    EXPECT_TRUE(failed_with(array->status_at(failed_loc), FailureReason::THERMO_FAILURE));
    EXPECT_TRUE(is_ok(array->status_at(array->flat_index(0, 0))));

    // The first failure represents the array and says how many entries failed
    const SolveStatus overall = array->status();
    EXPECT_TRUE(failed_with(overall, FailureReason::THERMO_FAILURE, "2 of 4"));
    EXPECT_NE(overall.stage.find("(1, 0)"), std::string::npos) << overall.describe();

    // Failed entries are NaN in every getter; good entries are burnt gas
    const Eigen::ArrayXXd T = array->temperature();
    EXPECT_TRUE(std::isnan(T(1, 0)));
    EXPECT_TRUE(std::isnan(T(1, 1)));
    EXPECT_GT(T(0, 0), 2500.0);
    EXPECT_GT(T(0, 1), 2500.0);
    for (const Eigen::ArrayXXd& values : {array->pressure(), array->enthalpy_mass(), array->entropy_mass(),
                                          array->mean_molecular_weight(), array->internal_energy_mole()}) {
        EXPECT_TRUE(std::isnan(values(1, 0)));
        EXPECT_TRUE(std::isfinite(values(0, 0)));
    }
}

TEST_F(ThermoArrayFailureTests, retryFailedLeavesGoodEntriesAloneAndCountsRemaining) {
    array->TP(temperatures, pressures);
    array->equilibrate("HP");
    const Eigen::ArrayXXd before = array->temperature();

    // Retrying cannot fix an entry outside the data range, and must not redo the good ones
    EXPECT_EQ(array->retry_failed("HP"), 2);
    const Eigen::ArrayXXd after = array->temperature();
    EXPECT_EQ(after(0, 0), before(0, 0));
    EXPECT_EQ(after(0, 1), before(0, 1));
    EXPECT_TRUE(std::isnan(after(1, 0)));
}

TEST_F(ThermoArrayFailureTests, settersAndSetStateResetStatuses) {
    array->TP(temperatures, pressures);
    array->equilibrate("HP");
    ASSERT_EQ(array->num_failed(), 2);

    // Writing a state to a failed entry makes it a valid entry again
    const int failed_loc = array->flat_index(1, 0);
    array->set_state(failed_loc, array->get_state(array->flat_index(0, 0)));
    EXPECT_TRUE(is_ok(array->status_at(failed_loc)));
    EXPECT_EQ(array->num_failed(), 1);
    EXPECT_TRUE(std::isfinite(array->temperature()(1, 0)));

    // A bulk setter resets every entry
    array->TP(temperatures, pressures);
    EXPECT_EQ(array->num_failed(), 0);
    EXPECT_TRUE(array->temperature().isFinite().all());
}

TEST_F(ThermoArrayFailureTests, setStatusMarksAnEntryFailedAndItsGettersNaN) {
    array->TP(temperatures, pressures);
    const int loc = array->flat_index(0, 1);
    array->set_status(loc, SolveStatus::failure(FailureReason::NOT_CONVERGED, "test failure"));

    EXPECT_EQ(array->num_failed(), 1);
    EXPECT_TRUE(failed_with(array->status_at(loc), FailureReason::NOT_CONVERGED, "test failure"));
    EXPECT_FALSE(array->ok()(0, 1));
    EXPECT_TRUE(array->ok()(0, 0));
    EXPECT_TRUE(std::isnan(array->temperature()(0, 1)));
    EXPECT_TRUE(failed_with(array->status(), FailureReason::NOT_CONVERGED, "test failure"));
    EXPECT_THROW(array->set_status(array->size(), SolveStatus{}), std::out_of_range);
    EXPECT_THROW(array->status_at(-1), std::out_of_range);
}

TEST_F(ThermoArrayFailureTests, bulkSetterThrowsOnAStateTheBackendRejects) {
    Eigen::ArrayXd negative(2);
    negative << 300.0, -5.0;
    EXPECT_TRUE(throws_solve_error([&] { array->TP(negative, pressures); }, FailureReason::THERMO_FAILURE));
}

TEST_F(ThermoArrayFailureTests, getStateOfAFailedEntryThrowsItsStatus) {
    array->TP(temperatures, pressures);
    array->equilibrate("HP");
    const int failed_loc = array->flat_index(1, 0);
    ASSERT_FALSE(array->status_at(failed_loc).ok());
    const std::vector<double> good_state = array->get_state(array->flat_index(0, 0));

    // The failed entry holds the state it had before the failed solve: there is no state to return.
    try {
        array->get_state(failed_loc);
        FAIL() << "expected a SolveError";
    } catch (const SolveError& error) {
        EXPECT_EQ(error.reason(), array->status_at(failed_loc).reason);
        EXPECT_EQ(error.status().message, array->status_at(failed_loc).message);
    }
    // Writing a state makes the entry readable again
    array->set_state(failed_loc, good_state);
    EXPECT_EQ(array->get_state(failed_loc), good_state);
}

TEST_F(ThermoArrayFailureTests, retryFailedStillReadsTheStateOfAFailedEntry) {
    array->TP(temperatures, pressures);
    array->equilibrate("HP");
    ASSERT_EQ(array->num_failed(), 2);
    // retry_failed starts each retry from the state the entry holds, which get_state refuses to give
    int still_failed = -1;
    EXPECT_NO_THROW(still_failed = array->retry_failed("HP", "vcs"));
    EXPECT_EQ(still_failed, 2);
    EXPECT_TRUE(failed_with(array->status_at(array->flat_index(1, 0)), FailureReason::THERMO_FAILURE));
}

TEST_F(ThermoArrayFailureTests, entryAccessorsRejectALocationOutsideTheArray) {
    array->TP(temperatures, pressures);
    const std::vector<double> state = array->get_state(0);
    for (int loc : {-1, array->size(), array->size() + 7}) {
        EXPECT_THROW(array->get_state(loc), std::out_of_range) << loc;
        EXPECT_THROW(array->set_state(loc, state), std::out_of_range) << loc;
        EXPECT_THROW(array->get_condensed_moles(loc), std::out_of_range) << loc;
        EXPECT_THROW(array->set_condensed_moles(loc, {}), std::out_of_range) << loc;
    }
}

TEST_F(ThermoArrayFailureTests, equilibrateRejectsAnUnknownPropertyPairOrSolver) {
    array->TP(temperatures, pressures);
    EXPECT_THROW(array->equilibrate("PT"), std::invalid_argument);
    EXPECT_THROW(array->equilibrate("TP", "gibs"), std::invalid_argument);
    EXPECT_THROW(array->retry_failed("TP", "nope"), std::invalid_argument);
    EXPECT_THROW(array->retry_failed("TT"), std::invalid_argument);
    // Nothing was solved or recorded
    EXPECT_EQ(array->num_failed(), 0);
}

// A bulk setter that throws leaves a consistent array: the entries written before the rejected
// one are ok and hold the state the setter gave them, the rejected one and those after it keep
// their previous state and status, and later reads agree with the stored entries.
TEST_F(ThermoArrayFailureTests, bulkSetterThatThrowsKeepsTheArrayConsistent) {
    Eigen::ArrayXd first_temperatures(2);
    first_temperatures << 300.0, 400.0;
    array->TP(first_temperatures, pressures);
    const std::vector<double> state_of_entry_1 = array->get_state(1);
    array->set_status(2, SolveStatus::failure(FailureReason::NOT_CONVERGED, "earlier failure"));

    // Enthalpy of this mixture at 600 K and 1 bar, from an independent Cantera phase
    auto reference = Cantera::newSolution("h2o2.yaml", "ohmech");
    reference->thermo()->setMoleFractionsByName("H2:2, O2:1");
    reference->thermo()->setState_TP(600.0, pressures[0]);
    const double enthalpy_at_600_K = reference->thermo()->enthalpy_mass();
    Eigen::ArrayXd enthalpies(2);
    // The second enthalpy is far beyond what any temperature in the data range gives
    enthalpies << enthalpy_at_600_K, 1.0e14;

    EXPECT_TRUE(throws_solve_error([&] { array->HP(enthalpies, pressures); },
                                   FailureReason::THERMO_FAILURE));

    // Entry 0 was written by the first, successful step of HP (600 K, not the 300 K it held
    // before and not whatever temperature the rejected iterate reached), and ok.
    EXPECT_TRUE(is_ok(array->status_at(0)));
    EXPECT_NEAR(array->temperature()(0, 0), 600.0, 1e-6);
    EXPECT_NEAR(array->enthalpy_mass()(0, 0), enthalpy_at_600_K, 1e-6 * std::abs(enthalpy_at_600_K));
    const std::vector<double> state_of_entry_0 = array->get_state(0);
    // A second read gives the same entry: reading does not disturb it
    EXPECT_EQ(array->get_state(0), state_of_entry_0);
    EXPECT_NEAR(array->temperature()(0, 0), 600.0, 1e-6);

    // The rejected entry and those after it keep their previous state and status
    EXPECT_TRUE(is_ok(array->status_at(1)));
    EXPECT_EQ(array->get_state(1), state_of_entry_1);
    EXPECT_NEAR(array->temperature()(1, 0), 400.0, 1e-9);
    EXPECT_TRUE(failed_with(array->status_at(2), FailureReason::NOT_CONVERGED, "earlier failure"));
    EXPECT_EQ(array->num_failed(), 1);
    // Entry 2 is still failed, so its state cannot be read; entry 3 never changed
    EXPECT_TRUE(throws_solve_error([&] { array->get_state(2); }, FailureReason::NOT_CONVERGED));
    EXPECT_NEAR(array->temperature()(1, 1), 400.0, 1e-9);
}

TEST(ThermoArrayThreeDimensional, okSelectsTheSlice) {
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    Goddard::ThermoArray array(sln, std::vector<long>{2, 2, 2});
    array.set_status(array.flat_index(1, 0, 1), SolveStatus::failure(FailureReason::NO_SOLUTION, "x"));
    EXPECT_TRUE(array.ok(0).all());
    const auto slice = array.ok(1);
    EXPECT_FALSE(slice(1, 0));
    EXPECT_EQ(slice.count(), 3);
    EXPECT_EQ(array.num_failed(), 1);
}
