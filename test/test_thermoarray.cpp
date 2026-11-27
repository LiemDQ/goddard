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

TEST(ThermoArrayConstruction, createThermoArrayWithShape) {

    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    std::vector<long> shape = {5,5,5};
    auto array = std::make_unique<Goddard::ThermoArray>(Goddard::Gas(sln), shape);

    ASSERT_EQ(array->shape(), shape);
    EXPECT_EQ(array->size(), 125);
    EXPECT_EQ(array->ndim(), 3);

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
        array = std::make_unique<Goddard::ThermoArray>(Goddard::Gas(sln), shape);
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
    
    ASSERT_EQ(array->size(), densities.size()*temperatures.size()) << "ThermoArray must have same size as thermodynamic states.";
    ASSERT_GE(array->get_state(0).size(), 2) << "State vector should have at least 2 entries";

    for (int i = 0; i < temperatures.size(); i++) {
        for (int j = 0; j < densities.size(); j++){
            std::vector<double> state = array->get_state(array->flat_index(i, j));
            EXPECT_DOUBLE_EQ(state[0], temperatures(i)) << "i = " << i << ", j = " << j;
            EXPECT_DOUBLE_EQ(state[1], densities(j)) << "i = " << i << ", j = " << j; 
        }
    }
}

TEST_F(ThermoArray2DTests, broadcastTP) {
    
    
    ASSERT_NO_THROW(array->TP(temperatures, pressures));
    
    ASSERT_EQ(array->size(), pressures.size()*temperatures.size()) << "ThermoArray must have same size as thermodynamic states.";
    ASSERT_GE(array->get_state(0).size(), 2) << "State vector should have at least 2 entries";

    auto MW = sln->thermo()->meanMolecularWeight();

    for (int i = 0; i < temperatures.size(); i++) {
        for (int j = 0; j < pressures.size(); j++){
            std::vector<double> state = array->get_state(array->flat_index(i, j));
            EXPECT_DOUBLE_EQ(state[0], temperatures(i)) << "i = " << i << ", j = " << j;
            
            //the state vector stores pressure implicitly as density, so needs to be calculated.
            double density = state[1];
            double calculated_pressure = Goddard::ideal_gas_D_to_P(density, temperatures(i), MW);
            EXPECT_DOUBLE_EQ(calculated_pressure, pressures(j)) << "i = " << i << ", j = " << j; 
        }
    }
}


TEST_F(ThermoArray2DTests, equilibrationHP) {
    ASSERT_NO_THROW(array->TP(temperatures, pressures));
    
    array->equilibrate("HP", "gibbs");

    /**
     * TODO: test results of equilibration.
     */

    ASSERT_EQ(array->size(), pressures.size()*temperatures.size()) << "ThermoArray must have same size as thermodynamic states.";
    ASSERT_GE(array->get_state(0).size(), 2) << "State vector should have at least 2 entries";

    std::shared_ptr<Cantera::Solution> ref_sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    std::shared_ptr<Cantera::ThermoPhase> ref_thermo = ref_sln->thermo();
    ref_thermo->setMoleFractionsByName("H2:1, N2:1, O2:1, AR:0.1");

    // auto MW = sln->thermo()->meanMolecularWeight();
    std::vector<double> prev_state(ref_thermo->stateSize());
    ref_thermo->saveState(prev_state);

    std::vector<double> ref_state(ref_thermo->stateSize());


    for (int i = 0; i < temperatures.size(); i++) {
        for (int j = 0; j < pressures.size(); j++){
            std::vector<double> state = array->get_state(array->flat_index(i, j));

            ref_thermo->setState_TP(temperatures(i), pressures(j));
            ref_thermo->equilibrate("HP", "gibbs", Goddard::DEFAULT_RELTOL);
            ref_thermo->saveState(ref_state);
            
            for (size_t k = 0; k < state.size(); k++){
                EXPECT_NEAR(state[k], ref_state[k], Goddard::max_fp_error(ref_state[k])) <<
                    "State vector should be identical. Indices (i,j,k): (" << i << "," << j << "," << k <<")" ;
            }

            ref_thermo->restoreState(prev_state);
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
        array = std::make_unique<Goddard::ThermoArray>(Goddard::Gas(sln), shape);
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
        << "ThermoArray must have same size as thermodynamic states.";
    ASSERT_GE(array->get_state(0).size(), 2) 
        << "State vector should have at least 2 entries";
    
    
    std::vector<double> ref_state(ref_thermo->stateSize());

    for (int k = 0; k < compositions.rows(); k++) {

        Eigen::ArrayXd composition = compositions.row(k);
        ref_thermo->setMoleFractions(composition.data());

        for (int i = 0; i < temperatures.size(); i++) {
            for (int j = 0; j < pressures.size(); j++){
                std::vector<double> state = array->get_state(array->flat_index(i, j, k));

                ref_thermo->setTemperature(temperatures(i));
                ref_thermo->setPressure(pressures(j));
                ref_thermo->saveState(ref_state);

                for (size_t l = 0; l < state.size(); l++) {
                    EXPECT_DOUBLE_EQ(state[l], ref_state[l])
                        << "State vectors should be identical. "
                        << "Indices (i,j,k,l): (" << i << "," << j << "," << k << "," << l << ")";
                }
            }
        }
    }
}


TEST_F(ThermoArray3DTests, equilibrateTPX){
    ASSERT_NO_THROW(array->TPX(temperatures, pressures, compositions));
    ASSERT_NO_THROW(array->equilibrate("HP", "gibbs")) << "Composition:\n" << compositions;


    std::shared_ptr<Cantera::Solution> ref_sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    std::shared_ptr<Cantera::ThermoPhase> ref_thermo = ref_sln->thermo();

    
    ASSERT_EQ(array->size(), pressures.size()*temperatures.size()*compositions.rows()) 
        << "ThermoArray must have same size as thermodynamic states.";
    ASSERT_GE(array->get_state(0).size(), 2) 
        << "State vector should have at least 2 entries";
    
    
    std::vector<double> ref_state(ref_thermo->stateSize());

    for (int k = 0; k < compositions.rows(); k++) {

        Eigen::ArrayXd composition = compositions.row(k);

        for (int i = 0; i < temperatures.size(); i++) {
            for (int j = 0; j < pressures.size(); j++){

                std::vector<double> state = array->get_state(array->flat_index(i, j, k));

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
            }
        }
    }
}

// ---- Storage order and property getters ----

TEST(ThermoArrayIndexing, flatIndexLastDimensionFastest) {
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    Goddard::ThermoArray array(Goddard::Gas(sln), std::vector<long>{2, 3, 4});

    EXPECT_EQ(array.flat_index(0, 0, 0), 0);
    EXPECT_EQ(array.flat_index(0, 0, 1), 1);
    EXPECT_EQ(array.flat_index(0, 1, 0), 4);
    EXPECT_EQ(array.flat_index(1, 0, 0), 12);
    EXPECT_EQ(array.flat_index(1, 2, 3), (1*3 + 2)*4 + 3);
    EXPECT_EQ(array.flat_index(std::vector<long>{1, 2, 3}), array.flat_index(1, 2, 3));
    EXPECT_THROW(array.flat_index(2, 0, 0), std::out_of_range);
    EXPECT_THROW(array.flat_index(0, 0, 4), std::out_of_range);
    EXPECT_THROW(array.flat_index(std::vector<long>{1, 2}), std::invalid_argument);

    Goddard::ThermoArray array_2d(Goddard::Gas(sln), std::vector<long>{2, 3});
    EXPECT_EQ(array_2d.flat_index(1, 2), 5);
    EXPECT_THROW(array_2d.flat_index(0, 0, 1), std::out_of_range);
}

TEST(ThermoArrayIndexing, constructorRejectsInvalidShapes) {
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    EXPECT_THROW(Goddard::ThermoArray(Goddard::Gas(sln), std::vector<long>{}), std::invalid_argument);
    EXPECT_THROW(Goddard::ThermoArray(Goddard::Gas(sln), std::vector<long>{1, 1, 1, 1}), std::invalid_argument);
    EXPECT_THROW(Goddard::ThermoArray(Goddard::Gas(sln), std::vector<long>{2, -1}), std::invalid_argument);
}

TEST(ThermoArrayIndexing, reshapeKeepsEntriesInStorageOrder) {
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    Goddard::ThermoArray array(Goddard::Gas(sln), std::vector<long>{2, 3});
    EXPECT_THROW(array.reshape({4, 2}), std::invalid_argument);
    array.reshape({3, 2});
    EXPECT_EQ(array.shape(), (std::vector<long>{3, 2}));
    EXPECT_EQ(array.size(), 6);
    EXPECT_EQ(array.flat_index(2, 1), 5);
}

TEST(ThermoArrayIndexing, setterDimensionMismatchThrows) {
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    Goddard::ThermoArray array(Goddard::Gas(sln), std::vector<long>{2, 3, 1});
    Eigen::ArrayXd temperatures(2);
    temperatures << 500.0, 600.0;
    Eigen::ArrayXd pressures(3);
    pressures << 1e5, 2e5, 3e5;
    EXPECT_THROW(array.TP(temperatures, pressures), std::length_error);

    // A 2-D array of the wrong extents is rejected too
    Goddard::ThermoArray array_2d(Goddard::Gas(sln), std::vector<long>{3, 2});
    EXPECT_THROW(array_2d.TP(temperatures, pressures), std::length_error);
}

TEST(ThermoArrayGetters, oneDimensionalArray) {
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    sln->thermo()->setState_TPX(300.0, 1e5, "H2:2, O2:1");
    Goddard::ThermoArray array(Goddard::Gas(sln), std::vector<long>{3});

    auto scratch = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    std::vector<double> state(scratch->stateSize());
    for (int loc = 0; loc < 3; loc++) {
        scratch->setState_TPX(1000.0 + 100.0 * loc, 1e5, "H2:2, O2:1");
        scratch->saveState(state);
        array.set_state(loc, state);
    }

    Eigen::ArrayXd T = array.temperature();
    ASSERT_EQ(T.size(), 3);
    for (int loc = 0; loc < 3; loc++) {
        EXPECT_DOUBLE_EQ(T(loc), 1000.0 + 100.0 * loc);
    }
}

TEST_F(ThermoArray2DTests, gettersReturnEntriesInStorageOrder) {
    ASSERT_NO_THROW(array->TP(temperatures, pressures));

    Eigen::ArrayXd T = array->temperature();
    Eigen::ArrayXd P = array->pressure();
    Eigen::ArrayXd H = array->enthalpy_mass();
    ASSERT_EQ(T.size(), temperatures.size() * pressures.size());
    ASSERT_EQ(P.size(), T.size());
    ASSERT_EQ(H.size(), T.size());

    auto ref_thermo = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    ref_thermo->setMoleFractionsByName("H2:1, N2:1, O2:1, AR:0.1");
    for (int i = 0; i < temperatures.size(); i++) {
        for (int j = 0; j < pressures.size(); j++) {
            const int loc = array->flat_index(i, j);
            ref_thermo->setState_TP(temperatures(i), pressures(j));
            EXPECT_DOUBLE_EQ(T(loc), temperatures(i)) << "i = " << i << ", j = " << j;
            EXPECT_NEAR(P(loc), pressures(j), 1e-9 * pressures(j)) << "i = " << i << ", j = " << j;
            EXPECT_NEAR(H(loc), ref_thermo->enthalpy_mass(), Goddard::max_fp_error(ref_thermo->enthalpy_mass()))
                << "i = " << i << ", j = " << j;
        }
    }
}

TEST_F(ThermoArray3DTests, gettersReturnEntriesInStorageOrder) {
    ASSERT_NO_THROW(array->TPX(temperatures, pressures, compositions));

    const Eigen::ArrayXd T = array->temperature();
    const Eigen::ArrayXd P = array->pressure();
    const Eigen::ArrayXd MW = array->molecular_weight();
    ASSERT_EQ(T.size(), array->size());

    auto ref_thermo = Cantera::newSolution("h2o2.yaml", "ohmech")->thermo();
    for (int k = 0; k < compositions.rows(); k++) {
        Eigen::ArrayXd composition = compositions.row(k);
        ref_thermo->setMoleFractions(composition.data());
        const double ref_molecular_weight = ref_thermo->meanMolecularWeight();

        for (int i = 0; i < temperatures.size(); i++) {
            for (int j = 0; j < pressures.size(); j++) {
                const int loc = array->flat_index(i, j, k);
                EXPECT_DOUBLE_EQ(T(loc), temperatures(i)) << "i,j,k = " << i << "," << j << "," << k;
                EXPECT_NEAR(P(loc), pressures(j), 1e-9 * pressures(j)) << "i,j,k = " << i << "," << j << "," << k;
                EXPECT_NEAR(MW(loc), ref_molecular_weight, 1e-12 * ref_molecular_weight)
                    << "i,j,k = " << i << "," << j << "," << k;
            }
        }
    }
}

TEST_F(ThermoArray3DTests, compositionGettersHaveOneRowPerEntry) {
    ASSERT_NO_THROW(array->TPX(temperatures, pressures, compositions));

    const std::vector<std::string> names = array->species_names();
    ASSERT_EQ(names.size(), NUM_H2O2_SPECIES);
    const Eigen::ArrayXXd x = array->mole_fractions();
    ASSERT_EQ(x.rows(), array->size());
    ASSERT_EQ(x.cols(), NUM_H2O2_SPECIES);
    for (int k = 0; k < compositions.rows(); k++) {
        const Eigen::ArrayXd expected = compositions.row(k) / compositions.row(k).sum();
        const int loc = array->flat_index(1, 2, k);
        for (size_t s = 0; s < NUM_H2O2_SPECIES; s++) {
            EXPECT_NEAR(x(loc, s), expected(s), 1e-12) << "k = " << k << ", species " << names[s];
        }
        EXPECT_NEAR(array->mass_fractions().row(loc).sum(), 1.0, 1e-12);
    }
}

TEST_F(ThermoArray2DTests, gasAtAndEvaluateAgreeWithGetters) {
    ASSERT_NO_THROW(array->TP(temperatures, pressures));

    const int loc = array->flat_index(2, 3);
    Goddard::Gas gas = array->gas_at(loc);
    EXPECT_DOUBLE_EQ(gas.temperature(), temperatures(2));
    EXPECT_NEAR(gas.pressure(), pressures(3), 1e-9 * pressures(3));

    const Eigen::ArrayXd cp = array->evaluate([](const Goddard::Gas& g) { return g.cp_mass(); });
    ASSERT_EQ(cp.size(), array->size());
    EXPECT_DOUBLE_EQ(cp(loc), gas.cp_mass());
    EXPECT_DOUBLE_EQ(cp(loc), array->cp_mass()(loc));
}


// ---- State isolation ----

TEST_F(ThermoArray2DTests, lastEntryCorrectImmediatelyAfterEquilibrate) {
    ASSERT_NO_THROW(array->TP(temperatures, pressures));
    array->equilibrate("HP", "gibbs");

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
    Goddard::ThermoArray array(Goddard::Gas(sln), std::vector<long>{1, 1});

    Eigen::ArrayXd temperatures(1);
    temperatures << 1000.0;
    Eigen::ArrayXd pressures(1);
    pressures << 1e5;
    array.TP(temperatures, pressures);
    array.equilibrate("HP", "gibbs");
    const std::vector<double> burnt = array.get_state(0);
    EXPECT_GT(burnt[0], 2500.0);

    sln->thermo()->setState_TP(300.0, 1e5);
    EXPECT_EQ(array.get_state(0), burnt);
    EXPECT_DOUBLE_EQ(array.temperature()(0), burnt[0]);
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
    array->equilibrate("HP", "gibbs");
    array->temperature();
    array->get_state(array->size() - 1);

    std::vector<double> after(sln->thermo()->stateSize());
    sln->thermo()->saveState(after);
    EXPECT_EQ(after, original);
}

TEST_F(ThermoArray2DTests, copyHasIndependentStates) {
    ASSERT_NO_THROW(array->TP(temperatures, pressures));
    Goddard::ThermoArray copy(*array);
    const std::vector<double> original = array->get_state(0);

    copy.TP(temperatures + 100.0, pressures);
    EXPECT_EQ(array->get_state(0), original);
    EXPECT_DOUBLE_EQ(copy.temperature()(0), temperatures(0) + 100.0);
    EXPECT_DOUBLE_EQ(array->temperature()(0), temperatures(0));
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
        array = std::make_unique<Goddard::ThermoArray>(Goddard::Gas(sln), std::vector<long>{2, 2});
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
    ASSERT_EQ(ok.size(), 4);
    EXPECT_TRUE(ok(array->flat_index(0, 0)));
    EXPECT_TRUE(ok(array->flat_index(0, 1)));
    EXPECT_FALSE(ok(array->flat_index(1, 0)));
    EXPECT_FALSE(ok(array->flat_index(1, 1)));

    const int failed_loc = array->flat_index(1, 0);
    EXPECT_TRUE(failed_with(array->status_at(failed_loc), FailureReason::THERMO_FAILURE));
    EXPECT_TRUE(is_ok(array->status_at(array->flat_index(0, 0))));

    // The first failure represents the array and says how many entries failed
    const SolveStatus overall = array->status();
    EXPECT_TRUE(failed_with(overall, FailureReason::THERMO_FAILURE, "2 of 4"));
    EXPECT_NE(overall.stage.find("(1, 0)"), std::string::npos) << overall.describe();

    // Failed entries are NaN in every getter; good entries are burnt gas
    const Eigen::ArrayXd T = array->temperature();
    EXPECT_TRUE(std::isnan(T(array->flat_index(1, 0))));
    EXPECT_TRUE(std::isnan(T(array->flat_index(1, 1))));
    EXPECT_GT(T(array->flat_index(0, 0)), 2500.0);
    EXPECT_GT(T(array->flat_index(0, 1)), 2500.0);
    for (const Eigen::ArrayXd& values : {array->pressure(), array->enthalpy_mass(), array->entropy_mass(),
                                         array->molecular_weight(), array->internal_energy_mass()}) {
        EXPECT_TRUE(std::isnan(values(array->flat_index(1, 0))));
        EXPECT_TRUE(std::isfinite(values(array->flat_index(0, 0))));
    }
}

TEST_F(ThermoArrayFailureTests, retryFailedLeavesGoodEntriesAloneAndCountsRemaining) {
    array->TP(temperatures, pressures);
    array->equilibrate("HP");
    const Eigen::ArrayXd before = array->temperature();

    // Retrying cannot fix an entry outside the data range, and must not redo the good ones
    EXPECT_EQ(array->retry_failed("HP"), 2);
    const Eigen::ArrayXd after = array->temperature();
    EXPECT_EQ(after(array->flat_index(0, 0)), before(array->flat_index(0, 0)));
    EXPECT_EQ(after(array->flat_index(0, 1)), before(array->flat_index(0, 1)));
    EXPECT_TRUE(std::isnan(after(array->flat_index(1, 0))));
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
    EXPECT_TRUE(std::isfinite(array->temperature()(failed_loc)));

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
    EXPECT_FALSE(array->ok()(loc));
    EXPECT_TRUE(array->ok()(array->flat_index(0, 0)));
    EXPECT_TRUE(std::isnan(array->temperature()(loc)));
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
    const std::vector<double> state_of_entry_2 = array->get_state(array->flat_index(1, 0));
    array->set_status(array->flat_index(1, 1), SolveStatus::failure(FailureReason::NOT_CONVERGED, "earlier failure"));

    // Enthalpy of this mixture at 600 K, from an independent Cantera phase
    auto reference = Cantera::newSolution("h2o2.yaml", "ohmech");
    reference->thermo()->setMoleFractionsByName("H2:2, O2:1");
    reference->thermo()->setState_TP(600.0, pressures[0]);
    const double enthalpy_at_600_K = reference->thermo()->enthalpy_mass();
    Eigen::ArrayXd enthalpies(2);
    // The second enthalpy is far beyond what any temperature in the data range gives
    enthalpies << enthalpy_at_600_K, 1.0e14;

    EXPECT_TRUE(throws_solve_error([&] { array->HP(enthalpies, pressures); },
                                   FailureReason::THERMO_FAILURE));

    // The entries of the first row were written by successful steps of HP (600 K, not the 300 K
    // they held before and not whatever temperature the rejected iterate reached), and are ok.
    for (int j = 0; j < 2; j++) {
        const int loc = array->flat_index(0, j);
        EXPECT_TRUE(is_ok(array->status_at(loc)));
        EXPECT_NEAR(array->temperature()(loc), 600.0, 1e-6);
        EXPECT_NEAR(array->enthalpy_mass()(loc), enthalpy_at_600_K, 1e-6 * std::abs(enthalpy_at_600_K));
    }
    const std::vector<double> state_of_entry_0 = array->get_state(0);
    // A second read gives the same entry: reading does not disturb it
    EXPECT_EQ(array->get_state(0), state_of_entry_0);
    EXPECT_NEAR(array->temperature()(0), 600.0, 1e-6);

    // The rejected entry and those after it keep their previous state and status
    const int rejected = array->flat_index(1, 0);
    EXPECT_TRUE(is_ok(array->status_at(rejected)));
    EXPECT_EQ(array->get_state(rejected), state_of_entry_2);
    EXPECT_NEAR(array->temperature()(rejected), 400.0, 1e-9);
    const int later = array->flat_index(1, 1);
    EXPECT_TRUE(failed_with(array->status_at(later), FailureReason::NOT_CONVERGED, "earlier failure"));
    EXPECT_EQ(array->num_failed(), 1);
    // The later entry is still failed, so its state cannot be read
    EXPECT_TRUE(throws_solve_error([&] { array->get_state(later); }, FailureReason::NOT_CONVERGED));
}

TEST(ThermoArrayThreeDimensional, okIsInStorageOrder) {
    auto sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    Goddard::ThermoArray array(Goddard::Gas(sln), std::vector<long>{2, 2, 2});
    array.set_status(array.flat_index(1, 0, 1), SolveStatus::failure(FailureReason::NO_SOLUTION, "x"));
    const auto ok = array.ok();
    ASSERT_EQ(ok.size(), 8);
    EXPECT_FALSE(ok(array.flat_index(1, 0, 1)));
    EXPECT_EQ(ok.count(), 7);
    EXPECT_EQ(array.num_failed(), 1);
}
