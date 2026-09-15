#include "goddard/thermoarray.hpp"
#include "goddard/thermo.hpp"
#include "goddard/utils.hpp"
#include "goddard/numerics.hpp"

#include "cantera/core.h"
#include "gtest/gtest.h"
#include "eigen3/Eigen/Dense"
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
