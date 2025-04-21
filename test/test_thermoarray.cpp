#include "goddard/thermoarray.hpp"
#include "goddard/thermo.hpp"
#include "goddard/utils.hpp"

#include "cantera/core.h"
#include "gtest/gtest.h"
#include "eigen3/Eigen/Dense"
#include <memory>
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
    EXPECT_EQ(array->size(), array->solutionarray()->size()) 
        << "Array should have the same size as its underlying `SolutionArray`";
    EXPECT_EQ(array->solution(), sln) 
        << "Array Solution pointer should point to the underlying Solution used to construct it.";

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
    EXPECT_EQ(array->size(), array->solutionarray()->size()) 
        << "Array should have the same size as its underlying `SolutionArray`";
    EXPECT_EQ(array->solution(), sln) 
        << "Array Solution poii=nter should point to the underlying Solution used to construct it.";

    EXPECT_EQ(array->solution()->thermo()->nSpecies(), NUM_H2O2_SPECIES)
     << "ohmech should have 10 species";
    
}

TEST_F(ThermoArray2DTests, broadcastTD) {
    ASSERT_NO_THROW(array->TD(temperatures, densities));
    auto states = array->solutionarray();
    
    int loc = 0;
    
    ASSERT_EQ(states->size(), densities.size()*temperatures.size()) << "SolutionArray must have same size as thermodynamic states.";
    ASSERT_GE(states->getState(0).size(), 2) << "State vector should have at least 2 entries";

    for (int j = 0; j < densities.size(); j++) {
        for (int i = 0; i < temperatures.size(); i++){
            std::vector<double> state = states->getState(loc);
            EXPECT_DOUBLE_EQ(state[0], temperatures(i)) << "i = " << i << ", j = " << j;
            EXPECT_DOUBLE_EQ(state[1], densities(j)) << "i = " << i << ", j = " << j; 
            loc++;
        }
    }
}

TEST_F(ThermoArray2DTests, broadcastTP) {
    
    
    ASSERT_NO_THROW(array->TP(temperatures, pressures));
    auto states = array->solutionarray();
    
    int loc = 0;
    
    ASSERT_EQ(states->size(), pressures.size()*temperatures.size()) << "SolutionArray must have same size as thermodynamic states.";
    ASSERT_GE(states->getState(0).size(), 2) << "State vector should have at least 2 entries";

    auto MW = sln->thermo()->meanMolecularWeight();

    for (int j = 0; j < pressures.size(); j++) {
        for (int i = 0; i < temperatures.size(); i++){
            std::vector<double> state = states->getState(loc);
            EXPECT_DOUBLE_EQ(state[0], temperatures(i)) << "i = " << i << ", j = " << j;
            
            //the state vector stores pressure implicitly as density, so needs to be calculated.
            double density = state[1];
            double calculated_pressure = ideal_gas_D_to_P(density, temperatures(i), MW);
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

    auto states = array->solutionarray();
    
    int loc = 0;
    
    ASSERT_EQ(states->size(), pressures.size()*temperatures.size()) << "SolutionArray must have same size as thermodynamic states.";
    ASSERT_GE(states->getState(0).size(), 2) << "State vector should have at least 2 entries";

    std::shared_ptr<Cantera::Solution> ref_sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    std::shared_ptr<Cantera::ThermoPhase> ref_thermo = ref_sln->thermo();
    ref_thermo->setMoleFractionsByName("H2:1, N2:1, O2:1, AR:0.1");

    // auto MW = sln->thermo()->meanMolecularWeight();
    std::vector<double> prev_state(ref_thermo->stateSize());
    ref_thermo->saveState(prev_state);

    std::vector<double> ref_state(ref_thermo->stateSize());


    for (int j = 0; j < pressures.size(); j++) {
        for (int i = 0; i < temperatures.size(); i++){
            std::vector<double> state = states->getState(loc);

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

        compositions = Eigen::ArrayXXd(NUM_H2O2_SPECIES, NUM_COMPOSITIONS);
        for (size_t i = 0; i < NUM_COMPOSITIONS; i++) {
            double o2 = static_cast<double>(NUM_COMPOSITIONS - i);
            //"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"
            compositions.row(i) << static_cast<double>(i), 1e-10, 1e-10, o2, 1e-10, 0.1, 1e-10, 0.5, 1.0;
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
    auto states = array->solutionarray();

    std::shared_ptr<Cantera::Solution> ref_sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    std::shared_ptr<Cantera::ThermoPhase> ref_thermo = ref_sln->thermo();

    
    ASSERT_EQ(states->size(), pressures.size()*temperatures.size()*compositions.rows()) 
        << "SolutionArray must have same size as thermodynamic states.";
    ASSERT_GE(states->getState(0).size(), 2) 
        << "State vector should have at least 2 entries";
    
    
    int loc = 0;
    std::vector<double> ref_state(ref_thermo->stateSize());

    for (int k = 0; k < compositions.rows(); k++) {
        
        auto composition = compositions.row(k);
        ref_thermo->setMoleFractions(composition.data());

        for (int j = 0; j < pressures.size(); j++) {

            
            for (int i = 0; i < temperatures.size(); i++){
                std::vector<double> state = states->getState(loc);
                
                ref_thermo->setTemperature(temperatures(i));
                ref_thermo->setPressure(pressures(j));
                ref_thermo->saveState(ref_state);

                for (size_t l = 0; l < state.size(); l++) {
                    EXPECT_DOUBLE_EQ(state[l], ref_state[l]) 
                        << "State vectors should be identical." 
                        << "Indices (i,j,k,l): (" << i << "," << j << "," << k << "," << l << ")";
                }
                loc++;
            }
        }
    }
}


TEST_F(ThermoArray3DTests, equilibrateTPX){
    ASSERT_NO_THROW(array->TPX(temperatures, pressures, compositions));
    array->equilibrate("HP", "gibbs", Goddard::DEFAULT_RELTOL);

    auto states = array->solutionarray();

    std::shared_ptr<Cantera::Solution> ref_sln = Cantera::newSolution("h2o2.yaml", "ohmech");
    std::shared_ptr<Cantera::ThermoPhase> ref_thermo = ref_sln->thermo();

    
    ASSERT_EQ(states->size(), pressures.size()*temperatures.size()*compositions.rows()) 
        << "SolutionArray must have same size as thermodynamic states.";
    ASSERT_GE(states->getState(0).size(), 2) 
        << "State vector should have at least 2 entries";
    
    
    int loc = 0;
    std::vector<double> ref_state(ref_thermo->stateSize());

    for (int k = 0; k < compositions.rows(); k++) {
        
        auto composition = compositions.row(k);
        
        for (int j = 0; j < pressures.size(); j++) {       
            for (int i = 0; i < temperatures.size(); i++){
                
                std::vector<double> state = states->getState(loc);
                
                ref_thermo->setMoleFractions(composition.data());
                ref_thermo->setTemperature(temperatures(i));
                ref_thermo->setPressure(pressures(j));
                ref_thermo->equilibrate("HP", "gibbs", Goddard::DEFAULT_RELTOL);
                ref_thermo->saveState(ref_state);

                for (size_t l = 0; l < state.size(); l++) {
                    EXPECT_NEAR(state[l], ref_state[l], Goddard::max_fp_error(ref_state[l])) 
                        << "State vectors should be identical." 
                        << "Indices (i,j,k,l): (" << i << "," << j << "," << k << "," << l << ")";
                }
                loc++;
            }
        }
    }
}