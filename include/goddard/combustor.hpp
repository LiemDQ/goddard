#pragma once
#include "cantera/core.h"
#include "cantera/equil/MultiPhase.h"

#include "eigen3/Eigen/Dense"
#include "goddard/mixture_ratio.hpp"
#include "goddard/thermo_states.hpp"
#include <memory>
#include <vector>
#include <utility>

namespace Goddard {

class Combustor {
    public:
    //this may lead to lots of unnecessary copies
    Combustor(std::shared_ptr<Cantera::Solution> fuel, std::shared_ptr<Cantera::Solution> oxidizer, Eigen::ArrayXd&& pressures, MixtureRatios&& mr);


    ThermoArray solve();

    protected:
    std::vector<Eigen::ArrayXd> generate_mole_fraction_tensor();
    
    std::shared_ptr<Cantera::Solution> fuel;
    std::shared_ptr<Cantera::Solution> oxidizer;
    Cantera::MultiPhase feed;
    //TODO: add a "reactant" solution object

    MixtureRatios mr;
    Eigen::ArrayXd feed_temperatures;
    Eigen::ArrayXd pressures;

};

} //namespace Goddard