#pragma once
#include "cantera/core.h"
#include "cantera/equil/MultiPhase.h"

#include "eigen3/Eigen/Dense"
#include "goddard/mixture_ratio.hpp"
#include "goddard/thermoarray.hpp"
#include "goddard/utils.hpp"
#include <memory>
#include <vector>
#include <utility>

namespace Goddard {

enum class CombustorType {
    INFINITE_AREA,
    FINITE_MASS_FLUX,
    FINITE_CONTRACTION_RATIO
};

struct CombustionOptions {
    CombustorType combustor_type;
    bool include_transport = false;
    bool include_ionized = false;
    double abstol = DEFAULT_ABSTOL;
    double reltol = DEFAULT_RELTOL;
    double trace_conc = DEFAULT_TRACE_CONCENTRATION;
};

/**
 * @brief Handles combustion reactions.
 */
class Combustor {
    public:
    //this may lead to lots of unnecessary copies
    Combustor(std::shared_ptr<Cantera::Solution> fuel, std::shared_ptr<Cantera::Solution> oxidizer, Eigen::ArrayXd&& temperatures, Eigen::ArrayXd&& pressures, MixtureRatios&& mr);


    ThermoArray solve(CombustionOptions options = {});

    protected:
    Eigen::ArrayXXd generate_mole_fraction_matrix();
    
    std::shared_ptr<Cantera::Solution> fuel;
    std::shared_ptr<Cantera::Solution> oxidizer;
    Cantera::MultiPhase feed;
    //TODO: add a "reactant" solution object
    std::shared_ptr<Cantera::Solution> products;

    MixtureRatios mr;
    Eigen::ArrayXd temperatures;
    Eigen::ArrayXd pressures;

};

} //namespace Goddard