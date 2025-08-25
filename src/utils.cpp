#include "goddard/utils.hpp"
#include "cantera/base/AnyMap.h"
#include "cantera/base/global.h"
#include "cantera/core.h"
#include <vector>

using Cantera::AnyMap;

namespace Goddard {

std::vector<double> save_thermo_state(const Cantera::ThermoPhase& thermo) {
    std::vector<double> state(thermo.stateSize());
    thermo.saveState(state);
    return state;
}

std::shared_ptr<Cantera::Solution> create_mixture_solution(Cantera::Solution& sol1, Cantera::Solution& sol2){
    auto thermo1 = sol1.thermo();
    auto thermo2 = sol2.thermo();

    AnyMap params = thermo1->parameters();
    AnyMap params2 = thermo2->parameters();
    params.update(params2);

    
    std::shared_ptr<Cantera::Solution> new_solution = Cantera::newSolution(params);
    return new_solution;
}

std::shared_ptr<Cantera::Solution> copy_solution(Cantera::Solution& sln) {
    std::vector<double> state(sln.thermo()->stateSize());
    sln.thermo()->saveState(state);

    std::shared_ptr<Cantera::Solution> p_new_sln = Cantera::newSolution(sln.source(), sln.name());
    p_new_sln->thermo()->restoreState(state);
    //TODO: is this appropriate?
    p_new_sln->setTransport(sln.transport());
    p_new_sln->setKinetics(sln.kinetics());

    return p_new_sln;
}

// std::shared_ptr<Cantera::Solution> speciate(const std::string& source, const std::vector<std::string>& elements) {
    
// }

double check_reltol(double reltol) {
    if (reltol <= 0) {
        Cantera::warn_user("check_reltol", "reltol is zero or negative. Setting reltol to default value: {0}", DEFAULT_RELTOL);
        return DEFAULT_RELTOL;
    }
    else return reltol;
}

double check_abstol(double abstol) {
    if (abstol <= 0) {
        Cantera::warn_user("check_abstol", "abstol is zero or negative. Setting abstol to default value: {0}", DEFAULT_ABSTOL);
        return DEFAULT_ABSTOL;
    }
    else return abstol;
}



}