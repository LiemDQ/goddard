#include "goddard/utils.hpp"
#include "cantera/base/AnyMap.h"
#include "cantera/base/global.h"
#include "cantera/core.h"

using Cantera::AnyMap;

namespace Goddard {

std::shared_ptr<Cantera::Solution> create_mixture_solution(Cantera::Solution& sol1, Cantera::Solution& sol2){
    auto thermo1 = sol1.thermo();
    auto thermo2 = sol2.thermo();

    AnyMap params;

    AnyMap params = thermo1->parameters();
    AnyMap params2 = thermo2->parameters();
    params.update(params2);
}

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