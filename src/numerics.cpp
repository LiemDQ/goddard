#include "goddard/numerics.hpp"
#include "cantera/core.h"

namespace Goddard {
    
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