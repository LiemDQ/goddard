#pragma once
#include <cmath>

namespace Goddard {

constexpr double DEFAULT_RELTOL = 1e-6;
constexpr double DEFAULT_ABSTOL = 1e-9;
constexpr double DEFAULT_TRACE_CONCENTRATION_THRESHOLD = 1e-5;

double check_reltol(double reltol); 
double check_abstol(double abstol); 


constexpr double max_fp_error(double val, double reltol = DEFAULT_RELTOL, double abstol = DEFAULT_ABSTOL) {
    double rel = std::abs(val*reltol);
    return rel >= abstol ? rel : abstol;
}

}