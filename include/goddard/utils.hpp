#pragma once
#include "cantera/core.h"
#include <memory>
#include <cmath>
#include <algorithm>

namespace Goddard {

std::shared_ptr<Cantera::Solution> create_mixture_solution(Cantera::Solution& sol1, Cantera::Solution& sol2);

std::shared_ptr<Cantera::Solution> speciate(const std::string& infile, const std::string& name);

constexpr double DEFAULT_RELTOL = 1e-6;
constexpr double DEFAULT_ABSTOL = 1e-9;
constexpr double DEFAULT_TRACE_CONCENTRATION_THRESHOLD = 1e-5;

double check_reltol(double reltol); 
double check_abstol(double abstol); 


constexpr double max_fp_error(double val, double reltol = DEFAULT_RELTOL, double abstol = DEFAULT_ABSTOL) {
    return std::max(abs(val*reltol), abstol);
}

}