#pragma once
#include "cantera/core.h"
#include "eigen3/Eigen/Dense"
#include <memory>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
namespace Goddard {

std::vector<double> save_thermo_state(const Cantera::ThermoPhase& sol);

std::shared_ptr<Cantera::Solution> create_mixture_solution(Cantera::Solution& sol1, Cantera::Solution& sol2);

std::shared_ptr<Cantera::Solution> copy_solution(Cantera::Solution& sln);

Eigen::ArrayXd vector_to_eigenarray(std::vector<double>& vec);

}