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

std::vector<double> build_state_TP(Cantera::Solution& sln, double T, double P, std::string& composition);

std::shared_ptr<Cantera::Solution> create_mixture_solution(Cantera::Solution& sol1, Cantera::Solution& sol2);

std::shared_ptr<Cantera::Solution> copy_solution(Cantera::Solution& sln);

Eigen::ArrayXd vector_to_eigenarray(std::vector<double>& vec);

/**
 * @brief Get molar mass from a state vector by restoring state and computing mean MW.
 */
double molar_mass_from_composition(Cantera::ThermoPhase& thermo, const std::vector<double>& state);

}