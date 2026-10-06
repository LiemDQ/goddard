#pragma once
#include "cantera/core.h"
#include "cantera/base/AnyMap.h"
#include "eigen3/Eigen/Dense"
#include <memory>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
namespace Goddard {

std::vector<double> save_thermo_state(const Cantera::ThermoPhase& sol);

Cantera::AnyMap load_root_node(const std::string& infile);

Eigen::ArrayXd vector_to_eigenarray(std::vector<double>& vec);

}