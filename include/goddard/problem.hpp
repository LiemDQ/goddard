#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <optional>
#include "cantera/core.h"
#include "goddard/case_options.hpp"

namespace Goddard {

struct RocketCaseParameters {
    std::string name;
    std::string problem_type;
    RocketCaseOptions options;
    std::vector<double> expansion_ratios;
    std::vector<double> pressure_ratios;
};

struct ChemicalParameters {
    std::string thermo_file;
    std::vector<std::string> species;
    std::vector<double> fuel_states;
    std::vector<double> oxidizer_states;
    std::vector<double> OF_ratios;
    std::vector<double> phi_ratios;
    std::vector<double> fuel_weight_percentages;
    std::vector<double> valance_equivalences;
};

struct RocketProblemInputs {
    std::vector<RocketCaseParameters> problem_cases;
    ChemicalParameters chemical_params;
    std::shared_ptr<Cantera::Solution> solution;
    bool trace_cutoff = 1e-6;
    bool include_transport = false;
};

struct RocketState {
    std::string name;
    std::vector<double> ct_state;
    double pressure_ratio;
    double area_ratio;
    double dlv_dlp_t;
    double dlv_dlt_p;
    double gamma_s;
    double speed_of_sound;
};

struct StateInfo {
    double pressure;
    double temperature;
    double density;
    double enthalpy;
    double internal_energy;
    double gibbs;
    double entropy;
    double molecular_weight;
    double dlv_dlp_t;
    double dlv_dlt_p;
    double gamma_s;
    double speed_of_sound;
    std::unordered_map<std::string, double> composition;
};

struct RocketPerformance {
    double pressure_ratio;
    double area_ratio;
    double mach_number;
    double cstar;
    double CF;
    double isp;
    double ivac;
};

struct RocketProblemCaseResult {
    std::string problem_type;
    RocketState inlet_state;
    std::vector<RocketState> nozzle_states;
    std::vector<RocketPerformance> nozzle_performance;
};

class RocketProblemResults {
    RocketProblemInputs params;
    std::unordered_map<std::string, RocketProblemCaseResult> cases;

    public: 
    StateInfo get_state_info(const std::string& case_name, const std::string& state_name);
};


} //namespace Goddard