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
    CombustorOptions combustor_options;
    NozzleOptions nozzle_options;
    
};

struct ChemicalParameters {
    std::string thermo_file;
    std::vector<std::string> species;
    std::vector<double> cantera_fuel_state;
    std::vector<double> cantera_oxidizer_state;
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
    bool include_ionized_species = false;
};

struct RocketState {
    std::string name;
    std::vector<double> cantera_state;
    double pressure_ratio;
    double area_ratio;
    double dlv_dlp_t;
    double dlv_dlt_p;
    double gamma_s;
    double speed_of_sound;
};

struct ThermoStateInfo {
    double pressure;
    double temperature;
    double density;
    double enthalpy;
    double internal_energy;
    double gibbs;
    double entropy;
    double molecular_weight;
    double gamma_s;
    double dlV_dlP_T;
    double dlV_dlT_P;
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

RocketPerformance performance_from_state() {

}
class RocketProblemBuilder {


};

class RocketProblem : public std::enable_shared_from_this<RocketProblem> {

    public:
    static std::shared_ptr<RocketProblem> create(const std::string& thermo_file);

    RocketProblemResults solve();
    
    inline std::shared_ptr<Cantera::Solution> solution() { return m_solution; }
    inline std::shared_ptr<Cantera::ThermoPhase> thermo() { return m_solution->thermo(); }
    bool include_transport = false;
    bool include_ionized_species = false;
    double m_trace_cutoff = 1e-6;
    
    private:
    RocketProblem(const std::string& thermo_file);
    std::vector<RocketCaseParameters> m_problem_cases;
    ChemicalParameters m_chemical_params;
    std::shared_ptr<Cantera::Solution> m_solution;

};


struct RocketProblemCaseResult {
    std::string problem_type;
    ThermoArray inlet_states;
    NozzleChemistryType chemistry;
    std::vector<NozzleResults> nozzle_states;
    // std::vector<RocketPerformance> nozzle_performance;
};

class RocketProblemResults {
    
    public: 
    RocketProblemResults(std::shared_ptr<RocketProblem>&& prob, std::unordered_map<std::string, RocketProblemCaseResult>&& cases);
    std::vector<ThermoStateInfo> extract_thermo_info(const std::string& case_name, std::size_t index);

    private:
    std::shared_ptr<RocketProblem> m_prob;
    std::unordered_map<std::string, RocketProblemCaseResult> m_cases;

    inline std::shared_ptr<Cantera::ThermoPhase> thermo() {m_prob->thermo();}
};



} //namespace Goddard