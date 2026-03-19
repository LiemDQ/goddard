#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>
#include <optional>
#include "cantera/core.h"
#include "goddard/combustor.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/thermo.hpp"

namespace Goddard {

struct RocketCaseParameters {
    std::string name;
    std::string problem_type;
    CombustorOptions combustor_options;
    NozzleOptions nozzle_options;
};

struct Speciation { //TODO: implement speciation functionality
    std::unordered_set<std::string> species;
    std::unordered_set<std::string> elements;
    int max_species = 10;
};

struct ChemicalParameters {
    std::string thermo_file;
    std::unordered_set<std::string> species;
    ThermodynamicState cantera_fuel_state;
    ThermodynamicState cantera_oxidizer_state;
    std::vector<double> OF_ratios;
    std::vector<double> phi_ratios;
    std::vector<double> fuel_weight_percentages;
    std::vector<double> valance_equivalences;
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

// Forward declaration
class RocketProblemResults;

// TODO: implement this function
// RocketPerformance performance_from_state();


class RocketProblem {

    public:
    
    RocketProblem(const ChemicalParameters& chem_params,
        const std::vector<RocketCaseParameters>& cases, 
        const std::string& name = "", 
        bool transport = false,
        bool ionized_species = false,
        double trace = 1e-6);
    
    RocketProblemResults solve();
    
    inline std::shared_ptr<Cantera::Solution> solution() { return m_sln; }
    inline std::shared_ptr<Cantera::ThermoPhase> thermo() { return m_sln->thermo(); }
    bool include_transport = false; //TODO: implement transport functionality
    bool include_ionized_species = false; //TODO: implement ionized species
    double trace_cutoff = 1e-6;
    std::vector<RocketCaseParameters> problem_cases;
    ChemicalParameters chemical_params;
    
    private:
    std::shared_ptr<Cantera::Solution> m_sln;

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
    RocketProblemResults(std::unordered_map<std::string, RocketProblemCaseResult>&& case_results, std::shared_ptr<Cantera::Solution> sln);

    /**
     * Extract all thermo states for a case at a given O/F index.
     * Returns [inlet, throat, exit1, exit2, ...] in order.
     */
    std::vector<ThermoStateInfo> extract_thermo_info(const std::string& case_name, std::size_t index);

    /**
     * Get chamber/inlet state for a specific case and O/F index.
     * @param case_name Name of the case
     * @param of_index Index into the O/F ratio array
     * @return Chamber state, or nullopt if not found
     */
    std::optional<ThermoStateInfo> get_chamber_state(const std::string& case_name, std::size_t of_index);

    /**
     * Get throat state for a specific case and O/F index.
     * @param case_name Name of the case
     * @param of_index Index into the O/F ratio array
     * @return Throat state, or nullopt if not found
     */
    std::optional<ThermoStateInfo> get_throat_state(const std::string& case_name, std::size_t of_index);

    /**
     * Get all exit states for a specific case and O/F index.
     * @param case_name Name of the case
     * @param of_index Index into the O/F ratio array
     * @return Vector of exit states (may be empty if no exits)
     */
    std::vector<ThermoStateInfo> get_exit_states(const std::string& case_name, std::size_t of_index);

    /**
     * Calculate rocket performance metrics from chamber, throat, and exit states.
     * @param chamber Chamber/inlet state
     * @param throat Throat state
     * @param exit Exit state
     * @return Performance metrics (cstar, CF, Isp, Ivac, etc.)
     */
    static RocketPerformance calculate_performance(
        const ThermoStateInfo& chamber,
        const ThermoStateInfo& throat,
        const ThermoStateInfo& exit);

    std::string report(const std::string& case_name = "") const;

    std::unordered_map<std::string, RocketProblemCaseResult> cases;

    private:
    std::shared_ptr<Cantera::Solution> m_sln;

    inline std::shared_ptr<Cantera::ThermoPhase> thermo() {return m_sln->thermo();}
};



} //namespace Goddard