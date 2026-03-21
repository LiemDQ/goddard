#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>

#include "cantera/core.h"
#include "goddard/nozzle.hpp"
#include "goddard/thermo.hpp"
#include "goddard/thermoarray.hpp"

namespace Goddard {


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
    ThermoArray inlet_states;
    NozzleChemistryType chemistry;
    std::vector<NozzleResults> nozzle_states;
    std::vector<double> OF_ratios;
    std::vector<double> pressures;
    std::vector<double> expansion_ratios;
    ExpansionType expansion_type;
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

    std::string report(const std::string& case_name = "");

    std::unordered_map<std::string, RocketProblemCaseResult> cases;

    private:
    std::shared_ptr<Cantera::Solution> m_sln;

    inline std::shared_ptr<Cantera::ThermoPhase> thermo() {return m_sln->thermo();}
};


} // namespace Goddard