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

enum class StationType {
    CHAMBER,
    THROAT,
    EXIT
};

struct RocketStation {
    std::string case_name;
    StationType type;
    std::size_t of_index;        // index into the case's OF_ratios vector
    std::size_t pressure_index;  // index into the case's pressures vector
    std::size_t expansion_index; // index into expansion_ratios; only meaningful for EXIT
    double area_ratio;           // 0.0 for CHAMBER, 1.0 for THROAT, >1 for EXIT
    ThermodynamicState thermo;
    bool converged;
};

// Internal type: used by RocketProblem::solve() to pass results to RocketProblemResults.
// Not part of the public API.
struct RocketProblemCaseResult {
    std::string problem_type;
    ThermoArray inlet_states;
    GasChemistry chemistry;
    std::vector<NozzleResults> nozzle_states;
    std::vector<double> OF_ratios;
    std::vector<double> pressures;
    std::vector<double> expansion_ratios;
    ExpansionType expansion_type;
};

class RocketProblemResults {

public:
    RocketProblemResults(std::unordered_map<std::string, RocketProblemCaseResult>&& case_results,
                         std::shared_ptr<Cantera::Solution> sln);

    // Direct access to the flat station list.
    const std::vector<RocketStation>& stations() const { return m_stations; }

    // All stations of a given type, optionally filtered to a single case.
    // If case_name is empty and there is exactly one case, that case is used.
    std::vector<RocketStation> stations_of_type(StationType type, const std::string& case_name = "") const;

    // Single-station convenience accessors.
    // If case_name is empty and there is exactly one case, that case is used.
    const RocketStation& chamber(std::size_t of_index = 0, const std::string& case_name = "") const;
    const RocketStation& throat(std::size_t of_index = 0, const std::string& case_name = "") const;
    std::vector<RocketStation> exits(std::size_t of_index = 0, const std::string& case_name = "") const;

    // Compute rocket performance for a given operating point and exit station.
    // If case_name is empty and there is exactly one case, that case is used.
    RocketPerformance performance(std::size_t of_index = 0, std::size_t exit_index = 0,
                                  const std::string& case_name = "") const;

    // List all case names.
    std::vector<std::string> case_names() const;

    // The OF ratios used for a case.
    // If case_name is empty and there is exactly one case, that case is used.
    const std::vector<double>& of_ratios(const std::string& case_name = "") const;

    // Generate a CEA-style formatted text report.
    // If case_name is empty, all cases are reported in sorted order.
    std::string report(const std::string& case_name = "") const;

    // Compute performance metrics from individual thermo states.
    static RocketPerformance calculate_performance(
        const ThermodynamicState& chamber,
        const ThermodynamicState& throat,
        const ThermodynamicState& exit);

private:
    std::vector<RocketStation> m_stations;

    struct CaseMeta {
        std::vector<double> of_ratios;
        std::vector<double> pressures;
        std::vector<double> expansion_ratios;
        GasChemistry chemistry;
        ExpansionType expansion_type;
    };
    std::unordered_map<std::string, CaseMeta> m_case_meta;
    std::shared_ptr<Cantera::Solution> m_sln;

    inline std::shared_ptr<Cantera::ThermoPhase> thermo() { return m_sln->thermo(); }

    // Resolve case_name: if empty, returns the single case name; throws if ambiguous.
    std::string resolve_case(const std::string& case_name) const;
};


} // namespace Goddard
