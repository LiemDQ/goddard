#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include <memory>

#include "cantera/core.h"
#include "goddard/combustor.hpp"
#include "goddard/gas.hpp"
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
    CHAMBER,        ///< Chamber state. For a finite-area combustor, the injector face.
    STAGNATION,     ///< Finite-area combustor only: stagnation state "inf" at P_inf.
    COMBUSTION_END, ///< Finite-area combustor only: end of the constant-area chamber.
    THROAT,
    EXIT
};

struct RocketStation {
    std::string case_name;
    StationType type;
    std::size_t of_index;        // index into the case's OF_ratios vector
    std::size_t pressure_index;  // index into the case's pressures vector
    std::size_t expansion_index; // index into expansion_ratios; only meaningful for EXIT
    double area_ratio;           // 0.0 for CHAMBER/STAGNATION, A_c/A_t for COMBUSTION_END, 1.0 for THROAT, >1 for EXIT
    /**
     * Mixture state of the station. Its `stagnation_enthalpy` [J/kg] is the enthalpy of the
     * stagnation state of the station's operating point (see `RocketProblemResults::stagnation`),
     * which the adiabatic expansion conserves.
     */
    ThermodynamicState thermo;
};

// Internal type: used by RocketProblem::solve() to pass results to RocketProblemResults.
// Not part of the public API.
struct RocketProblemCaseResult {
    ThermoArray inlet_states;
    GasChemistry chemistry;
    std::vector<NozzleResults> nozzle_states;
    std::vector<double> OF_ratios;
    std::vector<double> pressures;
    std::vector<double> expansion_ratios;
    ExpansionType expansion_type;
    CombustionProcess process;
    CombustorType combustor_type = CombustorType::INFINITE_AREA;
    /** Finite-area chambers, indexed like `nozzle_states`; empty for `INFINITE_AREA`. */
    std::vector<FiniteAreaChamber> finite_area_chambers;
    /** Freezing station of a FROZEN nozzle; see `NozzleOptions::frozen_NFZ`. */
    int frozen_NFZ = 0;
};

class RocketProblemResults {

public:
    /**
     * Build the flat station list from the results of every case.
     *
     * @param case_results Per-case combustion and nozzle states, consumed by this call.
     * @param gas Product mixture the states were computed with. It carries the candidate
     *            condensed species, so each station's state vector is read back with the
     *            condensed amounts it was saved with.
     */
    RocketProblemResults(std::unordered_map<std::string, RocketProblemCaseResult>&& case_results,
                         Gas gas);

    // Direct access to the flat station list.
    const std::vector<RocketStation>& stations() const { return m_stations; }

    // All stations of a given type, optionally filtered to a single case.
    // If case_name is empty and there is exactly one case, that case is used.
    std::vector<RocketStation> stations_of_type(StationType type, const std::string& case_name = "") const;

    /**
     * Chamber station of one operating point. For a finite-area combustor, the injector face.
     *
     * @param of_index Index into the case's mixture ratios [-].
     * @param pressure_index Index into the case's chamber pressures [-].
     * @param case_name Case to read. May be empty when the results hold a single case.
     * @throws std::runtime_error if the case or the station does not exist.
     */
    const RocketStation& chamber(std::size_t of_index = 0, std::size_t pressure_index = 0,
                                 const std::string& case_name = "") const;

    /**
     * Throat station of one operating point.
     *
     * @param of_index Index into the case's mixture ratios [-].
     * @param pressure_index Index into the case's chamber pressures [-].
     * @param case_name Case to read. May be empty when the results hold a single case.
     * @throws std::runtime_error if the case or the station does not exist.
     */
    const RocketStation& throat(std::size_t of_index = 0, std::size_t pressure_index = 0,
                                const std::string& case_name = "") const;

    /**
     * Exit stations of one operating point, ordered by expansion index.
     *
     * @param of_index Index into the case's mixture ratios [-].
     * @param pressure_index Index into the case's chamber pressures [-].
     * @param case_name Case to read. May be empty when the results hold a single case.
     * @throws std::runtime_error if the case does not exist.
     */
    std::vector<RocketStation> exits(std::size_t of_index = 0, std::size_t pressure_index = 0,
                                     const std::string& case_name = "") const;

    /**
     * Stagnation state that the nozzle expands from: the "inf" state of a finite-area combustor,
     * or the chamber state of an infinite-area combustor.
     *
     * @param of_index Index into the case's mixture ratios [-].
     * @param pressure_index Index into the case's chamber pressures [-].
     * @param case_name Case to read. May be empty when the results hold a single case.
     * @throws std::runtime_error if the case or the station does not exist.
     */
    const RocketStation& stagnation(std::size_t of_index = 0, std::size_t pressure_index = 0,
                                    const std::string& case_name = "") const;

    /**
     * Combustion-end station of a finite-area combustor.
     *
     * @param of_index Index into the case's mixture ratios [-].
     * @param pressure_index Index into the case's chamber pressures [-].
     * @param case_name Case to read. May be empty when the results hold a single case.
     * @throws std::runtime_error if the case uses an infinite-area combustor, or if the case or
     * the station does not exist.
     */
    const RocketStation& combustion_end(std::size_t of_index = 0, std::size_t pressure_index = 0,
                                        const std::string& case_name = "") const;

    /**
     * Rocket performance of one operating point at one exit station, from its stagnation,
     * throat and exit states (see `calculate_performance`).
     *
     * @param of_index Index into the case's mixture ratios [-].
     * @param pressure_index Index into the case's chamber pressures [-].
     * @param exit_index Index into the exit stations of the operating point [-].
     * @param case_name Case to read. May be empty when the results hold a single case.
     * @throws std::runtime_error if the case or a station does not exist, or if `exit_index`
     * is out of range.
     */
    RocketPerformance performance(std::size_t of_index = 0, std::size_t pressure_index = 0,
                                  std::size_t exit_index = 0,
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
        CombustionProcess process;
        CombustorType combustor_type;
        /**
         * Mass flux mdot/Ac [kg/(m^2 s)] of each operating point, row-major over
         * (of_index, pressure_index); empty for `CombustorType::INFINITE_AREA`. Not derivable
         * from a station's `ThermodynamicState`, so it is captured here for the report.
         */
        std::vector<double> mass_flux;
        /** Freezing station of a FROZEN nozzle; see `NozzleOptions::frozen_NFZ`. */
        int frozen_NFZ = 0;
    };
    std::unordered_map<std::string, CaseMeta> m_case_meta;
    /** Product mixture used to read the stored station states back. */
    Gas m_gas;

    // Resolve case_name: if empty, returns the single case name; throws if ambiguous.
    std::string resolve_case(const std::string& case_name) const;
};


} // namespace Goddard
