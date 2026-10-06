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


/** Standard gravity [m/s^2], for specific impulse in seconds. */
constexpr double STANDARD_GRAVITY = 9.80665;

/** Rocket performance of one exit, as CEA reports it. */
struct RocketPerformance {
    /** Chamber (or stagnation) pressure over exit pressure, Pc/Pe [-]. */
    double pressure_ratio;
    /** Exit area over throat area, Ae/At [-]. */
    double area_ratio;
    /** Exit Mach number [-]. */
    double mach_number;
    /** Characteristic velocity c* [m/s]. */
    double cstar;
    /** Thrust coefficient at optimum expansion [-]. */
    double CF;
    /** Specific impulse at optimum expansion, as an effective exhaust velocity [m/s]. */
    double isp;
    /** Vacuum specific impulse, as an effective exhaust velocity [m/s]. */
    double ivac;

    /** Specific impulse at optimum expansion in seconds, `isp / STANDARD_GRAVITY` [s]. */
    double isp_s() const { return isp / STANDARD_GRAVITY; }
    /** Vacuum specific impulse in seconds, `ivac / STANDARD_GRAVITY` [s]. */
    double ivac_s() const { return ivac / STANDARD_GRAVITY; }
};

/** Kind of a station of a rocket problem. */
enum class StationType {
    /// Chamber state. For a finite-area combustor, the injector face.
    CHAMBER,
    /// Finite-area combustor only: stagnation state "inf" at P_inf.
    STAGNATION,
    /// Finite-area combustor only: end of the constant-area chamber.
    COMBUSTION_END,
    /// Sonic throat.
    THROAT,
    /// Nozzle exit station, one per expansion ratio.
    EXIT
};

/** One station of one operating point of a rocket problem. */
struct RocketStation {
    /// Name of the case the station belongs to.
    std::string case_name;
    /// Kind of station.
    StationType type;
    /// Index into the case's mixture ratios (`ChemicalParameters::OF_ratios`) [-].
    std::size_t of_index;
    /// Index into the case's pressures (`CombustorOptions::pressures`) [-].
    std::size_t pressure_index;
    /// Index into the case's `NozzleOptions::expansion_ratios` [-]; only meaningful for EXIT.
    std::size_t expansion_index;
    /**
     * A/A_t [-]: 0 for CHAMBER and STAGNATION, A_c/A_t for COMBUSTION_END, 1 for THROAT, Ae/At
     * for EXIT (also for pressure-ratio exits).
     */
    double area_ratio;
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
    /** How `OF_ratios` are interpreted; see `CombustorOptions::mixture_type`. */
    MixtureRatioType mixture_type = MixtureRatioType::OF_RATIO;
};

/**
 * Results of a `RocketProblem`: a flat list of stations over every case, mixture ratio and
 * pressure, with accessors by operating point and a CEA-style report.
 */
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

    /** Every station of every case, as one flat list. */
    const std::vector<RocketStation>& stations() const { return m_stations; }

    /**
     * All stations of one type in one case.
     *
     * @param type Station type.
     * @param case_name Case to read. May be empty when the results hold a single case.
     * @throws std::runtime_error if the case does not exist, or if `case_name` is empty and the
     * results hold several cases.
     */
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

    /** Names of all cases, in no particular order. */
    std::vector<std::string> case_names() const;

    /**
     * Mixture ratios of a case [-], interpreted through its `CombustorOptions::mixture_type`.
     *
     * @param case_name Case to read. May be empty when the results hold a single case.
     */
    const std::vector<double>& of_ratios(const std::string& case_name = "") const;

    /**
     * CEA-style formatted text report.
     *
     * @param case_name Case to report. Empty reports every case, sorted by name.
     */
    std::string report(const std::string& case_name = "") const;

    /**
     * Rocket performance from the stagnation, throat and exit states.
     *
     * c* = P_0 / (rho_t a_t), with a_t the throat speed of sound; the exit velocity follows from
     * the enthalpy drop h_0 - h_e, the area ratio from mass conservation, and CF = u_e / c*.
     *
     * @param chamber Stagnation state the nozzle expands from (the chamber of an infinite-area
     * combustor).
     * @param throat Throat state.
     * @param exit Exit state.
     */
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
        MixtureRatioType mixture_type = MixtureRatioType::OF_RATIO;
    };
    std::unordered_map<std::string, CaseMeta> m_case_meta;
    /** Product mixture used to read the stored station states back. */
    Gas m_gas;

    // Resolve case_name: if empty, returns the single case name; throws if ambiguous.
    std::string resolve_case(const std::string& case_name) const;
};


} // namespace Goddard
