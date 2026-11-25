#pragma once
#include <map>
#include <optional>
#include <string>
#include <vector>

#include "goddard/combustor.hpp"
#include "goddard/gas.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/status.hpp"
#include "goddard/thermo.hpp"

namespace Goddard {

/** Standard gravity [m/s^2], for specific impulse in seconds. */
constexpr double STANDARD_GRAVITY = 9.80665;

/**
 * Rocket performance of one exit, as CEA reports it.
 *
 * When a station it is computed from failed (`status.ok()` is false) every value is NaN.
 */
struct RocketPerformance {
    /** Stagnation pressure over exit pressure, P_0/Pe [-]. */
    double pressure_ratio = NaN;
    /** Exit area over throat area, Ae/At [-], from the mass fluxes at the throat and the exit. */
    double area_ratio = NaN;
    /** Exit Mach number [-]. */
    double mach_number = NaN;
    /** Characteristic velocity c* [m/s]. */
    double cstar = NaN;
    /** Thrust coefficient at optimum expansion [-]. */
    double CF = NaN;
    /** Specific impulse at optimum expansion, as an effective exhaust velocity [m/s]. */
    double isp = NaN;
    /** Vacuum specific impulse, as an effective exhaust velocity [m/s]. */
    double ivac = NaN;
    /**
     * Ok when the stagnation, throat and exit stations were all solved; otherwise the failure
     * of the first of them that was not.
     */
    SolveStatus status;

    /** Specific impulse at optimum expansion in seconds, `isp / STANDARD_GRAVITY` [s]. */
    double isp_s() const { return isp / STANDARD_GRAVITY; }
    /** Vacuum specific impulse in seconds, `ivac / STANDARD_GRAVITY` [s]. */
    double ivac_s() const { return ivac / STANDARD_GRAVITY; }
};

/**
 * Rocket performance from the stagnation, throat and exit states.
 *
 * c* = P_0 / (rho_t a_t), with a_t the throat speed of sound; the exit velocity follows from
 * the enthalpy drop h_0 - h_e, the area ratio from mass conservation, and CF = u_e / c*.
 *
 * @param stagnation Stagnation state the nozzle expands from (the chamber of an infinite-area
 * combustor).
 * @param throat Throat state.
 * @param exit Exit state.
 */
RocketPerformance rocket_performance(const ThermodynamicState& stagnation,
                                     const ThermodynamicState& throat,
                                     const ThermodynamicState& exit);

/**
 * One operating point of a rocket: the reactants and how they burn, and the nozzle they
 * expand through.
 *
 * A sweep is a list of points; `RocketSolver::solve` returns one `RocketResult` per point, in
 * the same order. Fields that every point of a sweep shares are copied into each point.
 */
struct RocketPoint {
    /**
     * Fuel stream: temperature [K], pressure [Pa] and mole fractions. The species are those of
     * the reactant file of the solver, or product species when the solver has no reactant file.
     */
    PhaseSpecification fuel;
    /** Oxidizer stream; see `fuel`. */
    PhaseSpecification oxidizer;
    /** Mixture ratio [-], interpreted by `combustor.mixture_type` (O/F by default). */
    double mixture_ratio = NaN;
    /**
     * Pressure [Pa]: the chamber pressure of an infinite-area chamber, the injector-face pressure
     * of a finite-area chamber, or the initial reactant pressure for `CombustionProcess::ISOCHORIC`
     * combustion.
     */
    double chamber_pressure = NaN;
    /** Chamber model, mixture ratio type and combustion process. */
    CombustorOptions combustor;
    /** Nozzle chemistry and freezing station. */
    NozzleOptions nozzle;
    /** How `expansion_ratios` are interpreted. */
    ExpansionType expansion_type = ExpansionType::SUPERSONIC_AREA_RATIO;
    /**
     * Exit stations, in the order they are solved: area ratios A/A_t [-], or pressure ratios
     * P_0/P [-] (stagnation over station pressure, as CEA's `pi/p`). May be empty.
     *
     * In Python, assign a whole list: item assignment and append act on a copy.
     */
    std::vector<double> expansion_ratios;
    /** Optional name of the point, copied to the result and printed in the report. */
    std::string label;
};

/**
 * Result of one `RocketPoint`.
 *
 * For an infinite-area chamber `chamber`, `stagnation` and `combustion_end` are the same state:
 * the gas is at rest in the whole chamber. For a finite-area chamber `chamber` is the injector
 * face, `combustion_end` the end of the constant-area section and `stagnation` the state "inf"
 * the nozzle expands from.
 *
 * Every station carries the stagnation enthalpy of the point in `thermo.stagnation_enthalpy`.
 * A station that was not solved has a failed status and NaN values; the stations downstream of
 * it are `NOT_COMPUTED`.
 */
struct RocketResult {
    /** The point that was solved. */
    RocketPoint point;
    /** Chamber station; the injector face of a finite-area chamber. */
    NozzleStation chamber;
    /** Stagnation state the nozzle expands from. */
    NozzleStation stagnation;
    /** End of the combustion chamber. */
    NozzleStation combustion_end;
    /** Sonic throat. */
    ThroatCondition throat;
    /** One station per `point.expansion_ratios` entry, in order. */
    std::vector<NozzleStation> exits;
    /** Performance at each exit, indexed like `exits`. */
    std::vector<RocketPerformance> performance;
    /**
     * Contraction ratio A_c/A_t [-]: solved for a `FINITE_MASS_FLUX` chamber, infinite for an
     * infinite-area one. NaN when the chamber was not solved.
     */
    double contraction_ratio = NaN;
    /**
     * Mass flux through the chamber, mdot/A_c [kg/(m^2 s)]: derived for a
     * `FINITE_CONTRACTION_RATIO` chamber, 0 for an infinite-area one. NaN when the chamber was
     * not solved.
     */
    double mass_flux = NaN;
    /** Ok when every station was solved; otherwise the status of the first that was not. */
    SolveStatus status;
};

/**
 * CEA-style rocket problem: combustion followed by an isentropic nozzle expansion, for a list of
 * operating points.
 *
 * The product `Gas` fixes the species of the problem: its gas-phase species and its candidate
 * condensed species are the possible combustion products.
 *
 * Reactants are given in one of two ways:
 * - with a reactant file, the fuel and oxidizer compositions name species of that file (e.g.
 *   `H2(L)`, `RP-1`), which need not be product species; element amounts and enthalpies are
 *   transferred to the products by element name.
 * - without one, the compositions name product species and are blended in the product phase.
 *
 * Every point is checked before any is solved, so a bad argument raises before the sweep
 * starts. A point whose solve fails is reported in its result and the sweep goes on.
 */
class RocketSolver {
public:
    /**
     * @param products Product gas. The solver works on its own copy.
     * @param reactant_file Cantera YAML file holding the reactant species, e.g.
     * `nasa9_reactants.yaml`; empty when the reactants are product species.
     */
    explicit RocketSolver(const Gas& products, std::string reactant_file = "");

    /**
     * Solve one point.
     *
     * @throws std::invalid_argument for an invalid point (see `validate`).
     * @throws NotImplementedError for an unimplemented option.
     */
    RocketResult solve(const RocketPoint& point);

    /**
     * Solve a list of points.
     *
     * @return One result per point, in the order of `points`.
     * @throws std::invalid_argument if any point is invalid, before anything is solved.
     * @throws NotImplementedError for an unimplemented option.
     */
    std::vector<RocketResult> solve(const std::vector<RocketPoint>& points);

    /**
     * Check a point without solving it.
     *
     * @throws std::invalid_argument for a reactant without composition or with a nonpositive
     * temperature or pressure, a nonpositive chamber pressure, a mixture ratio out of range, an
     * invalid finite-area chamber, or an exit station or freezing station the nozzle rejects.
     * @throws NotImplementedError for an unimplemented option, such as `PHI_RATIO` with a
     * reactant file or a finite-area chamber frozen at the combustion end.
     */
    void validate(const RocketPoint& point) const;

    /** Not implemented: `solve` raises `NotImplementedError` if true. */
    bool include_transport = false;
    /** Not implemented: `solve` raises `NotImplementedError` if true. */
    bool include_ionized_species = false;
    /** Not implemented: `solve` raises `NotImplementedError` if non-zero. */
    double trace_cutoff = 0.0;

private:
    Gas m_products;
    std::string m_reactant_file;
    /** Reactant phases of `m_reactant_file`, by species list; built on first use. */
    mutable std::map<std::vector<std::string>, Gas> m_reactant_phases;
    /** Combustor of the most recent point, and the reactant streams it was built for. */
    std::optional<Combustor> m_combustor;
    PhaseSpecification m_combustor_fuel;
    PhaseSpecification m_combustor_oxidizer;
    /** Nozzle of the most recent point, and the options it was built with. */
    std::optional<Nozzle> m_nozzle;
    NozzleOptions m_nozzle_options;

    /** Reject the solver-level options that are not implemented. */
    void check_implemented() const;
    /** Solve a point that `validate` accepted. */
    RocketResult solve_validated(const RocketPoint& point);
    /** Combustion state of a point: a one-entry array whose status reports a failed solve. */
    ThermoArray combust(const RocketPoint& point);
    /** Nozzle with the given options, reusing the one of the previous point when they match. */
    Nozzle& nozzle_for(const NozzleOptions& options);
    /**
     * Phase of `m_reactant_file` holding the species of `composition`.
     *
     * @throws std::invalid_argument if a species is not in the reactant file.
     */
    const Gas& reactant_phase(const Composition& composition) const;
    /** Reactant stream `Gas` of `m_reactant_file` at the state of `stream`. */
    Gas reactant_gas(const PhaseSpecification& stream) const;
};

/**
 * CEA-style text report of one result.
 *
 * A point whose chamber or throat was not solved is reported as a short notice with the
 * failure; one with failed exits keeps its table, with NaN in the failed columns and the
 * failures listed below it.
 */
std::string rocket_report(const RocketResult& result);

/** CEA-style text report of a list of results, one page per result, in order. */
std::string rocket_report(const std::vector<RocketResult>& results);

} // namespace Goddard
