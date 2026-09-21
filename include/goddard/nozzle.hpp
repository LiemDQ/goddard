#pragma once
#include "cantera/core.h"
#include <memory>
#include <vector>
#include <optional>


#include "goddard/thermoarray.hpp"
#include "goddard/chemistry.hpp"
#include "goddard/combustor.hpp"
#include "goddard/gas.hpp"
#include "goddard/profile.hpp"

namespace Goddard {

enum class ExpansionType {
    SUPERSONIC_AREA_RATIO,
    SUBSONIC_AREA_RATIO,
    PRESSURE_RATIO
};


struct NozzleOptions {
    GasChemistry chemistry = GasChemistry::EQUILIBRIUM;
    ExpansionType expansion_type = ExpansionType::SUPERSONIC_AREA_RATIO;
    std::vector<double> expansion_ratios;
    /**
     * Freezing station for FROZEN chemistry [-]: 0 is the chamber, 1 the throat, and 2 onward the
     * expansion stations in the order they are solved. The flow is in equilibrium up to and
     * including this station and keeps its composition downstream of it. With 0 the chamber
     * derivatives are frozen as well.
     */
    int frozen_NFZ = 0;
    SolverOptions solver;
    double gamma = 1.4; // used only for PERFECT_GAS

    // KineticNozzle-specific (ignored by Nozzle)
    double dt_max = 1e-6;
    double dx_max = 1e-3;
    int max_steps = 100000;
};


struct ThroatCondition {
    bool converged;
    double speed_of_sound;
    double H_stagnation;
    double P_inlet;
    double S_inlet;
    double gamma_s;
    double dlV_dlP_T;
    double dlV_dlT_P;
    std::vector<double> state;
    /**
     * True when the throat sits exactly at a condensed phase transition, so its temperature is
     * pinned at the transition temperature [K] and both polymorphs coexist. The equilibrium
     * specific heat is then infinite and `gamma_s` is -1 / `dlV_dlP_T`; `dlV_dlT_P` is infinite.
     * Always false for a mixture without condensed phases.
     */
    bool pinned_transition = false;
};

struct NozzleStation {
    bool converged;
    double gamma_s;
    double dlV_dlP_T;
    double dlV_dlT_P;
    std::vector<double> state;
    /** True when the station sits at a condensed phase transition; see `ThroatCondition`. */
    bool pinned_transition = false;
};

struct NozzleResults {
    /**
     * Nozzle inlet (chamber) station. Its derivatives are in the chemistry of station 0:
     * frozen only for FROZEN chemistry with `frozen_NFZ` == 0, otherwise equilibrium.
     */
    NozzleStation inlet;
    ThroatCondition throat;
    std::vector<NozzleStation> expansions;
};

/**
 * Converged chamber of a finite-area combustor.
 *
 * The chamber has constant area A_c, and propellant enters it with no axial momentum. The flow
 * from the combustion end onward is the isentropic expansion from a hypothetical stagnation
 * state "inf" with the injector enthalpy and pressure P_inf < P_inj.
 */
struct FiniteAreaChamber {
    /** Injector-face station; always in equilibrium. */
    NozzleStation injector;
    /** Stagnation station "inf": equilibrium at the injector enthalpy and P_inf. */
    NozzleStation stagnation;
    /** Combustion-end station, the subsonic point at area ratio `contraction_ratio`. */
    NozzleStation combustion_end;
    /** Throat solved from the stagnation state; `throat.P_inlet` is P_inf [Pa]. */
    ThroatCondition throat;
    /** Injector-face pressure P_inj [Pa]. */
    double injector_pressure;
    /** Stagnation pressure P_inf [Pa]. */
    double stagnation_pressure;
    /** Contraction ratio A_c/A_t [-]; solved for in mass-flux mode. */
    double contraction_ratio;
    /** Mass flux through the chamber, mdot/A_c [kg/(m^2 s)]; derived in contraction-ratio mode. */
    double mass_flux;
    /** Number of momentum-balance iterations [-]. */
    int iterations;
};

class Nozzle {
    public:
    Nozzle(const Gas& gas, NozzleOptions options = {});
    Nozzle(const Gas& gas, std::vector<double> state, NozzleOptions options = {});

    NozzleResults solve(ExpansionType expansion_type, double ratio = 1.0);
    NozzleResults solve(ExpansionType expansion_type, const std::vector<double>& ratios);
    NozzleResults solve(const NozzleProfile& profile, int num_stations = 50);
    ThroatCondition solve_throat_conditions(double abstol = 4e-4);

    /**
     * Solve expansion stations from an already solved throat.
     *
     * @param throat_condition Throat of the current inlet state.
     * @param expansion_type How `ratios` are interpreted.
     * @param ratios Area ratios A/A_t [-] or pressure ratios P_inlet/P [-].
     * @return One station per ratio, in order.
     */
    std::vector<NozzleStation> solve_stations(const ThroatCondition& throat_condition,
        ExpansionType expansion_type, const std::vector<double>& ratios);

    /**
     * Solve expansion stations downstream of a finite-area combustor.
     *
     * @param chamber Chamber returned by `solve_finite_area_chamber`.
     * @param expansion_type How `ratios` are interpreted.
     * @param ratios Area ratios A/A_t [-] or pressure ratios P_inj/P [-]. Pressure ratios are
     * taken from the injector pressure, as CEA does, not from the stagnation pressure.
     * @return One station per ratio, in order.
     */
    std::vector<NozzleStation> solve_stations(const FiniteAreaChamber& chamber,
        ExpansionType expansion_type, const std::vector<double>& ratios);

    /**
     * Solve the chamber of a finite-area combustor.
     *
     * Iterates on the stagnation pressure P_inf until the chamber momentum balance
     * P_inj = P_c + rho_c u_c^2 holds, where c is the combustion end. The chamber is always in
     * equilibrium. On return the nozzle inlet state is the stagnation state, so
     * `solve_stations(result.throat, ...)` continues the expansion, and the gas holds the
     * stagnation state.
     *
     * The combustion-end velocity follows from the enthalpy drop h_inj - h_c. At very large
     * contraction ratios (Mach numbers of order 1e-3 and below) that drop is within the station
     * solve's error, and a warning says the combustion-end velocity and Mach number are not
     * resolved. Pressures and the momentum balance are unaffected.
     *
     * @param injector_state Equilibrium state at the injector face (HP at P_inj).
     * @param type `FINITE_CONTRACTION_RATIO` or `FINITE_MASS_FLUX`.
     * @param value Contraction ratio A_c/A_t [-] or mass flux mdot/A_c [kg/(m^2 s)].
     * @param reltol Tolerance on |1 - P_inj,calc / P_inj| [-].
     * @return Converged chamber.
     *
     * @throws std::invalid_argument if `type` is not a finite-area type, if the contraction
     * ratio is not above 1, or if the mass flux thermally chokes the chamber: it exceeds the
     * perfect-gas limit by more than 2 %, or A_c/A_t reaches 1 in the real-gas iteration.
     * @throws NotImplementedError for FROZEN chemistry with `frozen_NFZ` == 0.
     * @throws ConvergenceError if the momentum balance does not converge.
     */
    FiniteAreaChamber solve_finite_area_chamber(const std::vector<double>& injector_state,
        CombustorType type, double value, double reltol = 1e-6);

    double get_gamma_s();

    void reset_state();
    inline void set_inlet_state(const std::vector<double>& state) { 
        inlet_state = state;
        m_gas.restore_state(inlet_state);
        m_gas.set_current_state_as_reference();
    }
    inline std::vector<double> get_inlet_state() const { return inlet_state; }

    std::vector<double> inlet_state;

    private:
    Gas m_gas;
    NozzleOptions m_opts;

    /*
     * Station solvers. `station` is the index of the station being solved (see
     * `NozzleOptions::frozen_NFZ`), and `frozen_state` the state whose composition the station
     * keeps if it is frozen. An equilibrium station ignores it and returns its own state, which
     * is then the frozen state of the stations after it.
     */
    NozzleStation solve_supersonic_area_expansion(const ThroatCondition& throat_condition,
        double expansion_ratio, int station, const std::vector<double>& frozen_state,
        double abstol = 4.5e-5);
    NozzleStation solve_subsonic_area_expansion(const ThroatCondition& throat_condition,
        double expansion_ratio, int station, const std::vector<double>& frozen_state,
        double abstol = 4.5e-5);
    NozzleStation solve_pressure_ratio(const ThroatCondition& throat_condition,
        double pressure_ratio, int station, const std::vector<double>& frozen_state,
        double abstol = 0.5e-5);

    NozzleStation iterate_area_expansion(
        const ThroatCondition& throat_condition,
        double expansion_ratio, double pressure_ratio_guess, int station,
        const std::vector<double>& frozen_state, double abstol);

    /**
     * Temperature [K] of a frozen station: the isentrope at fixed composition.
     *
     * With condensed products the amounts of the condensed species are frozen too, so the
     * temperature cannot leave the data range of any species that is present.
     *
     * @param throat_condition Throat state the expansion starts from.
     * @param pressure_ratio Chamber pressure divided by the station pressure [-].
     * @param T_guess Starting temperature [K].
     * @param composition Gas-phase mole fractions held fixed [-].
     * @param station Index of the station being solved, for error messages [-].
     * @param abstol Convergence tolerance on d(log T) [-].
     * @return Station temperature [K].
     *
     * @throws FmtError if a present condensed species leaves its temperature range, as CEA
     * reports for a frozen expansion carried too far.
     */
    double iterate_temperature(
        const ThroatCondition& throat_condition,
        double pressure_ratio, double T_guess,
        const std::vector<double>& composition,
        int station,
        double abstol = 0.5e-5);

    /** True when station `station` is in equilibrium, from `NozzleOptions::frozen_NFZ`. */
    bool is_equilibrium_station(int station) const;

    /**
     * Set the gas chemistry for station `station`.
     *
     * @return True for an equilibrium station.
     */
    bool set_station_chemistry(int station);

    /**
     * Station at `state` with derivatives in the chemistry of station `station`. Leaves the gas
     * at `state`.
     */
    NozzleStation station_at_state(const std::vector<double>& state, int station);

    void throw_invalid_expansion_ratio(double expansion_ratio, double min_ratio) const;
};

} //namespace Goddard
