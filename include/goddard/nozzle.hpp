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

/** How the ratios of nozzle stations are interpreted. */
enum class ExpansionType {
    /// Area ratio A/A_t [-] downstream of the throat (supersonic branch).
    SUPERSONIC_AREA_RATIO,
    /// Area ratio A/A_t [-] upstream of the throat (subsonic branch).
    SUBSONIC_AREA_RATIO,
    /// Pressure ratio P_inlet/P [-], chamber over station pressure.
    PRESSURE_RATIO
};


/** Settings of a nozzle expansion. */
struct NozzleOptions {
    /**
     * Chemistry of the expansion: EQUILIBRIUM, or FROZEN downstream of station `frozen_NFZ`.
     * `Nozzle` rejects KINETIC (use `KineticNozzle`) and does not implement PERFECT_GAS.
     */
    GasChemistry chemistry = GasChemistry::EQUILIBRIUM;
    /// How `expansion_ratios` are interpreted. Read by `RocketProblem`.
    ExpansionType expansion_type = ExpansionType::SUPERSONIC_AREA_RATIO;
    /**
     * Station ratios, interpreted by `expansion_type`: area ratios A/A_t [-], or pressure ratios
     * P_inlet/P [-] (chamber over station pressure, as CEA's `pi/p`). Read by `RocketProblem`;
     * the `Nozzle::solve` overloads take their ratios as an argument.
     *
     * In Python, assign a whole list: item assignment and append act on a copy.
     */
    std::vector<double> expansion_ratios;
    /**
     * Freezing station for FROZEN chemistry [-], 0-based: 0 is the chamber (the default, and
     * CEA's default freezing point), 1 the throat, and 2 onward the expansion stations in the
     * order they are solved. The flow is in equilibrium up to and including this station and
     * keeps its composition downstream of it. With 0 the chamber derivatives are frozen as well.
     *
     * For an infinite-area combustor CEA's `nfz` is `frozen_NFZ + 1`. For a finite-area
     * combustor the count starts at the combustion end, so `nfz` is `frozen_NFZ + 3`; 1 is the
     * throat in both cases.
     */
    int frozen_NFZ = 0;
};


/**
 * Sonic throat of an isentropic expansion from the nozzle inlet state. Thermodynamic
 * derivatives are in the chemistry of the throat station (see `NozzleOptions::frozen_NFZ`).
 */
struct ThroatCondition {
    /// Speed of sound at the throat [m/s].
    double speed_of_sound;
    /// Stagnation enthalpy [J/kg]: the specific enthalpy of the inlet state.
    double H_stagnation;
    /// Inlet (chamber or stagnation) pressure [Pa].
    double P_inlet;
    /// Inlet specific entropy [J/(kg.K)], held constant through the expansion.
    double S_inlet;
    /// Isentropic exponent -(d log P / d log V)_s at the throat [-].
    double gamma_s;
    /// (d log V / d log P)_T at the throat [-].
    double dlV_dlP_T;
    /// (d log V / d log T)_P at the throat [-].
    double dlV_dlT_P;
    /// Raw state vector of the throat, for `Gas::restore_state`.
    std::vector<double> state;
    /**
     * True when the throat sits exactly at a condensed phase transition, so its temperature is
     * pinned at the transition temperature [K] and both polymorphs coexist. The equilibrium
     * specific heat is then infinite and `gamma_s` is -1 / `dlV_dlP_T`; `dlV_dlT_P` is infinite.
     * Always false for a mixture without condensed phases.
     */
    bool pinned_transition = false;
    /** Mixture state at the throat; see `NozzleStation::thermo`. */
    ThermodynamicState thermo;
    /** Flow velocity at the throat [m/s]. */
    double velocity = 0.0;
    /** Mach number at the throat [-]: 1 to within the throat solve's tolerance. */
    double mach = 0.0;
    /** Area ratio A/A_t [-]: 1 by definition. */
    double area_ratio = 1.0;
};

/**
 * One solved nozzle station.
 */
struct NozzleStation {
    /**
     * Mixture state at the station. `gamma_s`, `dlV_dlP_T`, `dlV_dlT_P`, `speed_of_sound` and
     * `pinned_transition` are in the chemistry of the station: frozen downstream of the freezing
     * station, equilibrium otherwise. `pinned_transition` is true when the station sits exactly at
     * a condensed phase transition (see `ThroatCondition::pinned_transition`).
     */
    ThermodynamicState thermo;
    /** Flow velocity [m/s], from the enthalpy drop below the stagnation enthalpy. */
    double velocity = 0.0;
    /** Mach number [-], from `velocity` and `thermo.speed_of_sound`. */
    double mach = 0.0;
    /**
     * Area ratio A/A_t [-], from the mass flux relative to the throat. 0 for a station where the
     * flow is at rest (the chamber of an infinite-area combustor, the injector face and the
     * stagnation state of a finite-area one).
     */
    double area_ratio = 0.0;
    /** Raw state vector, for `Gas::restore_state`. */
    std::vector<double> state;
};

/** Result of a `Nozzle::solve`: inlet, throat and the requested stations. */
struct NozzleResults {
    /**
     * Nozzle inlet (chamber) station. Its derivatives are in the chemistry of station 0:
     * frozen only for FROZEN chemistry with `frozen_NFZ` == 0, otherwise equilibrium.
     */
    NozzleStation inlet;
    /// Sonic throat.
    ThroatCondition throat;
    /// One station per requested ratio or profile point, in order.
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

/**
 * One-dimensional isentropic nozzle expansion with equilibrium or frozen chemistry.
 *
 * Solves the sonic throat from the inlet (chamber) state, then stations given by area ratio,
 * pressure ratio or a wall profile.
 */
class Nozzle {
    public:
    /**
     * @param gas gas at the inlet (chamber) state. The solver works on its own copy of `gas`;
     * the caller's `Gas` is not modified.
     * @param options chemistry, freezing station and stations to solve.
     * @throws std::invalid_argument for KINETIC chemistry.
     * @throws NotImplementedError for PERFECT_GAS chemistry.
     */
    Nozzle(const Gas& gas, NozzleOptions options = {});
    /**
     * @param gas product gas; only its phase definition is used. The solver works on its own copy.
     * @param state inlet (chamber) state vector, as returned by `Gas::save_state`.
     * @param options chemistry, freezing station and stations to solve.
     * @throws std::invalid_argument for KINETIC chemistry.
     * @throws NotImplementedError for PERFECT_GAS chemistry.
     */
    Nozzle(const Gas& gas, std::vector<double> state, NozzleOptions options = {});

    /**
     * Solve the throat and one station.
     *
     * @param expansion_type How `ratio` is interpreted.
     * @param ratio Area ratio A/A_t [-] or pressure ratio P_inlet/P [-].
     */
    NozzleResults solve(ExpansionType expansion_type, double ratio = 1.0);
    /**
     * Solve the throat and one station per ratio.
     *
     * @param expansion_type How `ratios` are interpreted.
     * @param ratios Area ratios A/A_t [-] or pressure ratios P_inlet/P [-].
     */
    NozzleResults solve(ExpansionType expansion_type, const std::vector<double>& ratios);
    /**
     * Solve the throat and supersonic stations along a wall profile.
     *
     * The first point of the profile (`x_min`) is taken as the throat, and the stations are
     * spaced evenly in x up to the last point. The area is that of an axisymmetric nozzle.
     *
     * @param profile Diverging wall contour starting at the throat.
     * @param num_stations Number of stations [-].
     */
    NozzleResults solve(const NozzleProfile& profile, int num_stations = 50);
    /**
     * Solve the sonic throat by isentropic expansion from the inlet state.
     *
     * @param abstol Tolerance on |1 - 1/M^2| at the throat [-].
     * @throws ConvergenceError if the throat iteration does not converge in 5 iterations.
     */
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

    /** Isentropic exponent of the solver's gas at its current state [-]. */
    double get_gamma_s();

    /** Restore the solver's gas to the inlet state. */
    void reset_state();
    /** Set the inlet (chamber) state the next solve expands from, and make it the reference. */
    inline void set_inlet_state(const std::vector<double>& state) {
        m_inlet_state = state;
        m_gas.restore_state(m_inlet_state);
        m_gas.set_current_state_as_reference();
    }
    /** Inlet (chamber) state vector. */
    inline std::vector<double> get_inlet_state() const { return m_inlet_state; }

    private:
    std::vector<double> m_inlet_state;
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
     * @throws std::runtime_error (FmtError) if a present condensed species leaves its temperature range, as CEA
     * reports for a frozen expansion carried too far.
     */
    double iterate_temperature(
        const ThroatCondition& throat_condition,
        double pressure_ratio, double T_guess,
        const std::vector<double>& composition,
        int station,
        double abstol = 0.5e-5);

    /**
     * Station at the current gas state, with derivatives `props` already solved for it.
     *
     * @param H_stagnation Stagnation enthalpy the velocity is measured from [J/kg].
     * @param area_per_mdot_throat Throat area per mass flow rate [m^2 s/kg]; 0 when there is no
     * throat yet, which leaves `area_ratio` at 0.
     */
    NozzleStation station_at_current_state(const ExpansionProperties& props, double H_stagnation,
        double area_per_mdot_throat);

    /** Throat area per mass flow rate [m^2 s/kg], from the throat state and velocity. */
    static double throat_area_per_mdot(const ThroatCondition& throat);

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
