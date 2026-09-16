#pragma once
#include "cantera/core.h"
#include <memory>
#include <vector>
#include <optional>


#include "goddard/thermoarray.hpp"
#include "goddard/chemistry.hpp"
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
    ThroatCondition throat;
    std::vector<NozzleStation> expansions;
};

class Nozzle {
    public:
    Nozzle(const Gas& gas, NozzleOptions options = {});
    Nozzle(const Gas& gas, std::vector<double> state, NozzleOptions options = {});

    NozzleResults solve(ExpansionType expansion_type, double ratio = 1.0);
    NozzleResults solve(ExpansionType expansion_type, const std::vector<double>& ratios);
    NozzleResults solve(const NozzleProfile& profile, int num_stations = 50);
    ThroatCondition solve_throat_conditions(double abstol = 4e-4);

    double get_gamma_s();

    void reset_state();
    inline void set_inlet_state(const std::vector<double>& state) { 
        inlet_state = state;
        m_gas.restore_state(inlet_state);
        m_gas.set_current_state_as_reference();
        m_current_station = 0;
    }
    inline std::vector<double> get_inlet_state() const { return inlet_state; }

    std::vector<double> inlet_state;

    private:
    Gas m_gas;
    NozzleOptions m_opts;
    int m_current_station = 0;

    NozzleStation solve_supersonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol = 4.5e-5);
    NozzleStation solve_subsonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol = 4.5e-5);
    NozzleStation solve_pressure_ratio(const ThroatCondition& throat_condition, double pressure_ratio, double abstol = 0.5e-5);

    NozzleStation iterate_area_expansion(
        const ThroatCondition& throat_condition,
        double expansion_ratio, double pressure_ratio_guess, double abstol);

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
        double abstol = 0.5e-5);

    bool determine_equilibrium_condition();

    void throw_invalid_expansion_ratio(double expansion_ratio, double min_ratio) const;
};

} //namespace Goddard
