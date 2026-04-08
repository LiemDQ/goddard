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
    unsigned int frozen_NFZ = 0;
    SolverOptions solver;
    double gamma; // used only for PERFECT_GAS

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
};

struct NozzleStation {
    bool converged;
    double gamma_s;
    double dlV_dlP_T;
    double dlV_dlT_P;
    std::vector<double> state;
};

struct NozzleResults {
    ThroatCondition throat;
    std::vector<NozzleStation> expansions;
};

class Nozzle {
    public:
    Nozzle(const Gas& gas, GasChemistry chemistry);
    Nozzle(const Gas& gas, GasChemistry chemistry, std::vector<double> state);

    NozzleResults solve(ExpansionType expansion_type, double ratio = 1.0);
    NozzleResults solve(ExpansionType expansion_type, const std::vector<double>& ratios);
    NozzleResults solve(const NozzleProfile& profile, int num_stations = 50);
    ThroatCondition solve_throat_conditions(double abstol = 4e-4);

    double get_gamma_s();

    void reset_state();
    inline void set_inlet_state(const std::vector<double>& state) { inlet_state = state; }
    inline std::vector<double> get_inlet_state() const { return inlet_state; }

    std::vector<double> inlet_state;

    private:
    Gas m_gas;
    NozzleOptions m_opts;

    void solve_chemistry();
    NozzleStation solve_supersonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol = 4.5e-5);
    NozzleStation solve_subsonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol = 4.5e-5);
    NozzleStation solve_pressure_ratio(const ThroatCondition& throat_condition, double pressure_ratio, double abstol = 0.5e-5);

    NozzleStation iterate_area_expansion(
        std::shared_ptr<Cantera::ThermoPhase>& gas_thermo,
        const ThroatCondition& throat_condition,
        double expansion_ratio, double pressure_ratio_guess, double abstol);

    double iterate_temperature(
        std::shared_ptr<Cantera::ThermoPhase>& gas_thermo,
        const ThroatCondition& throat_condition,
        double pressure_ratio, double T_guess,
        const std::vector<double>& composition,
        double abstol = 0.5e-5);
};

} //namespace Goddard
