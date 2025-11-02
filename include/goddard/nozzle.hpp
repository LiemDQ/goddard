#pragma once 
#include "cantera/core.h"
#include "eigen3/Eigen/Dense"
#include <memory>
#include <vector>


#include "goddard/thermoarray.hpp"
#include "goddard/case_options.hpp"

namespace Goddard {

// struct NozzleOptions {
//     ReactionType reaction_type;
//     ExpansionType expansion_type;
//     double reltol = 1e-3;
//     double abstol = 1e-6;
// };

struct ThroatCondition {
    bool converged;
    double H_stagnation;
    double P_inlet;
    double S_inlet;
    double gamma_s;
    double dlV_dlP_T;
    double dlV_dlT_P;
    std::vector<double> state;
};

struct NozzleResult {
    bool converged;
    double gamma_s;
    double dlV_dlP_T;
    double dlV_dlT_P;
    std::vector<double> state;
};

struct NozzleResults {
    ThroatCondition throat;
    std::vector<NozzleResult> expansions;
};
    
class NozzleBase {
    public:
    NozzleBase(Cantera::Solution& gas): m_gas(gas.shared_from_this()), m_inlet_state(gas.thermo()->stateSize()) {
        m_gas->thermo()->saveState(m_inlet_state);
    }

    NozzleBase(Cantera::Solution& gas, std::vector<double> state): m_gas(gas.shared_from_this()), m_inlet_state(state) {}

    NozzleResult solve(ExpansionType expansion_type, double ratio = 1.0);
    NozzleResults solve(ExpansionType expansion_type, std::vector<double> ratios);

    void reset_state();
    inline void set_inlet_state(const std::vector<double>& state) {m_inlet_state = state;}
    inline std::vector<double> get_inlet_state() const {return m_inlet_state;}

    /**
     * @brief Calculate chemical reaction-adjusted specific heat ratio from a thermodynamic state. 
     * @warning Can modify the underlying state!
     */
    virtual double get_gamma_s(Cantera::ThermoPhase& state) = 0;

    protected:
    std::shared_ptr<Cantera::Solution> m_gas;
    std::vector<double> m_inlet_state;
    ThroatCondition solve_throat_conditions(double abstol = 4e-4);
    virtual NozzleResult solve_supersonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol=4.5e-5) = 0;
    virtual NozzleResult solve_subsonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol=4.5e-5) = 0;
    virtual NozzleResult solve_pressure_ratio(const ThroatCondition& throat_condition, double pressure_ratio, double abstol=0.5e-5) = 0;
};

class EquilibriumNozzle : public NozzleBase {
    public:
    using NozzleBase::NozzleBase;
    double get_gamma_s(Cantera::ThermoPhase& state) override;

    protected:
    NozzleResult solve_supersonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol=4.5e-5) override;
    NozzleResult solve_subsonic_area_expansion(const ThroatCondition& throat_condition, double pressure_ratio, double abstol=4.5e-5) override;
    NozzleResult solve_pressure_ratio(const ThroatCondition& throat_condition, double pressure_ratio, double abstol=0.5e-5) override;

    NozzleResult iterate_area_expansion(
        std::shared_ptr<Cantera::ThermoPhase>& gas_thermo, 
        const ThroatCondition& throat_condition, 
        double expansion_ratio, double pressure_ratio_guess, double abstol);
};


class FrozenNozzle final : public NozzleBase {
    public:
    using NozzleBase::NozzleBase;
    double get_gamma_s(Cantera::ThermoPhase& state) override;

    protected:
    NozzleResult solve_supersonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol) override;
    NozzleResult solve_subsonic_area_expansion(const ThroatCondition& throat_condition, double pressure_ratio, double abstol) override;
    NozzleResult solve_pressure_ratio(const ThroatCondition& throat_condition, double pressure_ratio, double abstol) override;

    NozzleResult iterate_area_expansion(
        std::shared_ptr<Cantera::ThermoPhase>& gas_thermo, 
        const ThroatCondition& throat_condition, 
        double expansion_ratio, double pressure_ratio_guess, double abstol);

    double iterate_temperature(
        std::shared_ptr<Cantera::ThermoPhase>& gas_thermo,
        const ThroatCondition& throat_condition, 
        double pressure_ratio, double T_guess, double abstol=0.5e-5);
    
};

} //namespace Goddard