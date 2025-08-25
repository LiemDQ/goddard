#pragma once 
#include "cantera/core.h"
#include "eigen3/Eigen/Dense"
#include <memory>

#include "goddard/thermoarray.hpp"

namespace Goddard {

enum class ReactionType { //NOTE: should function call just depend on the nozzle type or via inheritance?
    FROZEN,
    EQUILIBRIUM,
    // KINETIC,
};

enum class ExpansionType {
    SUPERSONIC_AREA_RATIO,
    SUBSONIC_AREA_RATIO,
    PRESSURE_RATIO
};

struct NozzleConditions {
    Eigen::ArrayXd expansion_ratios;
    ExpansionType expansion_type;
    ReactionType reaction_type;
};

struct NozzleOptions {
    ReactionType reaction_type;
    ExpansionType expansion_type;
    double reltol = 1e-3;
    double abstol = 1e-6;
};

struct ThroatCondition {
    bool converged;
    double H_stagnation;
    double P_inlet;
    double S_inlet;
    std::vector<double> state;
};

struct ThroatResults {
    bool converged;
    Eigen::ArrayXd H_stagnation;
    Eigen::ArrayXd P_inlet;
    Eigen::ArrayXd S_inlet;
};

struct NozzleResult {
    bool converged;
    std::vector<double> state;
};

    
class NozzleBase {
    public:
    NozzleBase(Cantera::Solution& gas): m_gas(gas.shared_from_this()), m_initial_state(gas.thermo()->stateSize()) {
        m_gas->thermo()->saveState(m_initial_state);
    }

    NozzleResult solve(ExpansionType expansion_type = ExpansionType::PRESSURE_RATIO, double expansion_ratio = 0, double pressure_ratio = 0);
    void reset_state();

    /**
     * @warning Can modify the underlying state!
     */
    virtual double gamma(Cantera::ThermoPhase& state) = 0;

    protected:
    std::shared_ptr<Cantera::Solution> m_gas;
    std::vector<double> m_initial_state;
    ThroatCondition solve_throat_conditions(double abstol = 4e-4);
    virtual NozzleResult solve_supersonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol=4.5e-5) = 0;
    virtual NozzleResult solve_subsonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol=4.5e-5) = 0;
    virtual NozzleResult solve_pressure_ratio(const ThroatCondition& throat_condition, double pressure_ratio, double abstol=0.5e-5) = 0;
};

class EquilibriumNozzle : NozzleBase {
    public:    
    double gamma(Cantera::ThermoPhase& state) override;

    protected:
    NozzleResult solve_supersonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol=4.5e-5) override;
    NozzleResult solve_subsonic_area_expansion(const ThroatCondition& throat_condition, double pressure_ratio, double abstol=4.5e-5) override;
    NozzleResult solve_pressure_ratio(const ThroatCondition& throat_condition, double pressure_ratio, double abstol=0.5e-5) override;

    NozzleResult iterate_area_expansion(
        std::shared_ptr<Cantera::ThermoPhase>& gas_thermo, 
        const ThroatCondition& throat_condition, 
        double expansion_ratio, double pressure_ratio_guess, double abstol);
};


class FrozenNozzle : NozzleBase {
    public:    
    double gamma(Cantera::ThermoPhase& state) override;

    protected:
    NozzleResult solve_supersonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol) override;
    NozzleResult solve_subsonic_area_expansion(const ThroatCondition& throat_condition, double pressure_ratio, double abstol) override;
    NozzleResult solve_pressure_ratio(const ThroatCondition& throat_condition, double pressure_ratio, double abstol) override;

};

} //namespace Goddard