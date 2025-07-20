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

struct ThroatResult {
    bool converged;
    double H_stagnation;
    double P_inlet;
    double S_inlet;
};

struct ThroatResults {
    bool converged;
    Eigen::ArrayXd H_stagnation;
    Eigen::ArrayXd P_inlet;
    Eigen::ArrayXd S_inlet;
};

struct NozzleResult {
    bool converged;
    std::shared_ptr<Cantera::Solution> exit_gas;
};



/**
 * Get throat conditions.
 * 
 * @param frozen Whether an equilibrium or frozen gas assumption is used. 
 */
ThroatResults solve_throat_conditions(
    std::shared_ptr<Cantera::Solution> inlet_gas, 
    bool frozen = true, 
    double abstol = 4e-4
);
    
NozzleResult solve_equilibrium_supersonic_area_expansion(
    Cantera::Solution& throat_gas, 
    const ThroatResult& throat_result, 
    double expansion_ratio,
    double abstol = 4e-5
);

NozzleResult solve_equilibrium_pressure_ratio(
    Cantera::Solution& throat_gas,
    const ThroatResult& throat_result,
    double pressure_ratio,
    double abstol = 0.5e-5
);


    
class NozzleBase {
    public:
    NozzleBase(Cantera::Solution& gas): m_gas(gas.shared_from_this()), m_initial_state(gas.thermo()->stateSize()) {
        m_gas->thermo()->saveState(m_initial_state);
    }

    NozzleResult solve(const NozzleConditions& conditions);
    virtual NozzleResult solve_supersonic_area_expansion(double expansion_ratio, double abstol=4e-5) = 0;
    virtual NozzleResult solve_subsonic_area_expansion() = 0;
    virtual NozzleResult solve_pressure_ratio(double pressure_ratio, double abstol=0.5e-5) = 0;
    virtual ThroatResult solve_throat_conditions(double abstol = 4e-4) = 0;
    virtual double gamma() = 0;

    protected:
    std::shared_ptr<Cantera::Solution> m_gas;
    std::vector<double> m_initial_state;

};

class EquilibriumNozzle : NozzleBase {
    public:    
    NozzleResult solve_supersonic_area_expansion(double expansion_ratio, double abstol=4e-5) override;
    NozzleResult solve_subsonic_area_expansion() override;
    NozzleResult solve_pressure_ratio(double pressure_ratio, double abstol=0.5e-5) override;
    ThroatResult solve_throat_conditions(double abstol = 4e-4) override;
    protected:
    double gamma() override;

};

} //namespace Goddard