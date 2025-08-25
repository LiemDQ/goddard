#include "goddard/nozzle.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/error.hpp"
#include "goddard/utils.hpp"
#include "goddard/nozzle_utils.hpp"

#include "cantera/core.h"
#include <cmath>
#include <exception>
#include <vector>

namespace Goddard {

NozzleResult NozzleBase::solve(ExpansionType expansion_type, double expansion_ratio, double pressure_ratio) {
    
    ThroatCondition throat_condition = solve_throat_conditions();
    switch (expansion_type) {
        case ExpansionType::SUPERSONIC_AREA_RATIO:
            return this->solve_supersonic_area_expansion(throat_condition, expansion_ratio);
        case ExpansionType::SUBSONIC_AREA_RATIO:
            return this->solve_subsonic_area_expansion(throat_condition, expansion_ratio);
        default:
            return this->solve_pressure_ratio(throat_condition, expansion_ratio);
    }
}

void NozzleBase::reset_state(){
    m_gas->thermo()->restoreState(m_initial_state);
}

ThroatCondition NozzleBase::solve_throat_conditions(double abstol) {
    auto gas_state = m_gas->thermo();

    double P_inlet = gas_state->pressure();
    double S_inlet = gas_state->entropy_mass();
    double H_inlet = gas_state->enthalpy_mass();
    
    std::vector<double> X_inlet(gas_state->nSpecies());
    gas_state->getMoleFractions(X_inlet.data());

    EquilibriumProperties thermo_props = get_thermo_equilibrium_properties(*m_gas);

    double gamma_s = thermo_props.gamma_s;
    double P_throat = P_inlet / std::pow((gamma_s+1)/2,gamma_s/(gamma_s-1));
    
    const int max_iters = 5;
    int iter = 0;
    double Mach = 1.0; //throat mach number is 1 by definition
    double residual = 1.0;

    while (residual > abstol){
        if (iter >= max_iters){
            //show error message or throw exception
            return {};
        }
        P_throat = P_throat * (1 + gamma_s * Mach* Mach)/(1+ gamma_s);
        
        gas_state->setState_SP(S_inlet,P_throat);

        //if equilibrium conditions are selected, the composition must reach chemical
        //equilibrium in the throat.
        gamma_s = get_gamma_s(*gas_state);
        
        double velocity = gas_isenthalpic_velocity(*gas_state, H_inlet);
        double sonic_velocity = gas_sonic_velocity(*gas_state, gamma_s);

        Mach = velocity/sonic_velocity;
        residual = std::abs(1.0 - 1.0/(Mach * Mach));

        iter++;
    }

    return {true, H_inlet, P_inlet, S_inlet, save_gas_state(*gas_state)};
}

double EquilibriumNozzle::get_gamma_s(Cantera::ThermoPhase& state) {
    state.equilibrate("SP", "gibbs");
     //WARNING: if m_gas has a different ThermoPhase than `state` this will result in incorrect behavior!
    auto props = get_thermo_equilibrium_properties(*m_gas);
    return props.gamma_s;
}

NozzleResult EquilibriumNozzle::solve_subsonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol) {
    std::shared_ptr<Cantera::ThermoPhase> gas_thermo = m_gas->thermo();
    gas_thermo->restoreState(throat_condition.state);

    double ln_pressure_ratio = 0;
    double throat_pressure_ratio = throat_condition.P_inlet/gas_thermo->pressure();
    double ln_throat_ratio = std::log(throat_pressure_ratio);
    double ln_Ae_At = std::log(expansion_ratio);
    
    if (expansion_ratio >= 1.09) {
        ln_pressure_ratio = ln_throat_ratio/(expansion_ratio + 10.587*std::pow(ln_Ae_At, 3)+9.454*ln_Ae_At);
    }
    else if (expansion_ratio > 1.0001) {
        ln_pressure_ratio = 0.9* ln_throat_ratio / (expansion_ratio + 10.587*std::pow(ln_Ae_At, 3)+9.454*ln_Ae_At);
    }
    else {
        //invalid expansion ratio
        return {false, {}};
    }
    return iterate_area_expansion(gas_thermo, throat_condition, expansion_ratio, std::exp(ln_pressure_ratio), abstol);
}

NozzleResult EquilibriumNozzle::solve_supersonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol) {
    std::shared_ptr<Cantera::ThermoPhase> gas_thermo = m_gas->thermo();
    gas_thermo->restoreState(throat_condition.state);
    
    EquilibriumProperties equilibrium_props = get_thermo_equilibrium_properties(*m_gas);
    double gamma_s = equilibrium_props.gamma_s;

    double ln_pressure_ratio = 0;
    if (expansion_ratio >= 2) {
        ln_pressure_ratio = gamma_s+1.4*std::log(expansion_ratio);
    } 
    else if (expansion_ratio > 1.0001) {
        double throat_pressure_ratio = throat_condition.P_inlet/gas_thermo->pressure();
        double ln_throat_ratio = std::log(throat_pressure_ratio);
        double ln_Ae_At = std::log(expansion_ratio);
        ln_pressure_ratio = ln_throat_ratio + std::sqrt(3.294*ln_Ae_At*ln_Ae_At+1.535*ln_Ae_At);
    } else {
        //invalid expansion ratio
        return {false, {}};
    }
    
    return iterate_area_expansion(gas_thermo, throat_condition, expansion_ratio, std::exp(ln_pressure_ratio), abstol);
}

NozzleResult EquilibriumNozzle::iterate_area_expansion(
    std::shared_ptr<Cantera::ThermoPhase>& gas_thermo,
    const ThroatCondition& throat_condition, 
    double expansion_ratio, double pressure_ratio_guess, double abstol) {

    double pressure_ratio = pressure_ratio_guess;
    double P_exit = throat_condition.P_inlet/pressure_ratio;
    double velocity = gas_isenthalpic_velocity(*gas_thermo, throat_condition.H_stagnation);
    const double A_mdot_thrt = area_per_mdot(*gas_thermo, velocity);

    double gamma_s = get_gamma_s(*gas_thermo);
    double Ae_At = gas_thermo->temperature() / (gas_thermo->pressure() * velocity * gas_thermo->meanMolecularWeight())/A_mdot_thrt;

    int iters = 0;
    int max_iter = 10;
    double residual = 1.0;
    double sonic_velocity = 0.0;
    double dlogp_dlogA = 0.0;

    while (std::abs(residual) > abstol) {
        iters++;
        if (iters >= max_iter) {
            //show error message or throw exception
            return {false, {}};
        }
        velocity = gas_isenthalpic_velocity(*gas_thermo, throat_condition.H_stagnation);
        sonic_velocity = gas_sonic_velocity(*gas_thermo, gamma_s);
        Ae_At = area_per_mdot(*gas_thermo, velocity)/A_mdot_thrt;

        dlogp_dlogA = gamma_s * velocity * velocity / (velocity*velocity - sonic_velocity*sonic_velocity);
        residual = dlogp_dlogA * (std::log(expansion_ratio) - std::log(Ae_At));
        double log_pinf_pe = std::log(pressure_ratio) + residual;

        pressure_ratio = std::exp(log_pinf_pe);
        P_exit = throat_condition.P_inlet / pressure_ratio;

        gas_thermo->setState_SP(throat_condition.S_inlet, P_exit);
        gamma_s = get_gamma_s(*gas_thermo);
    }
    return {true, save_gas_state(*gas_thermo)};
}


NozzleResult EquilibriumNozzle::solve_pressure_ratio(const ThroatCondition& throat_condition, double pressure_ratio, double abstol) {
    std::shared_ptr<Cantera::ThermoPhase> gas_thermo = m_gas->thermo();
    gas_thermo->restoreState(throat_condition.state);

    double P_exit = throat_condition.P_inlet/pressure_ratio;
    gas_thermo->setState_SP(throat_condition.S_inlet, P_exit);
    gas_thermo->equilibrate("SP", "gibbs");
    
    //pressure ratio for equilibrium nozzle does not require iteration
    return {true, save_gas_state(*gas_thermo)};
}

double FrozenNozzle::get_gamma_s(Cantera::ThermoPhase& state) {
    //gamma_s = gamma. See CEA Part I Section 6.5.3.
    return state.cp_mass()/state.cv_mass();
}

NozzleResult FrozenNozzle::solve_supersonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol) {
    std::shared_ptr<Cantera::ThermoPhase> gas_thermo = m_gas->thermo();
    gas_thermo->restoreState(throat_condition.state);
    
    double gamma_s = get_gamma_s(*gas_thermo);

    double ln_pressure_ratio = 0;
    if (expansion_ratio >= 2) {
        ln_pressure_ratio = gamma_s+1.4*std::log(expansion_ratio);
    } 
    else if (expansion_ratio > 1.0001) {
        double throat_pressure_ratio = throat_condition.P_inlet/gas_thermo->pressure();
        double ln_throat_ratio = std::log(throat_pressure_ratio);
        double ln_Ae_At = std::log(expansion_ratio);
        ln_pressure_ratio = ln_throat_ratio + std::sqrt(3.294*ln_Ae_At*ln_Ae_At+1.535*ln_Ae_At);
    } else {
        //invalid expansion ratio
        return {false, {}};
    }

    return iterate_area_expansion(gas_thermo, throat_condition, expansion_ratio, std::exp(ln_pressure_ratio), abstol);
}

NozzleResult FrozenNozzle::solve_subsonic_area_expansion(const ThroatCondition& throat_condition, double expansion_ratio, double abstol) {
    std::shared_ptr<Cantera::ThermoPhase> gas_thermo = m_gas->thermo();
    gas_thermo->restoreState(throat_condition.state);

    double ln_pressure_ratio = 0;
    double throat_pressure_ratio = throat_condition.P_inlet/gas_thermo->pressure();
    double ln_throat_ratio = std::log(throat_pressure_ratio);
    double ln_Ae_At = std::log(expansion_ratio);
    
    if (expansion_ratio >= 1.09) {
        ln_pressure_ratio = ln_throat_ratio/(expansion_ratio + 10.587*std::pow(ln_Ae_At, 3)+9.454*ln_Ae_At);
    }
    else if (expansion_ratio > 1.0001) {
        ln_pressure_ratio = 0.9* ln_throat_ratio / (expansion_ratio + 10.587*std::pow(ln_Ae_At, 3)+9.454*ln_Ae_At);
    }
    else {
        //invalid expansion ratio
        return {false, {}};
    }
    return iterate_area_expansion(gas_thermo, throat_condition, expansion_ratio, std::exp(ln_pressure_ratio), abstol);
}

NozzleResult FrozenNozzle::solve_pressure_ratio(const ThroatCondition& throat_condition, double pressure_ratio, double abstol = 0.5e-5) {
    std::shared_ptr<Cantera::ThermoPhase> gas_thermo = m_gas->thermo();
    gas_thermo->restoreState(throat_condition.state);
    
    //To solve the gas state, we need to iterate to find the exit temperature
    //initial guess
    double P_exit = iterate_temperature(gas_thermo, throat_condition, pressure_ratio, gas_thermo->temperature(), abstol);
    if (P_exit < 0) {
        //iteration failed to find a solution
        return {false, {}};
    }
    else {
        return {true, save_gas_state(*gas_thermo)};
    }
}

NozzleResult FrozenNozzle::iterate_area_expansion(
    std::shared_ptr<Cantera::ThermoPhase>& gas_thermo, 
    const ThroatCondition& throat_condition, 
    double expansion_ratio, double pressure_ratio_guess, double abstol) {

    
    double pressure_ratio = pressure_ratio_guess;
    double P_exit = throat_condition.P_inlet/pressure_ratio;
    double velocity = gas_isenthalpic_velocity(*gas_thermo, throat_condition.H_stagnation);
    const double A_mdot_thrt = area_per_mdot(*gas_thermo, velocity);
    
    double Ae_At = gas_thermo->temperature() / (gas_thermo->pressure() * velocity * gas_thermo->meanMolecularWeight())/A_mdot_thrt;
    double T_exit = gas_thermo->temperature();

    int iters = 0;
    int max_iter = 10;
    double residual = 1.0;
    double sonic_velocity = 0.0;
    double dlogp_dlogA = 0.0;
    
    while (std::abs(residual) > abstol) {
        iters++;
        if (iters >= max_iter) {
            //show error message or throw exception
            return {false, {}};
        }
        T_exit = iterate_temperature(gas_thermo, throat_condition, pressure_ratio, T_exit);
        if (T_exit < 0) {
            return {false, {}};
        }

        double gamma_s = get_gamma_s(*gas_thermo);
        velocity = gas_isenthalpic_velocity(*gas_thermo, throat_condition.H_stagnation);
        sonic_velocity = gas_sonic_velocity(*gas_thermo, gamma_s);
        Ae_At = area_per_mdot(*gas_thermo, velocity)/A_mdot_thrt;

        dlogp_dlogA = gamma_s * velocity * velocity / (velocity*velocity - sonic_velocity*sonic_velocity);
        residual = dlogp_dlogA * (std::log(expansion_ratio) - std::log(Ae_At));
        double log_pinf_pe = std::log(pressure_ratio) + residual;

        pressure_ratio = std::exp(log_pinf_pe);
    }
    return {true, save_gas_state(*gas_thermo)};

}


double FrozenNozzle::iterate_temperature(
    std::shared_ptr<Cantera::ThermoPhase>& gas_thermo,
    const ThroatCondition& throat_condition, 
    double pressure_ratio, double T_guess, double abstol) {
    
    //To solve the gas state, we need to iterate to find the exit temperature
    //initial guess
    double T_exit = T_guess;

    gas_thermo->setState_TP(throat_condition.P_inlet/pressure_ratio, T_exit);
    gas_thermo->equilibrate("TPX");

    double Cp = gas_thermo->cp_mass(); //TODO: unsure if we can just use this or whether it needs to be equilibrium Cp
    double dlnT = (throat_condition.S_inlet - gas_thermo->entropy_mass())/Cp;

    int maxiter = 8;
    int iters = 0;

    while (std::abs(dlnT) > abstol) {
        iters++;
        if (iters >= maxiter){
            return -1;
        }
        
        double lnT_exit = std::log(T_exit) + dlnT;
        T_exit = std::exp(lnT_exit);

        gas_thermo->setState_TP(T_exit, gas_thermo->pressure());
        gas_thermo->equilibrate("TPX");
        Cp = gas_thermo->cp_mass();
        dlnT = (throat_condition.S_inlet - gas_thermo->entropy_mass())/Cp;
    }

    return T_exit;
}

} //namespace Goddard