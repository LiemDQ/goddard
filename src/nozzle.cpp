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

NozzleResult NozzleBase::solve(const NozzleConditions& conditions) {
    
    switch (conditions.expansion_type) {
        case ExpansionType::SUPERSONIC_AREA_RATIO:
            return this->solve_supersonic_area_expansion();
        case ExpansionType::SUBSONIC_AREA_RATIO:
            return this->solve_subsonic_area_expansion();
        default:
            return this->solve_pressure_ratio();
    }
}

ThroatResult solve_throat_conditions(
    Cantera::Solution& inlet_gas, 
    bool frozen, 
    double abstol
) {

    auto inlet_thermo = inlet_gas.thermo();

    double P_inlet = inlet_thermo->pressure();
    double S_inlet = inlet_thermo->entropy_mass();
    double H_inlet = inlet_thermo->enthalpy_mass();
    
    std::vector<double> X_inlet(inlet_thermo->nSpecies());
    inlet_thermo->getMoleFractions(X_inlet.data());

    EquilibriumProperties thermo_props = get_thermo_equilibrium_properties(inlet_gas);

    double gamma_s = thermo_props.gamma_s;
    double P_throat = P_inlet / std::pow((gamma_s+1)/2,gamma_s/(gamma_s-1));
    
    const int max_iters = 5;
    int iter = 0;
    double Mach = 1.0; //throat mach number is 1 by definition
    double residual = 1.0;

    std::shared_ptr<Cantera::Solution> throat_gas = copy_solution(inlet_gas);
    auto throat_state = throat_gas->thermo();

    while (residual > abstol){
        if (iter >= max_iters){
            //show error message or throw exception
            return {};
        }
        P_throat = P_throat * (1 + gamma_s * Mach* Mach)/(1+ gamma_s);
        
        throat_state->setState_SP(S_inlet,P_throat);

        //if equilibrium conditions are selected, the composition must reach chemical
        //equilibrium in the throat. 
        if (!frozen){
            throat_state->equilibrate("SP");
            thermo_props = get_thermo_equilibrium_properties(inlet_gas);
            gamma_s = thermo_props.gamma_s;
        }
        
        double velocity = gas_isenthalpic_velocity(*inlet_thermo, H_inlet);
        double sonic_velocity = gas_sonic_velocity(*inlet_thermo, gamma_s);

        Mach = velocity/sonic_velocity;
        residual = std::abs(1.0 - 1.0/(Mach * Mach));

        iter++;
    }

    return {true, H_inlet, P_inlet, S_inlet};
}


NozzleResult solve_equilibrium_supersonic_area_expansion(
    Cantera::Solution& throat_gas, 
    const ThroatResult& result, 
    double expansion_ratio,
    double abstol
) {
    auto thrt_thermo = throat_gas.thermo();

    double velocity = gas_isenthalpic_velocity(*thrt_thermo, result.H_stagnation);
    
    //throat area/mdot is constant
    const double A_mdot_thrt = area_per_mdot(*thrt_thermo, velocity);

    
    EquilibriumProperties equilibrium_props = get_thermo_equilibrium_properties(throat_gas);
    double gamma_s = equilibrium_props.gamma_s;
    
    //initial guess for area ratios > 2
    double pressure_ratio = std::exp(gamma_s+1.4*std::log10(expansion_ratio));
    double P_exit = result.P_inlet/pressure_ratio;

    std::shared_ptr<Cantera::Solution> exit_gas = copy_solution(throat_gas);
    auto exit_thermo = exit_gas->thermo();
    exit_thermo->setState_SP(result.S_inlet, P_exit);
    exit_thermo->equilibrate("SP");
    double Ae_At = exit_thermo->temperature() / (exit_thermo->pressure() * velocity * exit_thermo->meanMolecularWeight())/A_mdot_thrt;

    int iters = 0;
    int max_iter = 10;
    double residual = 1.0;
    double sonic_velocity = 0.0;
    double dlogp_dlogA = 0.0;

    while (std::abs(residual) > abstol) {
        iters++;
        if (iters >= max_iter) {
            //show error message or throw exception
            return {false, nullptr};
        }
        EquilibriumProperties eqprops = get_thermo_equilibrium_properties(*exit_gas);
        velocity = gas_isenthalpic_velocity(*exit_thermo, result.H_stagnation);
        sonic_velocity = gas_sonic_velocity(*exit_thermo, eqprops.gamma_s);
        Ae_At = area_per_mdot(*exit_thermo, velocity)/A_mdot_thrt;

        dlogp_dlogA = eqprops.gamma_s * velocity * velocity / (velocity*velocity - sonic_velocity*sonic_velocity);
        residual = dlogp_dlogA * (std::log(expansion_ratio) - std::log(Ae_At));
        double log_pinf_pe = std::log(pressure_ratio) + residual;

        pressure_ratio = std::exp(log_pinf_pe);
        P_exit = result.P_inlet / pressure_ratio;

        exit_thermo->setState_SP(result.S_inlet, P_exit);
        exit_thermo->equilibrate("SP");
    }
    return {true, exit_gas};
}

NozzleResult solve_equilibrium_pressure_ratio(
    Cantera::Solution& throat_gas,
    const ThroatResult& throat_result,
    double pressure_ratio,
    double abstol
) {
    auto throat_thermo = throat_gas.thermo();
    
    //To solve the gas state, we need to iterate to find the exit temperature
    //initial guess
    double T_exit = throat_thermo->temperature();

    std::shared_ptr<Cantera::Solution> exit_gas = copy_solution(throat_gas);
    auto exit_thermo = exit_gas->thermo();

    exit_thermo->setPressure(throat_result.P_inlet/pressure_ratio);
    exit_thermo->setTemperature(T_exit);
    exit_thermo->equilibrate("TP");

    EquilibriumProperties eqprops = get_thermo_equilibrium_properties(*exit_gas);
    double dlnT = (throat_result.S_inlet - exit_thermo->entropy_mass())/eqprops.spec_heat_p;

    int maxiter = 8;
    int iters = 0;

    while (std::abs(dlnT) > abstol) {
        iters++;
        if (iters >= maxiter){
            return {false, nullptr};
        }
        
        double lnT_exit = std::log(T_exit) + dlnT;
        T_exit = std::exp(lnT_exit);

        exit_thermo->setState_TP(T_exit, exit_thermo->pressure());
        exit_thermo->equilibrate("TP");
        eqprops = get_thermo_equilibrium_properties(*exit_gas);
        dlnT = (throat_result.S_inlet - exit_thermo->entropy_mass())/eqprops.spec_heat_p;
    }

    return {true, exit_gas};
}


} //namespace Goddard


