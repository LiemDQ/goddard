#include <cmath>
#include <vector>
#include <stdexcept>
#include <iostream>
#include "eigen3/Eigen/Dense"
#include "cantera/core.h"
#include "goddard/error.hpp"
#include "goddard/shocks.hpp"
#include "goddard/gas_dynamics.hpp"

namespace Goddard {

ShockResult normal_shock(double mach, double gamma) {

    ShockResult result;
    // check if input values are valid
    if (mach < 1.0 || gamma < 1.0){
        result.valid = false;
        result.static_pressure_ratio = -1.0;
        result.static_temperature_ratio = -1.0;
        return result;
    }
    
    //mach numbers
    double M1 = mach;
    double M2 = std::sqrt( (M1*M1 * (gamma - 1) + 2)/(2 *  gamma * M1 * M1 - (gamma - 1)));
    result.mach_in = M1;
    result.mach_out = M2;
    
    //dynamic pressures
    result.static_pressure_ratio = (2 * gamma * M1 * M1)/(gamma + 1) - (gamma - 1)/(gamma + 1);

    //static temperatures
    result.static_temperature_ratio = ((1 + (gamma - 1)/2*M1*M1)
        * (2*gamma/(gamma -1) * M1*M1 - 1)
        / (M1*M1 * (2 * gamma /(gamma - 1) + (gamma -1)/2)));

    //stagnation pressures
    result.total_pressure_ratio = std::pow(((gamma +1)/2 * M1 * M1)/(1 + (gamma-1)/2* M1 * M1), gamma/(gamma-1))
        * std::pow(1.0/(2*gamma/(gamma+1)*M1*M1 - (gamma-1)/(gamma+1)), 1/(gamma-1));

    result.valid = true;
    return result;
}

double normal_shock_control_factor(int iter) {
    double cf;
    if (iter > 20) {
        cf = 0.04879016;
    }
    else if (iter > 12) {
        cf = 0.09531018;
    }
    else if (iter > 4) {
        cf = 0.22314355;
    }
    else {
        cf = 0.40546511;
    }
    return cf;
}

ShockResult normal_shock(Cantera::ThermoPhase& thermo, double mach) {
    ShockResult result;
    result.mach_in = mach;

    if (mach < 1.0) {
        result.valid = false;
        result.static_pressure_ratio = -1.0;
        result.static_temperature_ratio = -1.0;
        return result;
    }
  
    const double divR = 1.0/Cantera::GasConstant;

    const double h1 = thermo.enthalpy_mass();
    const double P1 = thermo.pressure();
    const double T1 = thermo.temperature();
    const double rho1 = thermo.density();
    const double mw1 = thermo.meanMolecularWeight();

    std::vector<double> state1(thermo.stateSize());
    thermo.saveState(state1);
    std::vector<double> state2(state1);
    
    // TODO: change this based on whether equilibrium chemistry is used
    const double gamma1 = thermo.cp_mass()/thermo.cv_mass();
    const double u1 = gas_sonic_velocity(thermo, gamma1)*mach;
    const double h_stag = gas_stagnation_enthalpy(thermo, u1);
    const double P_stag1 = gas_stagnation_pressure(thermo, u1);
    // initial guesses
    double P2_P1 = (2*gamma1*mach*mach-gamma1+1)/(gamma1+1);
    // NOTE: for equilibrium we need to perform an isobaric equilibrium solve for temperature ratio
    double T2_T1 = P2_P1 * (2/(mach*mach) + gamma1-1)/(gamma1+1);
    double P2 = P2_P1 * P1;
    double T2 = T2_T1 * T1;

    const double dP_coeff = mw1*u1*u1*divR/T1; // MW1*u^2/(R*T1)
    const double dh_coeff = u1*u1*divR; //u1^2/R

    double logP2_P1 = log(P2_P1);
    double logT2_T1 = log(T2_T1);

    thermo.setState_TP(T2, P2);
    
    int k = 0;
    int max_iters = 100;
    double residual = 100.0;
    double control_factor_coeff = normal_shock_control_factor(k);
    const double abstol = 5e-5;

    // Use Newton's method to solve for shock conditions.
    // See NASA RP-1311 Part I, section 7. 
    while (residual >= abstol) {
        if (k > max_iters) 
            throw ConvergenceError("Normal shock properties failed to converge.", k, abstol, residual);
        std::cout << "Iter: " << k << "\n";
        std::cout << "P2/P1: " << P2_P1 << ", T2/T1: " << T2/T1 << "\n";
        double mw2 = thermo.meanMolecularWeight();
        double cp2 = thermo.cp_mass(); 
        double h2 = thermo.enthalpy_mass();
        double rho2 = thermo.density();
        double rho1_rho2 = rho1/rho2;
        double rho1_rho2_sq = rho1_rho2*rho1_rho2;
        
        //volumetric derivatives
        double dlogV_dlogT_P = 1.0;
        double dlogV_dlogP_T = -1.0;

        double dh_coeff_rho = dh_coeff*rho1_rho2_sq;

        // partial derivatives
        double dP_dlogP2P1 = -rho1_rho2 * dP_coeff * dlogV_dlogP_T - P2_P1;
        double dP_dlogT2T1 = -rho1_rho2 * dP_coeff * dlogV_dlogT_P; 
        double dh_dlogP2P1 = -dh_coeff_rho * dlogV_dlogP_T + T2/mw2 * (dlogV_dlogT_P - 1);
        double dh_dlogT2T1 = -dh_coeff_rho * dlogV_dlogT_P - T2*cp2*divR;
        
        double P2P1_minus_Pstar = P2_P1 - 1 + dP_coeff*(rho1_rho2 - 1);
        double h2_minus_hstar_R = (h2-h1)*divR - 0.5*dh_coeff*(1- rho1_rho2_sq);

       
        // directly solve system of equations
        double dlogT2_T1 = (P2P1_minus_Pstar - dP_dlogP2P1/dh_dlogP2P1 * h2_minus_hstar_R)/(dP_dlogT2T1 - dP_dlogP2P1/dh_dlogP2P1 * dh_dlogT2T1);
        double dlogP2_P1 = (h2_minus_hstar_R - dh_dlogT2T1*dlogT2_T1)/dh_dlogP2P1;

        double abs_dlogP2_P1 = abs(dlogP2_P1);
        double abs_dlogT2_T1 = abs(dlogT2_T1);
        
        control_factor_coeff = normal_shock_control_factor(k);
        double control_factor = std::min(control_factor_coeff/abs_dlogP2_P1, control_factor_coeff/abs_dlogT2_T1);
        control_factor = std::min(control_factor, 1.0);
        
        logP2_P1 += control_factor * dlogP2_P1;
        logT2_T1 += control_factor * dlogT2_T1;

        P2_P1 = exp(logP2_P1);
        P2 = P2_P1*P1;
        T2 = exp(logT2_T1)*T1;
        thermo.setState_TP(T2, P2);

        residual = std::max(abs_dlogP2_P1, abs_dlogT2_T1);
        k++;
    }
    
    double gamma2 = thermo.cp_mass()/thermo.cv_mass();
    double u2 = gas_isenthalpic_velocity(thermo, h_stag);
    double P_stag2 = gas_stagnation_pressure(thermo, u2);

    result.valid = true;
    result.mach_out = u2/gas_sonic_velocity(thermo, gamma2);
    result.static_pressure_ratio = P2_P1;
    result.static_temperature_ratio = T2/T1;
    result.total_pressure_ratio = P_stag2/P_stag1;

    return result;
}

std::pair<double, double> oblique_shock_wave_angle(double mach, double deflection_angle, double gamma) {
    double tan_theta = tan(deflection_angle);
    double M2 = mach*mach;
    double stagnation = stagnation_factor(mach, gamma);
    double lambda = std::sqrt((M2-1)*(M2-1) - 3*stagnation*(1 + (gamma+1)/2.0 * M2)*tan_theta*tan_theta);
    double chi = (pow(M2-1, 3.0)
        - 9*stagnation*(stagnation+(gamma+1)/4.0*M2*M2)
        *tan_theta*tan_theta)
        / pow(lambda, 3.0);

    double weak_beta_cos = cos((4*M_PI + acos(chi))/3);
    double strong_beta_cos = cos((acos(chi))/3);

    double weak_beta = atan((M2-1 + 2*lambda*weak_beta_cos)
        / (3 * stagnation * tan_theta));

    double strong_beta = atan((M2-1 + 2*lambda*strong_beta_cos)
        / (3 * stagnation * tan_theta));

    return {weak_beta, strong_beta};
}

double oblique_shock_deflection_angle(double mach, double wave_angle, double gamma) {
    double sin_beta = sin(wave_angle);
    double cos_2beta = cos(2*wave_angle);
    double cot_beta = 1.0/tan(wave_angle);
    return atan(2.0*cot_beta
        *(mach*mach*sin_beta*sin_beta -1)
        /((mach*mach )*(gamma + cos_2beta)+2));   
}

double oblique_shock_max_deflection(double /*mach*/, double /*gamma*/) {
    throw NotImplementedError("Max deflection is not implemented.");
}


ObliqueShockResult oblique_shock_from_deflection(
    double mach, double deflection_angle, double gamma, bool weak)
{
    ObliqueShockResult result;
    
    auto [weak_beta, strong_beta] = oblique_shock_wave_angle(mach, deflection_angle, gamma);

    if (weak) {
        result.beta = weak_beta;
    }
    else {
        result.beta = strong_beta;
    }
    double mach_n1 = mach * sin(result.beta);
    result.theta = deflection_angle;
    result.shock = normal_shock(mach_n1, gamma);
    result.mach_in = mach;
    result.mach_out = result.shock.mach_out/sin(result.beta - result.theta);
    result.valid = result.shock.valid;
    return result;
}

ObliqueShockResult oblique_shock_from_wave_angle(
    double mach, double wave_angle, double gamma)
{
    ObliqueShockResult result;
    double mach_n1 = mach * sin(wave_angle);

    result.shock = normal_shock(mach_n1, gamma);
    result.theta = oblique_shock_deflection_angle(mach, wave_angle, gamma);
    result.beta = wave_angle;
    result.mach_in = mach;
    result.mach_out = result.shock.mach_out/sin(wave_angle - result.theta);
    result.valid = result.shock.valid;
    return result;
}


ObliqueShockResult oblique_shock_from_wave_angle(
    Cantera::ThermoPhase& gas, double mach, double wave_angle)
{
    ObliqueShockResult result;
    double mach_n1 = mach * sin(wave_angle);
    double gamma = gas.cp_mass()/gas.cv_mass();
    double u1 = gas_sonic_velocity(gas, gamma)*mach_n1;
    
    result.shock = normal_shock(gas, mach_n1);
    gamma = gas.cp_mass()/gas.cv_mass();

    double u2 = gas_sonic_velocity(gas, gamma)*result.shock.mach_out;

    result.theta = wave_angle - atan(u2/u1 * tan(wave_angle));
    result.beta = wave_angle;
    result.mach_in = mach;
    result.mach_out = result.shock.mach_out/sin(wave_angle - result.theta);
    result.valid = result.shock.valid;
    return result;
}


ObliqueShockResult oblique_shock_from_deflection(
    Cantera::ThermoPhase& gas, double mach, double deflection_angle, bool weak) 
{
    ObliqueShockResult result;
    const double gamma1 = gas.cp_mass()/gas.cv_mass();

    // this is the point where tan(beta-theta)/tan(beta) is maximized
    // 
    const double beta_peak = deflection_angle/2 + M_PI/4; 

    std::vector<double> state1(gas.stateSize());
    gas.saveState(state1);

    //initial guess
    auto [weak_beta, strong_beta] = oblique_shock_wave_angle(mach, deflection_angle, gamma1);
    
    double left_bound;
    double right_bound;
    double beta;
    
    if (weak) {
        left_bound = 0.0;
        right_bound = beta_peak;
        beta = weak_beta;
    }
    else {
        left_bound = beta_peak;
        right_bound = M_PI/2;
        beta = strong_beta;
    }

    ObliqueShockResult left_result = oblique_shock_from_wave_angle(gas, mach, left_bound);
    gas.restoreState(state1);
    ObliqueShockResult right_result = oblique_shock_from_wave_angle(gas, mach, right_bound);

    double left_residual = left_result.theta - deflection_angle;
    double right_residual = right_result.theta - deflection_angle;
    
    if ((left_residual)*(right_residual) >= 0) {
        throw std::runtime_error("Bisection method error: product of bounds should be negative.");
    }

    // bisection method to find root
    // the objective function is theta - theta_{calculated} where the second term 
    // is calculated from an oblique shock with the current beta angle.
    // In the regime bracketed by beta_peak, there is at least one root. 
    double residual = 1.0;
    double abstol = 1e-8;
    int k = 0;
    int max_iters = 100;
    while (abs(residual) > abstol) {
        if (k > max_iters)
            throw ConvergenceError("Wave angle failed to converge.", k, abstol, residual);
        
        gas.restoreState(state1);
        result = oblique_shock_from_wave_angle(gas, mach, beta);
        residual = result.theta - deflection_angle;
        if (residual*left_residual <= 0) {
            right_bound = beta;
            right_residual = residual;
        }
        else if (residual*right_residual <= 0) {
            left_bound = beta;
            left_residual = residual;            
        }
        else {
            throw std::runtime_error("Bisection method failed: Bounds are not of opposite sign.");
        }

        beta = (left_bound + right_bound)/2.0;
        k++;        
    }
    result.theta = deflection_angle;
    return result;
}


} //namespace Goddard