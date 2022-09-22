'''
Utility functions for various nozzle conditions. Note that the classical compressible 
flow equations cannot be directly applied, because they assume that the heat capacity 
ratio is constant, which holds for frozen flow but not for equilibrium flow.
'''

import numpy as np
from scipy.optimize import brentq
import cantera as ct 
from thermo import get_thermo_derivatives, get_thermo_properties, to_si
from pint import Quantity

def calculate_throat_conditions(gas: ct.Solution, p_chamber: Quantity, gamma_chamber):
    '''
    Assumptions: 
    - Flow is isentropic
    '''
    P_throat = P_chamber / np.power((gamma_chamber + 1)/2., gamma_chamber/(gamma_chamber-1))
    # CEA defaults
    max_iter_throat = 5
    tolerance_throat = 0.4e-4

    gamma_s = gamma_chamber
    
    M = 1.0 #throat mach is 1.0 by definition
    num_iter = 0
    residual = 1
    while residual > tolerance_throat:
        num_iter += 1
        if num_iter == max_iter_throat:
            # TODO: show error message
            break

        P_throat = P_throat * (1 + gamma_s * M**2)/(1 + gamma_s)
        gas.SPX = entropy_chamber, to_si(P_throat), mole_fractions_chamber
        gas.equilibrate('SP')
        derivs = get_thermo_derivatives(gas)
        dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(
            gas, derivs[0], derivs[1], derivs[2]
        )

        velocity = _get_velocity(gas, enthalpy_chamber)
        sonic_velocity = _get_sonic_velocity(gas, gamma_s)
        M = velocity/sonic_velocity

        residual = np.abs(1.0 - 1/M**2)
    
    return gas

def calculate_exit_conditions(gas_throat: ct.Solution, area_ratio, p_chamber, gamma_s, velocity, filename = 'nasa_gas.yaml'):
    A_mdot_thr = gas_throat.T / (gas_throat.P * velocity * gas_throat.mean_molecular_weight) #remains constant
    #TODO: some filename passing to ensure that the exit gas uses the same file as the throat gas
    gas_exit = ct.Solution(thermo=gas_throat.thermo_model, species=gas_throat.species())
    gas_exit.SPX = gas_throat.SPX

    
    pressure_ratio = np.exp(gamma_s + 1.4 * np.log(area_ratio)) #initial guess
    p_exit = to_si(p_chamber)/pressure_ratio

    gas_exit.SP = entropy_chamber, p_exit
    gas_exit.equilibrate('SP')

    Ae_At = gas_exit.T / (gas_exit.P * velocity * gas_exit.mean_molecular_weight)/A_mdot_thr
    num_iter = 0
    max_iter_exit = 10
    tolerance = 4e-5

    residual = 1
    while np.abs(residual) > tolerance:
        num_iter += 1
        if num_iter == max_iter_exit:
            #TODO: error message if nozzle fails to converge
            break

        derivs = get_thermo_derivatives(gas_exit)
        dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(
            gas_exit, derivs[0], derivs[1], derivs[2]
        )
        velocity = _get_velocity(gas_exit, enthalpy_chamber)
        sonic_velocity = _get_sonic_velocity(gas_exit, gamma_s)

        Ae_At = gas_exit.T / (gas_exit.P * velocity * gas_exit.mean_molecular_weight) / A_mdot_thr
        dlogp_dlogA = gamma_s * velocity** 2 / (velocity**2 - sonic_velocity**2)
        residual = dlogp_dlogA * (np.log(area_ratio) - np.log(Ae_At))
        log_pinf_pe = np.log(pressure_ratio) + residual

        pressure_ratio = np.exp(log_pinf_pe)
        p_exit = to_si(p_chamber) /  pressure_ratio

        gas_exit.SP = entropy_chamber, p_exit
        gas_exit.equilibrate('SP')
    
    return gas_exit


def _get_velocity(gas: ct.Solution, enthalpy):
    return np.sqrt(2*(enthalpy - gas.enthalpy_mass))

def _get_sonic_velocity(gas: ct.Solution, gamma):
    return np.sqrt(ct.gas_constant,* gas.T * gamma / gas.mean_molecular_weight)


def calculate_c_star(gamma, temperature, molecular_weight):
    '''
    Calculate the C* value for a given combustion temperature and set of gas properties. C* is also known as the characteristic velocity,
    and is a measure of the energy present in the fluid as a result of combustion. It is used as a metric for engine performance independent
    of nozzle expansion.
    '''
    return (
        np.sqrt(ct.gas_constant * temperature / (molecular_weight * gamma)) *
        np.power(2 / (gamma + 1), -(gamma + 1) / (2*(gamma - 1)))
        )

def calculate_mach_from_area_ratio_subsonic(gamma, ratio):
    '''
    For a given area ratio, calculate the mach in the smaller area given the mach number in the larger area.
    There are two possible solutions for a given ratio (one subsonic, one supersonic). 
    This function always returns the subsonic solution.
    '''
    return brentq(lambda mach: calculate_area_ratio_from_mach(gamma, mach)-ratio, 0, 1.0)
    

def calculate_mach_from_area_ratio_supersonic(gamma, ratio):
    '''
    For a given area ratio, calculate the mach in the smaller area given the mach number in the larger area.
    There are two possible solutions for a given ratio (one subsonic, one supersonic). 
    This function always returns the supersonic solution.
    '''
    upper_mach_limit = 10,000,000.0
    return brentq(lambda mach: calculate_area_ratio_from_mach(gamma, mach)-ratio, 1.0, upper_mach_limit)

def calculate_area_ratio_from_mach(gamma, M):
    '''
    Area ratio (w.r.t throat) for a given mach number and isentropic expansion factor.
    '''
    return 1/M * ((1+ (gamma-1)/2*M**2) /((gamma+1)/2))**((gamma+1)/(2*(gamma-1)))

def calculate_isp():
    pass

def calculate_isp_vac():
    pass

def calculate_pressure_ratio():
    pass

def calculate_expansion_ratio():
    pass
