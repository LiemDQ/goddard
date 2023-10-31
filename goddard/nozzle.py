'''
Utility functions for various nozzle conditions. Note that the classical compressible 
flow equations cannot be directly applied, because they assume that the heat capacity 
ratio is constant, which holds for frozen flow but not for equilibrium flow.
'''

import numpy as np
from scipy.optimize import brentq
import cantera as ct 
from thermo import get_thermo_derivatives, get_thermo_properties, speed_of_sound
from utils import to_si, copy_solution
from pint import Quantity
from abc import ABC, abstractmethod, abstractstaticmethod
from typing import Optional

class ExitConditions:
    subsonic_ratio: np.array
    supersonic_ratio: np.array
    pressure_ratio: np.array
    
    def __init__(self, supersonic_ratio=None, subsonic_ratio=None, pressure_ratio=None) -> None:
        self.supersonic_ratio = supersonic_ratio
        self.subsonic_ratio = subsonic_ratio
        self.pressure_ratio = pressure_ratio
        if not self.supersonic_ratio and not subsonic_ratio and not pressure_ratio:
            raise ValueError("No nozzle exit conditions are specified.")

class Nozzle(ABC):
    
    def __init__(self, inlet_gas: ct.Solution, exit_conditions: ExitConditions, gamma_s: float | None = None) -> None:
        super().__init__()
        self.inlet = inlet_gas
        self.exit_conditions = exit_conditions
        self.inlet_pressure = inlet_gas.P
        self.inlet_gamma = gamma_s if gamma_s else self._get_gamma()

        #update ct solution values
        self.throat = self.get_throat_conditions(copy_solution(inlet_gas))
        self.sonic_velocity = speed_of_sound(self.throat, self.inlet_gamma)
        self.exit = self.get_exit_conditions()
    
    def get_cstar(self):
        '''
        Calculate the C* value for a given combustion temperature and set of gas properties. C* is also known as the characteristic velocity,
        and is a measure of the energy present in the fluid as a result of combustion. It is used as a metric for engine performance independent
        of nozzle expansion.
        '''
        gamma = self.inlet_gamma
        mw = self.inlet.mean_molecular_weight
        temperature = self.inlet.T
        return get_cstar(gamma, temperature, mw)
    

    @abstractmethod
    def _get_gamma(self):
        pass

    @abstractmethod    
    def get_throat_conditions(self, gas: ct.Solution) -> ct.Solution:
        pass
    
    @abstractmethod
    def get_exit_conditions(self) -> ct.Solution:
        pass

    def update(self):
        self.throat = self.get_throat_conditions(copy_solution(self.inlet))
        self.sonic_velocity = speed_of_sound(self.throat, self.inlet_gamma)
        self.inlet_pressure = self.inlet.P
        self.exit = self.get_exit_conditions()


class EquilibriumNozzle(Nozzle):

    def _get_gamma(self):
        return get_thermo_properties(self.inlet)
        
    def get_throat_conditions(self, gamma_chamber: float) -> ct.Solution:
        '''
        Assumptions: 
        - Flow is isentropic
        - Equilibrium flow
        '''
        gas = self.inlet
        P_throat = self.inlet_pressure / np.power((gamma_chamber + 1)/2., gamma_chamber/(gamma_chamber-1))
        # CEA defaults
        max_iter_throat = 5
        tolerance_throat = 0.4e-4

        gamma_s = gamma_chamber
        
        M = 1.0 #throat mach is 1.0 by definition
        num_iter = 0
        residual = 1
        # print(f"P throat: {P_throat}")
        while residual > tolerance_throat:
            num_iter += 1
            if num_iter == max_iter_throat:
                # TODO: show error message
                break

            P_throat = P_throat * (1 + gamma_s * M**2)/(1 + gamma_s)
            # print(f"P_throat: {to_si(P_throat)}")
            gas.SPX = self.inlet.s, P_throat, self.inlet.X
            gas.equilibrate('SP')
            derivs = get_thermo_derivatives(gas)
            dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(gas)

            velocity = _get_velocity(gas, self.inlet.h)
            sonic_velocity = _get_sonic_velocity(gas, gamma_s)
            # print(f"Velocity: {velocity}, sonic velocity: {sonic_velocity}")
            M = velocity/sonic_velocity

            residual = np.abs(1.0 - 1/M**2)
        
        return gas

    def get_exit_conditions(self) -> ct.Solution:
        '''
        Calculate exit conditions, based on expansion area ratio.
        '''
        A_mdot_thr = self.throat.T / (self.throat.P * velocity * self.throat.mean_molecular_weight) #remains constant
        #TODO: some filename passing to ensure that the exit gas uses the same file as the throat gas
        gas_exit = ct.Solution(thermo=self.throat.thermo_model, species=self.throat.species())
        gas_exit.SPX = self.throat.SPX

        #TODO: initial guess only valid for area ratios > 2
        pressure_ratio = np.exp(gamma_s + 1.4 * np.log(area_ratio)) #initial guess 
        p_exit = self.inlet_pressure/pressure_ratio

        gas_exit.SP = self.inlet.s, p_exit
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
            dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(gas_exit)
            velocity = _get_velocity(gas_exit, self.inlet.h)
            sonic_velocity = _get_sonic_velocity(gas_exit, gamma_s)

            Ae_At = gas_exit.T / (gas_exit.P * velocity * gas_exit.mean_molecular_weight) / A_mdot_thr
            dlogp_dlogA = gamma_s * velocity** 2 / (velocity**2 - sonic_velocity**2)
            residual = dlogp_dlogA * (np.log(area_ratio) - np.log(Ae_At))
            log_pinf_pe = np.log(pressure_ratio) + residual

            pressure_ratio = np.exp(log_pinf_pe)
            p_exit = to_si(self.inlet_pressure) /  pressure_ratio

            gas_exit.SP = self.inlet.s, p_exit
            gas_exit.equilibrate('SP')
        
        return gas_exit

class FrozenNozzle(Nozzle):
    def __init__(self, inlet_gas: ct.Solution, exit_conditions: ExitConditions, gamma_s: float | None = None, NFZ: int = 1) -> None:
        super().__init__(inlet_gas, exit_conditions, gamma_s)
        if NFZ < 1 or NFZ > 5 or not isinstance(NFZ, int):
            raise ValueError("Invalid NFZ value: must be integer between 1 and 5. See NASA RP1311 Part II (Users Manual), p. 18-19")
        self.NFZ = NFZ
    
    def get_throat_conditions(self, gas: ct.Solution) -> ct.Solution:
        """
        Throat conditions for frozen flow. See Gordon & McBride, 1994, Section 6.5.2.
        """
        gas = self.inlet
        P_throat = self.inlet_pressure / np.power((self.inlet_gamma + 1)/2., self.inlet_gamma/(self.inlet_gamma-1))
        # CEA defaults
        max_iter_throat = 5
        tolerance_throat = 0.4e-4

        gamma_s = self.inlet_gamma
        
        M = 1.0 # throat mach is 1.0 by definition
        num_iter = 0
        residual = 1
        while residual > tolerance_throat:
            num_iter += 1
            if num_iter == max_iter_throat:
                # TODO: show error message
                break

            P_throat = P_throat * (1 + gamma_s * M**2)/(1 + gamma_s)
            gas.SPX = self.inlet.s, P_throat, self.inlet.X
            # unlike with the equilibrium calculation, the gas mixture is never equilibrated
            # and thus gamma never changes

            # gas.equilibrate('SP')
            # derivs = get_thermo_derivatives(gas)
            velocity = _get_velocity(gas, self.inlet.h)
            sonic_velocity = speed_of_sound(gas, gamma_s)
            M = velocity/sonic_velocity

            residual = np.abs(1.0 - 1/M**2)
        
        return gas
    
    def get_exit_conditions(self) -> ct.Solution:
        """
        See Gordon & McBride, 1994, Section 6.5.
        """
        A_mdot_thr = self.throat.T / (self.throat.P * velocity * self.throat.mean_molecular_weight) #remains constant
        # we iterate to obtain the nozzle exit temperature 
        T_e = self.throat.T #initial guess
        gas_exit = copy_solution(self.throat)

        gamma_s = gas_exit.cp/gas_exit.cv
        
        pressure_ratio = np.exp(gamma_s + 1.4 * np.log(area_ratio)) #initial guess 
        # print(f"Area ratio: {area_ratio}")
        # print(f"P_chamber: {p_chamber}")
        # print(f"Pressure ratio: {pressure_ratio}")
        p_exit = self.inlet_pressure/pressure_ratio
        # print(f"P_exit: {p_exit}")
        gas_exit.TPX = T_e, p_exit, gas_exit.X

        gamma_s = gas_exit.cp/gas_exit.cv

        dlnT = (self.inlet.s - gas_exit.s)/gas_exit.cp

        dlogT_tolerance = 0.5e-4
        pressure_tolerance = 0.4e-4
        maxiter = 10
        maxiter_temp = 8
        n = 0

        residual = 1
        while np.abs(residual) > pressure_tolerance:
            n += 1
            if n == maxiter:
                #TODO: error message if nozzle fails to converge
                break
            # print(f"Iteration: {n}")
            m = 0
            dlnT = 1000
            while np.abs(dlnT) >= dlogT_tolerance:
                m += 1
                if m >= maxiter_temp:
                    break
                gas_exit.TPX = T_e, p_exit, gas_exit.X
                dlnT = (self.inlet.s - gas_exit.s)/gas_exit.cp
                T_e = np.exp(np.log(T_e) + dlnT)
                # print(f"T_exit: {T_e}")

            gamma_s = gas_exit.cp/gas_exit.cv #gamma_s = gamma for frozen flow
            velocity = _get_velocity(gas_exit, self.inlet.h)
            sonic_velocity = _get_sonic_velocity(gas_exit, gamma_s)

            Ae_At = gas_exit.T / (gas_exit.P * velocity * gas_exit.mean_molecular_weight) / A_mdot_thr
            dlogp_dlogA = gamma_s * velocity** 2 / (velocity**2 - sonic_velocity**2)
            residual = dlogp_dlogA * (np.log(area_ratio) - np.log(Ae_At))
            log_pinf_pe = np.log(pressure_ratio) + residual


            pressure_ratio = np.exp(log_pinf_pe)
            p_exit = self.inlet_pressure /  pressure_ratio

        print(f"number of iterations: {n}")

        return gas_exit

def get_exit_conditions(gas_throat: ct.Mixture, area_ratio: float, P_chamber: float, gas_chamber: ct.Solution, gamma_s, velocity):
    '''
    Calculate exit conditions, based on expansion area ratio.
    '''
    A_mdot_thr = gas_throat.T / (gas_throat.P * velocity * gas_throat.mean_molecular_weight) #remains constant
    #TODO: some filename passing to ensure that the exit gas uses the same file as the throat gas
    gas_exit = ct.Solution(thermo=gas_throat.thermo_model, species=gas_throat.species())
    gas_exit.SPX = gas_throat.SPX

    #TODO: initial guess only valid for area ratios > 2
    pressure_ratio = np.exp(gamma_s + 1.4 * np.log(area_ratio)) #initial guess 
    p_exit = P_chamber/pressure_ratio

    gas_exit.SP = gas_chamber.s, p_exit
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
        dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(gas_exit)
        velocity = _get_velocity(gas_exit, gas_chamber.h)
        sonic_velocity = _get_sonic_velocity(gas_exit, gamma_s)

        Ae_At = gas_exit.T / (gas_exit.P * velocity * gas_exit.mean_molecular_weight) / A_mdot_thr
        dlogp_dlogA = gamma_s * velocity** 2 / (velocity**2 - sonic_velocity**2)
        residual = dlogp_dlogA * (np.log(area_ratio) - np.log(Ae_At))
        log_pinf_pe = np.log(pressure_ratio) + residual

        pressure_ratio = np.exp(log_pinf_pe)
        p_exit = to_si(P_chamber) /  pressure_ratio

        gas_exit.SP = gas_chamber.s, p_exit
        gas_exit.equilibrate('SP')
    
    return gas_exit

def get_exit_conditions_pressure_ratio(gas_throat: ct.Mixture, pressure_ratio):
    pass

def get_throat_conditions_frozen(gas: ct.Mixture, P_chamber: float, gas_chamber: ct.Mixture, gamma_chamber: float):
    """
    Throat conditions for frozen flow. See Gordon & McBride, 1994, Section 6.5.2.
    """
    P_throat = P_chamber / np.power((gamma_chamber + 1)/2., gamma_chamber/(gamma_chamber-1))
    # CEA defaults
    max_iter_throat = 5
    tolerance_throat = 0.4e-4

    gamma_s = gamma_chamber
    
    M = 1.0 # throat mach is 1.0 by definition
    num_iter = 0
    residual = 1
    while residual > tolerance_throat:
        num_iter += 1
        if num_iter == max_iter_throat:
            # TODO: show error message
            break

        P_throat = P_throat * (1 + gamma_s * M**2)/(1 + gamma_s)
        gas.SPX = gas_chamber.s, P_throat, gas_chamber.X
        # unlike with the equilibrium calculation, the gas mixture is never equilibrated
        # and thus gamma never changes

        # gas.equilibrate('SP')
        # derivs = get_thermo_derivatives(gas)
        velocity = _get_velocity(gas, gas_chamber.h)
        sonic_velocity = _get_sonic_velocity(gas, gamma_s)
        M = velocity/sonic_velocity

        residual = np.abs(1.0 - 1/M**2)
    
    return gas



def get_exit_conditions_frozen(gas_throat: ct.Mixture, area_ratio, p_chamber, gas_chamber: ct.Solution, gamma_s, velocity):
    """
    See Gordon & McBride, 1994, Section 6.5.
    """
    A_mdot_thr = gas_throat.T / (gas_throat.P * velocity * gas_throat.mean_molecular_weight) #remains constant
    # we iterate to obtain the nozzle exit temperature 
    T_e = gas_throat.T #initial guess
    gas_exit = copy_solution(gas_throat)

    gamma_s = gas_exit.cp/gas_exit.cv
    
    pressure_ratio = np.exp(gamma_s + 1.4 * np.log(area_ratio)) #initial guess 
    # print(f"Area ratio: {area_ratio}")
    # print(f"P_chamber: {p_chamber}")
    # print(f"Pressure ratio: {pressure_ratio}")
    p_exit = to_si(p_chamber)/pressure_ratio
    # print(f"P_exit: {p_exit}")
    gas_exit.TPX = T_e, p_exit, gas_exit.X

    gamma_s = gas_exit.cp/gas_exit.cv

    dlnT = (gas_chamber.s - gas_exit.s)/gas_exit.cp

    dlogT_tolerance = 0.5e-4
    pressure_tolerance = 0.4e-4
    maxiter = 10
    maxiter_temp = 8
    n = 0

    residual = 1
    while np.abs(residual) > pressure_tolerance:
        n += 1
        if n == maxiter:
            #TODO: error message if nozzle fails to converge
            break
        # print(f"Iteration: {n}")
        m = 0
        dlnT = 1000
        while np.abs(dlnT) >= dlogT_tolerance:
            m += 1
            if m >= maxiter_temp:
                break
            gas_exit.TPX = T_e, p_exit, gas_exit.X
            dlnT = (gas_chamber.s - gas_exit.s)/gas_exit.cp
            T_e = np.exp(np.log(T_e) + dlnT)
            # print(f"T_exit: {T_e}")

        gamma_s = gas_exit.cp/gas_exit.cv #gamma_s = gamma for frozen flow
        velocity = _get_velocity(gas_exit, gas_chamber.h)
        sonic_velocity = _get_sonic_velocity(gas_exit, gamma_s)

        Ae_At = gas_exit.T / (gas_exit.P * velocity * gas_exit.mean_molecular_weight) / A_mdot_thr
        dlogp_dlogA = gamma_s * velocity** 2 / (velocity**2 - sonic_velocity**2)
        residual = dlogp_dlogA * (np.log(area_ratio) - np.log(Ae_At))
        log_pinf_pe = np.log(pressure_ratio) + residual

        
        # print(f"Velocity: {velocity}")
        # print(f"Gamma: {gamma_s}")
        # print(f"Sonic velocity: {sonic_velocity}")
        # print(f"dlogp_dlogA: {dlogp_dlogA}")

        # print(f"Ae_at: {Ae_At}")
        # print(f"Pressure residual: {residual}")

        pressure_ratio = np.exp(log_pinf_pe)
        p_exit = to_si(p_chamber) /  pressure_ratio

    

    print(f"number of iterations: {n}")

    return gas_exit
    
def get_exit_conditions_pressure_ratio_frozen(gas_throat: ct.Mixture, pressure_ratio, p_chamber, gas_chamber: ct.Solution, velocity):
    # we iterate to obtain the nozzle exit temperature 
    
    T_e = gas_chamber.T #initial guess
    gas_exit = copy_solution(gas_throat)
    gas_exit.P = p_chamber / pressure_ratio

    dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(gas_exit)
    dlnT = (gas_chamber.s - gas_exit.s)/cp

    dlogT_residual = 0.5e-4
    maxiter = 8
    n = 0
    while np.abs(dlnT) >= dlogT_residual:
        n += 1
        if n >= maxiter:
            break
        
        lnT_e = np.log(T_e) + dlnT
        T_e = np.exp(lnT_e)
        gas_exit.TP = T_e, gas_exit.P

        _,_,cp,_ = get_thermo_properties(gas_exit)
        dlnT = (gas_chamber.s - gas_exit.s)/cp
    
    return gas_exit

def _estimate_pressure_ratio():
    pass

def _get_velocity(gas: ct.Mixture, stagnation_enthalpy):
    '''
    Velocity in isentropic supersonic flow can be found from the difference in enthalpy
    between the gas in motion and the stagnation enthalpy. 
    '''
    return np.sqrt(2*(stagnation_enthalpy - gas.enthalpy_mass))

def _get_sonic_velocity(gas: ct.Solution, gamma):
    return np.sqrt(ct.gas_constant * gas.T * gamma / gas.mean_molecular_weight)


def get_cstar(gamma, temperature, molecular_weight):
    '''
    Calculate the C* value for a given combustion temperature and set of gas properties. C* is also known as the characteristic velocity,
    and is a measure of the energy present in the fluid as a result of combustion. It is used as a metric for engine performance independent
    of nozzle expansion.
    '''
    return (
        np.sqrt(ct.gas_constant * temperature / (molecular_weight * gamma)) *
        np.power(2 / (gamma + 1), -(gamma + 1) / (2*(gamma - 1)))
        )

def get_mach_from_subsonic_area_ratio(gamma, ratio):
    '''
    For a given area ratio, calculate the mach in the smaller area given the mach number in the larger area.
    There are two possible solutions for a given ratio (one subsonic, one supersonic). 
    This function always returns the subsonic solution.
    '''
    return brentq(lambda mach: get_area_ratio_from_mach_num(gamma, mach)-ratio, 0, 1.0)
    

def get_mach_from_supersonic_area_ratio(gamma, ratio):
    '''
    For a given area ratio, calculate the mach in the smaller area given the mach number in the larger area.
    There are two possible solutions for a given ratio (one subsonic, one supersonic). 
    This function always returns the supersonic solution.
    '''

    #ideally, should be inf, but this would not converge for certain ill-conditioned edge cases.
    upper_mach_limit = 10,000,000.0 
    return brentq(lambda mach: get_area_ratio_from_mach_num(gamma, mach)-ratio, 1.0, upper_mach_limit)

def get_area_ratio_from_mach_num(gamma, M):
    '''
    Area ratio (w.r.t throat) for a given mach number and isentropic expansion factor.
    '''
    return 1/M * ((1+ (gamma-1)/2*M**2) /((gamma+1)/2))**((gamma+1)/(2*(gamma-1)))

def get_isp(gas, gamma, enthalpy):
    '''
    Calculate specific impulse given an exhaust gas thermodynamic state and combustion chamber enthalpy. 
    '''
    derivs = get_thermo_derivatives(gas)
    velocity = _get_velocity(gas, enthalpy)
    cstar = get_cstar(gamma, gas.T, gas.mean_molecular_weight)
    CF = velocity/cstar
    g0 = 9.80655 #gravitational acceleration
    return velocity, velocity/g0
    

def get_ivac(gas, isp):
    '''
    Ivac also includes the thrust from pressure forces. This assumes isp is provided in (m/s)
    and not s.
    '''
    ivac = isp + gas.T * ct.gas_constant / (isp * gas.mean_molecular_weight)
    g0 = 9.80655 #gravitational acceleration

    return ivac, ivac/g0