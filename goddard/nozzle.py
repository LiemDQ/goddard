'''
Utility functions for various nozzle conditions. Note that the classical compressible 
flow equations cannot be directly applied, because they assume that the heat capacity 
ratio is constant, which holds for frozen flow but not for equilibrium flow.
'''

from abc import ABC, abstractmethod
from typing import Type
from dataclasses import dataclass

import numpy as np
import cantera as ct 

from thermo import get_thermo_derivatives, get_thermo_properties, get_speed_of_sound
from isp import get_cstar, get_velocity
from utils import to_si, copy_ct_solution, args_to_np_array
from error import SolverError
class ExpansionConditions:
    
    def __init__(self, supersonic_ratio=None, subsonic_ratio=None, pressure_ratio=None) -> None:
        self.supersonic_ratio = supersonic_ratio
        self.subsonic_ratio = subsonic_ratio
        self.pressure_ratio = pressure_ratio
        if not self.supersonic_ratio and not subsonic_ratio and not pressure_ratio:
            raise ValueError("No nozzle exit conditions are specified.")

def subsonic_ratio(value, *args):
    values = args_to_np_array(value, args)
    return ExpansionConditions(subsonic_ratio=values)

def supersonic_ratio(value, *args):
    values = args_to_np_array(value, args)
    return ExpansionConditions(supersonic_ratio=values)

def pressure_ratio(value, *args):
    values = args_to_np_array(value, args)
    return ExpansionConditions(pressure_ratio=values)

class NozzleBase(ABC):

    def __init__(self, inlet_gas: ct.Solution, inlet_states: ct.SolutionArray, expansion_conditions: ExpansionConditions, gamma_s: np.array | None = None) -> None:

        super().__init__()
        self.inlet = inlet_gas
        self.inlet_states = inlet_states 
        self.expansion_conditions = expansion_conditions
        self.inlet_pressures = inlet_states.P
        self.inlet_gammas = gamma_s if gamma_s else self._get_gamma()

        #update ct solution values
        self.throat_states = self.get_throat_conditions()
        self.sonic_velocity = get_speed_of_sound(self.throat_states, self.inlet_gammas)
        self.exit_states = self.get_exit_conditions()

    def get_exit_conditions(self) -> ct.SolutionArray:
        """Solve the nozzle equations. 

        Returns:
            ct.SolutionArray: A solution array of size n x m, where n is the number of inlet conditions 
        and m is the number of expansion conditions specified. If the input is an array of dimension N, 
        then the output is of dimension N+1. 
        """
        if self.expansion_conditions.supersonic_ratio:
            return self._get_exit_conditions_supr()
        elif self.expansion_conditions.pressure_ratio:
            return self._get_exit_conditions_pr()
        else:
            raise NotImplementedError("Subsonic expansion is not implemented.")
    
    def get_cstar(self):
        '''
        Calculate the C* value for a given combustion temperature and set of gas properties. C* is also known as the characteristic velocity,
        and is a measure of the energy present in the fluid as a result of combustion. It is used as a metric for engine performance independent
        of nozzle expansion.
        '''
        gamma = self.inlet_gammas
        mw = self.inlet.mean_molecular_weight
        temperature = self.inlet.T
        return get_cstar(gamma, temperature, mw)
        
    def update(self):
        self.throat_states = self.get_throat_conditions(copy_ct_solution(self.inlet))
        self.sonic_velocity = get_speed_of_sound(self.throat_states, self.inlet_gammas)
        self.exit_states = self.get_exit_conditions()
        
    @abstractmethod
    def _get_gamma(self) -> np.array:
        pass

    @abstractmethod    
    def get_throat_conditions(self) -> ct.SolutionArray:
        pass
    
    @abstractmethod
    def _get_exit_conditions_pr(self) -> ct.SolutionArray:
        pass
    
    @abstractmethod
    def _get_exit_conditions_supr(self) -> ct.SolutionArray:
        pass
    
    @abstractmethod
    def _get_exit_conditions_subr(self) -> ct.SolutionArray:
        pass
    

class EquilibriumNozzle(NozzleBase):

    def _get_gamma(self):
        return get_thermo_properties(self.inlet)
        
    def get_throat_conditions(self) -> ct.SolutionArray:
        '''
        Assumptions: 
        - Flow is isentropic
        - Equilibrium flow
        
        Generate the throat conditions for inlet gases. 
        '''
        throat_states = ct.SolutionArray(self.inlet, self.inlet_states.shape)
        #Throat pressure 
        P_throat = self.inlet_pressures / np.power((self.inlet_gammas + 1)/2., self.inlet_gammas/(self.inlet_gammas-1))
        # CEA defaults
        max_iter_throat = 5
        tolerance_throat = 0.4e-4

        gamma_s = self.inlet_gammas
        
        M = 1.0 #throat mach is 1.0 by definition
        num_iter = 0
        residual = np.ones(throat_states.shape)
        # print(f"P throat: {P_throat}")
        while not np.all(residual < tolerance_throat):
            num_iter += 1
            if num_iter == max_iter_throat:
                # TODO: show error message
                break

            P_throat = P_throat * (1 + gamma_s * M**2)/(1 + gamma_s)
            # print(f"P_throat: {to_si(P_throat)}")
            throat_states.SPX = self.inlet_states.s, P_throat, self.inlet_states.X
            throat_states.equilibrate('SP')
            dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(throat_states)

            velocity = get_velocity(throat_states, self.inlet.h)
            sonic_velocity = get_speed_of_sound(throat_states, gamma_s)
            M = velocity/sonic_velocity

            residual = np.abs(1.0 - 1/M**2)
        
        return throat_states

    def _get_exit_conditions_supr(self) -> ct.SolutionArray:
        '''
        Calculate exit conditions, based on expansion area ratio.
        '''
        A_mdot_thr = self.throat_states.T / (self.throat_states.P * velocity * self.throat_states.mean_molecular_weight) #remains constant
        
        area_ratio = self.expansion_conditions.supersonic_ratio
        
        exit_states = ct.SolutionArray(self.inlet, (*self.inlet_states.shape, len(area_ratio)))
        
        exit_states.SPX = self.throat_states.SPX
        

        #TODO: initial guess only valid for area ratios > 2
        pressure_ratio = np.exp(gamma_s + 1.4 * np.log(self.expansion_conditions.supersonic_ratio)) #initial guess 
        p_exit = self.inlet_pressure/pressure_ratio

        exit_states.SP = self.inlet.s, p_exit
        exit_states.equilibrate('SP')

        Ae_At = exit_states.T / (exit_states.P * velocity * exit_states.mean_molecular_weight)/A_mdot_thr
        num_iter = 0
        max_iter_exit = 10
        tolerance = 4e-5

        residual = 1
        while np.abs(residual) > tolerance:
            num_iter += 1
            if num_iter == max_iter_exit:
                #TODO: error message if nozzle fails to converge
                break

            derivs = get_thermo_derivatives(exit_states)
            dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(exit_states)
            velocity = get_velocity(exit_states, self.inlet.h)
            sonic_velocity = get_speed_of_sound(exit_states, gamma_s)

            Ae_At = exit_states.T / (exit_states.P * velocity * exit_states.mean_molecular_weight) / A_mdot_thr
            dlogp_dlogA = gamma_s * velocity** 2 / (velocity**2 - sonic_velocity**2)
            residual = dlogp_dlogA * (np.log(self.expansion_conditions.supersonic_ratio) - np.log(Ae_At))
            log_pinf_pe = np.log(pressure_ratio) + residual

            pressure_ratio = np.exp(log_pinf_pe)
            p_exit = self.inlet_pressures /  pressure_ratio

            exit_states.SP = self.inlet.s, p_exit
            exit_states.equilibrate('SP')
        
        return exit_states
    
    def _get_exit_conditions_pr(self) -> ct.SolutionArray:
         # we iterate to obtain the nozzle exit temperature 
        
        T_e = self.inlet_states.T #initial guess
        pressure_ratio = self.expansion_conditions.pressure_ratio
        
        exit_shape = (*self.inlet_states.shape, len(pressure_ratio))
        exit_states = ct.SolutionArray(self.inlet, exit_shape)
        exit_states.P = self.inlet_pressure / pressure_ratio
        exit_states.T = T_e
        exit_states.equilibrate('TP')

        dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(exit_states)
        dlnT = (self.inlet.s - exit_states.s)/cp

        dlogT_residual = 0.5e-4
        maxiter = 8
        n = 0
        while not np.all(np.abs(dlnT) < dlogT_residual):
            n += 1
            if n >= maxiter:
                raise SolverError("Max iterations for exit conditions reached.")
            
            lnT_e = np.log(T_e) + dlnT
            T_e = np.exp(lnT_e)
            
            exit_states.TP = T_e, exit_states.P
            exit_states.equilibrate('TP')
            dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(exit_states)
            dlnT = (self.inlet.s - exit_states.s)/cp
        
        return exit_states
    
    def _get_exit_conditions_subr(self) -> ct.SolutionArray:
        raise NotImplementedError("Subsonic expansion is not implemented for equilibrium nozzles.")

class FrozenNozzle(NozzleBase):
    def __init__(self, inlet_gas: ct.Solution, exit_conditions: ExpansionConditions, gamma_s: float | None = None, NFZ: int = 1) -> None:
        super().__init__(inlet_gas, exit_conditions, gamma_s)
        if NFZ < 1 or NFZ > 5 or not isinstance(NFZ, int):
            raise ValueError("Invalid NFZ value: must be integer between 1 and 5. See NASA RP1311 Part II (Users Manual), p. 18-19.")
        self.NFZ = NFZ
    
    def get_throat_conditions(self, throat_states: ct.Solution) -> ct.Solution:
        """
        Throat conditions for frozen flow. See Gordon & McBride, 1994, Section 6.5.2.
        """
        throat_states = ct.SolutionArray(self.inlet, self.inlet_states.shape)
        P_throat = self.inlet_pressures / np.power((self.inlet_gammas + 1)/2., self.inlet_gammas/(self.inlet_gammas-1))
        # CEA defaults
        max_iter_throat = 5
        tolerance_throat = 0.4e-4

        gamma_s = self.inlet_gammas
        
        M = 1.0 # throat mach is 1.0 by definition
        num_iter = 0
        residual = np.ones(throat_states.shape)
        while not np.all(residual < tolerance_throat):
            num_iter += 1
            if num_iter == max_iter_throat:
                # TODO: show error message
                break

            P_throat = P_throat * (1 + gamma_s * M**2)/(1 + gamma_s)
            throat_states.SPX = self.inlet_states.s, P_throat, self.inlet_states.X
            
            #if NFZ = 1, then the frozen composition is at the combustion point.
            #otherwise frozen composition is at the throat, or downstream of it. 
            if self.NFZ > 1: 
                throat_states.equilibrate('SP')
                dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(throat_states)


            velocity = get_velocity(throat_states, self.inlet_states.h)
            sonic_velocity = get_speed_of_sound(throat_states, gamma_s)
            M = velocity/sonic_velocity

            residual = np.abs(1.0 - 1/M**2)
        
        return throat_states
    
    def _get_exit_conditions_supr(self) -> ct.SolutionArray:
        """
        See Gordon & McBride, 1994, Section 6.5.
        """
        area_ratio = self.expansion_conditions.supersonic_ratio
        A_mdot_thr = self.throat_states.T / (self.throat_states.P * velocity * self.throat_states.mean_molecular_weight) #remains constant
        exit_shape = (*self.inlet_states.shape, len(area_ratio)) #we take all of the inlet states x the number of area ratios to calculate
        # we iterate to obtain the nozzle exit temperature 
        T_e = self.throat_states.T #initial guess
        exit_states = ct.SolutionArray(self.inlet, exit_shape)

        gamma_s = exit_states.cp/exit_states.cv
        
        pressure_ratio = np.exp(gamma_s + 1.4 * np.log(area_ratio)) #initial guess 

        p_exit = self.inlet_pressure/pressure_ratio #TODO: probably mismatched dimensions
        # print(f"P_exit: {p_exit}")
        exit_states.TPX = T_e, p_exit, self.throat_states.X

        gamma_s = exit_states.cp/exit_states.cv

        dlnT = (self.inlet.s - exit_states.s)/exit_states.cp

        dlogT_tolerance = 0.5e-4
        pressure_tolerance = 0.4e-4
        maxiter = 10
        maxiter_temp = 8
        n = 0

        residual = np.ones(exit_shape)
        while not np.all(np.abs(residual) < pressure_tolerance):
            n += 1
            if n == maxiter:
                #TODO: error message if nozzle fails to converge
                break
            # print(f"Iteration: {n}")
            m = 0
            dlnT = np.ones(exit_states.shape) #arbitrary initial value
            while not np.all(np.abs(dlnT) >= dlogT_tolerance):
                m += 1
                if m >= maxiter_temp:
                    break
                exit_states.TPX = T_e, p_exit, exit_states.X
                dlnT = (self.inlet.s - exit_states.s)/exit_states.cp
                T_e = np.exp(np.log(T_e) + dlnT)
                # print(f"T_exit: {T_e}")

            gamma_s = exit_states.cp/exit_states.cv #gamma_s = gamma for frozen flow
            velocity = get_velocity(exit_states, self.inlet.h)
            sonic_velocity = get_speed_of_sound(exit_states, gamma_s)

            Ae_At = exit_states.T / (exit_states.P * velocity * exit_states.mean_molecular_weight) / A_mdot_thr
            dlogp_dlogA = gamma_s * velocity**2 / (velocity**2 - sonic_velocity**2)
            residual = dlogp_dlogA * (np.log(area_ratio) - np.log(Ae_At))
            log_pinf_pe = np.log(pressure_ratio) + residual


            pressure_ratio = np.exp(log_pinf_pe)
            p_exit = self.inlet_pressures /  pressure_ratio

        print(f"number of iterations: {n}")

        return exit_states
    
    def _get_exit_conditions_pr(self) -> ct.SolutionArray:
        # we iterate to obtain the nozzle exit temperature 
        
        T_e = self.inlet_states.T #initial guess
        pressure_ratio = self.expansion_conditions.pressure_ratio
        
        exit_shape = (*self.inlet_states.shape, len(pressure_ratio))
        exit_states = ct.SolutionArray(self.inlet, exit_shape)
        exit_states.P = self.inlet_pressure / pressure_ratio

        # dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(exit_states)
        dlnT = (self.inlet.s - exit_states.s)/exit_states.cp

        dlogT_residual = 0.5e-4
        maxiter = 8
        n = 0
        while not np.all(np.abs(dlnT) < dlogT_residual):
            n += 1
            if n >= maxiter:
                break
            
            lnT_e = np.log(T_e) + dlnT
            T_e = np.exp(lnT_e)
            
            exit_states.TP = T_e, exit_states.P
            dlnT = (self.inlet.s - exit_states.s)/exit_states.cp
        
        return exit_states
    
    def _get_exit_conditions_subr(self) -> ct.SolutionArray:
        raise NotImplementedError("Subsonic expansion is not implemented for frozen nozzles.")

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
            raise SolverError("Max iterations for exit conditions reached.")
            break

        derivs = get_thermo_derivatives(gas_exit)
        dlogV_dlogT_P, dlogV_dlogP_T, cp, gamma_s = get_thermo_properties(gas_exit)
        velocity = get_velocity(gas_exit, gas_chamber.h)
        sonic_velocity = get_speed_of_sound(gas_exit, gamma_s)

        Ae_At = gas_exit.T / (gas_exit.P * velocity * gas_exit.mean_molecular_weight) / A_mdot_thr
        dlogp_dlogA = gamma_s * velocity** 2 / (velocity**2 - sonic_velocity**2)
        residual = dlogp_dlogA * (np.log(area_ratio) - np.log(Ae_At))
        log_pinf_pe = np.log(pressure_ratio) + residual

        pressure_ratio = np.exp(log_pinf_pe)
        p_exit = to_si(P_chamber) /  pressure_ratio

        gas_exit.SP = gas_chamber.s, p_exit
        gas_exit.equilibrate('SP')
    
    return gas_exit


def _get_nozzle_constructor(nozzle_class: Type[NozzleBase], **kwargs):
    return lambda inlet_gas, exit_conditions, gamma=None: nozzle_class(inlet_gas, exit_conditions, gamma, **kwargs)

# API functions
def frozen_nozzle(NFZ: int = 1):
    return _get_nozzle_constructor(FrozenNozzle, NFZ=NFZ)

def equilibrium_nozzle():
    return _get_nozzle_constructor(EquilibriumNozzle)

def kinetic_nozzle():
    raise NotImplementedError("Kinetic nozzles are not implemented.")

@dataclass
class NozzleOutputs:
    """
    Organizes all the calculations from the nozzle calculations not directly 
    contained in the Cantera objects (cstar, isp, ivac, gamma) into a convenient 
    helper class.
    """
    isp: np.ndarray
    ivac: np.ndarray
    cstar: np.ndarray
    gamma: np.ndarray
    
    def append(self, isp, ivac, cstar, gamma) -> None:
        """Append 

        Args:
            isp (bool): _description_
            ivac (_type_): _description_
            cstar (_type_): _description_
            gamma (_type_): _description_
        """
        np.append(self.isp, isp)
        np.append(self.ivac, ivac)
        np.append(self.cstar, cstar)
        np.append(self.gamma, gamma)