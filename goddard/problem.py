from abc import ABC, abstractmethod
from pint import Quantity
from dataclasses import dataclass
import numpy as np

import cantera as ct

from . import utils
from .combustion import CombustorBase
from .nozzle import NozzleBase, ExpansionConditions
from .result import Result

class ProblemBase(ABC):
    def __init__(self, dfiles = None, **kwargs) -> None:
        self.files = dfiles
        self.info = kwargs
        self.inputs = None
    
    @abstractmethod
    def report(self) -> None:
        pass
    
    @abstractmethod
    def solve(self):
        pass
    
    @abstractmethod
    def _generate_input_dict(self):
        pass

class CombustionProblem(ProblemBase):
    """
    
    """
    def __init__(self, fuel: str | ct.Solution, oxidizer: str | ct.Solution, pressure: float | Quantity, mixture_ratio, exit_conditions: ExpansionConditions, combustor, **kwargs) -> None:
        super().__init__(**kwargs)
        pressure_SI = map(utils.normalize_input_units_to_si, pressure)
        pressure_array = np.fromiter(pressure_SI)
        #TODO: accept cantera objects directly as well?
        self.fuel = utils.generate_ct_solution_from_text_input(fuel, pressure_array[0]) 
        self.oxidizer = utils.generate_ct_solution_from_text_input(oxidizer, pressure_array[0])
        self.mixture_ratio = mixture_ratio(self.fuel, self.oxidizer)
        self.exit_conditions = exit_conditions
        self.combustor = combustor(self.fuel, self.oxidizer, self.mixture_ratio)
        
    def solve(self):
        outputs = self.combustor.solve()
        #possibly calculate additional parameters before returning values
        return super().solve()

class NozzleProblem(ProblemBase):
    def __init__(self, gas, nozzle, pressure, expansion_conditions, **kwargs) -> None:
        super().__init__(**kwargs)
        
        self.inlet = utils.copy_ct_solution(gas)
        self.inlet_states = ct.SolutionArray(self.inlet, len(pressure))
        self.inlet_states.P = pressure
        #TODO: accept nozzle objects directly rather than just a nozzle factory function
        self._nozzle_builder = nozzle
        self.expansion_conditions = expansion_conditions
    
    def solve(self):
        nozzle = self._nozzle_builder(self.inlet, self.inlet_states, self.expansion_conditions)
        output_states = nozzle.get_exit_conditions()
        return Result(output_states)
        
class RocketProblem(CombustionProblem):
    def __init__(self, 
                 fuel: str | ct.Solution, 
                 oxidizer: str | ct.Solution, 
                 pressure: float | Quantity, 
                 mixture_ratio, 
                 exit_conditions: ExpansionConditions, 
                 combustor, 
                 nozzle, 
                 **kwargs) -> None:
        super().__init__(fuel, oxidizer, pressure, mixture_ratio, exit_conditions, combustor, **kwargs)
        self._nozzle_builder = nozzle
        
    def solve(self) -> Result:
        combustor_states = self.combustor.solve()
        self.nozzle = self._nozzle_builder(self.combustor.reactants, combustor_states, self.exit_conditions)
        output = self.nozzle.get_exit_conditions()
        return Result(output)
        
class KineticProblem(ProblemBase):
    def __init__(self, **kwargs) -> None:
        raise NotImplementedError("KineticProblem is not implemented.")
