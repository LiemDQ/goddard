from abc import ABC, abstractmethod
from pint import Quantity
from dataclasses import dataclass

import cantera as ct


from . import utils
from .combustion import Combustor
from .nozzle import Nozzle, ExitConditions
from .result import Result

class Problem:
    
    def __init__(self, **kwargs) -> None:
        self.info = kwargs
    
    def report() -> None:
        pass
    
    @abstractmethod
    def solve():
        pass

class CombustionProblem(Problem):
    pass

class NozzleProblem(Problem):
    def __init__(self, gas, nozzle, pressure, exit_conditions, **kwargs) -> None:
        super().__init__(**kwargs)
        
        self.inlet = utils.copy_ct_solution(gas)
        self.inlet.P = pressure
        self.nozzle = nozzle(gas, exit_conditions)
        
        
class RocketProblem(Problem):
    def __init__(self, fuel: str, oxidizer: str, pressure: float | Quantity, mixture_ratio, exit_conditions: ExitConditions, combustor, nozzle, **kwargs) -> None:
        super().__init__(**kwargs)
        pressure_SI = utils.normalize_input_units_to_si(pressure)
        #TODO: accept cantera objects directly as well?
        self.fuel = utils.generate_ct_solution_from_text_input(fuel, pressure_SI) 
        self.oxidizer = utils.generate_ct_solution_from_text_input(oxidizer, pressure_SI)
        self.mixture_ratio = mixture_ratio(self.fuel, self.oxidizer)
        self.exit_conditions = exit_conditions
        self.combustor = combustor(self.fuel, self.oxidizer, self.mixture_ratio)
        self._nozzle_builder = nozzle
        
    def solve(self) -> Result:
        combustor_gas = self.combustor.solve()
        self.nozzle = self._nozzle_builder(combustor_gas, self.exit_conditions)
        output = self.nozzle.get_exit_conditions()
        
class KineticProblem(Problem):
    pass

class Solution:
    pass

@dataclass
class SolverOptions:
    ionized: bool = False
    transport: bool = False
    use_mole_frac: bool = False
    
def solve(problem: Problem, options: SolverOptions):
    pass
