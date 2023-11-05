from abc import ABC, abstractmethod
from pint import Quantity
from dataclasses import dataclass
import numpy as np

import cantera as ct

from . import utils
from .combustion import Combustor
from .nozzle import Nozzle, ExpansionConditions
from .result import Result

class Problem:
    def __init__(self, data_files = None, **kwargs) -> None:
        self.files = data_files
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
    def __init__(self, fuel: str | ct.Solution, oxidizer: str | ct.Solution, pressure: float | Quantity, mixture_ratio, exit_conditions: ExpansionConditions, combustor, nozzle, **kwargs) -> None:
        super().__init__(**kwargs)
        pressure_SI = map(utils.normalize_input_units_to_si, pressure)
        pressure_array = np.fromiter(pressure_SI)
        #TODO: accept cantera objects directly as well?
        self.fuel = utils.generate_ct_solution_from_text_input(fuel, pressure_array[0]) 
        self.oxidizer = utils.generate_ct_solution_from_text_input(oxidizer, pressure_array[0])
        self.mixture_ratio = mixture_ratio(self.fuel, self.oxidizer)
        self.exit_conditions = exit_conditions
        self.combustor = combustor(self.fuel, self.oxidizer, self.mixture_ratio)
        self._nozzle_builder = nozzle
        
    def solve(self) -> Result:
        combustor_states = self.combustor.solve()
        self.nozzle = self._nozzle_builder(self.combustor.reactants, combustor_states, self.exit_conditions)
        output = self.nozzle.get_exit_conditions()
        return Result(output)
        
class KineticProblem(Problem):
    pass

@dataclass
class SolverOptions:
    ionized: bool = False
    transport: bool = False
    mole_basis: bool = False
    verbose: bool = False
    debug_output: bool = False
    
def solve(problem: Problem, options: SolverOptions = SolverOptions()):
    pass
