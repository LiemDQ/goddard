from abc import ABC, abstractmethod
from combustion import 
class Problem:
    def __init__(self, species, *args, **kwargs) -> None:
        self.species = species
        self.info = kwargs
    
    def report() -> None:
        pass


class Solver(ABC):
    def __init__(self) -> None:
        super().__init__()
    
    @abstractmethod
    def report():
        pass

    @abstractmethod
    def solve(problem: Problem):
        pass


class CombustionSolver(Solver):
    def solve(problem: Problem):
        
        return 

class NozzleSolver(Solver):
    pass

class RocketSolver(CombustionSolver):
    def solve(problem: Problem):
        chamber_conditions = super().solve(problem)
        return None

class KineticSolver(Solver):
    """
    Not implemented yet!
    """
    pass

class Solution:
    pass


def solve(problem: Problem, solver: Solver) -> Solution:
    pass
