import numpy as np

class Inputs:
    def __init__(self, vals: list[float] = [], names: list[str] = [] ) -> None:
        self.vals = vals
        self.names = names
        assert len(self.vals) == len(self.names), "There must be the same number of names as input conditions"
    
    def append(self, val: np.array, name: str):
        self.names.append(name)
        self.vals.append(val)
        assert len(self.vals) == len(self.names), "There must be the same number of names as input conditions"
        
    def __getitem__(self, idx: int) -> (str, np.array):
        return self.name[idx], self.vals[idx]