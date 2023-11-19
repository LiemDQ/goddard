import numpy as np

class Inputs:
    """Class representing a set of numerical inputs to a problem, provided in the order they are given.
    """
    def __init__(self, vals: list[np.ndarray] = [], names: list[str] = [] ) -> None:
        """_summary_

        Args:
            vals (list[float], optional): _description_. Defaults to [].
            names (list[str], optional): _description_. Defaults to [].
        """
        self.vals = vals
        self.names = names
        assert len(self.vals) == len(self.names), "There must be the same number of names as input conditions"
    
    def append(self, val: np.array, name: str):
        self.names.append(name)
        self.vals.append(val)
        assert len(self.vals) == len(self.names), "There must be the same number of names as input conditions"
        
    def __getitem__(self, idx: int) -> (str, np.ndarray):
        return self.name[idx], self.vals[idx]