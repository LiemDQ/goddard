from typing import Any
import cantera as ct
import pandas as pd

from problem import ProblemBase
from inputs import Inputs

STANDARD_OUTPUT_COLUMNS = ["P", "T", "s", "D", "u", "g", "h", "cp_mass", "Y", "mean_molecular_weight"]

def ct_SolutionArray_to_df(slnarr: ct.SolutionArray, inputs: Inputs, cols: list[str] = STANDARD_OUTPUT_COLUMNS) -> pd.DataFrame:
    """Converts a cantera `SolutionArray` to a pandas `DataFrame`, with indices based
    on the `Inputs` object. Works for up to 3D arrays (unlike the native cantera functions, 
    which only work for 1D arrays).
    

    Args:
        slnarr (ct.SolutionArray): The `SolutionArray` to be converted to a `DataFrame`. 
        inputs (Inputs): An input object representing the inputs to the system. This is used to create indexes for the results. 
        cols (list[str]): a list of strings indicating columns in the output `DataFrame`. These strings must correspond to properties
            from a Cantera `Thermophase` object. See the official Cantera documentation for more details. 

    Raises:
        ValueError: _description_

    Returns:
        pd.DataFrame: `DataFrame` containing the properties of the `SolutionArray` listed above, with the input
        values used as indices. 
    """
    
    
    #TODO: make sure this works as expected in the 1D case
    idx = pd.Multiindex.from_product(inputs.vals, names=inputs.names)
    df = pd.DataFrame(index=idx)
    
    if slnarr.ndim == 1:
        df = slnarr.to_pandas(cols=cols)
        
    elif slnarr.ndim == 2:
        nrows = slnarr.shape[0] 
        
        for i in range(nrows):
            slndf = slnarr[i, :].to_pandas(cols=cols)
            slndf[inputs.names[0]] = inputs.vals[0]
            slndf[inputs.names[1]] = inputs.vals[1][i] #second index val will be the same for all entries in the row
            df = pd.concat([df, slndf])
    
    elif slnarr.ndim == 3:
        nrows,ncols,_ = slnarr.shape
        idx = pd.MultiIndex.from_product(inputs.vals, names=inputs.names)
        df = pd.DataFrame(index=idx)
        for i in range(nrows):
            for j in range(ncols):
                #TODO: not sure if this iteration ordering is correct
                slndf = slnarr[j, i, :].to_pandas(cols=cols)
                slndf[inputs.names[0]] = inputs.vals[0]
                slndf[inputs.names[1]] = inputs.vals[1][j]
                slndf[inputs.names[2]] = inputs.vals[2][i]
                df = pd.concat([df, slndf])
    
    else:
        raise ValueError("More than 3 dimensions is not supported for SolutionArrays.")
    
    return df


class Result:
    
    df: pd.DataFrame
    def __init__(self, problem: ProblemBase, output_states: ct.SolutionArray) -> None:
        self.output_solutions = output_states
        self.problem = problem
        self.inputs = problem.inputs
        self.dfs = [ct_SolutionArray_to_df(output_states, problem.inputs)]    
    
    def __getitem__(self, idx):
        pass
    
    def __call__(self, *args: Any, **kwds: Any) -> Any:
        pass
    
    def report(self):
        print(self.df)