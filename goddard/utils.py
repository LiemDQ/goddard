import pint
import cantera as ct

#Earth gravitational acceleration
g0 = 9.8067 

def to_si(quant: pint.Quantity):
    '''
    Convert pint Quantity to magnitude in base SI units.
    '''
    return quant.to_base_units().magnitude

def extract_reaction_species(
    fuel: ct.Solution, 
    oxidizer: ct.Solution, 
    filename: str = "nasa_gas.yaml"
    ) -> ct.Solution:
    '''
    Extract candidate reactants, based on the elements contained in the fuel and oxidizer.
    '''
    elements = [*fuel.element_names, *oxidizer.element_names]
    full_species = ct.Species.list_from_file(filename)

    species = [S for S in full_species if all(x in elements for x in S.composition)]
    
    return ct.Solution(thermo='IdealGas', species=species)

def copy_solution(solution: ct.Solution) -> ct.Solution:
    result = ct.Solution(thermo=solution.thermo_model, species=solution.species())
    result.SPX = solution.SPX
    return result
