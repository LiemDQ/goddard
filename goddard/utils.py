import pint
import cantera as ct
import numpy as np

#Earth gravitational acceleration
g0 = 9.8067 
#TODO: determine if additional finangling is needed to get the file path to play nicely
reactant_files = ["../data/nasa_reactants.yaml", "nasa_gas.yaml", "nasa_condensed.yaml"] 
product_files = ["nasa_gas.yaml"]

def to_si(quant: pint.Quantity):
    '''
    Convert pint Quantity to magnitude in base SI units.
    '''
    return quant.to_base_units().magnitude

def species_list_to_dict(species_list):
    return {species.name:species for species in species_list}

def process_text_input(text: str):
    """
    String format: CH3OH(L):0.4|C2H5OH(L):0.6|mass

    """
    input_species = text.split("|")
    composition_type = "mass"
    last = input_species[-1]
    if last == "mass" or last == "mole":
        composition_type = last
        input_species.pop()
    species_dict = {}
    for species in input_species:
        s = species.split(":")
        name = s[0]
        if len(s) > 2:
            raise ValueError("Encountered multiple ':' following a species name. Is the input string correctly formatted?")
        elif len(s) == 2:
            composition = s[1]
        else:
            composition = None
        species_dict[name] = composition
    
    return composition_type, species_dict
    

def normalize_compositions(species_dict: dict[str, float]):
    total = sum(species_dict.values())
    species_dict.update((k, v/total) for k,v in species_dict.items())
    #TODO: handle cases where some of the compositions are None
    return species_dict

def extract_reaction_species(
    fuel: ct.Solution, 
    oxidizer: ct.Solution, 
    filenames: str = product_files
    ) -> ct.Solution:
    '''
    Extract candidate reactants, based on the elements contained in the fuel and oxidizer.
    
    '''
    elements = [*fuel.element_names, *oxidizer.element_names]
    full_species = []
    for filename in filenames:
      full_species += ct.Species.list_from_file(filename)

    species = [S for S in full_species if all(x in elements for x in S.composition)]
    
    return ct.Solution(thermo='IdealGas', species=species)

def copy_solution(solution: ct.Solution) -> ct.Solution:
    result = ct.Solution(thermo=solution.thermo_model, species=solution.species())
    result.SPX = solution.SPX
    return result
