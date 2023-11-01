import pint
import cantera as ct
import numpy as np
import itertools

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

def normalize_input_units_to_si(value: float | pint.Quantity) -> float:
    """
    Convert input value to base SI units. This is used to ensure all downstream calculations are in SI base units.
    """
    if isinstance(value, float):
        return value
    else:
        return value.to_base_units().magnitude

def species_list_to_dict(species_list):
    return {species.name:species for species in species_list}

def process_text_input(text: str) :
    """
    String format: species name:basis fraction|species 2 name:basis fraction|...|TX (temperature X)|basis (mass or molar)
    String format: CH3OH(L):0.4|C2H5OH(L):0.6|T293.0|mass
    """
    input_species = text.split("|")
    composition_type = "mass"
    temperature = 273.15 #default is 0 degC
    if len(input_species) > 1: #potentially extract temperature and composition type data
        last = input_species[-1]
        second_last = input_species[-2]

        if last == "mass" or last == "molar": #case 1: composition type is specified
            composition_type = last
            input_species.pop()
            if second_last[0] == "T" and second_last[1].isdigit():
                temperature = float(second_last[1:])
                input_species.pop()
        elif last[0] == "T" and last[1].isdigit(): # case 2: composition type is not specified, temperature is specified
            temperature = float(last[1:])
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
            if len(input_species) > 1:
                raise ValueError(f"Composition of multi-species mixture was not specified properly. Expected a value for {name}, but none was found.")
            composition = 1.0 #if composition is not specified, then it is the balance of the remaining compositions
        species_dict[name] = composition
    
    return species_dict, temperature, composition_type
    

def normalize_compositions(species: dict[str, float]):
    total = sum(species.values())
    species.update((k, v/total) for k,v in species.items())
    
    return species

def generate_ct_solution(input_species_dict: dict[str, float], temperature: float, pressure: float, composition_type: str, reactant_file = reactant_files[0]) -> ct.Solution:
    species_dict = species_list_to_dict(ct.Species.list_from_file(reactant_file))
    
    input_species = [species_dict[name] for name in input_species_dict.keys()]

    solution = ct.Solution(thermo="IdealGas", species = input_species)
    #Could add optional function parameter to support non-ideal gas models
    #TODO: add check on the valid temperature range for the thermodynamic correlations, and give an error if the 
    #temperature is out of bounds
    if composition_type == "mass":
        solution.TPY = temperature, pressure, input_species_dict
    elif composition_type == "mole":
        solution.TPX = temperature, pressure, input_species_dict
    else:
        raise ValueError("Composition type must be one of either 'mass' or 'mole'.")
    
    
    return solution

def generate_ct_solution_from_text_input(text: str, pressure: float) -> ct.Solution:
    species, temperature, comp_type = process_text_input(text)
    species = normalize_compositions(species)

    return generate_ct_solution(species, temperature, pressure, comp_type)

def set_species_file(filename: str):
    reactant_files.append(filename)

def extract_reaction_species(fuel: ct.Solution, oxidizer: ct.Solution, filenames: str = product_files) -> ct.Solution:
    '''
    Extract candidate reactants, based on the elements contained in the fuel and oxidizer.
    
    '''
    elements = [*fuel.element_names, *oxidizer.element_names]
    full_species = []
    for filename in filenames:
      full_species += ct.Species.list_from_file(filename)

    species = [S for S in full_species if all(x in elements for x in S.composition)]
    
    return ct.Solution(thermo='IdealGas', species=species)

def copy_ct_solution(solution: ct.Solution) -> ct.Solution:
    result = ct.Solution(thermo=solution.thermo_model, species=solution.species())
    result.SPX = solution.SPX
    return result

def mixed_flatten(iterable):
    """Flattens an iterable containing mixed iterables and scalars."""
    for item in iterable:
        try:
            iterator = iter(item)
        except TypeError: #item is scalar
            yield item
        else: #item is an iterable
            for subitem in item: yield subitem

def args_to_np_array(value, *args) -> np.array:
    """
    Convert arguments into a numpy array for further numerical processing.
    
    Cases that can be handled:
    1. Value is a scalar, args are a series of scalars
    2. Value is an iterable, args are a series of scalars
    3. Value is an iterable, args are a series of iterables
    4. Value is a numpy array, args are a series of scalars
    5. Value is a numpy array, args are a series of iterables
    """                
    res = np.fromiter(mixed_flatten(args))
    
    if isinstance(value, np.array):
        return np.concatenate((value, res))
    else:
        try:
            iterator = iter(value)
        except TypeError: #value is a scalar
            return np.concatenate((np.array([value]), res))
        else:
            return np.concatenate((np.fromiter(value), res))
            