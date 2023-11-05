
import cantera as ct

reactant_files = ["../data/nasa_reactants.yaml", "nasa_gas.yaml", "nasa_condensed.yaml"] 
product_files = ["nasa_gas.yaml"]

def add_source(source, priority=9):
    #TODO: implement adding data sources
    pass

def set_species_file(filename: str):
    reactant_files.append(filename)
    

def species_list_to_dict(species_list):
    return {species.name:species for species in species_list}

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

def generate_ct_solution_from_file(input_species_dict: dict[str, float], temperature: float, pressure: float, composition_type: str, reactant_file = reactant_files[0]) -> ct.Solution:
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

species_aliases = {
    'kerosene':'Jet-A',
    'RP-1':'Jet-A',
    'RP-1(L)':'Jet-A(L)',
    #TODO: implement paraffin
}