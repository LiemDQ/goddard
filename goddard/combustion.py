import utils
import thermo
import nozzle as nz
import data
from enum import Enum, auto
from dataclasses import dataclass
import cantera as ct
import numpy as np
from abc import ABC, abstractmethod

class MixtureRatio:
    def __init__(self, M_fuel, M_oxidizer, OF, phi, equiv, fuel_percent) -> None:
        self.OF_ratio = OF
        self.phi = phi
        self.req_ratio = equiv
        self.fuel_percent = fuel_percent
        self.M_fuel = M_fuel
        self.M_ox = M_oxidizer
        self.molar_ratio = self._OF_to_molar_ratio(OF)
        self.size = len(self.OF_ratio)
    
    @property
    def M_fuel(self):
        return self._M_fuel
    
    @property
    def M_ox(self):
        return self._M_ox
    
    @property
    def phi(self):
        return self._phi
    
    @property
    def req_ratio(self):
        return self._req_ratio
    
    @property
    def OF_ratio(self):
        return self._OF_ratio
    
    @property
    def fuel_percent(self):
        return self._fuel_percent
    
    @property.setter
    def phi(self):
        pass 
    #TODO: I am not sure if it is even possible to convert between all these different ratios without stoichiometric data
    
    def _OF_to_molar_ratio(self, OF):
        return OF / (self.M_ox / self.M_fuel)
    
    def fuel_mole_frac(self):
        return 1 - self.molar_ratio/(1+self.molar_ratio)
    
    def ox_mole_frac(self):
        return self.molar_ratio/(1+self.molar_ratio)

def _create_mixture_ratio(**kwargs):
    return lambda fuel, oxidizer: MixtureRatio(fuel.mean_molecular_weight, oxidizer.mean_molecular_weight)

def OF_ratio(value, *args):
    values = utils.args_to_np_array(value, args)
    return _create_mixture_ratio(OF=values)

def phi_ratio(value, *args):
    values = utils.args_to_np_array(value, args)
    return _create_mixture_ratio(phi=values)

def equiv_ratio(value, *args):
    values = utils.args_to_np_array(value, args)
    return _create_mixture_ratio(equiv=values)

def fuel_pct(value, *args):
    values = utils.args_to_np_array(value, args)
    return _create_mixture_ratio(fuel_percent=values)


class CombustorBase(ABC):
    def __init__(self, fuel: ct.Solution, oxidizer: ct.Solution, pressures: np.array, mr: MixtureRatio) -> None:
        super().__init__()
        self.fuel = fuel
        self._fuel_output = utils.copy_ct_solution(fuel)
        self.oxidizer = oxidizer
        self._oxidizer_output = utils.copy_ct_solution(oxidizer)
        self.reactants = data.extract_reaction_species(self._fuel_output, self._oxidizer_output)
        self.pressures = pressures
        self.mixture_ratio = mr
        
    @abstractmethod
    def solve(self) -> ct.SolutionArray:
        pass
    
    def _generate_mole_frac_matrix(self) -> np.array:
        """
        Create a 2D array where each entry is the set of mole fractions reflecting a given mixture ratio.
        Used as an input to SolutionArray.
        """
        remove_L = lambda s: s.replace('(L)', '')
        
        mole_ratios = self.mixture_ratio.OF_ratio / (self.oxidizer.mean_molecular_weight / self.fuel.mean_molecular_weight)
        moles_ox = mole_ratios / ( 1 + mole_ratios)
        moles_f = 1 - moles_ox
        
        mole_frac_matrix = np.zeros((len(self.pressures), len(self.reactants.species_names)))
        
        for i,(f, o) in enumerate(zip(moles_f, moles_ox)):
            fidx = [self.reactants.species_index(remove_L(species)) for species in self.fuel.species_names]
            oidx = [self.reactants.species_index(remove_L(species)) for species in self.oxidizer.species_names]
            for idx in fidx:
                mole_frac_matrix[i][idx] = f
            for idx in oidx:
                mole_frac_matrix[i][idx] = o
        
        return mole_frac_matrix[:, np.newaxis, :]
        
class FiniteAreaCombustor(CombustorBase):
    def __init__(self, fuel: ct.Solution, oxidizer: ct.Solution, mr: MixtureRatio, mass_flux=None, contraction_ratio = None) -> None:
        super().__init__(fuel, oxidizer, mr)
        

class InfiniteAreaCombustor(CombustorBase):
    def solve(self):
        combustion_states = ct.SolutionArray(self.reactants, (len(self.pressures), self.mixture_ratio.size))
        mole_fracs = self._generate_mole_frac_matrix()
        combustion_states.TPX = self.reactants.T, self.pressures, mole_fracs
         
        # solve for combustion chamber composition
        combustion_states.equilibrate('HP',solver="gibbs")
        return combustion_states

def finite_area_combustor(mass_flux=None, contraction_ratio=None):
    pass

def infinite_area_combustor():
    return lambda fuel, oxidizer, pressures, mr: InfiniteAreaCombustor(fuel, oxidizer, pressures, mr)