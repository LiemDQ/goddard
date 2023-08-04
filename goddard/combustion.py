import utils
import thermo
import nozzle as nz
from enum import Enum, auto
from dataclasses import dataclass
import cantera as ct
import numpy as np
from collections.abc import Iterable

class NozzleType(Enum):
    EQ = auto()
    FROZEN = auto()
    KINETIC = auto()

class CombustorArgs:
    def __init__(self, CR = 0, flux_ratio = 0, use_CR = True, nozzle_type = NozzleType.EQ):
        """
        CR: contraction ratio
        """
        self.CR = CR
        self.flux_ratio = flux_ratio
        self.use_CR = use_CR
        self.nozzle_type = nozzle_type
    
    def __repr__(self) -> str:
        #TODO proper implementation of repr
        return "Combustor"

class CombustionResult:
    def __init__(
        self,
        nozzle_type: NozzleType,
        of_ratio,
        cstar,
        fuel,
        oxidizer,
        conditions,
        ):
        self.nozzle_type = nozzle_type
        self.of_ratio = of_ratio
        self.cstar = cstar
        self.equiv_ratio = 0
        self.phi_ratio = 0
        self.fuel = fuel
        self.oxidizer = oxidizer
        self.nozzle_conditions = conditions

    def full_output(self):
        output = ""
        # section 1: inputs
        return output

    def get_isp(self):
        return [point.isp for point in self.nozzle_conditions]
    
    def get_ivac(self):
        return [point.ivac for point in self.nozzle_conditions]

@dataclass
class CombustionPoint:
    M: float
    CF: float
    isp: float
    ivac: float
    expansion_ratio: float
    pressure_ratio: float
    exhaust_gas: ct.Solution

class NozzleConditions:
    pass

class CombustionAnalysis:
    def __init__(
        self,
        Pcs,
        fuel: ct.Solution,
        oxidizer: ct.Solution,
        MRs,
        expRs,
        combustor: CombustorArgs
        ):
        self.fuel = utils.copy_solution(fuel) #avoids mutating state outside of the analysis object
        self.oxidizer = utils.copy_solution(oxidizer)
        

        self.chamber_pressure = np.array(Pcs)
        self.mixture_ratio = np.array(MRs)
        self.exit_conditions = np.array(expRs)

        self.combustor = combustor
    
    def run(self) -> CombustionResult:
        
        #normalize to molar basis
        #TODO: this assumes that mixture ratio is a single value as opposed to an iterable.
        molar_ratio = self.mixture_ratio / (self.oxidizer.mean_molecular_weight / self.fuel.mean_molecular_weight)
        print(f"Molar ratio: {molar_ratio}")
        moles_ox = molar_ratio / (1 + molar_ratio)
        print(f"Moles ox: {moles_ox}")
        moles_f = 1 - moles_ox

        # Generate results matrix
        # nesting order: Pc -> MR -> ExpR

        #
        gas_chamber = utils.extract_reaction_species(self.fuel, self.oxidizer)

        mixture = ct.Mixture([(self.fuel, moles_f), (self.oxidizer, moles_ox), (gas_chamber, 0.0)])
        
        # solve for combustion chamber composition
        mixture.equilibrate('HP',solver="gibbs")
        # print("CHAMBER CONDITIONS: ")
        # gas_chamber()

        _,_,_,gamma = thermo.get_thermo_properties(gas_chamber)
        cstar = nz.get_cstar(gamma, gas_chamber.T, gas_chamber.mean_molecular_weight)
        
        # print("THROAT CONDITIONS: ")
        if self.combustor.nozzle_type == NozzleType.FROZEN:
            gas_throat = utils.copy_solution(gas_chamber)
            gas_throat = nz.get_throat_conditions_frozen(gas_throat, self.chamber_pressure, gas_chamber, gamma)
            # gas_throat()
            sonic_velocity = thermo.speed_of_sound(gas_throat, gamma)


            gas_exit = utils.copy_solution(gas_throat)
            gas_exit = nz.get_exit_conditions_frozen(gas_exit, self.exit_conditions, self.chamber_pressure, gas_chamber, gamma, sonic_velocity)
            # gas_exit()
            
        elif self.combustor.nozzle_type == NozzleType.EQ:
            gas_throat = utils.copy_solution(gas_chamber)
            gas_throat = nz.get_throat_conditions(gas_throat, self.chamber_pressure, gas_chamber, gamma)
            # gas_throat()
            sonic_velocity = thermo.speed_of_sound(gas_throat, gamma)

            gas_exit = utils.copy_solution(gas_throat)
            gas_exit = nz.get_exit_conditions(gas_exit, self.exit_conditions, self.chamber_pressure, gas_chamber, gamma, sonic_velocity)
            # gas_exit()
        else:
            raise NotImplementedError
        isp = nz.get_isp(gas_exit, gas_chamber.h)
        
        isp_vac = nz.get_ivac(gas_exit, isp[0])
        
        print("EXIT CONDITIONS: ")
        gas_exit()

        point = CombustionPoint(
            0,
            0,
            isp,
            isp_vac,
            self.exit_conditions,
            0,
            gas_exit
        )

        return CombustionResult(
            NozzleType.EQ,
            self.mixture_ratio,
            cstar,
            self.fuel,
            self.oxidizer,
            [point])


        

    def report(self, detailed = False):
        pass


