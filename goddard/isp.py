from combustion import *
import cantera as ct
from utils import g0
from scipy import brentq

def get_cstar(gamma, temperature, molecular_weight):
    '''
    Calculate the C* value for a given combustion temperature and set of gas properties. C* is also known as the characteristic velocity,
    and is a measure of the energy present in the fluid as a result of combustion. It is used as a metric for engine performance independent
    of nozzle expansion.
    '''
    return (
        np.sqrt(ct.gas_constant * temperature / (molecular_weight * gamma)) *
        np.power(2 / (gamma + 1), -(gamma + 1) / (2*(gamma - 1)))
        )

def get_mach_from_subsonic_area_ratio(gamma, ratio):
    '''
    For a given area ratio, calculate the mach in the smaller area given the mach number in the larger area.
    There are two possible solutions for a given ratio (one subsonic, one supersonic). 
    This function always returns the subsonic solution.
    '''
    return brentq(lambda mach: get_area_ratio_from_mach_num(gamma, mach)-ratio, 0, 1.0)
    

def get_mach_from_supersonic_area_ratio(gamma, ratio):
    '''
    For a given area ratio, calculate the mach in the smaller area given the mach number in the larger area.
    There are two possible solutions for a given ratio (one subsonic, one supersonic). 
    This function always returns the supersonic solution.
    '''

    #ideally, should be inf, but this would not converge for certain ill-conditioned edge cases.
    upper_mach_limit = 10,000,000.0 
    return brentq(lambda mach: get_area_ratio_from_mach_num(gamma, mach)-ratio, 1.0, upper_mach_limit)

def get_area_ratio_from_mach_num(gamma, M):
    '''
    Area ratio (w.r.t throat) for a given mach number and isentropic expansion factor.
    '''
    return 1/M * ((1+ (gamma-1)/2*M**2) /((gamma+1)/2))**((gamma+1)/(2*(gamma-1)))

def get_velocity(gas: ct.Solution, stagnation_enthalpy):
    '''
    Velocity in isentropic supersonic flow can be found from the difference in enthalpy
    between the gas in motion and the stagnation enthalpy. 
    '''
    return np.sqrt(2*(stagnation_enthalpy - gas.enthalpy_mass))

def get_isp(gas, gamma, enthalpy):
    '''
    Calculate specific impulse given an exhaust gas thermodynamic state and combustion chamber enthalpy. 
    '''
    velocity = get_velocity(gas, enthalpy)
    cstar = get_cstar(gamma, gas.T, gas.mean_molecular_weight)
    CF = velocity/cstar
    g0 = 9.80655 #gravitational acceleration
    return velocity, velocity/g0
    

def get_ivac(gas, isp):
    '''
    Ivac also includes the thrust from pressure forces. This assumes isp is provided in (m/s)
    and not s.
    '''
    ivac = isp + gas.T * ct.gas_constant / (isp * gas.mean_molecular_weight)
    g0 = 9.80655 #gravitational acceleration

    return ivac, ivac/g0