import numpy as np
import cantera as ct 

def calculate_c_star(gamma, temperature, molecular_weight):
    '''Calculate the C* value for a given combustion temperature and set of gas properties. C* is also known as the characteristic velocity,
    and is a measure of the energy present in the fluid as a result of combustion. It is used as a metric for engine performance independent
    of nozzle expansion.
    '''
    return (
        np.sqrt(ct.gas_constant * temperature / (molecular_weight * gamma)) *
        np.power(2 / (gamma + 1), -(gamma + 1) / (2*(gamma - 1)))
        )

def calculate_mach_number(gamma, expansion_ratio):
    pass

def calculate_isp():
    pass

def calculate_isp_vac():
    pass

def calculate_pressure_ratio():
    pass

def calculate_expansion_ratio():
    pass
