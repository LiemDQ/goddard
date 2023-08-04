from combustion import *
import cantera as ct
from utils import g0

def get_cstar(fuel, oxidizer, pressure, of_ratio, expansion_ratio):
    '''
    Calculate characteristic combustion velocity. 
    '''
    pass

def get_isp(fuel, oxidizer, pressure, OF_ratio, expansion_ratio, frozen=False, return_as_seconds=True):
    '''
    Calculate sea-level specific impulse.
    '''
    result = _get_isp(fuel, oxidizer, pressure, OF_ratio, expansion_ratio, frozen)
    isp = result.get_isp()

    if not return_as_seconds:
        return isp
    else:
        return list(map( lambda i: i[1], isp))



def get_ivac(fuel, oxidizer, pressure, OF_ratio, expansion_ratio, frozen=False, return_as_seconds=True):
    '''
    Calculate vacuum specific impulse. 
    '''
    result = _get_isp(fuel, oxidizer, pressure, OF_ratio, expansion_ratio, frozen)

    ivac = result.get_ivac()

    if not return_as_seconds:
        return ivac
    else:
        return list(map(lambda i: i[1], ivac))

def _get_isp(fuel, oxidizer, pressure, OF_ratio, expansion_ratio, frozen=False):
    if frozen:
        nt = NozzleType.FROZEN
    else:
        nt = NozzleType.EQ

    analysis = CombustionAnalysis(pressure, fuel, oxidizer, OF_ratio, expansion_ratio, CombustorArgs(nozzle_type=nt))
    result = analysis.run()

    return result

