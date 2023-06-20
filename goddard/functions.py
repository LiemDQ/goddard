from combustion import *
import cantera as ct
from utils import g0

def get_cstar(fuel, oxidizer, pressure, of_ratio, expansion_ratio):
    '''
    Calculate characteristic combustion velocity. 
    '''
    pass

def calc_isp_sl(fuel, oxidizer, pressure, of_ratio, expansion_ratio, frozen=False, return_as_seconds=True):
    '''
    Calculate sea-level specific impulse.
    '''
    result = _get_isp(fuel, oxidizer, pressure, of_ratio, expansion_ratio, frozen)
    isp = result.get_isp()

    if not return_as_seconds:
        return isp
    else:
        return list(map( lambda i: i[1], isp))



def calc_isp_vac(fuel, oxidizer, pressure, of_ratio, expansion_ratio, frozen=False, return_as_seconds=True):
    '''
    Calculate vacuum specific impulse. 
    '''
    result = _get_isp(fuel, oxidizer, pressure, of_ratio, expansion_ratio, frozen)

    ivac = result.get_ivac()

    if not return_as_seconds:
        return ivac
    else:
        return list(map(lambda i: i[1], ivac))

def _get_isp(fuel, oxidizer, pressure, of_ratio, expansion_ratio, frozen=False):
    if frozen:
        nt = NozzleType.FROZEN
    else:
        nt = NozzleType.EQ

    analysis = CombustionAnalysis(pressure, fuel, oxidizer, of_ratio, expansion_ratio, CombustorArgs(is_frozen=nt))
    result = analysis.run()

    return result


def calc_exhaust_temperature(fuel, oxidizer, pressure, of_ratio, expansion_ratio):
    pass

def calc_analyze(fuel, oxidizer, pressure, of_ratio, expansion_ratio):
    pass