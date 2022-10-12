from combustion import *
import cantera as ct
from utils import g0

def calc_cstar(fuel, oxidizer, pressure, of_ratio, expansion_ratio):
    '''
    Calculate characteristic combustion velocity. 
    '''
    pass

def calc_isp_sl(fuel, oxidizer, pressure, of_ratio, expansion_ratio, frozen=False, return_as_seconds=True):
    '''
    Calculate sea-level specific impulse.
    '''
    analysis = CombustionAnalysis(pressure, fuel, oxidizer, of_ratio, expansion_ratio, CombustorArgs())
    result = analysis.run()

    isp = result.get_isp()

    if not return_as_seconds:
        return isp
    else:
        return list(map( lambda i: i[1], isp))



def calc_isp_vac(fuel, oxidizer, pressure, of_ratio, expansion_ratio, frozen=False, return_as_seconds=True):
    '''
    Calculate vacuum specific impulse. 
    '''
    pass

def calc_exhaust_temperature(fuel, oxidizer, pressure, of_ratio, expansion_ratio):
    pass

def calc_analyze(fuel, oxidizer, pressure, of_ratio, expansion_ratio):
    pass