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

def get_isp(fuel: str, oxidizer: str, pressure: float, OF_ratio: float, expansion_ratio: float, frozen=False, output_in_seconds=True):
    fuel = utils.generate_solution_from_text_input(fuel, utils.to_si(pressure))
    oxidizer = utils.generate_solution_from_text_input(oxidizer, utils.to_si(pressure))

    result = _get_isp(fuel, oxidizer, pressure, OF_ratio, expansion_ratio, frozen)
    isp = result.get_isp()

    if not output_in_seconds:
        return isp
    else:
        return list(map(lambda i: i[1], isp)) #why are we doing this?


def get_isp_vac(fuel: str, oxidizer: str, pressure: float, OF_ratio: float, expansion_ratio: float, frozen=False, output_in_seconds=True):
    fuel = utils.generate_solution_from_text_input(fuel, utils.to_si(pressure))
    oxidizer = utils.generate_solution_from_text_input(oxidizer, utils.to_si(pressure))

    result = _get_isp(fuel, oxidizer, pressure, OF_ratio, expansion_ratio, frozen)
    ivac = result.get_ivac()

    if not output_in_seconds:
        return ivac
    else:
        return list(map(lambda i: i[1], ivac)) #why are we doing this?


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