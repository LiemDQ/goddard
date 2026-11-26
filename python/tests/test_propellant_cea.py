"""Compare RocketSolver against NASA CEA (the `cea` Python module) on common propellants.

Infinite-area rocket problems, each solved with an equilibrium and a frozen (throat) nozzle
over a grid of mixture ratios (lean to rich) and chamber pressures:

* LCH4/LOX and gaseous CH4/O2
* N2O4/MMH, storable hypergolic
* RP-1/LOX, at engine and gas-generator mixture ratios
* H2O2 monopropellant, 90 % and 98 %
* N2H4 monopropellant
* AP/Al/HTPB composite solid propellant
* N2O4/UDMH and N2O4/Aerozine-50
* LF2/LH2

Both codes run on the same NASA9 data and are given CEA's own product list for the reactants,
condensed species included, so the tolerances are those of the RP-1311 example 8 comparison in
`test_condensed_cea.py`, not the cross-database ones of `test_cea_comparison.py`.

CEA's rocket solver is not a reliable oracle everywhere (see "Known defects of the oracle" in
`docs/contributing/cea_oracle.md`): its stations are checked before they are compared, and the
equilibrium exits are also checked against CEA's plain TP equilibrium, which does not depend on
its rocket solver.

Tests are skipped automatically if the `cea` package is not installed.
"""
from dataclasses import dataclass, replace
from importlib.metadata import version

import numpy as np
import pytest

cea = pytest.importorskip("cea")

import goddard
from conftest import assert_close_rel
from test_condensed_cea import (BAR, NASA9_CONDENSED, NASA9_GAS, NASA9_REACTANTS,
                                run_isolated, split_product_names)

pytestmark = pytest.mark.slow


@dataclass(frozen=True)
class PropellantCase:
    """A propellant at one operating point. Each stream maps reactant names to mass fractions; a
    monopropellant is the same stream as fuel and oxidizer."""
    name: str
    fuel: dict
    fuel_temperature: float          # [K]
    oxidizer: dict
    oxidizer_temperature: float      # [K]
    of_ratio: float
    chamber_pressure_bar: float
    area_ratios: tuple
    frozen: bool = True              # also solve with a frozen nozzle

    def at(self, of_ratio, chamber_pressure_bar):
        """The same propellant at another operating point, named after it."""
        return replace(self, name=f"{self.name}-of{of_ratio:g}-pc{chamber_pressure_bar:g}",
                       of_ratio=of_ratio, chamber_pressure_bar=chamber_pressure_bar)


# Reactants in neither database, defined by formula and heat of formation as in CEA's input.
# Goddard has no API for them (see the coverage report), so the test appends them to a copy of
# the reactant file as assigned-enthalpy species, the form the converter gives CEA's own.
CUSTOM_REACTANTS = {
    # A representative HTPB-like binder, not a particular product: both codes get the same input
    "HTPB": {"formula": {"C": 7.075, "H": 10.65, "O": 0.223, "N": 0.063},
             "enthalpy_kj_per_mol": -5.81, "temperature": 298.15},
}

H2O2_90 = {"H2O2(L)": 0.9, "H2O(L)": 0.1}
A50 = {"C2H8N2(L),UDMH": 0.5, "N2H4(L)": 0.5}
AP_AL_HTPB = {"NH4CLO4(I)": 0.68, "AL(cr)": 0.18, "HTPB": 0.14}
H2O2_98 = {"H2O2(L)": 0.98, "H2O(L)": 0.02}

# Each propellant on a grid of mixture ratios (lean to rich, the nominal one among them) times
# low, nominal and high chamber pressure.
NOMINAL_CASES = [
    PropellantCase("lch4_lox", {"CH4(L)": 1.0}, 111.64, {"O2(L)": 1.0}, 90.17, 3.4, 100.0,
                   (20.0, 50.0)),
    PropellantCase("ch4_o2_gas", {"CH4": 1.0}, 298.15, {"O2": 1.0}, 298.15, 3.4, 100.0,
                   (20.0, 50.0)),
    PropellantCase("nto_mmh", {"CH6N2(L)": 1.0}, 298.15, {"N2O4(L)": 1.0}, 298.15, 2.0, 10.0,
                   (40.0, 100.0)),
    PropellantCase("rp1_lox", {"RP-1": 1.0}, 298.15, {"O2(L)": 1.0}, 90.17, 2.6, 70.0,
                   (16.0, 40.0)),
    # Gas-generator mixture ratio: a cool, very rich chamber that deposits graphite
    PropellantCase("rp1_lox_rich", {"RP-1": 1.0}, 298.15, {"O2(L)": 1.0}, 90.17, 0.3, 70.0,
                   (4.0, 10.0)),
    # Hydrogen peroxide monopropellant, 90 % and 98 % by mass in water: one stream, a cool chamber,
    # water that can condense in the nozzle
    PropellantCase("h2o2_90", H2O2_90, 298.15, H2O2_90, 298.15, 1.0, 10.0, (20.0, 50.0)),
    PropellantCase("h2o2_98", H2O2_98, 298.15, H2O2_98, 298.15, 1.0, 10.0, (20.0, 50.0)),
    # Hydrazine monopropellant. Equilibrium decomposes the ammonia completely, which a catalyst
    # bed does not; the comparison is of the equilibrium model. The chamber is at 870-950 K, and
    # a frozen expansion beyond AR 10 falls below the 200 K end of the H2 data.
    PropellantCase("n2h4", {"N2H4(L)": 1.0}, 298.15, {"N2H4(L)": 1.0}, 298.15, 1.0, 10.0,
                   (4.0, 8.0)),
    # Aluminized composite solid propellant: one stream, alumina liquid in the chamber that
    # solidifies in the nozzle. Frozen at the throat, the liquid alumina leaves its data range
    # (2327 K) between AR 2 and AR 10, and both codes stop there; frozen expansions are compared
    # in the near-throat case only.
    PropellantCase("ap_al_htpb", AP_AL_HTPB, 298.15, AP_AL_HTPB, 298.15, 1.0, 70.0, (10.0, 40.0),
                   frozen=False),
    PropellantCase("ap_al_htpb_near", AP_AL_HTPB, 298.15, AP_AL_HTPB, 298.15, 1.0, 70.0,
                   (1.5, 2.0)),
    # Storables: NTO with UDMH, and with Aerozine-50 (a two-species fuel stream)
    PropellantCase("nto_udmh", {"C2H8N2(L),UDMH": 1.0}, 298.15, {"N2O4(L)": 1.0}, 298.15, 2.6, 20.0,
                   (40.0, 100.0)),
    PropellantCase("nto_a50", A50, 298.15, {"N2O4(L)": 1.0}, 298.15, 2.0, 20.0, (40.0, 100.0)),
    # Fluorine: HF products and a chamber near 4500 K
    PropellantCase("lf2_lh2", {"H2(L)": 1.0}, 20.27, {"F2(L)": 1.0}, 85.02, 12.0, 70.0,
                   (40.0, 100.0)),
]
SWEEP_OF_RATIOS = {
    "lch4_lox": (2.0, 2.8, 3.4, 4.0, 5.0),
    "ch4_o2_gas": (2.0, 2.8, 3.4, 4.0, 5.0),
    "nto_mmh": (1.2, 1.65, 2.0, 2.5, 3.0),
    "rp1_lox": (1.6, 2.0, 2.6, 3.0, 3.6),
    "rp1_lox_rich": (0.25, 0.3, 0.35, 0.45, 0.6),
    # A monopropellant: the "mixture ratio" of a stream with itself
    "h2o2_90": (1.0,),
    "h2o2_98": (1.0,),
    "n2h4": (1.0,),
    "ap_al_htpb": (1.0,),
    "nto_udmh": (1.2, 1.6, 2.0, 2.6, 3.2),
    "nto_a50": (1.2, 1.6, 2.0, 2.6, 3.2),
    "ap_al_htpb_near": (1.0,),
    # Stoichiometric is 19
    "lf2_lh2": (6.0, 9.0, 12.0, 15.0, 19.0),
}
SWEEP_PRESSURES_BAR = (5.0, 200.0)   # with the nominal pressure

PROPELLANT_CASES = [
    case.at(of_ratio, pressure)
    for case in NOMINAL_CASES
    for pressure in (SWEEP_PRESSURES_BAR[0], case.chamber_pressure_bar, SWEEP_PRESSURES_BAR[1])
    for of_ratio in SWEEP_OF_RATIOS[case.name]
]

CHEMISTRIES = ["equilibrium", "frozen"]

# Cantera 3.x's "gibbs" solver (MultiPhaseEquil) cannot converge when the equilibrium amount of a
# species underflows a double: here heavy hydrocarbons in cool (900-1150 K), oxygen-rich exits. Its
# trace-species test is never applied (fixed on Cantera main, to be released as 4.0, where these
# cases pass).
GIBBS_TRACE_SPECIES_CASES = ("nto_mmh-of3-", "nto_a50-of2.6-", "nto_a50-of3.2-")
GIBBS_TRACE_SPECIES_BUG = pytest.mark.xfail(
    int(version("cantera").split(".")[0]) < 4, strict=True,
    reason="Cantera < 4 MultiPhaseEquil fails when a trace species underflows")


# Cantera 3.x's "gibbs" loses element mass in the AP/Al/HTPB chamber at 5 bar: carbon changes by
# 1.3e-7 relative, past Goddard's 1e-7 check, from every starting guess. On Cantera main
# ("Fix element conservation for trace minor species in MultiPhaseEquil::step") the case passes.
GIBBS_ELEMENT_CONSERVATION_BUG = pytest.mark.xfail(
    int(version("cantera").split(".")[0]) < 4, strict=True,
    reason="Cantera < 4 MultiPhaseEquil does not conserve carbon to 1e-7 here")


def _marks(case, chemistry):
    if chemistry == "equilibrium" and case.name.startswith(GIBBS_TRACE_SPECIES_CASES):
        return [GIBBS_TRACE_SPECIES_BUG]
    if case.name.startswith("ap_al_htpb") and case.name.endswith("-pc5"):  # in the chamber
        return [GIBBS_ELEMENT_CONSERVATION_BUG]
    return []


# Exits where CEA's rocket solver fails, so they cannot be compared: (name, chemistry, exit index).
# See docs/contributing/cea_oracle.md. Equilibrium exits are still checked against CEA's TP
# equilibrium (`test_propellant_exit_entropy`).
KNOWN_CEA_ROCKET_FAILURES = {
    # Rich RP-1/LOX at 5 bar, equilibrium: "Too many singular matrices" in the expansion with
    # graphite; the area-ratio stations come back NaN while `converged` is True.
    ("rp1_lox_rich-of0.25-pc5", "equilibrium", 0), ("rp1_lox_rich-of0.25-pc5", "equilibrium", 1),
    ("rp1_lox_rich-of0.3-pc5", "equilibrium", 0), ("rp1_lox_rich-of0.3-pc5", "equilibrium", 1),
    ("rp1_lox_rich-of0.35-pc5", "equilibrium", 1),
    # Rich RP-1/LOX at 5 bar, frozen: CEA stops ("Frozen calculations stopped: temperature is
    # more than 50 K outside the range of a condensed species") before AR 10, where Goddard is at
    # 466-509 K, although graphite has data down to 200 K.
    *((f"rp1_lox_rich-of{of:g}-pc5", "frozen", 1) for of in (0.25, 0.3, 0.35, 0.45, 0.6)),
}


def _cea_species(reactant_names):
    """Reactant names for CEA, custom ones (`CUSTOM_REACTANTS`) as `cea.Reactant` definitions."""
    species = []
    for name in reactant_names:
        custom = CUSTOM_REACTANTS.get(name)
        species.append(name if custom is None else cea.Reactant(
            name=name, formula=custom["formula"], enthalpy=custom["enthalpy_kj_per_mol"],
            enthalpy_units="kJ/mol", temperature=custom["temperature"]))
    return species


def _cea_reactants(reactant_names):
    """CEA mixture of the reactants, custom ones included."""
    return cea.Mixture(_cea_species(reactant_names))


def _cea_product_names(reactant_names):
    """Species CEA itself would use for these reactants."""
    return list(cea.Mixture(_cea_species(reactant_names), products_from_reactants=True).species_names)


def _cea_weights_and_enthalpy(reactant_names, oxidizer_weights, fuel_weights, of_ratio,
                              temperatures):
    """Normalized reactant mass fractions and the mixture enthalpy [J/kg] at `temperatures`."""
    reactants = _cea_reactants(reactant_names)
    weights = reactants.of_ratio_to_weights(np.asarray(oxidizer_weights),
                                            np.asarray(fuel_weights), of_ratio)
    weights = weights / weights.sum()
    enthalpy = reactants.calc_property(cea.ENTHALPY, weights, np.asarray(temperatures))
    return list(weights), float(enthalpy)


def _cea_moles(reactant_names, weights):
    """Mole fractions of a stream from its reactant mass fractions."""
    moles = _cea_reactants(reactant_names).weights_to_moles(np.asarray(weights, dtype=float))
    return [float(value) for value in moles / moles.sum()]


def _cea_rocket(reactant_names, product_names, weights, chamber_pressure_bar, hc, pi_p, supar,
                n_frz=None):
    """Solve an infinite-area-combustor rocket problem with `cea.RocketSolver`."""
    reactants = _cea_reactants(reactant_names)
    products = cea.Mixture(product_names)
    solver = cea.RocketSolver(products, reactants=reactants)
    solution = cea.RocketSolution(solver)

    options = {"n_frz": n_frz} if n_frz is not None else {}
    if supar:  # pycea rejects an empty list
        options["supar"] = list(supar)
    solver.solve(solution, np.asarray(weights), chamber_pressure_bar, list(pi_p),
                 iac=True, hc=hc, **options)
    return {
        "S": [float(value) for value in solution.entropy],   # [kJ/(kg K)]
        "T": [float(value) for value in solution.T],
        "P": [float(value) for value in solution.P],
        "M": [float(value) for value in solution.M],
        "gamma_s": [float(value) for value in solution.gamma_s],
        "Mach": [float(value) for value in solution.Mach],
        "ae_at": [float(value) for value in solution.ae_at],
        "c_star": [float(value) for value in solution.c_star],
        "Isp": [float(value) for value in solution.Isp],
        "Ivac": [float(value) for value in solution.Isp_vacuum],
        "CF": [float(value) for value in solution.coefficient_of_thrust],
        "num_pts": int(solution.num_pts),
        "converged": bool(solution.converged),
    }


def _cea_tp_entropy(reactant_names, product_names, weights, states):
    """Entropy [kJ/(kg K)] of CEA's TP equilibrium at each (T [K], P [bar]) of `states`."""
    reactants = _cea_reactants(reactant_names)
    products = cea.Mixture(product_names)
    solver = cea.EqSolver(products, reactants=reactants)
    solution = cea.EqSolution(solver)
    entropies = []
    for temperature, pressure in states:
        solver.solve(solution, cea.TP, temperature, pressure, np.asarray(weights))
        entropies.append(float(solution.entropy) if solution.converged else float("nan"))
    return entropies


# Infinite area: chamber, throat, the placeholder pressure-ratio station, then the area ratios
FIRST_AREA_RATIO_STATION = 3


def cea_exit_station(case, rocket, exit_index):
    """Index of the CEA station of exit `exit_index` of `case`.

    @raises ValueError if CEA has no station at that area ratio: a solve CEA stopped, or whose
    stations are NaN.
    """
    area_ratio = case.area_ratios[exit_index]
    for index in range(FIRST_AREA_RATIO_STATION, rocket["num_pts"]):
        if abs(rocket["ae_at"][index] / area_ratio - 1.0) <= 1e-3:
            return index
    raise ValueError(f"CEA has no station at area ratio {area_ratio:g}")


def cea_problem(case, rocket, exit_index):
    """Why CEA's station of an exit cannot be used as the reference, or an empty string.

    `converged` alone does not tell: CEA can report convergence for an expansion whose stations
    are NaN, and a solve it stops partway keeps valid stations before the stop.
    """
    ratio = case.area_ratios[exit_index]
    stopped = "" if rocket["converged"] else " (CEA reports the solve did not converge)"
    try:
        index = cea_exit_station(case, rocket, exit_index)
    except ValueError as error:
        return f"{error}{stopped}"
    values = [rocket[key][index] for key in ("T", "P", "Mach", "ae_at")]
    if not all(np.isfinite(values)):
        return f"CEA's station at area ratio {ratio:g} is not finite: T, P, Mach, ae_at = {values}"
    if not rocket["T"][index] < rocket["T"][1]:
        return (f"CEA's station at area ratio {ratio:g} is hotter than the throat: "
                f"{rocket['T'][index]:.6g} K against {rocket['T'][1]:.6g} K")
    return ""


def require_cea_station(case, chemistry, rocket, exit_index):
    """Fail unless CEA's station of the exit is usable; xfail a known CEA failure, which must
    fail."""
    problem = cea_problem(case, rocket, exit_index)
    if (case.name, chemistry, exit_index) in KNOWN_CEA_ROCKET_FAILURES:
        if problem:
            pytest.xfail(f"known CEA rocket solver failure: {problem}")
        pytest.fail("listed in KNOWN_CEA_ROCKET_FAILURES, but CEA solved it: remove it")
    assert not problem, f"CEA's rocket solve cannot serve as the reference: {problem}"


@pytest.fixture(scope="module")
def reactant_file(tmp_path_factory):
    """The NASA9 reactant file with `CUSTOM_REACTANTS` appended as assigned-enthalpy species.

    The converter writes the `species` section last, so a list item appended to the file is a
    species. Each one is a constant-cp species with cp0 = 0: h(T) = h_formation, as CEA takes it.
    """
    lines = []
    for name, custom in CUSTOM_REACTANTS.items():
        composition = ", ".join(f"{element}: {atoms!r}" for element, atoms in custom["formula"].items())
        lines += [f"- name: {name}",
                  f"  composition: {{{composition}}}",
                  f"  thermo: {{model: constant-cp, reference-pressure: 100000.0, "
                  f"T0: {custom['temperature']!r}, h0: {custom['enthalpy_kj_per_mol'] * 1e3!r} J/mol, "
                  f"s0: 0.0, cp0: 0.0}}"]
    path = tmp_path_factory.mktemp("reactants") / "nasa9_reactants_custom.yaml"
    with open(NASA9_REACTANTS) as source:
        text = source.read()
    path.write_text(text.rstrip("\n") + "\n" + "\n".join(lines) + "\n")
    for name in CUSTOM_REACTANTS:
        goddard.Gas(str(path), phase_name="reactants", species={name})   # loads as a species
    return str(path)


@pytest.fixture(scope="module",
                params=[pytest.param((case, chemistry), id=f"{case.name}-{chemistry}",
                                     marks=_marks(case, chemistry))
                        for case in PROPELLANT_CASES for chemistry in CHEMISTRIES
                        if case.frozen or chemistry == "equilibrium"])
def propellant_rocket(request, reactant_file):
    """One case solved by both codes: a RocketSolver with reactant streams, and cea.RocketSolver."""
    case, chemistry = request.param
    reactant_names = list(dict.fromkeys([*case.fuel, *case.oxidizer]))
    names = run_isolated(_cea_product_names, reactant_names)
    gas_names, condensed_names, missing = split_product_names(names)
    assert not missing, f"CEA products absent from Goddard's data: {missing}"

    temperatures = [case.fuel_temperature if name in case.fuel else case.oxidizer_temperature
                    for name in reactant_names]
    weights, enthalpy = run_isolated(_cea_weights_and_enthalpy, reactant_names,
                                     [case.oxidizer.get(name, 0.0) for name in reactant_names],
                                     [case.fuel.get(name, 0.0) for name in reactant_names],
                                     case.of_ratio, temperatures)
    # Goddard frozen_NFZ=1 (freeze at the throat) is CEA n_frz=2 for an infinite-area chamber.
    # CEA needs at least one pressure ratio; area-ratio stations are found by their ae_at. The
    # placeholder must stay warm: a frozen expansion that goes below CEA's data range stops, and
    # drops every station after it (Pc/P = 1000 did, for N2H4 and fuel-rich RP-1/LOX).
    rocket = run_isolated(_cea_rocket, reactant_names, names, weights, case.chamber_pressure_bar,
                          enthalpy / cea.R, [3.0], case.area_ratios,
                          n_frz=2 if chemistry == "frozen" else None)
    rocket["reactants"] = (reactant_names, names, weights)

    chamber_pressure = case.chamber_pressure_bar * BAR
    products = goddard.Gas(NASA9_GAS, "gas", species=set(gas_names),
                           condensed_file=NASA9_CONDENSED,
                           condensed_species=set(condensed_names))
    nozzle = (goddard.equilibrium_nozzle() if chemistry == "equilibrium"
              else goddard.frozen_nozzle(frozen_NFZ=1))
    fuel_moles = run_isolated(_cea_moles, list(case.fuel), list(case.fuel.values()))
    oxidizer_moles = run_isolated(_cea_moles, list(case.oxidizer), list(case.oxidizer.values()))
    point = goddard.RocketPoint(
        fuel=goddard.PhaseSpecification(case.fuel_temperature, chamber_pressure,
                                        dict(zip(case.fuel, fuel_moles))),
        oxidizer=goddard.PhaseSpecification(case.oxidizer_temperature, chamber_pressure,
                                            dict(zip(case.oxidizer, oxidizer_moles))),
        mixture_ratio=case.of_ratio,
        chamber_pressure=chamber_pressure,
        combustor=goddard.infinite_area_combustor(),
        nozzle=nozzle,
        expansion_ratios=list(case.area_ratios),
        label=case.name)

    result = goddard.RocketSolver(products, reactant_file).solve(point)
    return case, chemistry, result, rocket


def test_propellant_chamber(propellant_rocket):
    case, _, result, rocket = propellant_rocket
    assert result.status.ok, result.status.message
    chamber = result.chamber.thermo

    assert_close_rel(chamber.temperature, rocket["T"][0], 1e-4, "chamber temperature")
    assert_close_rel(chamber.molecular_weight, rocket["M"][0], 2e-4, "chamber M")
    assert_close_rel(chamber.gamma_s, rocket["gamma_s"][0], 3e-5, "chamber gamma_s")
    assert_close_rel(chamber.pressure, case.chamber_pressure_bar * BAR, 1e-6, "chamber pressure")


def test_propellant_throat(propellant_rocket):
    _, _, result, rocket = propellant_rocket
    assert result.status.ok, result.status.message
    throat = result.throat.thermo
    station = 1

    assert_close_rel(throat.pressure / BAR, rocket["P"][station], 1e-5, "throat pressure")
    assert_close_rel(throat.temperature, rocket["T"][station], 1e-4, "throat temperature")
    assert_close_rel(throat.molecular_weight, rocket["M"][station], 2e-4, "throat M")
    assert_close_rel(throat.gamma_s, rocket["gamma_s"][station], 3e-5, "throat gamma_s")
    # c* depends on the chamber and the throat only, so any exit station reports the same value.
    assert_close_rel(result.performance[0].cstar, rocket["c_star"][station], 5e-5, "c*")


@pytest.mark.parametrize("exit_index", [0, 1])
def test_propellant_exits(propellant_rocket, exit_index):
    case, chemistry, result, rocket = propellant_rocket
    assert result.status.ok, result.status.message
    require_cea_station(case, chemistry, rocket, exit_index)
    station = result.exits[exit_index].thermo
    index = cea_exit_station(case, rocket, exit_index)
    label = f"AR {case.area_ratios[exit_index]}"

    assert_close_rel(station.temperature, rocket["T"][index], 3e-4, f"{label} temperature")
    assert_close_rel(station.pressure / BAR, rocket["P"][index], 2e-4, f"{label} pressure")
    assert_close_rel(station.molecular_weight, rocket["M"][index], 2e-4, f"{label} M")
    assert_close_rel(station.gamma_s, rocket["gamma_s"][index], 5e-5, f"{label} gamma_s")

    performance = result.performance[exit_index]
    assert_close_rel(performance.area_ratio, rocket["ae_at"][index], 1e-4, f"{label} area_ratio")
    assert_close_rel(performance.cstar, rocket["c_star"][index], 5e-5, f"{label} c*")
    assert_close_rel(performance.isp, rocket["Isp"][index], 3e-5, f"{label} Isp")
    assert_close_rel(performance.ivac, rocket["Ivac"][index], 3e-5, f"{label} Ivac")
    assert_close_rel(performance.mach_number, rocket["Mach"][index], 1e-4, f"{label} Mach")
    assert_close_rel(performance.CF, rocket["CF"][index], 3e-5, f"{label} CF")


@pytest.mark.parametrize("exit_index", [0, 1])
def test_propellant_exit_entropy(propellant_rocket, exit_index):
    """An equilibrium exit is CEA's TP equilibrium at its temperature and pressure, with the
    chamber's entropy. CEA's TP solve does not go through its rocket solver, so this checks the
    exits also where that fails."""
    case, chemistry, result, rocket = propellant_rocket
    if chemistry != "equilibrium":
        pytest.skip("a frozen exit is not an equilibrium state")
    assert result.status.ok, result.status.message
    station = result.exits[exit_index].thermo
    reactant_names, product_names, weights = rocket["reactants"]
    [entropy] = run_isolated(_cea_tp_entropy, reactant_names, product_names, weights,
                             [(station.temperature, station.pressure / BAR)])
    if not np.isfinite(entropy):
        # CEA's TP solve, started afresh, can fail where its rocket solver, continuing from the
        # previous station, does not (cool H2-rich LF2/LH2 products)
        assert not cea_problem(case, rocket, exit_index), (
            "neither CEA's TP solve nor its rocket solver gives a reference for this exit")
        pytest.skip("CEA's TP solve does not converge here; test_propellant_exits compares it")
    label = f"AR {case.area_ratios[exit_index]}"
    if station.pinned_transition:
        # Pinned at a phase transition: the temperature is exact and cp is reported as 0
        assert_close_rel(entropy, rocket["S"][0], 1e-6, f"{label} entropy")
        return
    # Compared as the temperature error the entropy mismatch implies, dT/T = ds/cp: where water
    # condenses, a 1e-7 temperature difference moves the entropy by 1e-6. The tolerance allows
    # for the 1-2e-6 by which exit temperatures also differ from CEA's rocket stations (LCH4/LOX).
    temperature_error = (entropy - rocket["S"][0]) * 1e3 / station.cp
    assert abs(temperature_error) <= 5e-6, (
        f"{label}: entropy {entropy:.9g} kJ/(kg K) against the chamber's {rocket['S'][0]:.9g}, "
        f"a relative temperature error of {temperature_error:.3e}")
