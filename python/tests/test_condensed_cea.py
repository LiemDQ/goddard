"""Compare Goddard's condensed-phase equilibrium against NASA CEA (the `cea` Python module).

Both codes run on the same NASA9 thermodynamic data -- Goddard reads the YAML files converted
from `data/nasa9.dat`, CEA reads `thermo.lib` -- and are given the same product species, so the
tolerances here are much tighter than the cross-database ones of `test_cea_comparison.py`.

Two properties of pycea shape the code below:

* A Fortran `STOP 1` (an unknown species name, or a re-insertion at a transition temperature)
  kills the interpreter, so every pycea call runs in a forked child through `run_isolated`.
  Each child solves a whole group of points and returns plain data, so the fork cost is paid
  once per test rather than once per point.
* `Mixture.of_ratio_to_weights` returns unnormalized weights, which must be normalized before
  they are used as mass fractions.

Species names are spelled identically in both databases (including `H2O(cr)` for ice), so the
product list of one code transfers to the other unchanged.

Tests are skipped automatically if the `cea` package is not installed.
"""
import multiprocessing as mp
import os
from functools import lru_cache

import numpy as np
import pytest

cea = pytest.importorskip("cea")

import goddard
from goddard import Gas, GasChemistry, MixtureRatioType, reactant_gas
from conftest import find_data_dir, assert_close_abs, assert_close_rel

DATA_DIR = find_data_dir()
NASA9_GAS = os.path.join(DATA_DIR, "nasa9_gas.yaml")
NASA9_CONDENSED = os.path.join(DATA_DIR, "nasa9_condensed.yaml")
NASA9_REACTANTS = os.path.join(DATA_DIR, "nasa9_reactants.yaml")

BAR = 1.0e5


# ---------------------------------------------------------------------------
# Running pycea out of process
# ---------------------------------------------------------------------------

def run_isolated(fn, *args, timeout=300, **kwargs):
    """Run `fn` in a forked child so a Fortran STOP cannot kill the test session.

    Returns the function's result. A crash or an exception in the child fails the test.
    """
    context = mp.get_context("fork")
    queue = context.Queue()

    def target(queue):
        try:
            queue.put(("ok", fn(*args, **kwargs)))
        except BaseException as error:  # noqa: BLE001 - reported to the parent as a failure
            queue.put(("err", f"{type(error).__name__}: {error}"))

    process = context.Process(target=target, args=(queue,))
    process.start()
    try:
        status, value = queue.get(timeout=timeout)
    except Exception:  # noqa: BLE001 - an empty queue means the child died
        status, value = "crash", "the pycea child process died"
    process.join(1)
    if process.is_alive():
        process.kill()

    assert status == "ok", f"pycea call failed: {value}"
    return value


# The four functions below run inside the forked child and return plain data.

def _cea_product_names(reactant_names, omit=()):
    """Species CEA itself would use for these reactants."""
    products = cea.Mixture(reactant_names, products_from_reactants=True, omit=list(omit))
    return list(products.species_names)


def _cea_weights_and_enthalpy(reactant_names, oxidizer_weights, fuel_weights, of_ratio,
                              temperatures):
    """Normalized reactant mass fractions and the mixture enthalpy [J/kg] at `temperatures`."""
    reactants = cea.Mixture(reactant_names)
    weights = reactants.of_ratio_to_weights(np.asarray(oxidizer_weights),
                                            np.asarray(fuel_weights), of_ratio)
    weights = weights / weights.sum()
    enthalpy = reactants.calc_property(cea.ENTHALPY, weights, np.asarray(temperatures))
    return list(weights), float(enthalpy)


def _cea_reactant_enthalpy(reactant_names, weights, temperatures):
    reactants = cea.Mixture(reactant_names)
    return float(reactants.calc_property(cea.ENTHALPY, np.asarray(weights),
                                         np.asarray(temperatures)))


def _cea_equilibrium(reactant_names, product_names, kind, points, weights):
    """Solve `points` (pairs of the two held properties) with `cea.EqSolver`."""
    reactants = cea.Mixture(reactant_names)
    products = cea.Mixture(product_names)
    solver = cea.EqSolver(products, reactants=reactants)
    solution = cea.EqSolution(solver)

    results = []
    for first, second in points:
        solver.solve(solution, kind, first, second, np.asarray(weights))
        results.append({
            "T": float(solution.T),
            "P": float(solution.P),
            "M": float(solution.M),
            "gamma_s": float(solution.gamma_s),
            "cp_eq": float(solution.cp_eq),
            "X": {name: float(value) for name, value in dict(solution.mole_fractions).items()},
            "converged": bool(solution.converged),
        })
    return results


def _cea_rocket(reactant_names, product_names, weights, chamber_pressure_bar, pressure_ratios,
                hc, supar=None, **solver_kwargs):
    """Solve an infinite-area-combustor rocket problem with `cea.RocketSolver`."""
    reactants = cea.Mixture(reactant_names)
    products = cea.Mixture(product_names)
    solver = cea.RocketSolver(products, reactants=reactants, **solver_kwargs)
    solution = cea.RocketSolution(solver)

    extra = {"supar": list(supar)} if supar else {}
    solver.solve(solution, np.asarray(weights), chamber_pressure_bar, list(pressure_ratios),
                 iac=True, hc=hc, **extra)
    return {
        "T": [float(value) for value in solution.T],
        "P": [float(value) for value in solution.P],
        "M": [float(value) for value in solution.M],
        "gamma_s": [float(value) for value in solution.gamma_s],
        "X": {name: [float(v) for v in values]
              for name, values in dict(solution.mole_fractions).items()},
        "num_pts": int(solution.num_pts),
        "converged": bool(solution.converged),
    }


# ---------------------------------------------------------------------------
# Building the matching Goddard problem
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def _database_names():
    """Species names of the gas and condensed NASA9 databases."""
    import cantera as ct
    gas = frozenset(species.name for species in ct.Species.list_from_file(NASA9_GAS))
    condensed = frozenset(species.name
                          for species in ct.Species.list_from_file(NASA9_CONDENSED))
    return gas, condensed


def split_product_names(names):
    """Split CEA's product list into (gas, condensed, names absent from Goddard's data)."""
    gas_names, condensed_names = _database_names()
    return ([name for name in names if name in gas_names],
            [name for name in names if name in condensed_names],
            [name for name in names
             if name not in gas_names and name not in condensed_names])


def goddard_products(product_names, chemistry=GasChemistry.FROZEN):
    """Product Gas holding exactly the species CEA was given."""
    gas_names, condensed_names, missing = split_product_names(product_names)
    assert not missing, f"CEA products absent from Goddard's data: {missing}"
    return Gas(NASA9_GAS, phase_name="products", species=set(gas_names),
               condensed_file=NASA9_CONDENSED, condensed_species=set(condensed_names),
               chemistry=chemistry)


def element_moles_of(stream, products):
    """Element amounts of `stream` [kmol/kg] in the element order of `products`."""
    amounts = dict(zip(stream.element_names, stream.element_moles))
    return np.array([amounts.get(element, 0.0) for element in products.element_names])


def condensed_mole_fraction(gas, name):
    """Mole fraction of a condensed species among all species, as CEA reports it."""
    names = gas.condensed_species_names
    if name not in names:
        return 0.0
    return gas.condensed_moles[names.index(name)] * gas.mixture_molecular_weight


def present_condensed(gas, threshold=1e-6):
    """Names of the condensed species present in more than `threshold` mole fraction."""
    return {name for name in gas.condensed_species_names
            if condensed_mole_fraction(gas, name) > threshold}


def cea_condensed(result, index=None, threshold=1e-6):
    """Same, read off a pycea result dictionary."""
    _, condensed_names = _database_names()
    fractions = result["X"] if index is None else {
        name: values[index] for name, values in result["X"].items()}
    return {name for name, value in fractions.items()
            if name in condensed_names and value > threshold}


# ---------------------------------------------------------------------------
# (a) CH4/O2 with graphite: constant temperature and pressure
# ---------------------------------------------------------------------------

METHANE_OXYGEN_TP = [(T, P, OF)
                     for T in (1000.0, 1500.0, 2500.0)
                     for P in (1.0 * BAR, 50.0 * BAR)
                     for OF in (0.5, 1.0, 2.0)]


@pytest.fixture(scope="module")
def methane_oxygen():
    """Product Gas and CEA product list for CH4/O2."""
    names = run_isolated(_cea_product_names, ["CH4", "O2"])
    return goddard_products(names, chemistry=GasChemistry.EQUILIBRIUM), names


def methane_oxygen_state(products, of_ratio):
    """Element amounts [kmol/kg] and specific enthalpy [J/kg] of CH4/O2 at 298.15 K."""
    fuel_fraction = 1.0 / (1.0 + of_ratio)
    stream = reactant_gas(NASA9_REACTANTS,
                          {"CH4": fuel_fraction, "O2": 1.0 - fuel_fraction},
                          298.15, 101325.0, basis="mass")
    return element_moles_of(stream, products), stream.enthalpy_mass


@pytest.fixture(scope="module")
def methane_oxygen_tp_reference(methane_oxygen):
    """CEA's constant-T, constant-P solutions, one child for all the points."""
    _, names = methane_oxygen
    reference = {}
    for of_ratio in {OF for _, _, OF in METHANE_OXYGEN_TP}:
        weights, _ = run_isolated(_cea_weights_and_enthalpy, ["CH4", "O2"],
                                  [0.0, 1.0], [1.0, 0.0], of_ratio, [298.15, 298.15])
        points = [(T, P / BAR) for T, P, OF in METHANE_OXYGEN_TP if OF == of_ratio]
        results = run_isolated(_cea_equilibrium, ["CH4", "O2"], names, cea.TP, points, weights)
        for (T, pressure_bar), result in zip(points, results):
            reference[(T, pressure_bar * BAR, of_ratio)] = result
    return reference


@pytest.mark.parametrize("point", METHANE_OXYGEN_TP,
                         ids=[f"T{int(T)}_P{int(P / BAR)}_OF{OF}" for T, P, OF in
                              METHANE_OXYGEN_TP])
def test_methane_oxygen_tp(methane_oxygen, methane_oxygen_tp_reference, point):
    """Fuel-rich CH4/O2 deposits graphite; compare the amount and the mixture properties.

    Observed worst errors over the 18 points: graphite 4e-4 absolute, M 1.0e-3,
    gamma_s 2.6e-4, cp_eq 7.4e-3 (all relative except the graphite fraction).
    """
    temperature, pressure, of_ratio = point
    products, _ = methane_oxygen
    reference = methane_oxygen_tp_reference[point]

    elements, _ = methane_oxygen_state(products, of_ratio)
    products.set_element_moles(elements, temperature, pressure)
    products.equilibrate_TP(temperature, pressure)

    assert reference["converged"]
    assert_close_abs(condensed_mole_fraction(products, "C(gr)"),
                     reference["X"].get("C(gr)", 0.0), 2e-3, "X[C(gr)]")
    assert_close_rel(products.mean_molecular_weight, reference["M"], 2e-3, "M")
    assert_close_rel(products.gamma_s, reference["gamma_s"], 3e-3, "gamma_s")

    properties = goddard.equilibrium_properties(products)
    assert_close_rel(properties.spec_heat_p, reference["cp_eq"] * 1e3, 1e-2, "cp_eq")


# ---------------------------------------------------------------------------
# (b) CH4/O2 through the Combustor: constant enthalpy and pressure
# ---------------------------------------------------------------------------

METHANE_OXYGEN_HP = [(P, OF) for P in (1.0 * BAR, 20.0 * BAR) for OF in (1.0, 3.0)]


@pytest.fixture(scope="module")
def methane_oxygen_hp_reference(methane_oxygen):
    """CEA's constant-H, constant-P solutions for reactants at 298.15 K."""
    _, names = methane_oxygen
    reference = {}
    for of_ratio in {OF for _, OF in METHANE_OXYGEN_HP}:
        weights, enthalpy = run_isolated(_cea_weights_and_enthalpy, ["CH4", "O2"],
                                         [0.0, 1.0], [1.0, 0.0], of_ratio, [298.15, 298.15])
        pressures = [P for P, OF in METHANE_OXYGEN_HP if OF == of_ratio]
        results = run_isolated(_cea_equilibrium, ["CH4", "O2"], names, cea.HP,
                               [(enthalpy / cea.R, P / BAR) for P in pressures], weights)
        for pressure, result in zip(pressures, results):
            reference[(pressure, of_ratio)] = result
    return reference


@pytest.mark.parametrize("point", METHANE_OXYGEN_HP,
                         ids=[f"P{int(P / BAR)}_OF{OF}" for P, OF in METHANE_OXYGEN_HP])
def test_methane_oxygen_combustor(methane_oxygen, methane_oxygen_hp_reference, point):
    """Burn CH4 in O2 through the reactant-stream Combustor.

    Observed worst errors: flame temperature 8.0e-4 relative, graphite 2.2e-4 absolute.
    """
    pressure, of_ratio = point
    products, _ = methane_oxygen
    reference = methane_oxygen_hp_reference[point]

    fuel = reactant_gas(NASA9_REACTANTS, {"CH4": 1.0}, 298.15, pressure)
    oxidizer = reactant_gas(NASA9_REACTANTS, {"O2": 1.0}, 298.15, pressure)
    options = goddard.infinite_area_combustor([pressure])
    options.mixture_type = MixtureRatioType.OF_RATIO

    states = goddard.Combustor(products, fuel, oxidizer).solve(
        np.array([pressure]), np.array([of_ratio]), options)
    assert states.size() == 1
    assert states.num_condensed() == len(products.condensed_species_names)
    assert states.condensed_species_names() == products.condensed_species_names
    condensed_moles = states.get_condensed_moles(0)
    products.restore_state(states.get_state(0))
    assert condensed_moles == pytest.approx(products.condensed_moles)

    assert_close_rel(products.temperature, reference["T"], 2e-3, "chamber temperature")
    assert_close_abs(condensed_mole_fraction(products, "C(gr)"),
                     reference["X"].get("C(gr)", 0.0), 2e-3, "X[C(gr)]")
    assert_close_rel(products.pressure, pressure, 1e-6, "chamber pressure")


# ---------------------------------------------------------------------------
# (c) H2/O2 at 0.05 atm: the water dew point (RP-1311 example 14)
# ---------------------------------------------------------------------------

WATER_TEMPERATURES = (300.0, 270.0, 250.0)


@pytest.fixture(scope="module")
def hydrogen_oxygen():
    """Products, element amounts and CEA reference for 100 mol H2 : 60 mol O2 at 0.05 atm."""
    names = run_isolated(_cea_product_names, ["H2", "O2"])
    products = goddard_products(names)

    stream = reactant_gas(NASA9_REACTANTS, {"H2": 100.0, "O2": 60.0}, 298.15, 101325.0)
    elements = element_moles_of(stream, products)

    total_mass = 100 * 2.01588 + 60 * 31.9988
    weights = [100 * 2.01588 / total_mass, 60 * 31.9988 / total_mass]
    pressure = 0.05 * 101325.0

    temperatures = WATER_TEMPERATURES + (304.0,)
    results = run_isolated(_cea_equilibrium, ["H2", "O2"], names, cea.TP,
                           [(T, pressure / BAR) for T in temperatures], weights)
    reference = dict(zip(temperatures, results))
    return products, elements, pressure, reference


@pytest.mark.parametrize("temperature", WATER_TEMPERATURES)
def test_water_condensation(hydrogen_oxygen, temperature):
    """Liquid water above 273.15 K and ice below it, in the amounts CEA reports.

    Observed worst error: 9.5e-3 absolute at 300 K, where the liquid fraction moves by about
    0.25 per kelvin, so the station is unusually sensitive to the equilibrium temperature.
    """
    products, elements, pressure, reference = hydrogen_oxygen
    expected = reference[temperature]["X"]

    products.set_element_moles(elements, temperature, pressure)
    products.equilibrate_TP(temperature, pressure)

    for name in ("H2O(L)", "H2O(cr)"):
        assert_close_abs(condensed_mole_fraction(products, name), expected.get(name, 0.0),
                         1e-2, f"X[{name}] at {temperature} K")


def test_water_dew_point_neighbourhood(hydrogen_oxygen):
    """Within a degree of the dew point only convergence and the order of magnitude are checked."""
    products, elements, pressure, reference = hydrogen_oxygen

    products.set_element_moles(elements, 304.0, pressure)
    products.equilibrate_TP(304.0, pressure)

    liquid = condensed_mole_fraction(products, "H2O(L)")
    assert 0.1 < liquid < 0.5
    assert reference[304.0]["converged"]


# ---------------------------------------------------------------------------
# (d) Aluminized ammonium perchlorate: the Al2O3 melting transition
# ---------------------------------------------------------------------------

PROPELLANT_PRESSURES = (34.47 * BAR, 3.447 * BAR, 0.3447 * BAR)


@pytest.fixture(scope="module")
def propellant():
    """Products and CEA reference for AP/Al at two aluminium fractions."""
    reactant_names = ["NH4CLO4(I)", "AL(cr)"]
    names = run_isolated(_cea_product_names, reactant_names)
    products = goddard_products(names)

    reference = {}
    for aluminium in (0.20, 0.05):
        weights = [1.0 - aluminium, aluminium]
        enthalpy = run_isolated(_cea_reactant_enthalpy, reactant_names, weights, [298.15] * 2)
        results = run_isolated(_cea_equilibrium, reactant_names, names, cea.HP,
                               [(enthalpy / cea.R, P / BAR) for P in PROPELLANT_PRESSURES],
                               weights)
        reference[aluminium] = dict(zip(PROPELLANT_PRESSURES, results))
    return products, reference


@pytest.mark.parametrize("aluminium,expected_polymorph", [(0.20, "AL2O3(L)"), (0.05, "AL2O3(a)")])
def test_aluminized_perchlorate(propellant, aluminium, expected_polymorph):
    """Burn ammonium perchlorate with aluminium at three pressures.

    20 % aluminium burns above the melting point of alumina and leaves the liquid, 5 % stays
    below it and leaves the solid. Observed worst errors: temperature 5.6e-4 relative,
    alumina 3e-5 absolute, with the same polymorph selected as CEA everywhere.
    """
    products, reference = propellant
    composition = {"NH4CLO4(I)": 1.0 - aluminium, "AL(cr)": aluminium}
    stream = reactant_gas(NASA9_REACTANTS, composition, 298.15, 101325.0, basis="mass")
    elements = element_moles_of(stream, products)

    products.set_element_moles(elements, 3000.0, PROPELLANT_PRESSURES[0])
    for pressure in PROPELLANT_PRESSURES:
        expected = reference[aluminium][pressure]
        products.equilibrate_HP(stream.enthalpy_mass, pressure)

        label = f"{aluminium:.2f} Al at {pressure / BAR:.4f} bar"
        assert_close_rel(products.temperature, expected["T"], 2e-3, f"temperature, {label}")
        assert present_condensed(products) == cea_condensed(expected), label
        assert_close_abs(condensed_mole_fraction(products, expected_polymorph),
                         expected["X"][expected_polymorph], 2e-3,
                         f"X[{expected_polymorph}], {label}")
        assert_close_rel(products.enthalpy_mass, stream.enthalpy_mass, 1e-7, f"enthalpy, {label}")


# ---------------------------------------------------------------------------
# (e) RP-1311 example 13: N2H4/Be with H2O2, whose beryllia passes two transitions
# ---------------------------------------------------------------------------

BERYLLIUM_REACTANTS = ["N2H4(L)", "Be(a)", "H2O2(L)"]
BERYLLIUM_CHAMBER_PRESSURE = 206.8419 * BAR
BERYLLIUM_PRESSURE_RATIOS = (3.0, 10.0, 30.0, 300.0)


@pytest.fixture(scope="module")
def beryllium_rocket():
    """Products, reactant stream and the CEA rocket solution for example 13.

    The product list is CEA's own, which includes gaseous Be(OH)2, and both codes get it.
    """
    names = run_isolated(_cea_product_names, BERYLLIUM_REACTANTS)
    products = goddard_products(names)

    # 67 % fuel by mass, the fuel being 80 % N2H4(L) and 20 % Be(a).
    weights, enthalpy = run_isolated(_cea_weights_and_enthalpy, BERYLLIUM_REACTANTS,
                                     [0.0, 0.0, 1.0], [0.8, 0.2, 0.0], 33.0 / 67.0,
                                     [298.15] * 3)
    stream = reactant_gas(NASA9_REACTANTS, dict(zip(BERYLLIUM_REACTANTS, weights)),
                          298.15, 101325.0, basis="mass")

    rocket = run_isolated(_cea_rocket, BERYLLIUM_REACTANTS, names, weights,
                          BERYLLIUM_CHAMBER_PRESSURE / BAR, BERYLLIUM_PRESSURE_RATIOS,
                          enthalpy / cea.R, trace=1e-10, insert=["BeO(L)"])
    assert rocket["converged"]

    products.set_element_moles(element_moles_of(stream, products), 3000.0,
                              BERYLLIUM_CHAMBER_PRESSURE)
    products.equilibrate_HP(stream.enthalpy_mass, BERYLLIUM_CHAMBER_PRESSURE)
    # The chamber entropy drives every expansion station below. It is captured here because
    # the Gas is shared by the tests of this group, which leave it at their own station.
    return products, stream, rocket, products.entropy_mass


def test_beryllium_chamber(beryllium_rocket):
    """The chamber of example 13, at constant enthalpy and pressure.

    Observed errors: temperature 7e-5 relative, M 2e-5 relative, BeO(L) 4e-6 absolute.
    """
    products, stream, rocket, _ = beryllium_rocket

    assert_close_rel(products.temperature, rocket["T"][0], 2e-3, "chamber temperature")
    assert_close_rel(products.mean_molecular_weight, rocket["M"][0], 1e-3, "chamber M")
    assert_close_abs(condensed_mole_fraction(products, "BeO(L)"),
                     rocket["X"]["BeO(L)"][0], 2e-3, "chamber X[BeO(L)]")
    assert not products.at_phase_transition
    assert_close_rel(products.enthalpy_mass, stream.enthalpy_mass, 1e-9, "chamber enthalpy")


def test_beryllium_pinned_stations(beryllium_rocket):
    """Stations 1 and 2 sit exactly on the BeO(b)/BeO(L) melting point.

    The expansion entropy falls inside the latent heat there, so the temperature stops moving
    and the two polymorphs split to match. Observed worst error on the split: 4.3e-4.
    """
    products, _, rocket, entropy = beryllium_rocket

    for station in (1, 2):
        pressure = rocket["P"][station] * BAR
        products.equilibrate_SP(entropy, pressure)

        assert products.at_phase_transition, f"station {station}"
        assert products.pinned_polymorphs != (-1, -1)
        assert products.temperature == pytest.approx(2851.0, abs=1e-6)
        for name in ("BeO(L)", "BeO(b)"):
            assert_close_abs(condensed_mole_fraction(products, name),
                             rocket["X"][name][station], 3e-3,
                             f"X[{name}] at station {station}")
        assert_close_rel(products.entropy_mass, entropy, 1e-8, f"entropy at station {station}")

    # The pin can also be cleared and re-declared by name, as a caller restoring a split would.
    products.clear_phase_transition()
    assert not products.at_phase_transition
    products.set_phase_transition("BeO(b)", "BeO(L)")
    assert products.at_phase_transition


@pytest.mark.parametrize("station", [3, 4, 5])
def test_beryllium_expansion_stations(beryllium_rocket, station):
    """Below the melting point the expansion follows BeO(b) and then BeO(a).

    Observed worst error on the temperature: 5e-5 relative.
    """
    products, _, rocket, entropy = beryllium_rocket
    pressure = rocket["P"][station] * BAR

    products.equilibrate_SP(entropy, pressure)

    assert_close_rel(products.temperature, rocket["T"][station], 2e-3,
                     f"temperature at station {station}")
    assert not products.at_phase_transition
    assert present_condensed(products) == cea_condensed(rocket, index=station)


# ---------------------------------------------------------------------------
# (f) RP-1311 example 8: H2(L)/O2(L) through RocketProblem
# ---------------------------------------------------------------------------

def test_cryogenic_rocket_chamber():
    """Build example 8 as a RocketProblem with reactant streams and check the chamber.

    Only the chamber is compared: the nozzle stations carry condensed-phase fields that the
    results structs do not expose yet (work package D).
    """
    reactant_names = ["H2(L)", "O2(L)"]
    of_ratio = 5.55157
    chamber_pressure = 53.3172 * BAR

    names = run_isolated(_cea_product_names, reactant_names)
    gas_names, _, missing = split_product_names(names)
    assert not missing, missing

    weights, enthalpy = run_isolated(_cea_weights_and_enthalpy, reactant_names,
                                     [0.0, 1.0], [1.0, 0.0], of_ratio, [20.27, 90.17])
    rocket = run_isolated(_cea_rocket, reactant_names, names, weights,
                          chamber_pressure / BAR, [1000.0], enthalpy / cea.R, supar=[5.0])
    assert rocket["converged"]

    chemistry = goddard.ChemicalParameters()
    chemistry.thermo_file = NASA9_GAS
    chemistry.species = set(gas_names)
    chemistry.reactant_file = NASA9_REACTANTS
    chemistry.cantera_fuel_state = goddard.PhaseSpecification(20.27, chamber_pressure,
                                                              {"H2(L)": 1.0})
    chemistry.cantera_oxidizer_state = goddard.PhaseSpecification(90.17, chamber_pressure,
                                                                  {"O2(L)": 1.0})
    chemistry.mixture_ratio_type = MixtureRatioType.OF_RATIO
    chemistry.OF_ratios = [of_ratio]

    case = goddard.RocketCaseParameters()
    case.name = "ex8"
    case.problem_type = "rocket"
    case.combustor_options = goddard.infinite_area_combustor([chamber_pressure])
    case.nozzle_options = goddard.equilibrium_nozzle(5.0)

    results = goddard.RocketProblem(chemistry, [case], "gas").solve()
    chamber = results.chamber(0, "ex8").thermo

    assert_close_rel(chamber.temperature, rocket["T"][0], 2e-3, "chamber temperature")
    assert_close_rel(chamber.molecular_weight, rocket["M"][0], 1e-3, "chamber M")
    assert_close_rel(chamber.pressure, chamber_pressure, 1e-6, "chamber pressure")


# ---------------------------------------------------------------------------
# (g) Finite differences on the equilibrium properties with a condensed phase present
# ---------------------------------------------------------------------------

def test_equilibrium_cp_matches_finite_differences(methane_oxygen):
    """cp of the Gordon & McBride system equals dh/dT along the equilibrium path.

    The state carries graphite, so this exercises the condensed rows of the derivative system
    from Python. Observed error: below 1e-6 relative.
    """
    products, _ = methane_oxygen
    temperature, pressure = 1500.0, 1.0 * BAR

    elements, _ = methane_oxygen_state(products, 0.5)
    products.set_element_moles(elements, temperature, pressure)
    products.equilibrate_TP(temperature, pressure)
    assert products.has_condensed_phases

    spec_heat_p = goddard.equilibrium_properties(products).spec_heat_p

    step = 1e-4 * temperature
    products.equilibrate_TP(temperature + step, pressure)
    enthalpy_high = products.enthalpy_mass
    products.equilibrate_TP(temperature - step, pressure)
    enthalpy_low = products.enthalpy_mass

    assert_close_rel(spec_heat_p, (enthalpy_high - enthalpy_low) / (2.0 * step), 1e-4, "cp")
