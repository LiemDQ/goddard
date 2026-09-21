"""Compare Goddard's finite-area combustor (FAC) against NASA CEA (the `cea` Python module).

Both codes are given the same product species -- every gas species of `nasa9_gas.yaml`
whose only elements are H and O -- so the tolerances here follow the same-database
comparisons of `test_condensed_cea.py` rather than the looser cross-database ones of
`test_cea_comparison.py`.

The reference case is gaseous H2/O2 at 298.15 K, O/F 5.55157, P_inj 53.3172 bar --
RP-1311 example 8's propellants at a gaseous rather than cryogenic state, and burned
through a finite-area chamber instead of an infinite-area one.

Facts about `cea.RocketSolver` in FAC mode, verified directly against pycea before writing
these tests:

* The solution array holds points in the fixed order [injector, inf, comb end, throat,
  pi_p exits..., supar exits...] (`solution.num_pts` total), unlike the infinite-area order
  [chamber, throat, exits...].
* `pi_p` is P_inj/P in FAC mode (not P_chamber/P as for an infinite-area combustor).
* Two `pi_p` values are always requested here, not one. With a single `pi_p` value and
  `n_frz` set past the throat, the `supar` (area-ratio) stations silently disappear from
  the solution (`num_pts` drops and they are simply absent) -- a pycea/CEA quirk that two
  values avoids; every station is present with two or more `pi_p` values.
* Frozen station numbers: `n_frz=4` freezes at the throat (Goddard `frozen_NFZ=1`),
  `n_frz=5` freezes one station later, at the first `pi_p` exit (Goddard `frozen_NFZ=2`).

Tests are skipped automatically if the `cea` package is not installed.
"""
import multiprocessing as mp
import os
from functools import lru_cache

import numpy as np
import pytest

cea = pytest.importorskip("cea")

import goddard
from goddard import MixtureRatioType
from conftest import find_data_dir, assert_close_abs, assert_close_rel

DATA_DIR = find_data_dir()
NASA9_GAS = os.path.join(DATA_DIR, "nasa9_gas.yaml")

BAR = 1.0e5

REACTANT_NAMES = ["H2", "O2"]
REACTANT_TEMPERATURE = 298.15
OF_RATIO = 5.55157
INJECTOR_PRESSURE = 53.3172 * BAR
CONTRACTION_RATIO = 1.58
MASS_FLUX = 1333.9
PI_P = (3.0, 1000.0)
AREA_RATIOS = (10.0, 25.0)

# CEA FAC station indices for pi_p=(3.0, 1000.0), supar=(10.0, 25.0); verified empirically.
IDX_INJECTOR = 0
IDX_STAGNATION = 1
IDX_COMBUSTION_END = 2
IDX_THROAT = 3
IDX_PI_P = (4, 5)
IDX_SUPAR = (6, 7)


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


def _cea_fac_rocket(product_names, fac_kwargs):
    """Solve the H2/O2 FAC rocket problem with `cea.RocketSolver`, inside the forked child."""
    reactants = cea.Mixture(REACTANT_NAMES)
    products = cea.Mixture(product_names)
    solver = cea.RocketSolver(products, reactants=reactants)
    solution = cea.RocketSolution(solver)

    T_reactant = np.array([REACTANT_TEMPERATURE, REACTANT_TEMPERATURE])
    fuel_weights = np.array([1.0, 0.0])
    oxidant_weights = np.array([0.0, 1.0])
    weights = reactants.of_ratio_to_weights(oxidant_weights, fuel_weights, OF_RATIO)
    hc = reactants.calc_property(cea.ENTHALPY, weights, T_reactant) / cea.R

    solver.solve(solution, weights, INJECTOR_PRESSURE / BAR, list(PI_P),
                 supar=list(AREA_RATIOS), hc=hc, iac=False, **fac_kwargs)
    return {
        "T": [float(v) for v in solution.T],
        "P": [float(v) for v in solution.P],
        "M": [float(v) for v in solution.M],
        "Mach": [float(v) for v in solution.Mach],
        "gamma_s": [float(v) for v in solution.gamma_s],
        "density": [float(v) for v in solution.density],
        "ae_at": [float(v) for v in solution.ae_at],
        "c_star": [float(v) for v in solution.c_star],
        "Isp": [float(v) for v in solution.Isp],
        "Ivac": [float(v) for v in solution.Isp_vacuum],
        "CF": [float(v) for v in solution.coefficient_of_thrust],
        "Y": {name: [float(v) for v in values]
              for name, values in dict(solution.mass_fractions).items()},
        "num_pts": int(solution.num_pts),
        "converged": bool(solution.converged),
    }


# ---------------------------------------------------------------------------
# Building the matching Goddard problem
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def _hydrogen_oxygen_species():
    """Every gas species of nasa9_gas.yaml whose only elements are H and O."""
    import cantera as ct
    species_list = ct.Species.list_from_file(NASA9_GAS)
    return sorted(sp.name for sp in species_list if set(sp.composition) <= {"H", "O"})


PRODUCT_SPECIES = _hydrogen_oxygen_species()


def _build_problem(name, combustor_options, nozzle_options):
    """RocketProblem for gaseous H2/O2 through a finite-area combustor."""
    chemistry = goddard.ChemicalParameters()
    chemistry.thermo_file = NASA9_GAS
    chemistry.species = set(PRODUCT_SPECIES)
    chemistry.cantera_fuel_state = goddard.PhaseSpecification(
        REACTANT_TEMPERATURE, INJECTOR_PRESSURE, {"H2": 1.0})
    chemistry.cantera_oxidizer_state = goddard.PhaseSpecification(
        REACTANT_TEMPERATURE, INJECTOR_PRESSURE, {"O2": 1.0})
    chemistry.mixture_ratio_type = MixtureRatioType.OF_RATIO
    chemistry.OF_ratios = [OF_RATIO]

    case = goddard.RocketCaseParameters()
    case.name = name
    case.problem_type = "rocket"
    case.combustor_options = combustor_options
    case.nozzle_options = nozzle_options

    return goddard.RocketProblem(chemistry, [case], "gas")


# ---------------------------------------------------------------------------
# (a) Equilibrium chemistry: contraction-ratio and mass-flux modes
# ---------------------------------------------------------------------------

FAC_MODES = {
    "contraction_ratio": dict(
        combustor=lambda: goddard.finite_contraction_ratio_combustor(
            CONTRACTION_RATIO, [INJECTOR_PRESSURE]),
        cea_kwargs=dict(ac_at=CONTRACTION_RATIO),
    ),
    "mass_flux": dict(
        combustor=lambda: goddard.finite_mass_flux_combustor(
            MASS_FLUX, [INJECTOR_PRESSURE]),
        cea_kwargs=dict(mdot=MASS_FLUX),
    ),
}


@pytest.fixture(scope="module", params=list(FAC_MODES), ids=list(FAC_MODES))
def fac_rocket(request):
    """Equilibrium FAC rocket problem solved by both codes, for each FAC mode."""
    mode = FAC_MODES[request.param]
    rocket = run_isolated(_cea_fac_rocket, PRODUCT_SPECIES, mode["cea_kwargs"])
    assert rocket["converged"]

    problem = _build_problem("fac", mode["combustor"](), goddard.equilibrium_nozzle(*AREA_RATIOS))
    results = problem.solve()
    return request.param, results, rocket


def test_fac_injector(fac_rocket):
    """The injector-face station is `chamber()` for a finite-area combustor."""
    mode, results, rocket = fac_rocket
    chamber = results.chamber(0, "fac").thermo

    assert_close_rel(chamber.pressure / BAR, rocket["P"][IDX_INJECTOR], 1e-6, f"{mode} injector P")
    assert_close_rel(chamber.temperature, rocket["T"][IDX_INJECTOR], 2e-3, f"{mode} injector T")
    assert_close_rel(chamber.molecular_weight, rocket["M"][IDX_INJECTOR], 1e-3, f"{mode} injector M")


def test_fac_stagnation(fac_rocket):
    """The hypothetical stagnation state "inf" the nozzle expands from."""
    mode, results, rocket = fac_rocket
    stagnation = results.stagnation(0, "fac").thermo

    assert_close_rel(stagnation.pressure / BAR, rocket["P"][IDX_STAGNATION], 3e-3,
                     f"{mode} stagnation P")
    assert_close_rel(stagnation.temperature, rocket["T"][IDX_STAGNATION], 2e-3,
                     f"{mode} stagnation T")
    # Stagnation pressure is at most the injector pressure (some is lost to momentum).
    assert stagnation.pressure <= INJECTOR_PRESSURE


def test_fac_combustion_end(fac_rocket):
    """The subsonic end of the constant-area chamber, at area ratio contraction_ratio."""
    mode, results, rocket = fac_rocket
    combustion_end = results.combustion_end(0, "fac")
    thermo = combustion_end.thermo

    assert_close_rel(thermo.pressure / BAR, rocket["P"][IDX_COMBUSTION_END], 3e-3,
                     f"{mode} comb-end P")
    assert_close_rel(thermo.temperature, rocket["T"][IDX_COMBUSTION_END], 2e-3,
                     f"{mode} comb-end T")
    assert_close_rel(combustion_end.area_ratio, rocket["ae_at"][IDX_COMBUSTION_END], 3e-3,
                     f"{mode} comb-end Ac/At")


def test_fac_throat(fac_rocket):
    mode, results, rocket = fac_rocket
    throat = results.throat(0, "fac").thermo

    assert_close_rel(throat.pressure / BAR, rocket["P"][IDX_THROAT], 3e-3, f"{mode} throat P")
    assert_close_rel(throat.temperature, rocket["T"][IDX_THROAT], 2e-3, f"{mode} throat T")
    assert_close_rel(throat.molecular_weight, rocket["M"][IDX_THROAT], 1e-3, f"{mode} throat M")


@pytest.mark.parametrize("exit_index,area_ratio", list(enumerate(AREA_RATIOS)),
                         ids=[f"AR{int(ar)}" for ar in AREA_RATIOS])
def test_fac_exits(fac_rocket, exit_index, area_ratio):
    mode, results, rocket = fac_rocket
    cea_index = IDX_SUPAR[exit_index]
    exits = results.exits(0, "fac")
    station = exits[exit_index].thermo
    label = f"{mode} AR{area_ratio}"

    assert_close_rel(station.pressure / BAR, rocket["P"][cea_index], 3e-3, f"{label} P")
    assert_close_rel(station.temperature, rocket["T"][cea_index], 2e-3, f"{label} T")
    assert_close_rel(station.density, rocket["density"][cea_index], 3e-3, f"{label} density")
    for name in ("H2", "H2O", "O2"):
        assert_close_abs(station.composition.get(name, 0.0), rocket["Y"][name][cea_index],
                         5e-3, f"{label} Y[{name}]")

    performance = results.performance(0, exit_index, "fac")
    assert_close_rel(performance.area_ratio, rocket["ae_at"][cea_index], 3e-3,
                     f"{label} area_ratio")
    assert_close_rel(performance.cstar, rocket["c_star"][cea_index], 1e-2, f"{label} c*")
    assert_close_rel(performance.CF, rocket["CF"][cea_index], 1e-2, f"{label} CF")
    assert_close_rel(performance.isp, rocket["Isp"][cea_index], 1e-2, f"{label} Isp")
    assert_close_rel(performance.ivac, rocket["Ivac"][cea_index], 1e-2, f"{label} Ivac")


def test_fac_mass_flux_contraction_ratio():
    """In mass-flux mode, the contraction ratio is derived rather than given."""
    mode = FAC_MODES["mass_flux"]
    rocket = run_isolated(_cea_fac_rocket, PRODUCT_SPECIES, mode["cea_kwargs"])
    assert rocket["converged"]

    problem = _build_problem("fac_mdot", mode["combustor"](), goddard.equilibrium_nozzle(*AREA_RATIOS))
    results = problem.solve()

    combustion_end = results.combustion_end(0, "fac_mdot")
    assert_close_rel(combustion_end.area_ratio, rocket["ae_at"][IDX_COMBUSTION_END], 3e-3,
                     "derived Ac/At")


def test_fac_report_contains_comb_end(fac_rocket):
    _, results, _ = fac_rocket
    report = results.report()
    assert "COMB END" in report


# ---------------------------------------------------------------------------
# (b) Frozen chemistry: NFZ=1 (throat) and NFZ=2 (first exit) vs n_frz=4/5
# ---------------------------------------------------------------------------
# CEA orders its exits pi_p first, then supar, so n_frz=5 freezes at the first pi_p station. The
# Goddard problem therefore expands to the same pressure ratios, whose first exit is that station.

FROZEN_NFZ_TO_N_FRZ = {1: 4, 2: 5}


@pytest.fixture(scope="module", params=[1, 2], ids=["NFZ1", "NFZ2"])
def fac_frozen_rocket(request):
    frozen_NFZ = request.param
    rocket = run_isolated(_cea_fac_rocket, PRODUCT_SPECIES,
                          dict(ac_at=CONTRACTION_RATIO, n_frz=FROZEN_NFZ_TO_N_FRZ[frozen_NFZ]))
    assert rocket["converged"]

    combustor_options = goddard.finite_contraction_ratio_combustor(
        CONTRACTION_RATIO, [INJECTOR_PRESSURE])
    nozzle_options = goddard.pressure_ratio(*PI_P)
    nozzle_options.chemistry = goddard.GasChemistry.FROZEN
    nozzle_options.frozen_NFZ = frozen_NFZ
    results = _build_problem(f"fac_frozen_{frozen_NFZ}", combustor_options, nozzle_options).solve()
    return frozen_NFZ, results, rocket


def test_fac_frozen_throat_unaffected(fac_frozen_rocket):
    """Chamber, stagnation and throat don't depend on the freezing station."""
    frozen_NFZ, results, rocket = fac_frozen_rocket
    throat = results.throat(0, f"fac_frozen_{frozen_NFZ}").thermo

    assert_close_rel(throat.pressure / BAR, rocket["P"][IDX_THROAT], 3e-3, f"NFZ{frozen_NFZ} throat P")
    assert_close_rel(throat.temperature, rocket["T"][IDX_THROAT], 2e-3, f"NFZ{frozen_NFZ} throat T")


@pytest.mark.parametrize("exit_index,pressure_ratio", list(enumerate(PI_P)),
                         ids=[f"pi_p{int(ratio)}" for ratio in PI_P])
def test_fac_frozen_exits(fac_frozen_rocket, exit_index, pressure_ratio):
    frozen_NFZ, results, rocket = fac_frozen_rocket
    cea_index = IDX_PI_P[exit_index]
    exits = results.exits(0, f"fac_frozen_{frozen_NFZ}")
    station = exits[exit_index].thermo
    label = f"NFZ{frozen_NFZ} pi_p{pressure_ratio}"

    assert_close_rel(station.pressure / BAR, rocket["P"][cea_index], 3e-3, f"{label} P")
    assert_close_rel(station.temperature, rocket["T"][cea_index], 3e-3, f"{label} T")


# ---------------------------------------------------------------------------
# (c) PRESSURE_RATIO expansion: Goddard's ratios are P_inj/P for a finite-area combustor
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def fac_pressure_ratio_rocket():
    rocket = run_isolated(_cea_fac_rocket, PRODUCT_SPECIES, dict(ac_at=CONTRACTION_RATIO))
    assert rocket["converged"]

    combustor_options = goddard.finite_contraction_ratio_combustor(
        CONTRACTION_RATIO, [INJECTOR_PRESSURE])
    results = _build_problem("fac_pi_p", combustor_options, goddard.pressure_ratio(*PI_P)).solve()
    return results, rocket


@pytest.mark.parametrize("exit_index,pressure_ratio", list(enumerate(PI_P)),
                         ids=[f"pi_p{int(ratio)}" for ratio in PI_P])
def test_fac_pressure_ratio_exits(fac_pressure_ratio_rocket, exit_index, pressure_ratio):
    """A finite-area combustor's PRESSURE_RATIO expansion ratios are P_inj/P, like CEA's pi_p."""
    results, rocket = fac_pressure_ratio_rocket
    cea_index = IDX_PI_P[exit_index]
    exits = results.exits(0, "fac_pi_p")
    label = f"pi_p {pressure_ratio}"

    assert exits[exit_index].area_ratio == pressure_ratio  # the ratio requested, echoed back

    station = exits[exit_index].thermo
    assert_close_rel(station.pressure, INJECTOR_PRESSURE / pressure_ratio, 1e-3,
                     f"{label} P matches P_inj / ratio")
    assert_close_rel(station.pressure / BAR, rocket["P"][cea_index], 3e-3, f"{label} P")
    assert_close_rel(station.temperature, rocket["T"][cea_index], 2e-3, f"{label} T")

    # Performance is referenced to the stagnation state, so its pressure ratio is P_inf/P.
    stagnation_pressure = results.stagnation(0, "fac_pi_p").thermo.pressure
    performance = results.performance(0, exit_index, "fac_pi_p")
    assert_close_rel(performance.pressure_ratio,
                     pressure_ratio * stagnation_pressure / INJECTOR_PRESSURE, 1e-3,
                     f"{label} pressure_ratio")
