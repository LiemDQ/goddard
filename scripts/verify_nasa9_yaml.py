#!/usr/bin/env python3
"""Verification for the one-off ``nasa9.dat`` to Cantera YAML conversion (WP-E).

Run once, next to ``scripts/parse_nasa_data.py``, in an environment that has Cantera
and pycea:

    .pixi/envs/test/bin/python scripts/verify_nasa9_yaml.py

It checks four things:

1. Round trip. Every species in the three generated files is read back with
   ``ct.Species.list_from_file`` and compared field by field against the record parsed
   from ``data/nasa9.dat``: name, composition, temperature ranges and all nine NASA9
   coefficients per interval, with exact float equality. The one quantity that cannot
   be exact is the constant-cp ``h0``, which Cantera converts from J/mol to J/kmol;
   its deviation is reported in units in the last place.
2. pycea cross-check. Cantera and pycea read the same source data, so their enthalpy
   and heat capacity must agree to round-off once the two gas constants are divided
   out (Cantera 3.2 uses R = 8314.462618 J/kmol/K, CEA 2002 uses 8314.51), and once
   pycea's mass basis is converted with the molecular weight printed in nasa9.dat.
   The raw mass-basis difference is reported as well; it is dominated by those two
   constants and is about 1e-5.
3. Counts and loadability, including ``ct.Solution(nasa9_gas.yaml, "gas")``.
4. Spot checks on polymorph temperature ranges and the condensed equation of state.

pycea aborts the whole process with a Fortran ``STOP 1`` when it does not recognise a
species name, so every pycea call runs in a forked child process.
"""

from __future__ import annotations

import math
import multiprocessing
import os
import sys

import numpy as np
import cantera as ct

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(SCRIPT_DIR)
sys.path.insert(0, SCRIPT_DIR)

import parse_nasa_data as generator  # noqa: E402

SOURCE = os.path.join(REPO_ROOT, "data", "nasa9.dat")
DATA_DIR = os.path.join(REPO_ROOT, "data")
GAS_FILE = os.path.join(DATA_DIR, "nasa9_gas.yaml")
CONDENSED_FILE = os.path.join(DATA_DIR, "nasa9_condensed.yaml")
REACTANTS_FILE = os.path.join(DATA_DIR, "nasa9_reactants.yaml")


def section(title: str) -> None:
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def ulp_distance(actual: float, expected: float) -> float:
    if actual == expected:
        return 0.0
    scale = math.ulp(max(abs(actual), abs(expected)))
    return abs(actual - expected) / scale


# ---------------------------------------------------------------------------
# 0. Element symbols
# ---------------------------------------------------------------------------


def check_element_map() -> dict[str, str]:
    section("0. Element symbols")
    element_map = generator.canonical_element_symbols()

    title_case = {
        symbol: symbol[0].upper() + symbol[1:].lower()
        for symbol in generator.ALL_CEA_ELEMENT_SYMBOLS
    }
    disagreements = {s: (element_map[s], title_case[s]) for s in title_case if element_map[s] != title_case[s]}
    print(f"CEA element symbols in nasa9.dat: {len(element_map)}")
    print(f"Cantera-derived map vs. title-case fallback: {len(disagreements)} disagreements")

    unknown = []
    for cea_symbol, canonical in element_map.items():
        if cea_symbol in generator.UNSUPPORTED_ELEMENTS:
            continue
        try:
            ct.Element(canonical)
        except Exception:
            unknown.append((cea_symbol, canonical))
    print(f"symbols not in Cantera's element table: {unknown}")
    print(f"renamed by the conversion: "
          f"{sorted(s for s in element_map if element_map[s] != s)}")
    return element_map


# ---------------------------------------------------------------------------
# 1. Round trip
# ---------------------------------------------------------------------------


def expected_records() -> dict[str, list]:
    products, reactants = generator.parse_nasa9_file(SOURCE)

    def supported(record) -> bool:
        return not any(s in generator.UNSUPPORTED_ELEMENTS for s, _ in record.composition)

    products = [record for record in products if supported(record)]
    reactants = [record for record in reactants if supported(record)]
    products, _ = generator.merge_duplicate_names(products)
    reactants, _ = generator.merge_duplicate_names(reactants)

    gas = [record for record in products if record.is_gas]
    condensed = [record for record in products if not record.is_gas]
    return {
        GAS_FILE: gas,
        CONDENSED_FILE: condensed,
        REACTANTS_FILE: reactants + gas + condensed,
    }


def check_round_trip(element_map: dict[str, str]) -> None:
    section("1. Round trip: nasa9.dat vs. the generated YAML, exact float equality")
    total_species = 0
    total_values = 0
    mismatches: list[str] = []
    worst_h0_ulp = 0.0
    constant_cp_species = 0

    for path, records in expected_records().items():
        loaded = ct.Species.list_from_file(path)
        name = os.path.basename(path)
        if len(loaded) != len(records):
            mismatches.append(f"{name}: {len(loaded)} species loaded, {len(records)} expected")
            continue
        file_values = 0
        for record, species in zip(records, loaded):
            total_species += 1
            if species.name != record.name:
                mismatches.append(f"{name}: name {species.name!r} != {record.name!r}")
                continue
            file_values += 1

            expected_composition: dict[str, float] = {}
            for symbol, atoms in record.composition:
                canonical = element_map[symbol]
                expected_composition[canonical] = expected_composition.get(canonical, 0.0) + atoms
            if species.composition != expected_composition:
                mismatches.append(
                    f"{name}/{record.name}: composition {species.composition} "
                    f"!= {expected_composition}"
                )
            file_values += len(expected_composition)

            coefficients = list(species.thermo.coeffs)
            if record.is_constant_cp:
                constant_cp_species += 1
                if coefficients[0] != record.assigned_temperature:
                    mismatches.append(
                        f"{name}/{record.name}: T0 {coefficients[0]!r} "
                        f"!= {record.assigned_temperature!r}"
                    )
                if coefficients[2] != 0.0 or coefficients[3] != 0.0:
                    mismatches.append(f"{name}/{record.name}: s0/cp0 not zero")
                # h0 is the only value Cantera unit-converts (J/mol -> J/kmol).
                worst_h0_ulp = max(
                    worst_h0_ulp,
                    ulp_distance(coefficients[1], record.heat_of_formation * 1000.0),
                )
                file_values += 3
                continue

            num_regions = int(coefficients[0])
            if num_regions != len(record.intervals):
                mismatches.append(
                    f"{name}/{record.name}: {num_regions} intervals, "
                    f"{len(record.intervals)} expected"
                )
                continue
            for index, interval in enumerate(record.intervals):
                base = 1 + 11 * index
                if coefficients[base] != interval.T_min:
                    mismatches.append(
                        f"{name}/{record.name}[{index}]: T_min {coefficients[base]!r} "
                        f"!= {interval.T_min!r}"
                    )
                if coefficients[base + 1] != interval.T_max:
                    mismatches.append(
                        f"{name}/{record.name}[{index}]: T_max {coefficients[base + 1]!r} "
                        f"!= {interval.T_max!r}"
                    )
                for offset, expected in enumerate(interval.coefficients):
                    actual = coefficients[base + 2 + offset]
                    if actual != expected:
                        mismatches.append(
                            f"{name}/{record.name}[{index}] a{offset + 1}: "
                            f"{actual!r} != {expected!r}"
                        )
                file_values += 11
        print(f"{name:24s} {len(loaded):5d} species, {file_values:7d} values compared")
        total_values += file_values

    print()
    print(f"species compared:            {total_species}")
    print(f"values compared (exact ==):  {total_values}")
    print(f"mismatches:                  {len(mismatches)}")
    print(f"constant-cp species:         {constant_cp_species}")
    print(f"worst h0 deviation:          {worst_h0_ulp:g} ulp "
          f"(J/mol -> J/kmol conversion inside Cantera)")
    for problem in mismatches[:20]:
        print("  ! " + problem)


# ---------------------------------------------------------------------------
# 2. pycea cross-check
# ---------------------------------------------------------------------------

CROSS_CHECK_SPECIES = [
    # (name, class)
    ("H2", "gas"),
    ("O2", "gas"),
    ("N2", "gas"),
    ("CO", "gas"),
    ("CO2", "gas"),
    ("H2O", "gas"),
    ("OH", "gas"),
    ("NO", "gas"),
    ("H", "gas"),
    ("CH4", "gas"),
    ("C2H4", "gas"),
    ("NH3", "gas"),
    ("HCL", "gas"),
    ("ALCL2", "gas"),
    ("BeOH", "gas"),
    ("C(gr)", "condensed"),
    ("AL(cr)", "condensed"),
    ("AL(L)", "condensed"),
    ("AL2O3(a)", "condensed"),
    ("AL2O3(L)", "condensed"),
    ("BeO(a)", "condensed"),
    ("BeO(b)", "condensed"),
    ("BeO(L)", "condensed"),
    ("H2O(cr)", "condensed"),
    ("H2O(L)", "condensed"),
    ("MgO(cr)", "condensed"),
    ("SiO2(a-qz)", "condensed"),
    ("H2(L)", "reactant"),
    ("O2(L)", "reactant"),
    ("CH4(L)", "reactant"),
    ("N2H4(L)", "reactant"),
    ("RP-1", "reactant"),
    ("Jet-A(L)", "reactant"),
    ("NH4CLO4(I)", "reactant"),
]


def _cea_worker(name: str, temperatures: list[float], connection) -> None:
    import cea

    mixture = cea.Mixture([name])
    weights = np.array([1.0])
    enthalpies = []
    heat_capacities = []
    for temperature in temperatures:
        grid = np.array([float(temperature)])
        enthalpies.append(float(np.asarray(mixture.calc_property(cea.ENTHALPY, weights, grid)).ravel()[0]))
        heat_capacities.append(float(np.asarray(mixture.calc_property(cea.FROZEN_CP, weights, grid)).ravel()[0]))
    connection.send((enthalpies, heat_capacities))
    connection.close()


def cea_properties(name: str, temperatures: list[float]):
    """Enthalpy [J/kg] and frozen cp [J/kg/K] from pycea, in a forked child."""
    context = multiprocessing.get_context("fork")
    parent, child = context.Pipe(duplex=False)
    process = context.Process(target=_cea_worker, args=(name, temperatures, child))
    process.start()
    child.close()
    try:
        payload = parent.recv()
    except EOFError:
        payload = None
    process.join(60)
    return payload


def sample_temperatures(record) -> list[float]:
    if record.is_constant_cp:
        base = record.assigned_temperature
        return [base, base + 10.0, base + 50.0]
    T_min = record.intervals[0].T_min
    T_max = min(record.intervals[-1].T_max, 6000.0)
    span = T_max - T_min
    return [T_min + fraction * span for fraction in (0.17, 0.43, 0.81)]


def check_against_pycea(element_map: dict[str, str]) -> None:
    section("2. pycea cross-check (same source data)")
    import cea

    records = {}
    for group in expected_records().values():
        for record in group:
            records.setdefault(record.name, record)

    species_nodes = {}
    for path in (GAS_FILE, CONDENSED_FILE, REACTANTS_FILE):
        for species in ct.Species.list_from_file(path):
            species_nodes.setdefault(species.name, species)

    R_cantera = ct.gas_constant  # J/kmol/K
    R_cea = cea.R  # J/kmol/K
    print(f"R(Cantera) = {R_cantera!r} J/kmol/K, R(CEA) = {R_cea!r} J/kmol/K, "
          f"ratio - 1 = {R_cantera / R_cea - 1:.3e}")
    print()
    print("NASA9 species are compared as h/R and cp/R, which removes the two gas")
    print("constants; CEA assigned-enthalpy (constant-cp) species carry an absolute")
    print("enthalpy that R never multiplies, so they are compared per unit mass.")
    print()
    header = (f"{'species':14s} {'class':10s} {'T [K]':>9s} {'basis':>7s} "
              f"{'rel err h':>12s} {'rel err cp':>13s} {'raw rel err h':>14s}")
    print(header)
    print("-" * len(header))

    worst_h = 0.0
    worst_cp = 0.0
    worst_raw = 0.0
    failures = []
    for name, kind in CROSS_CHECK_SPECIES:
        record = records.get(name)
        species = species_nodes.get(name)
        if record is None or species is None:
            failures.append(f"{name}: not present in the generated files")
            continue
        temperatures = sample_temperatures(record)
        payload = cea_properties(name, temperatures)
        if payload is None:
            failures.append(f"{name}: pycea aborted (unknown name or out-of-range request)")
            continue
        enthalpies, heat_capacities = payload
        molecular_weight = record.molecular_weight  # kg/kmol, as printed in nasa9.dat
        scale = molecular_weight if record.is_constant_cp else R_cantera
        cea_scale = molecular_weight if record.is_constant_cp else R_cea
        basis = "h/M" if record.is_constant_cp else "h/R"
        for temperature, h_mass, cp_mass in zip(temperatures, enthalpies, heat_capacities):
            h_cantera = species.thermo.h(temperature) / scale
            cp_cantera = species.thermo.cp(temperature) / scale
            h_cea = h_mass if record.is_constant_cp else h_mass * molecular_weight / cea_scale
            cp_cea = cp_mass if record.is_constant_cp else cp_mass * molecular_weight / cea_scale
            error_h = abs(h_cantera - h_cea) / abs(h_cea) if h_cea else abs(h_cantera - h_cea)
            error_cp = abs(cp_cantera - cp_cea) / abs(cp_cea) if cp_cea else abs(cp_cantera - cp_cea)
            raw = species.thermo.h(temperature) / molecular_weight
            error_raw = abs(raw - h_mass) / abs(h_mass) if h_mass else abs(raw - h_mass)
            worst_h = max(worst_h, error_h)
            worst_cp = max(worst_cp, error_cp)
            worst_raw = max(worst_raw, error_raw)
            print(f"{name:14s} {kind:10s} {temperature:9.2f} {basis:>7s} "
                  f"{error_h:12.3e} {error_cp:13.3e} {error_raw:14.3e}")

    print()
    print(f"species compared:              {len(CROSS_CHECK_SPECIES) - len(failures)}")
    print(f"max relative error in h:       {worst_h:.3e}")
    print(f"max relative error in cp:      {worst_cp:.3e}")
    print(f"max raw mass-basis h error:    {worst_raw:.3e}  (gas constant + atomic weights)")
    for failure in failures:
        print("  ! " + failure)


# ---------------------------------------------------------------------------
# 3. Counts and loadability
# ---------------------------------------------------------------------------


def check_counts() -> None:
    section("3. Counts and loadability")
    for path in (GAS_FILE, CONDENSED_FILE, REACTANTS_FILE):
        loaded = ct.Species.list_from_file(path)
        size = os.path.getsize(path)
        print(f"{os.path.basename(path):24s} {len(loaded):5d} species  {size / 1024:8.1f} KiB")

    solution = ct.Solution(GAS_FILE, "gas")
    print(f"\nct.Solution('data/nasa9_gas.yaml', 'gas'): {solution.n_species} species, "
          f"{solution.n_elements} elements, T = {solution.T} K, P = {solution.P:.1f} Pa")
    print(f"elements: {' '.join(sorted(solution.element_names))}")

    reactant_phase = ct.Solution(
        thermo="ideal-gas",
        species=[
            species
            for species in ct.Species.list_from_file(REACTANTS_FILE)
            if species.name in ("H2(L)", "O2(L)", "RP-1", "AL(cr)", "CH4", "O2")
        ],
    )
    reactant_phase.TPX = 20.27, 101325.0, "H2(L):1"
    print(f"reactant phase from nasa9_reactants.yaml: {reactant_phase.n_species} species, "
          f"h(H2(L), 20.27 K) = {reactant_phase.enthalpy_mole:.6e} J/kmol")


# ---------------------------------------------------------------------------
# 4. Spot checks
# ---------------------------------------------------------------------------

SPOT_RANGES = [
    ("C(gr)", CONDENSED_FILE, [200.0, 600.0, 2000.0, 6000.0]),
    ("AL2O3(a)", CONDENSED_FILE, [200.0, 500.0, 1200.0, 2327.0]),
    ("AL2O3(L)", CONDENSED_FILE, [2327.0, 6000.0]),
    ("BeO(a)", CONDENSED_FILE, [200.0, 1000.0, 2373.0]),
    ("BeO(b)", CONDENSED_FILE, [2373.0, 2851.0]),
    ("BeO(L)", CONDENSED_FILE, [2851.0, 6000.0]),
    ("H2O(cr)", CONDENSED_FILE, [200.0, 273.15]),
    ("H2O(L)", CONDENSED_FILE, [273.15, 373.15, 600.0]),
]


def check_spot_values() -> None:
    section("4. Spot checks: polymorph ranges and the condensed equation of state")
    condensed = {species.name: species for species in ct.Species.list_from_file(CONDENSED_FILE)}

    for name, path, expected in SPOT_RANGES:
        species = {s.name: s for s in ct.Species.list_from_file(path)}[name]
        coefficients = list(species.thermo.coeffs)
        num_regions = int(coefficients[0])
        edges = [coefficients[1]]
        edges.extend(coefficients[1 + 11 * index + 1] for index in range(num_regions))
        status = "ok" if edges == expected else "MISMATCH"
        print(f"{name:12s} temperature-ranges = {edges}  {status}")

    print()
    print("polymorph continuity (same substance, adjacent ranges):")
    for lower, upper in [("AL2O3(a)", "AL2O3(L)"), ("BeO(a)", "BeO(b)"),
                         ("BeO(b)", "BeO(L)"), ("H2O(cr)", "H2O(L)")]:
        boundary = condensed[lower].thermo.max_temp
        status = "ok" if boundary == condensed[upper].thermo.min_temp else "MISMATCH"
        print(f"  {lower:10s} -> {upper:10s} at {boundary} K  {status}")

    print()
    missing_eos = []
    wrong_density = []
    for name, species in condensed.items():
        phase = ct.Solution(thermo="fixed-stoichiometry", species=[species])
        phase.TP = max(species.thermo.min_temp, 300.0), 101325.0
        if phase.density != generator.CONDENSED_DENSITY:
            wrong_density.append((name, phase.density))
    with open(CONDENSED_FILE) as handle:
        eos_lines = sum(
            1 for line in handle
            if line.strip() == "equation-of-state: {model: constant-volume, density: 1000000.0}"
        )
    print(f"condensed species with the constant-volume EOS line: {eos_lines} / {len(condensed)}")
    print(f"species whose phase density != {generator.CONDENSED_DENSITY:g}: {wrong_density}")
    print(f"species missing an EOS: {missing_eos}")

    reference = ct.Solution(thermo="fixed-stoichiometry", species=[condensed["AL2O3(L)"]])
    reference.TP = 3000.0, 70e5
    correction = (reference.P - ct.one_atm) / reference.density / (ct.gas_constant / reference.mean_molecular_weight * reference.T)
    print(f"(P - P0) V / (R T) for AL2O3(L) at 70 bar, 3000 K: {correction:.2e}")


def main() -> int:
    element_map = check_element_map()
    check_round_trip(element_map)
    check_counts()
    check_spot_values()
    check_against_pycea(element_map)
    return 0


if __name__ == "__main__":
    sys.exit(main())
