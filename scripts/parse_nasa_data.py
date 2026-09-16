#!/usr/bin/env python3
"""One-off converter from CEA's ``thermo.inp`` (``data/nasa9.dat``) to Cantera YAML.

THIS IS A ONE-OFF GENERATOR, NOT A BUILD PIPELINE. The YAML files it produced are
the data of record and are committed to the repository; ``data/nasa9.dat`` is kept
for provenance only. Nothing in the build, the tests or the Python package runs this
script. It is kept here so the conversion can be audited and, if the source data ever
changes, repeated by hand.

Generated files (written into ``data/``):

* ``nasa9_gas.yaml``        gaseous product species (CEA phase flag 0) plus an
                            ``ideal-gas`` phase named ``gas``
* ``nasa9_condensed.yaml``  condensed product species (CEA phase flag != 0), each
                            carrying ``equation-of-state: {model: constant-volume,
                            density: 1.0e6}``
* ``nasa9_reactants.yaml``  the REACTANTS section plus the gaseous and condensed
                            products, so that any CEA reactant name (``H2(L)``,
                            ``RP-1``, ``AL(cr)``, ``CH4``, ``O2``, ...) resolves from
                            one file; no equation-of-state entries

Provenance of the committed data:
    source      data/nasa9.dat (CEA 2002 thermo.inp, dated 9/09/04)
    commit      ef6511352aa98964f62fdb18360cd933f37890a3
    generated   2026-09-15
    command     python scripts/parse_nasa_data.py --input data/nasa9.dat --outdir data

Determinism: species are emitted in source-file order, ``yaml.safe_dump`` is called
with ``sort_keys=False``, and floats are written with ``repr`` (shortest round-trip
representation, exact to the last bit). Re-running the script on the same input
reproduces the committed files byte for byte.

Source record layout (CEA thermo.inp, fixed columns, 0-based Python slices):

    record 1  [0:18]   species name          [18:80]  comment / reference
    record 2  [0:2]    number of T intervals [3:9]    reference-date code
              [10:50]  formula, 5 x (A2 symbol + F6.2 atoms)
              [50:52]  phase flag (0 = gas, >0 = condensed phase index)
              [52:65]  molecular weight      [65:80]  heat of formation [J/mol]
    record 3  [0:11]   T_min [K]             [11:22]  T_max [K]
              [22]     number of coefficients (always 7 in this file)
              [23:63]  8 x F5.1 temperature exponents (always -2 -1 0 1 2 3 4 0)
              [65:80]  H(298.15) - H(0) [J/mol] (unused by Cantera)
    record 4  5 x D16.9   a1..a5
    record 5  2 x D16.9   a6 a7, 16 blank columns, 2 x D16.9  b1 b2

A record with zero intervals is a CEA "assigned enthalpy" reactant: record 3 carries
only the assigned temperature and the species has a constant (zero) heat capacity.
"""

from __future__ import annotations

import argparse
import os
import sys

try:
    import yaml
except ImportError:  # the parsing half of this module is usable without PyYAML
    yaml = None

# ---------------------------------------------------------------------------
# Fixed-column record layout
# ---------------------------------------------------------------------------

NAME_COLUMNS = (0, 18)
COMMENT_COLUMNS = (18, 80)
NUM_INTERVALS_COLUMNS = (0, 2)
FORMULA_COLUMNS = (10, 50)
PHASE_FLAG_COLUMNS = (50, 52)
MOLECULAR_WEIGHT_COLUMNS = (52, 65)
HEAT_OF_FORMATION_COLUMNS = (65, 80)
INTERVAL_T_MIN_COLUMNS = (0, 11)
INTERVAL_T_MAX_COLUMNS = (11, 22)
NUM_COEFFICIENTS_COLUMN = 22
EXPONENT_COLUMNS = (23, 63)

NUM_FORMULA_SLOTS = 5
FORMULA_SLOT_WIDTH = 8
COEFFICIENT_WIDTH = 16
EXPECTED_EXPONENTS = (-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 0.0)

# CEA element symbols that Cantera cannot represent. Cantera's element table has no
# atomic weight for radon (it has no stable isotope) and refuses to build a phase
# containing it, so the two species that use it -- Rn and Rn+, both monatomic gases of
# no interest to propulsion -- are dropped from the generated files.
UNSUPPORTED_ELEMENTS = frozenset({"RN"})

# Density [kg/m^3] assigned to every condensed species. Gordon & McBride neglect the
# volume of condensed phases entirely; Cantera's StoichSubstance needs *some* density,
# and a large one makes the (P - P0) * M / (rho * R * T) term it adds to the chemical
# potential negligible: 2.8e-5 for AL2O3(L) at 70 bar and 3000 K, i.e. 0.7 J/mol.
CONDENSED_DENSITY = 1.0e6


class Nasa9Interval:
    """One NASA9 temperature interval: [T_min, T_max] and the nine coefficients."""

    def __init__(self, T_min: float, T_max: float, coefficients: list[float]) -> None:
        self.T_min = T_min
        self.T_max = T_max
        self.coefficients = coefficients


class SpeciesRecord:
    """One species as it appears in nasa9.dat, before any Cantera-specific mapping."""

    def __init__(
        self,
        name: str,
        comment: str,
        composition: list[tuple[str, float]],
        phase_flag: int,
        molecular_weight: float,
        heat_of_formation: float,
        intervals: list[Nasa9Interval],
        assigned_temperature: float | None,
        section: str,
    ) -> None:
        self.name = name
        self.comment = comment
        self.composition = composition
        self.phase_flag = phase_flag
        self.molecular_weight = molecular_weight
        self.heat_of_formation = heat_of_formation
        self.intervals = intervals
        self.assigned_temperature = assigned_temperature
        self.section = section

    @property
    def is_gas(self) -> bool:
        return self.phase_flag == 0

    @property
    def is_constant_cp(self) -> bool:
        return self.assigned_temperature is not None

    def __repr__(self) -> str:
        return f"SpeciesRecord({self.name!r}, phase={self.phase_flag}, section={self.section})"


# ---------------------------------------------------------------------------
# Element symbols
# ---------------------------------------------------------------------------


def canonical_element_symbols() -> dict[str, str]:
    """Map CEA's upper-case element symbols onto Cantera's spelling.

    The table is derived from Cantera's own element list when ``cantera`` can be
    imported; otherwise the equivalent title-case rule is applied. The two agree for
    every symbol appearing in nasa9.dat -- ``scripts/verify_nasa9_yaml.py`` asserts
    this in an environment where Cantera is available.
    """
    symbols: set[str] = set()
    try:
        import cantera as ct

        for element_name in ct.Element.element_names:
            try:
                symbols.add(ct.Element(element_name).symbol)
            except Exception:
                # Elements without a stable isotope (technetium, promethium, ...)
                # carry no weight in Cantera's table and cannot be instantiated.
                continue
        for extra in ("D", "E"):  # deuterium and the electron pseudo-element
            try:
                symbols.add(ct.Element(extra).symbol)
            except Exception:
                continue
    except ImportError:
        pass

    table = {symbol.upper(): symbol for symbol in symbols}

    def canonical(cea_symbol: str) -> str:
        return table.get(cea_symbol, cea_symbol[0].upper() + cea_symbol[1:].lower())

    return {symbol: canonical(symbol) for symbol in ALL_CEA_ELEMENT_SYMBOLS}


# Every element symbol that occurs in data/nasa9.dat, in the order the file uses.
ALL_CEA_ELEMENT_SYMBOLS = (
    "AG AL AR B BA BE BR C CA CD CL CO CR CS CU D E F FE GA GE H HE HG I IN K KR LI "
    "MG MN MO N NA NB NE NI O P PB RB RN S SC SI SN SR TA TH TI U V W XE ZN ZR"
).split()


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------


def parse_fortran_double(field: str) -> float:
    """Parse one D16.9 field, e.g. ``' 1.009950160D+04'``."""
    text = field.strip()
    if not text:
        return 0.0
    return float(text.replace("D", "E").replace("d", "e"))


def parse_formula(formula_field: str) -> list[tuple[str, float]]:
    """Parse the 40-column formula field into (CEA symbol, atom count) pairs."""
    composition: list[tuple[str, float]] = []
    for slot in range(NUM_FORMULA_SLOTS):
        begin = slot * FORMULA_SLOT_WIDTH
        entry = formula_field[begin : begin + FORMULA_SLOT_WIDTH]
        symbol = entry[0:2].strip()
        if not symbol:
            continue
        atoms = float(entry[2:8])
        if atoms == 0.0:
            continue
        composition.append((symbol, atoms))
    return composition


def parse_interval(lines: list[str], index: int) -> tuple[Nasa9Interval, int]:
    """Parse one three-line interval block starting at ``index``."""
    header = lines[index]
    T_min = float(header[slice(*INTERVAL_T_MIN_COLUMNS)])
    T_max = float(header[slice(*INTERVAL_T_MAX_COLUMNS)])

    num_coefficients = int(header[NUM_COEFFICIENTS_COLUMN])
    if num_coefficients != 7:
        raise ValueError(f"line {index + 1}: expected 7 coefficients, got {num_coefficients}")

    exponent_field = header[slice(*EXPONENT_COLUMNS)]
    exponents = tuple(float(exponent_field[5 * n : 5 * n + 5]) for n in range(8))
    if exponents != EXPECTED_EXPONENTS:
        raise ValueError(f"line {index + 1}: unexpected temperature exponents {exponents}")

    # a1..a5
    coefficients = [
        parse_fortran_double(lines[index + 1][n * COEFFICIENT_WIDTH : (n + 1) * COEFFICIENT_WIDTH])
        for n in range(5)
    ]
    # a6, a7, (16 blank columns), b1, b2
    third_line = lines[index + 2]
    for n in (0, 1, 3, 4):
        coefficients.append(
            parse_fortran_double(third_line[n * COEFFICIENT_WIDTH : (n + 1) * COEFFICIENT_WIDTH])
        )

    return Nasa9Interval(T_min, T_max, coefficients), index + 3


def parse_species(lines: list[str], index: int, section: str) -> tuple[SpeciesRecord, int]:
    """Parse one species record starting at ``index``; return it and the next index."""
    name = lines[index][slice(*NAME_COLUMNS)].strip()
    comment = lines[index][slice(*COMMENT_COLUMNS)].strip()

    record = lines[index + 1]
    num_intervals = int(record[slice(*NUM_INTERVALS_COLUMNS)])
    composition = parse_formula(record[slice(*FORMULA_COLUMNS)])
    phase_flag = int(record[slice(*PHASE_FLAG_COLUMNS)])
    molecular_weight = float(record[slice(*MOLECULAR_WEIGHT_COLUMNS)])
    heat_of_formation = float(record[slice(*HEAT_OF_FORMATION_COLUMNS)])

    index += 2
    intervals: list[Nasa9Interval] = []
    assigned_temperature: float | None = None
    if num_intervals == 0:
        # CEA "assigned enthalpy" reactant: one line holding the assigned temperature.
        assigned_temperature = float(lines[index][slice(*INTERVAL_T_MIN_COLUMNS)])
        index += 1
    else:
        for _ in range(num_intervals):
            interval, index = parse_interval(lines, index)
            intervals.append(interval)

    return (
        SpeciesRecord(
            name,
            comment,
            composition,
            phase_flag,
            molecular_weight,
            heat_of_formation,
            intervals,
            assigned_temperature,
            section,
        ),
        index,
    )


def parse_nasa9_file(path: str) -> tuple[list[SpeciesRecord], list[SpeciesRecord]]:
    """Parse the whole file; return (product records, reactant records)."""
    with open(path, "r") as source:
        lines = source.read().split("\n")

    index = 0
    while not lines[index].startswith("thermo"):
        index += 1
    index += 2  # skip the "thermo" keyword and the global temperature-range line

    products: list[SpeciesRecord] = []
    reactants: list[SpeciesRecord] = []
    section = "products"
    while index < len(lines):
        line = lines[index]
        if line.startswith("END PRODUCTS"):
            section = "reactants"
            index += 1
            continue
        if line.startswith("END REACTANTS") or not line.strip():
            break
        species, index = parse_species(lines, index, section)
        (products if section == "products" else reactants).append(species)

    return products, reactants


# ---------------------------------------------------------------------------
# Records that share a CEA name
# ---------------------------------------------------------------------------


def merge_duplicate_names(records: list[SpeciesRecord]) -> tuple[list[SpeciesRecord], list[str]]:
    """Merge records that share a CEA name into one multi-interval species.

    Ten condensed substances appear more than once in the PRODUCTS section under the
    same name -- CEA splits them at a lambda or Curie transition and distinguishes the
    pieces by the phase-flag column alone (``Cr(cr)`` 200-311.5 K and 311.5-1000 K,
    ``Fe(a)``, ``Fe2O3(cr)``, ...). Cantera species names must be unique, and the
    pieces cover contiguous, non-overlapping temperature ranges, so they are
    concatenated into a single NASA9 species. Evaluating the result at any temperature
    picks exactly the same interval CEA would.
    """
    merged: list[SpeciesRecord] = []
    by_name: dict[str, SpeciesRecord] = {}
    notes: list[str] = []

    for record in records:
        first = by_name.get(record.name)
        if first is None:
            by_name[record.name] = record
            merged.append(record)
            continue

        if first.is_constant_cp or record.is_constant_cp:
            raise ValueError(f"cannot merge constant-cp duplicate {record.name!r}")
        if first.composition != record.composition:
            raise ValueError(f"duplicate {record.name!r} with different composition")
        if first.intervals[-1].T_max != record.intervals[0].T_min:
            raise ValueError(
                f"duplicate {record.name!r} ranges are not contiguous: "
                f"{first.intervals[-1].T_max} != {record.intervals[0].T_min}"
            )

        notes.append(
            f"{record.name}: merged phase-flag {record.phase_flag} record "
            f"({record.intervals[0].T_min:g}-{record.intervals[-1].T_max:g} K) into the "
            f"phase-flag {first.phase_flag} record"
        )
        first.intervals.extend(record.intervals)
        first.comment = f"{first.comment} | {record.comment}"

    return merged, notes


# ---------------------------------------------------------------------------
# Cantera YAML nodes
# ---------------------------------------------------------------------------


def format_float(value: float) -> str:
    """Shortest representation that round-trips exactly (Python's ``repr``)."""
    return repr(float(value))


def species_node(
    record: SpeciesRecord, element_map: dict[str, str], with_equation_of_state: bool
) -> dict:
    """Build the Cantera YAML node for one species."""
    composition = {}
    for symbol, atoms in record.composition:
        canonical = element_map[symbol]
        composition[canonical] = composition.get(canonical, 0.0) + atoms

    node: dict = {"name": record.name, "composition": composition}

    if record.is_constant_cp:
        # CEA assigned-enthalpy reactant: h(T) = h_formation, cp = 0, s = 0.
        node["thermo"] = {
            "model": "constant-cp",
            "T0": record.assigned_temperature,
            "h0": f"{format_float(record.heat_of_formation)} J/mol",
            "s0": 0.0,
            "cp0": 0.0,
        }
    else:
        temperature_ranges = [record.intervals[0].T_min]
        temperature_ranges.extend(interval.T_max for interval in record.intervals)
        node["thermo"] = {
            "model": "NASA9",
            "temperature-ranges": temperature_ranges,
            "data": [list(interval.coefficients) for interval in record.intervals],
        }

    if with_equation_of_state:
        node["equation-of-state"] = {"model": "constant-volume", "density": CONDENSED_DENSITY}

    if record.comment:
        node["note"] = record.comment

    return node


GAS_PHASE_NODE = {
    "name": "gas",
    "thermo": "ideal-gas",
    "species": "all",
    "skip-undeclared-elements": True,
    "state": {"T": 300.0, "P": 101325.0},
}


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------


def _configure_dumper() -> type:
    """A SafeDumper that writes floats with full ``repr`` precision."""

    class ExactFloatDumper(yaml.SafeDumper):
        pass

    def represent_exact_float(dumper, value):
        if value != value:
            text = ".nan"
        elif value == float("inf"):
            text = ".inf"
        elif value == float("-inf"):
            text = "-.inf"
        else:
            text = repr(value)
            if "." not in text and "e" not in text and "n" not in text:
                text += ".0"
            elif "." not in text and "e" in text:
                text = text.replace("e", ".0e", 1)
        return dumper.represent_scalar("tag:yaml.org,2002:float", text)

    ExactFloatDumper.add_representer(float, represent_exact_float)
    return ExactFloatDumper


def write_yaml(path: str, header_lines: list[str], document: dict) -> None:
    if yaml is None:
        raise RuntimeError("PyYAML is required to write the YAML files")
    text = yaml.dump(
        document,
        Dumper=_configure_dumper(),
        default_flow_style=None,
        sort_keys=False,
        width=1_000_000,
    )
    with open(path, "w") as output:
        for line in header_lines:
            output.write(f"# {line}\n" if line else "#\n")
        output.write("\n")
        output.write(text)


PROVENANCE = [
    "Autogenerated by scripts/parse_nasa_data.py from data/nasa9.dat",
    "(CEA 2002 thermo.inp, dated 9/09/04) at commit",
    "ef6511352aa98964f62fdb18360cd933f37890a3 on 2026-09-15.",
    "",
    "This is a one-off conversion, not a build step: this file is the data of",
    "record and is edited only by re-running the generator. Species names are",
    "exactly as CEA writes them; element symbols use Cantera's spelling.",
    "",
    "The only species dropped in the conversion are Rn and Rn+: Cantera's element",
    "table has no atomic weight for radon, which has no stable isotope.",
]

CONDENSED_HEADER = [
    "",
    "Condensed (solid and liquid) product species.",
    "",
    "Every species carries equation-of-state {model: constant-volume, density:",
    "1.0e6 kg/m^3}. Gordon & McBride neglect the volume of condensed phases",
    "altogether: the NASA9 polynomials below give the standard-state properties",
    "at 1 bar and no pressure correction is intended. Cantera's StoichSubstance",
    "does add (P - P0) * M / (rho * R * T) to the chemical potential, so the",
    "density is chosen large enough to make that term negligible: 2.8e-5 for",
    "AL2O3(L) at 70 bar and 3000 K, which is 0.7 J/mol against fit uncertainties",
    "of kilojoules per mole. The value is a numerical placeholder and must not be",
    "read as a physical density.",
]

REACTANTS_HEADER = [
    "",
    "Every species that can name a CEA reactant, in one file: the REACTANTS",
    "section (H2(L), O2(L), RP-1, NH4CLO4(I), ...) followed by the gaseous",
    "products (CH4, O2, ...) and then the condensed products (AL(cr), C(gr),",
    "AL2O3(a), ...), which CEA also accepts as propellant ingredients.",
    "",
    "Records with no temperature intervals are CEA assigned-enthalpy reactants and",
    "become constant-cp species with cp0 = 0, so h(T) = h_formation exactly as CEA",
    "evaluates them.",
    "",
    "No species here carries an equation-of-state entry: a reactant stream only",
    "contributes an enthalpy and its element amounts, and is evaluated inside an",
    "ideal-gas phase node. Use data/nasa9_condensed.yaml for equilibrium products.",
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Convert CEA's nasa9.dat into Cantera YAML (one-off generator)."
    )
    parser.add_argument("--input", default="data/nasa9.dat", help="path to nasa9.dat")
    parser.add_argument("--outdir", default="data", help="directory for the YAML files")
    arguments = parser.parse_args(argv)

    element_map = canonical_element_symbols()
    products, reactants = parse_nasa9_file(arguments.input)

    def uses_unsupported_element(record: SpeciesRecord) -> bool:
        return any(symbol in UNSUPPORTED_ELEMENTS for symbol, _ in record.composition)

    skipped = [record.name for record in products + reactants if uses_unsupported_element(record)]
    products = [record for record in products if not uses_unsupported_element(record)]
    reactants = [record for record in reactants if not uses_unsupported_element(record)]

    products, merge_notes = merge_duplicate_names(products)
    reactants, reactant_merge_notes = merge_duplicate_names(reactants)
    merge_notes.extend(reactant_merge_notes)

    gas_products = [record for record in products if record.is_gas]
    condensed_products = [record for record in products if not record.is_gas]

    gas_document = {
        "phases": [GAS_PHASE_NODE],
        "species": [species_node(record, element_map, False) for record in gas_products],
    }
    condensed_document = {
        "species": [species_node(record, element_map, True) for record in condensed_products]
    }
    reactant_document = {
        "species": [
            species_node(record, element_map, False)
            for record in reactants + gas_products + condensed_products
        ]
    }

    write_yaml(
        os.path.join(arguments.outdir, "nasa9_gas.yaml"),
        PROVENANCE + ["", "Gaseous product species (CEA phase flag 0)."],
        gas_document,
    )
    write_yaml(
        os.path.join(arguments.outdir, "nasa9_condensed.yaml"),
        PROVENANCE + CONDENSED_HEADER,
        condensed_document,
    )
    write_yaml(
        os.path.join(arguments.outdir, "nasa9_reactants.yaml"),
        PROVENANCE + REACTANTS_HEADER,
        reactant_document,
    )

    print(f"gas products:       {len(gas_products)}")
    print(f"condensed products: {len(condensed_products)}")
    print(f"reactants:          {len(reactants)}")
    print(
        f"reactants file:     "
        f"{len(reactants) + len(gas_products) + len(condensed_products)} species"
    )
    print(f"skipped:            {skipped}")
    for note in merge_notes:
        print(f"merged: {note}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
