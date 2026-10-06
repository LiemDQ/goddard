"""Smoke test for an installed goddard package.

Runs from any working directory without the source tree: the data files must come from the
package itself and be found by bare name.
"""
import os
import re
import sys

import goddard

PRODUCT_SPECIES = {"H2", "O2", "H2O", "OH", "H", "O"}
CHAMBER_PRESSURE = 20.0e5  # Pa
OF_RATIO = 6.0


def check_version():
    version = goddard.__version__
    assert re.fullmatch(r"\d+\.\d+\.\d+", version), version
    assert version == goddard._core.__version__
    print(f"goddard {version}")


def check_data_dir():
    data_dir = goddard.data_dir
    assert data_dir.is_dir(), data_dir
    assert (data_dir / "nasa9_gas.yaml").is_file(), sorted(os.listdir(data_dir))
    print(f"data directory: {data_dir}")


def check_gas_by_bare_name():
    gas = goddard.Gas("nasa9_gas.yaml", "gas", goddard.GasChemistry.EQUILIBRIUM, PRODUCT_SPECIES)
    gas.set_state_TP(300.0, 101325.0)
    assert abs(gas.temperature - 300.0) < 1e-9, gas.temperature


def check_rocket_chamber():
    chemistry = goddard.ChemicalParameters()
    chemistry.thermo_file = "nasa9_gas.yaml"
    chemistry.species = PRODUCT_SPECIES
    chemistry.cantera_fuel_state = goddard.PhaseSpecification(300.0, 101325.0, "H2:1")
    chemistry.cantera_oxidizer_state = goddard.PhaseSpecification(300.0, 101325.0, "O2:1")
    chemistry.mixture_ratio_type = goddard.MixtureRatioType.OF_RATIO
    chemistry.OF_ratios = [OF_RATIO]

    case = goddard.RocketCaseParameters()
    case.name = "smoke"
    case.problem_type = "rocket"
    case.combustor_options = goddard.infinite_area_combustor([CHAMBER_PRESSURE])
    case.nozzle_options = goddard.equilibrium_nozzle(10.0)

    results = goddard.RocketProblem(chemistry, [case], "gas").solve()
    chamber = results.chamber(0, case_name="smoke").thermo
    print(f"H2/O2 chamber temperature at O/F {OF_RATIO}: {chamber.temperature:.1f} K")
    assert 2500.0 < chamber.temperature < 4000.0, chamber.temperature
    assert abs(chamber.pressure - CHAMBER_PRESSURE) < 1e-6 * CHAMBER_PRESSURE, chamber.pressure


def main():
    check_version()
    check_data_dir()
    check_gas_by_bare_name()
    check_rocket_chamber()
    print("smoke test passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
