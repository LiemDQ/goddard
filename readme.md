# Goddard

Goddard is a software library for simulating chemically reactive compressible flows. It is primarily intended for analyzing flows with strongly coupled chemistry and thermodynamics, such as in hypersonic flows and rocket engines. 

At a basic level, Goddard provides similar functionality as NASA's [Chemical Equilibrium with Applications](https://www1.grc.nasa.gov/research-and-engineering/ceaweb/) (CEA), but with a modern API, more features, more extensibility and more flexibility in problem specification. Goddard is built on top of the [Cantera](https://cantera.org/) C++ API and provides very similar results as NASA CEA in applicable situations. However, if you prefer to use CEA, consider [RocketCEA](https://rocketcea.readthedocs.io/en/latest/) or [cea](https://nasa.github.io/cea/) for wrappers around the CEA Fortran code.

## Features

- **Combustion** — isobaric, recirculated, and isochoric combustion
- **Thermochemical physics** — chemical equilibrium and kinetics with self-consistent thermodynamics.
- **Compressible flow** — 1-D equilibrium, frozen, and kinetic nozzle expansions
- **Method of Characteristics** — automated nozzle contour design and analysis of 2D & axisymmetric supersonic flows.
- **Shock relations** — thermodynamically consistent normal, reflected, and oblique shocks
- **Thermodynamic data** — integrated with the NASA Glenn thermodynamic database, covering nearly 2,000 gaseous and condensed species.

## Quick start

### Installation

```bash
# Clone and install with pixi
git clone https://github.com/LiemDQ/goddard.git
cd goddard
pixi shell
pixi run configure
pixi run compile
pixi run pip-install
```

### Minimal example

```python
from goddard import Gas, gas_from_yaml

# Create a gas from a YAML thermodynamic data file
gas = gas_from_yaml("h2o2.yaml", phase_name="gas")
gas.set_state_TP(3000.0, 101325.0)

print(f"Temperature: {gas.temperature:.1f} K")
print(f"Speed of sound: {gas.speed_of_sound:.1f} m/s")
print(f"gamma_s: {gas.gamma_s:.4f}")
```

## Contact
Questions? Reach out by opening a discussion topic, or via email at [dq@liem.ca](mailto:dq@liem.ca). Keep in mind that I work on Goddard in my spare time, so updates may be sporadic.


## License
Goddard is subject to the terms of the GNU Lesser General Public License License, v. 2.1.

Copyright © 2026 Liem Dam-Quang
