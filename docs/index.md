# Goddard

Goddard is a C++/Python rocket engine simulation toolkit. It wraps [Cantera](https://cantera.org/) for chemical kinetics and thermodynamics and provides higher-level solvers for:

- **Combustion** -- isobaric combustion with infinite-area, finite mass-flux, and finite contraction-ratio modes
- **Nozzle flow** -- 1-D equilibrium, frozen, and kinetic nozzle expansions
- **Method of Characteristics** -- 2-D/axisymmetric supersonic flow fields for nozzle design and analysis
- **Shock relations** -- normal, reflected, and oblique shocks

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

## Documentation sections

| Section | Description |
|---------|-------------|
| [Theory](theory/thermodynamics.md) | Background on the thermodynamics and flow physics |
| [Tutorials](tutorials/getting_started.md) | Step-by-step guides |
| [Examples](examples/basic_rocket.ipynb) | Jupyter notebooks |
| [API Reference](api/gas.md) | Auto-generated Python API docs |
