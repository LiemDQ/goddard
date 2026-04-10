# Getting Started

!!! note "Placeholder"
    This tutorial will walk through a first calculation with Goddard.

## Prerequisites

- Python >= 3.11
- [Pixi](https://pixi.sh) package manager

## Installation

```bash
git clone https://github.com/LiemDQ/goddard.git
cd goddard
pixi shell
pixi run configure
pixi run compile
pixi run pip-install
```

## Your first calculation

```python
from goddard import gas_from_yaml

gas = gas_from_yaml("h2o2.yaml", phase_name="gas")
gas.set_state_TP(3000.0, 101325.0)
print(gas.temperature, gas.pressure)
```
