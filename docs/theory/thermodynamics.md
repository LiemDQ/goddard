# Thermodynamics

!!! note "Placeholder"
    This page is a placeholder. Theory content will be added as the documentation grows.

## Overview

Goddard uses Cantera as its thermodynamic backend. The `Gas` class provides a
chemistry-aware interface to Cantera's `Solution` object, adding compressible-flow
calculations such as the isentropic exponent $\gamma_s$, speed of sound, and
stagnation properties.

## Chemistry modes

Goddard supports several chemistry modes that control how derived thermodynamic
properties are computed:

- **PERFECT_GAS** -- uses the simple $c_p / c_v$ ratio
- **FROZEN** -- holds composition fixed during isentropic calculations
- **EQUILIBRIUM** -- allows composition to shift to maintain chemical equilibrium
- **KINETIC** -- finite-rate chemistry (used in kinetic nozzle calculations)
