# Thermodynamics

!!! note "Placeholder"
    This page is a placeholder. Theory content will be added as the documentation grows.

## Overview

Goddard uses Cantera as its thermodynamic backend. The `Gas` class wraps Cantera's
`Solution` object and adds compressible-flow calculations such as the isentropic
exponent $\gamma_s$, speed of sound, and stagnation properties.

## Frozen and equilibrium properties

A `Gas` holds a state, not a chemistry model. Derived properties depend on how the
composition responds when the state is perturbed, which `CompositionModel` selects:

- **FROZEN** -- the composition is held fixed: $\gamma_s = c_p / c_v$.
- **EQUILIBRIUM** -- the composition shifts to remain at chemical equilibrium.

In Python, `gamma_s` and `speed_of_sound` are the frozen values, and
`equilibrium_gamma_s` and `equilibrium_speed_of_sound` the equilibrium ones. Methods such
as `snapshot`, `expansion_properties` and `stagnation_pressure` take a `model` argument,
frozen by default. In C++ every one of them takes the model as a required argument.

## Solver chemistry

Flow solvers take a `GasChemistry` in their options or constructor:

- **PERFECT_GAS** -- constant ratio of specific heats
- **FROZEN** -- composition held fixed along the flow
- **EQUILIBRIUM** -- composition at chemical equilibrium at every state
- **KINETIC** -- finite-rate chemistry (used in kinetic nozzle calculations)
