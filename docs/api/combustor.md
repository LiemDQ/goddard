# Combustor

Isobaric combustion with support for infinite-area, finite mass-flux,
and finite contraction-ratio modes.

::: goddard.Combustor

::: goddard.CombustorOptions

## Results

`Combustor.solve` returns a `ThermoArray`: one product state per operating point, on a grid
of shape (mixture ratios, pressures), or (temperatures, pressures, mixture ratios) for the
overload that takes reactant temperatures. Every scalar property of `Gas` is an array of that
shape, the compositions add a species axis, and indexing gives a `Gas` at one entry:

```python
states = combustor.solve(300.0, 300.0, pressures, mixture_ratios)

states.temperature                  # shape (len(mixture_ratios), len(pressures))
states.equilibrium_gamma_s
states.mole_fractions[..., states.species_names.index("H2O")]
gas = states[i, j]                  # independent Gas in the state of entry (i, j)
states.snapshot((i, j), model=goddard.CompositionModel.EQUILIBRIUM)
states.evaluate(lambda gas: gas.cstar(goddard.CompositionModel.EQUILIBRIUM))
```

Properties are those of the whole mixture, condensed phases included, as `Gas` reports them.
A failed entry is NaN in every property array (see
[Handling failures](../tutorials/handling_failures.md)).

::: goddard.ThermoArray
