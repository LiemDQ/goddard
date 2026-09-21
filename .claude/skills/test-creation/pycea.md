# pycea as a test oracle

The `cea` package (NASA's modernized CEA, Python bindings over Fortran) is the reference for Goddard's
CEA comparison tests. It lives in the `test` pixi environment alongside Cantera 3.2. That environment has
no scipy. Source: <https://github.com/nasa/cea> (`source/rocket.f90` for the rocket solver).

## Version

`pixi.toml` pins `cea >= 3.3.1`. Earlier versions decide frozen versus equilibrium once per group of
points (throat, the whole `pi_p` list, the whole `supar` list) from the group's first index. Points after
`n_frz` inside the same list then stay in equilibrium. Example: FAC, `n_frz=5`, `pi_p=[3, 1000]`: 3.1.1 gives
the equilibrium 1254.0 K at `pi_p=1000`, 3.3.4 gives the frozen 1043.2 K. Fixed in 3.3.1 (CHANGELOG, "Fixed").

## Isolate every call

Unknown species names, and CEA's "re-insertion likely to cause singular matrix" (TP exactly at a
polymorph transition temperature), end in a Fortran `STOP 1` that kills the interpreter. Call the oracle in
a forked child: see `run_isolated` in `python/tests/test_fac_cea.py`. Run scripts with `python -u` so output
printed before a crash is not lost.

## Call signatures

- `Mixture.of_ratio_to_weights(oxidant_weights, fuel_weights, OF)` returns **unnormalized** weights.
  `calc_property` is per kg regardless, but normalize before deriving element amounts per kg.
- `EqSolver.solve(solution, cea.HP, h0 / cea.R, P_bar, weights)`. TP takes `(T, P_bar)`; UV takes
  `(u / cea.R, v)`.
- `RocketSolver.solve(solution, weights, P_bar, pi_p, supar=..., hc=h / cea.R, iac=..., n_frz=...,
  ac_at=... | mdot=...)`. At least one `pi_p` value is required.
- Species-level thermo: `cea.Mixture([name]).calc_property(cea.ENTHALPY, [1.0], [T])`, in J/kg.
- Long explicit product lists can end in "Too many singular matrices". Use
  `Mixture(reactants, products_from_reactants=True, omit=[...])` instead.

## Station order

Solution arrays are 0-based, in this order:

| Combustor | Order |
|---|---|
| Infinite area (`iac=True`) | chamber, throat, `pi_p` exits, `supar` exits |
| Finite area (`iac=False`) | injector, inf, combustion end, throat, `pi_p` exits, `supar` exits |

For finite area, `pi_p` is P_inj/P, not P_chamber/P. Composition arrays are `solution.mass_fractions[name]`
and `solution.mole_fractions[name]`.

## Frozen chemistry (`n_frz`)

- `n_frz` is a 1-based point number in the order above, so `n_frz` = array index + 1. The freeze point is
  itself an equilibrium point; later points use its composition.
- Goddard's `frozen_NFZ` counts from the last equilibrium station before the nozzle: infinite area
  `n_frz = frozen_NFZ + 1`, finite area `n_frz = frozen_NFZ + 3`. Both map `frozen_NFZ=1` to the throat.
- An `n_frz` outside `1..num_pts` disables frozen mode with only an info-level log message; the whole run
  is silently equilibrium.
- CEA2 scheduling rule, all versions: when freezing after the throat, `pi_p` values below the freeze
  point's pressure ratio and `supar` values at or below its Ae/At are omitted, again with only an info log.
  `num_pts` shrinks and later indices shift. Setting `n_frz` to the last `pi_p` point therefore drops every
  `supar` point.

## Naming and phase selection

- CEA `H2O(cr)` is Cantera `H2O(s)`; otherwise NASA names match (`AL2O3(a)`, `C(gr)`).
- RP-1311 example 13 needs `insert=["BeO(L)"]` for CEA itself to select phases correctly.

## Cantera side

Cantera's `Mixture` has no `enthalpy` or `entropy` attributes. Sum `phase.enthalpy_mole * phase_moles`
over its phases instead.
