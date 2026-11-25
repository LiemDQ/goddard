# Handling failures

Every solver in Goddard reports failure the same way.

- A solve that returns a result reports a numerical failure in the result's `status`. It does not raise, so a sweep keeps the points that did solve.
- The values that depend on a failed solve are `NaN`. A failed result cannot be mistaken for a solved one in arithmetic or in a plot.
- A function with no result to carry a status raises `goddard.SolveError`. These are the state-changing methods of `Gas` (`set_state_*`, `equilibrate_*`) and the few solves that return a single number.
- A bad argument raises `ValueError`, and an option that is not implemented raises `NotImplementedError`. These are mistakes in the call, not failures of the numerics, and they are raised before any solving starts.

## Checking a result

```python
import goddard as gd

shock = gd.ShockSolver(gas, gd.GasChemistry.FROZEN).normal_shock(0.8)

shock.status.ok          # False
shock.status.reason      # FailureReason.NO_SOLUTION
shock.status.message     # "no shock exists: the upstream normal velocity ..."
shock.mach_out           # nan
```

`status` is a [`SolveStatus`][goddard.SolveStatus]:

| Attribute | Meaning |
|---|---|
| `ok` | `True` when the solve succeeded |
| `reason` | A [`FailureReason`][goddard.FailureReason]: the kind of failure |
| `message` | What happened, with the values involved |
| `stage` | Which part of the solve failed, e.g. `"throat"` or `"station 3 (supersonic area ratio 40)"` |
| `iterations`, `tolerance`, `residual` | Data of the iteration that failed, when it reports them |
| `warnings` | Conditions that did not stop the solve but may affect the result |

`status.describe()` gives all of it on one line.

## What the reason tells you

| Reason | Meaning | What to do |
|---|---|---|
| `NOT_CONVERGED` | An iteration ran out of steps or stalled | Check the inputs named in the message; a smaller step or a looser tolerance may help |
| `NO_SOLUTION` | The inputs have no solution in this model: a detached oblique shock, a wave slower than the Chapman-Jouguet speed, condensed products in a gas-only solver | Treat the point as outside the model; retrying does not help |
| `OUT_OF_DATA_RANGE` | A state left the range of the thermodynamic data or of a table | Use data that covers the state, or stop the expansion earlier |
| `THERMO_FAILURE` | A thermodynamic state or chemical equilibrium could not be computed | Check that the state in the message is physical |
| `INCOMPLETE` | The solve stopped early; what was computed is valid | Raise the step limit named in the message |
| `NOT_COMPUTED` | The solve was skipped because a stage it depends on failed | Look at the stage named in the message |

## Sweeps

A result made of several parts has a status for the whole and one for each part. One failed part does not affect the others.

```python
results = gd.Nozzle(chamber).solve_ratios(
    gd.ExpansionType.SUPERSONIC_AREA_RATIO, [10.0, 40.0, 100.0])

results.status.ok                              # False if any station failed
[s.status.ok for s in results.expansions]      # one flag per station
[s.mach for s in results.expansions]           # nan where a station failed
```

A `RocketProblem` solves every operating point it can:

```python
results = problem.solve()

results.status.ok               # True when every station was solved
results.failed_stations()       # the stations that were not, each with its status
results.performance(0, exit_index=1).isp   # nan if that exit, its throat or its chamber failed
```

The grid a `Combustor` returns has a status per operating point:

```python
states = combustor.solve(300.0, 300.0, pressures, mixture_ratios)

states.ok()             # boolean array, same shape as states.temperature()
states.num_failed
states.status_at(states.flat_index(i, j))
```

A failed entry has no state: `states.get_state(loc)` raises `SolveError` with the entry's status instead of returning the state the entry held before the failed solve.

## Stopping at the first failure

A script that should not continue past a failed solve calls `raise_if_failed()`:

```python
results = problem.solve()
results.status.raise_if_failed()   # raises goddard.SolveError if anything failed
```

It raises [`ConvergenceError`][goddard.ConvergenceError] for `NOT_CONVERGED` and [`SolveError`][goddard.SolveError] for every other reason. `ConvergenceError` is a `SolveError`, and both carry the failure as `error.status`.

## Functions that raise

`Gas` methods change the state of the gas and return nothing, so a failure there raises:

```python
try:
    gas.equilibrate_HP(enthalpy, pressure)
except gd.SolveError as error:
    print(error.status.reason, error.status.message)
```

The same holds for solves that return a number, such as `ShockSolver.max_deflection`, and for the accessors that return a `Gas`: `ShockSolver.post_shock_state()` raises if the last solve failed, since there is then no state to return.

## Results that keep a partial solution

Two solvers keep what they computed before a failure:

- `KineticNozzle.solve` keeps the stations integrated up to the failing step. `status.reason` is `INCOMPLETE` if the step limit ran out and `NOT_CONVERGED` if the integration failed.
- `MocNozzle.solve` keeps the characteristic net marched so far. `MocResult.failure` adds the MoC-specific detail: an error code, the location and the marching pass. `MocResult.exit_mach`, `nozzle_length` and `area_ratio` are `NaN` and `exit_plane` is empty unless the solve succeeded, and `compute_thrust_coefficient` on a failed result returns `NaN` values with a `NOT_COMPUTED` status.
