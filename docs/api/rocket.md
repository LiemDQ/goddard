# Rocket

CEA-style rocket problems: combustion followed by an isentropic nozzle expansion, solved for a
list of operating points.

A `RocketPoint` holds everything one operating point needs: the reactants and their mixture
ratio, the chamber pressure, the combustor and nozzle options, and the exit stations. A
`RocketSolver` built on the product `Gas` solves one point or any iterable of them, and returns
one `RocketResult` per point, in order.

```python
import itertools
import goddard as gd

products = gd.Gas.from_elements("nasa9_gas.yaml", "products", ["H", "O"])
solver = gd.RocketSolver(products, reactant_file="nasa9_reactants.yaml")

lh2 = gd.PhaseSpecification(20.27, 101325.0, "H2(L):1")
lox = gd.PhaseSpecification(90.17, 101325.0, "O2(L):1")
points = [gd.RocketPoint(fuel=lh2, oxidizer=lox, mixture_ratio=of_ratio,
                         chamber_pressure=pressure, expansion_ratios=[40.0])
          for of_ratio, pressure in itertools.product([4.0, 5.0, 6.0], [50e5, 100e5])]

results = solver.solve(points)
results.performance("ivac_s")   # array of shape (6 points, 1 exit)
results[0].chamber.thermo.temperature
print(results.report())
```

::: goddard.RocketSolver

::: goddard.RocketPoint

::: goddard.RocketResult

::: goddard.RocketResults

::: goddard.RocketPerformance

::: goddard.rocket_performance

::: goddard.rocket_report
