"""H2/O2 rocket performance over a sweep of mixture ratios and chamber pressures."""
import itertools

import goddard as gd

bar = 1e5

# The product gas fixes the possible combustion products: here every species of the ohmech
# phase of h2o2.yaml.
products = gd.Gas("h2o2.yaml", "ohmech")
solver = gd.RocketSolver(products)

fuel = gd.PhaseSpecification(273.15, 101325.0, "H2:1")
oxidizer = gd.PhaseSpecification(273.15, 101325.0, "O2:1")

# One point per operating condition; the sweep can be any iterable of points.
points = [
    gd.RocketPoint(fuel=fuel, oxidizer=oxidizer, mixture_ratio=of_ratio,
                   chamber_pressure=chamber_pressure,
                   expansion_ratios=[2.0, 5.0, 10.0, 50.0],
                   label=f"O/F {of_ratio:g}, {chamber_pressure / bar:g} bar")
    for of_ratio, chamber_pressure in itertools.product([3.0, 6.0, 9.0], [10 * bar, 50 * bar])
]
results = solver.solve(points)

# One row per point, one column per exit.
print(results.performance("ivac_s"))
print(results.report())
