# Kinetic Nozzle

1-D supersonic nozzle with finite-rate chemistry via Cantera's
`IdealGasMoleReactor`. Produces spatially-resolved profiles of composition,
temperature, and per-species Damkohler numbers.

The wall profile is given in throat radii (radius 1 at the throat), as for the
method of characteristics, and `throat_radius` [m] sets the physical scale. The
reactor integration is in SI units. Station positions and `dx_max` are in throat
radii.

::: goddard.KineticNozzle

::: goddard.KineticNozzleResults

::: goddard.KineticNozzleStation
