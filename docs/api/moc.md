# Method of Characteristics

2-D/axisymmetric supersonic nozzle flow solver using the Method of
Characteristics. Supports design (minimum-length and Rao thrust-optimized) and
analysis modes.

Angles are in radians throughout, with one deliberate exception: the
[`NozzleProfile`](#goddard.NozzleProfile) contour generators take their shape
angles in degrees, and their parameters carry a `_deg` suffix to say so.

## Solver

::: goddard.MocNozzle

::: goddard.MocOptions

::: goddard.NozzleGeometry

## Geometry

::: goddard.NozzleProfile

## Results

::: goddard.MocResult

::: goddard.ExitPlane

::: goddard.ThrustCoefficient

::: goddard.compute_thrust_coefficient

## Flow field

::: goddard.CharacteristicNet

::: goddard.CharacteristicPoint

::: goddard.ChainMetadata

::: goddard.PointMembership

## Diagnostics

::: goddard.MocFailure

::: goddard.MocPassDiagnostics

::: goddard.MocCrossings

::: goddard.find_like_characteristic_crossings

::: goddard.validate_moc_options

## Plotting

::: goddard.plotting
