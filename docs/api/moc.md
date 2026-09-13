# Method of Characteristics

2-D/axisymmetric supersonic nozzle flow solver using the Method of
Characteristics. Supports design (minimum-length and Rao thrust-optimized) and
analysis modes.

Angles are in radians throughout, with one deliberate exception: the
[`NozzleProfile`](#goddard.NozzleProfile) contour generators take their shape
angles in degrees, and their parameters carry a `_deg` suffix to say so.

## Solver

::: goddard.MocNozzle

### Marching kernel

`solve()` advances the characteristic net with one of two internal kernels,
selected automatically from `MocOptions.mode` -- there is no user-facing switch.
`MocMode.DESIGN_MIN_LENGTH` marches characteristic chains: a front of chain
leading edges advances by pairing neighbours, and the contour is defined by the
characteristics it absorbs at the wall. `MocMode.ANALYSIS` (planar and
axisymmetric) and `MocMode.DESIGN_RAO` instead march reference-plane fronts:
every point of a front is prescribed and its two characteristics are traced back
to the previous front, so both characteristic families stay resolved at the same
density everywhere. In axisymmetric flow the chain-ladder's chain densities
diverge near the axis, and on a faceted contour its wall angle is quantized per
facet; the reference-plane march has neither problem.

::: goddard.MocFlowKind

::: goddard.MocMode

::: goddard.MocLogLevel

::: goddard.MocStartLine

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

::: goddard.CharacteristicFamily

::: goddard.ChainMetadata

::: goddard.ChainTermination

::: goddard.PointMembership

## Diagnostics

::: goddard.MocFailure

::: goddard.MocErrorCode

::: goddard.MocPassDiagnostics

::: goddard.MocStepLimiter

::: goddard.MocInitDiagnostics

::: goddard.MocFrontShear

::: goddard.summarize_front_shear

::: goddard.MocCrossings

::: goddard.find_like_characteristic_crossings

::: goddard.validate_moc_options

## Plotting

::: goddard.plotting
