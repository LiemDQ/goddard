# Method of Characteristics

2-D/axisymmetric supersonic nozzle flow solver using the Method of
Characteristics. Supports design (minimum-length and Rao thrust-optimized) and
analysis modes.

Angles are in radians throughout, with one deliberate exception: the
[`NozzleProfile`](#goddard.NozzleProfile) contour generators take their shape
angles in degrees, and their parameters carry a `_deg` suffix to say so.

## Solver

::: goddard.MocNozzle

### Marching schemes

`MocOptions.march_scheme` selects which kernel `solve()` uses to advance the
characteristic net. **DIRECT** is the original chain-pairing kernel: it advances a
front of chain leading edges by pairing neighbours, and is used for
`MocMode.DESIGN_MIN_LENGTH` and (by default) planar flow, since its contour or
mesh is defined by the characteristics it absorbs. **INVERSE** instead prescribes
every point of a reference-plane front and traces its two characteristics back to
the previous front, so both characteristic families stay resolved at the same
density everywhere; it is the default for axisymmetric `ANALYSIS` and
`DESIGN_RAO`, where the DIRECT kernel's chain densities otherwise diverge near the
axis. See `instructions/moc_fix/diagnosis.md` for the full failure analysis behind
that default (not published with these docs, but present in the repository).

::: goddard.MocMarchScheme

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

::: goddard.ChainMetadata

::: goddard.PointMembership

## Diagnostics

::: goddard.MocFailure

::: goddard.MocPassDiagnostics

::: goddard.MocStepLimiter

::: goddard.MocCrossings

::: goddard.find_like_characteristic_crossings

::: goddard.validate_moc_options

## Plotting

::: goddard.plotting
