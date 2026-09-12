#pragma once
#include "goddard/moc.hpp"
#include "goddard/characteristics.hpp"

namespace Goddard {

/**
 * A CharacteristicPoint together with the error code (if any) produced while
 * computing it. Unit-process solvers return this instead of a bare
 * CharacteristicPoint so a numerical failure is reported as data -- instead of
 * being laundered into the point's fields (e.g. a sentinel Mach number) or
 * thrown as an exception.
 */
struct PointResult {
    CharacteristicPoint point;
    MocErrorCode error = MocErrorCode::NONE;
};

/**
 * Validity of a computed point: finite fields, supersonic, non-negative Prandtl-Meyer angle,
 * and (by default) non-negative flow angle.
 *
 * The flow-angle requirement is what the DIRECT kernel uses to catch a march that has gone
 * wrong. It is optional because a negative flow angle is physically admissible in the
 * analysis of an arbitrary contour (any locally converging wall) and because, wherever the
 * exact angle is zero -- the uniform exit region of a minimum-length nozzle, the axis
 * neighbourhood -- a computed angle is discretization noise of either sign, which a
 * tolerance of the solver's abstol cannot absorb.
 *
 * @param pt  Point to validate.
 * @param tol Absolute tolerance for the nu/theta non-negativity checks (axis points
 *            legitimately pin theta to exactly 0.0, and interior points can carry
 *            small negative roundoff).
 * @return MocErrorCode::NONE if the point is valid, else the first violated
 *         condition in priority order: nu, theta, mach, finiteness.
 */
MocErrorCode check_point_validity(const CharacteristicPoint& pt, double tol,
                                  bool require_nonnegative_theta = true);

} // namespace Goddard
