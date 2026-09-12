#include <cmath>
#include "goddard/moc_unit_processes.hpp"

namespace Goddard {

MocErrorCode check_point_validity(const CharacteristicPoint& pt, double tol,
                                  bool require_nonnegative_theta) {
    if (pt.nu < -tol) return MocErrorCode::NEGATIVE_NU;
    if (require_nonnegative_theta && pt.theta < -tol) return MocErrorCode::NEGATIVE_THETA;
    if (pt.mach < 1.0) return MocErrorCode::SUBSONIC_MACH;
    if (!std::isfinite(pt.x) || !std::isfinite(pt.y) || !std::isfinite(pt.theta) ||
        !std::isfinite(pt.nu) || !std::isfinite(pt.mach) || !std::isfinite(pt.mu)) {
        return MocErrorCode::NONFINITE_VALUE;
    }
    return MocErrorCode::NONE;
}

} // namespace Goddard
