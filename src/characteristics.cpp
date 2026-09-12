#include <cmath>
#include <utility>
#include "goddard/characteristics.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/error.hpp"
namespace Goddard {


void CharacteristicPoint::update_Ks() {
    K_plus = theta - nu;
    K_minus = theta + nu;
}

std::string_view to_string(MocErrorCode code) {
    switch (code) {
        case MocErrorCode::NONE: return "NONE";
        case MocErrorCode::NEGATIVE_NU: return "NEGATIVE_NU";
        case MocErrorCode::NEGATIVE_THETA: return "NEGATIVE_THETA";
        case MocErrorCode::SUBSONIC_MACH: return "SUBSONIC_MACH";
        case MocErrorCode::NONFINITE_VALUE: return "NONFINITE_VALUE";
        case MocErrorCode::PM_INVERSION_FAILED: return "PM_INVERSION_FAILED";
        case MocErrorCode::TABLE_RANGE_EXCEEDED: return "TABLE_RANGE_EXCEEDED";
        case MocErrorCode::NON_DOWNSTREAM_POINT: return "NON_DOWNSTREAM_POINT";
        case MocErrorCode::WALL_QUERY_OUT_OF_BOUNDS: return "WALL_QUERY_OUT_OF_BOUNDS";
        case MocErrorCode::INITIALIZATION_FAILED: return "INITIALIZATION_FAILED";
        case MocErrorCode::MAX_ITERATIONS_REACHED: return "MAX_ITERATIONS_REACHED";
        default: return "UNKNOWN";
    }
}

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

void characteristic_isentropic_PT_from_parent(CharacteristicPoint& point, const CharacteristicPoint& parent){
    double parent_stagnation_factor = stagnation_factor(parent.mach, parent.gamma_s);
    double current_stagnation_factor = stagnation_factor(point.mach, point.gamma_s);
    double ratio = parent_stagnation_factor/current_stagnation_factor;

    double average_gamma = 0.5*(parent.gamma_s + point.gamma_s);
    
    point.temperature = parent.temperature * ratio;
    point.pressure = parent.pressure * pow(ratio, average_gamma / (average_gamma - 1.0));
}

std::pair<double, double> characteristic_intersection_with_angle(
    const CharacteristicPoint& p1, 
    const CharacteristicPoint& p2)
{
    double angle1 = p1.theta - p1.mu;
    double angle2 = p2.theta + p2.mu;

    return characteristic_intersection_with_angle(p1, p2, angle1, angle2);
}

std::pair<double, double> characteristic_intersection_with_angle(
    const CharacteristicPoint& p1, 
    const CharacteristicPoint& p2, 
    double angle1,
    double angle2)
{
    double x = (p1.x*tan(angle1) - p2.x*tan(angle2) + p2.y - p1.y)/(tan(angle1) - tan(angle2));
    double y = (x  - p2.x) * tan(angle2) + p2.y;

    return {x,y};
} 

// -- CharacteristicNet --


} // namespace Goddard