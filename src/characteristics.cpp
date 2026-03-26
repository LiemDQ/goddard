#include <cmath>
#include "goddard/characteristics.hpp"
#include "goddard/gas_dynamics.hpp"
namespace Goddard {

void characteristic_isentropic_PT_from_parent(CharacteristicPoint& point, const CharacteristicPoint& parent){
    double parent_stagnation_factor = stagnation_factor(parent.mach, parent.gamma_s);
    double current_stagnation_factor = stagnation_factor(point.mach, point.gamma_s);
    double ratio = parent_stagnation_factor/current_stagnation_factor;

    double average_gamma = 0.5*(parent.gamma_s + point.gamma_s);
    
    point.temperature = parent.temperature * ratio;
    point.pressure = parent.pressure * pow(ratio, average_gamma / (average_gamma - 1.0));
}

std::pair<double, double> characteristic_intersection_coordinates(
    const CharacteristicPoint& p1, 
    const CharacteristicPoint& p2, 
    double angle1,
    double angle2)
{
    double x = (p1.x*tan(angle1) - p2.x*tan(angle2) + p2.y - p1.y)/(tan(angle1) - tan(angle2));
    double y = (x  - p2.x) * tan(angle2) + p2.y;

    return {x,y};
} 
    
int CharacteristicNet::row_offset(int j) const {
    // row 0 starts at 0, row 1 at N, row 2 at (N-1), etc.
    return j * num_c_plus - j * (j -1) / 2;
}

} // namespace Goddard