#include "goddard/moc_initialization.hpp"

namespace Goddard {

std::vector<CharacteristicPoint> MocInitialization::initialize_sauer(const ThroatGeometry& geometry, const ThroatCondition& throat) {
    
    CharacteristicPoint throat_point;
    throat_point.y = 1.0;
    throat_point.x = 0.0;
    throat_point.theta = 0.0;

    double gamma = throat.gamma_s;
    double alpha = sauer_alpha(gamma, geometry);
    size_t num_points = static_cast<size_t>(m_options.num_characteristics);
    double dy = throat_point.y/num_points;
    double delt = delta();

    double sonic_x = (gamma + 1)*alpha/(2*(3+delt));

    std::vector<CharacteristicPoint> points(num_points);

    for (size_t i = 0; i < num_points; i++) {
        CharacteristicPoint pt = points[i];
        pt.y = dy * i;
        pt.x = (gamma + 1)*alpha/(2*(3+delta()))*(1.0 - pt.y*pt.y);
        pt.theta = 0.0; // by definition
        double u_prime = alpha*pt.x + (gamma+1)*alpha*alpha*pt.y*pt.y/(2*(1+delt));
        pt.mach = pm_table.interpolate_mach_from_V(u_prime);
        pt.nu = pm_table.interpolate_nu_from_mach(pt.mach);
        // TODO: update thermodynamic state for each point
    }


    return points;
}

} // namespace Goddard