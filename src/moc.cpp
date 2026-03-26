#include <cmath>
#include <algorithm>
#include "goddard/equilibrium.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/moc.hpp"
#include "goddard/prandtlmeyer.hpp"
#include "goddard/error.hpp"

namespace Goddard {
    
// -- NozzleProfile --

double NozzleProfile::slope_at(double x_query) const {
    size_t idx = find_index(x_query); // first index with value larger than x
    double x2 = x[idx];
    double x1 = x[idx-1];
    
    double y2 = y[idx];
    double y1 = y[idx-1];

    return (y2-y1)/(x2-x1);
}

double NozzleProfile::theta_at(double x_query) const {
    
    return atan(slope_at(x_query));
}

std::pair<double, double> NozzleProfile::at(size_t idx) const {
    return {x[idx], y[idx]};
}

void NozzleProfile::push_back(std::pair<double, double>&& coords) {
    x.push_back(coords.first);
    y.push_back(coords.second);
}

size_t NozzleProfile::size() const {
    return x.size();
}

NozzleProfile NozzleProfile::load_profile_csv(const std::string& /*filename*/) {
    throw NotImplementedError("Loading csv profiles is not implemented.");
}

void NozzleProfile::save_profile_csv(const std::string& /*filename*/) {
    throw NotImplementedError("Saving csv profiles is not implemented.");
}

size_t NozzleProfile::find_index(double x_query) const {
    // linear scan acceptable for small grids
    for (size_t i = 0; i < x.size(); i++){
        if (x[i] > x_query) return i;
    }
    // query beyond profile: return last segment
    return x.size() - 1;
}

// -- MocNozzle --
MocResult MocNozzle::solve() {
    m_messages.clear();
    m_theta_schedule.clear();

    CharacteristicNet net;
    if (m_options.mode == MocMode::DESIGN_MIN_LENGTH) {
        net.topology = NetTopology::TRIANGULAR;
    } 
    else {
        net.topology = NetTopology::WAVEFRONT;
    }

    
    std::vector<CharacteristicPoint> data_line;

    if (m_options.chemistry == MocChemistry::PERFECT_GAS) {
        // Pure algebraic path: no Cantera dependency
        m_L_ref = m_options.geometry.throat_radius;
        m_P_ref = 1.0; // dimensionless stagnation pressure
        m_T_ref = 1.0; // dimensionless stagnation temperature

        data_line = generate_initial_data_line_perfect_gas(
            m_options.num_characteristics);
    } 
    else {
        // Cantera-backed path for frozen/equilibrium chemistry
        std::unique_ptr<NozzleBase> nozzle;
        switch (m_options.chemistry) {
            case MocChemistry::FROZEN: {
                nozzle = std::make_unique<FrozenNozzle>(*m_gas);
                break;
            }
            case MocChemistry::EQUILIBRIUM: {
                nozzle = std::make_unique<EquilibriumNozzle>(*m_gas);
                break;
            }
            default:
                break;
        }
        auto throat = nozzle->solve_throat_conditions();

        m_L_ref = m_options.geometry.throat_radius;
        m_P_ref = throat.P_inlet;
        m_T_ref = m_gas->thermo()->temperature();

        data_line = generate_initial_data_line(throat, m_options.geometry, m_options.num_characteristics);
    }
    net.num_c_minus = data_line.size();
    net.num_c_plus = data_line.size();
    net.wavefronts.push_back(data_line);
    net.wall_x.push_back(0.0);
    net.wall_y.push_back(1.0);

    solve_kernel_region(net);
    solve_wall_region(net);

    MocResult result;
    result.net = net;
    result.messages = m_messages;
    result.converged = true; // no invalid points detected (TODO: check during solve)

    // Populate performance fields
    if (!net.wall_x.empty()) {
        result.nozzle_length = net.wall_x.back();
    }
    if (!net.wall_y.empty() && net.wall_y.front() > 0.0) {
        if (m_options.flow_type == MocFlowKind::PLANAR) {
            result.area_ratio = net.wall_y.back() / net.wall_y.front();
        } else {
            double y_ratio = net.wall_y.back() / net.wall_y.front();
            result.area_ratio = y_ratio * y_ratio;
        }
    }

    // Exit Mach: from the last wavefront's last point (outermost at exit plane)
    // For a min-length nozzle, the last wavefront has one point and the exit
    // plane Mach should be uniform.
    if (!net.wavefronts.empty() && !net.wavefronts.back().empty()) {
        result.exit_mach = net.wavefronts.back().back().mach;
    }

    // Build wall profile from computed wall coordinates
    result.profile.x = net.wall_x;
    result.profile.y = net.wall_y;

    return result;
}

std::vector<CharacteristicPoint> MocNozzle::generate_initial_data_line(
    const ThroatCondition& throat,
    const ThroatGeometry& /*geometry*/,
    size_t num_points)
{
    auto thermo = m_gas->thermo();
    thermo->restoreState(throat.state);

    std::vector<CharacteristicPoint> data_line;
    data_line.reserve(num_points);
    

    CharacteristicPoint sonic_point;
    sonic_point.temperature = thermo->temperature()/m_T_ref;
    sonic_point.pressure = thermo->pressure()/m_P_ref;
    sonic_point.cantera_state = throat.state;
    sonic_point.mach = 1.0; // by definition
    sonic_point.theta = 0.0; // by definition
    sonic_point.nu = 0.0;
    sonic_point.x = 0.0;
    sonic_point.y = 1.0; // dimensionless throat radius
    sonic_point.K_minus = 0.0;
    sonic_point.K_plus = 0.0;

    sonic_point.gamma_s = get_gamma_s();
    double sonic_stagnation_factor = stagnation_factor(sonic_point.mach, sonic_point.gamma_s);
    
    
    switch (m_options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            // the first theta should be small to minimize approximation error.
            double dtheta_initial = m_options.theta_max/(num_points*10);
            double dtheta = (m_options.theta_max-dtheta_initial)/(num_points - 1);
            CharacteristicPoint upstream_point;
            
            //for a minimum length nozzle, the sonic line is straight.
            //the initial data line is generated through a single Prandtl-Meyer expansion at the throat.
            for (size_t i = 0; i < num_points; i++){
                //TODO: this may only be valid for the algebraic case. Need to confirm.
                CharacteristicPoint expansion_point;
                expansion_point.x = sonic_point.x;
                expansion_point.y = sonic_point.y; // dimensionless throat radius
                expansion_point.theta = dtheta_initial + i*dtheta;
                // save theta for later use in wall solver
                m_theta_schedule.push_back(expansion_point.theta);
                expansion_point.nu = expansion_point.theta;
                expansion_point.K_plus = expansion_point.theta - expansion_point.nu;
                expansion_point.K_minus = expansion_point.theta + expansion_point.nu;
                expansion_point.mach = mach_from_prandtl_meyer(expansion_point.nu, m_options.gamma, 1.0);
                expansion_point.mu = asin(1.0/expansion_point.mach);
                expansion_point.gamma_s = m_options.gamma;
                
                double stagnation_ratio = sonic_stagnation_factor/stagnation_factor(expansion_point.mach, expansion_point.gamma_s);
                expansion_point.temperature = sonic_point.temperature*stagnation_ratio;
                expansion_point.pressure = sonic_point.pressure*std::pow(stagnation_ratio, expansion_point.gamma_s/(expansion_point.gamma_s-1.0));
                
                // special case: the first point is on the centerline
                if (i == 0) {
                    upstream_point = solve_axis_point(expansion_point);
                    data_line.push_back(upstream_point);
                } 
                else {
                    upstream_point = solve_interior_point(expansion_point, upstream_point);
                    data_line.push_back(upstream_point);
                }
            }
            break;
        }
        case MocMode::DESIGN_RAO: {
            throw NotImplementedError("Rao design mode initialization not implemented.");
        }
        case MocMode::ANALYSIS: {
            throw NotImplementedError("Analysis mode initialization not implemented.");
        }
    }
    return data_line;
}

std::vector<CharacteristicPoint> MocNozzle::generate_initial_data_line_perfect_gas(
    size_t num_points)
{
    double gamma = m_options.gamma;

    // Build or use theta schedule
    if (!m_options.theta_schedule.empty()) {
        m_theta_schedule = m_options.theta_schedule;
        num_points = m_theta_schedule.size();
    } else {
        // Auto-generate: small first step, then uniform spacing
        double dtheta_initial = m_options.theta_max / (num_points * 10);
        double dtheta = (m_options.theta_max - dtheta_initial) / (num_points - 1);
        m_theta_schedule.resize(num_points);
        for (size_t i = 0; i < num_points; i++) {
            m_theta_schedule[i] = dtheta_initial + i * dtheta;
        }
    }

    // Sonic point at throat (dimensionless)
    CharacteristicPoint sonic_point{};
    sonic_point.mach = 1.0;
    sonic_point.theta = 0.0;
    sonic_point.nu = 0.0;
    sonic_point.x = 0.0;
    sonic_point.y = 1.0;
    sonic_point.K_minus = 0.0;
    sonic_point.K_plus = 0.0;
    sonic_point.gamma_s = gamma;
    sonic_point.temperature = 1.0 / stagnation_factor(1.0, gamma); // T/T0 at M=1
    sonic_point.pressure = std::pow(sonic_point.temperature, gamma / (gamma - 1.0)); // p/p0 at M=1

    double sonic_stag = stagnation_factor(1.0, gamma);

    std::vector<CharacteristicPoint> data_line;
    data_line.reserve(num_points);

    CharacteristicPoint upstream_point;

    for (size_t i = 0; i < num_points; i++) {
        // Each C- characteristic from the centered expansion fan
        CharacteristicPoint expansion_point{};
        expansion_point.x = sonic_point.x;
        expansion_point.y = sonic_point.y;
        expansion_point.theta = m_theta_schedule[i];
        expansion_point.nu = expansion_point.theta; // centered fan: nu = theta
        expansion_point.K_plus = 0.0; // theta - nu = 0 for all fan rays
        expansion_point.K_minus = 2.0 * expansion_point.theta; // theta + nu
        expansion_point.mach = mach_from_prandtl_meyer(expansion_point.nu, gamma, 1.0);
        expansion_point.mu = asin(1.0 / expansion_point.mach);
        expansion_point.gamma_s = gamma;

        double stag_ratio = sonic_stag / stagnation_factor(expansion_point.mach, gamma);
        expansion_point.temperature = sonic_point.temperature * stag_ratio;
        expansion_point.pressure = sonic_point.pressure
            * std::pow(stag_ratio, gamma / (gamma - 1.0));

        if (i == 0) {
            upstream_point = solve_axis_point(expansion_point);
        } else {
            upstream_point = solve_interior_point(expansion_point, upstream_point);
        }
        data_line.push_back(upstream_point);
    }

    return data_line;
}

void MocNozzle::solve_kernel_region(CharacteristicNet& net) {
    switch (m_options.mode) {
        //TODO: I believe the logic below is generalizable to all cases as long as 
        //there is special handling for reflections in the expansion region.
        case MocMode::DESIGN_MIN_LENGTH: {
            // for a minimum length nozzle, the expansion region is infinitely small
            // so we are guaranteed to not have Mach wave reflections off the upper wall.
            // This means the number of characteristics is known upfront.
            size_t num_characteristics = static_cast<size_t>(net.num_c_plus);
            // the last wavefront at num_characteristics-1 should be skipped.
            for (size_t i = 0; i < num_characteristics-1; i++) {
                CharacteristicNet::Wavefront next_wavefront;
                auto& curr_wavefront = net.wavefronts[i];
                // skip the first node in a wavefront -- it is the centerline node which doesn't impact
                // downstream nodes.
                CharacteristicPoint prev_point;
                for (size_t j = 1; j < curr_wavefront.size(); j++) {
                    const CharacteristicPoint& parent_point = curr_wavefront[j];
                    // first non-centerline node impacts the downstream centerline node.
                    if (j == 1) {
                        prev_point = solve_axis_point(parent_point);
                    }
                    else {
                        prev_point = solve_interior_point(parent_point, prev_point);
                    }
                    next_wavefront.push_back(prev_point);
                }
                net.wavefronts.push_back(next_wavefront);
            }
            break;
        }
        case MocMode::DESIGN_RAO: {
            break;
        }
        case MocMode::ANALYSIS: {
            break;
        }
    }
}

void MocNozzle::solve_wall_region(CharacteristicNet& net) {
    switch (m_options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            CharacteristicPoint wall_point;
            // For a minimum length nozzle, the first wall point is located at the throat.
            // We only use the coordinates and angle to determine the position of the next wall point,
            // so the other fields can be left empty.
            wall_point.theta = m_options.theta_max;
            wall_point.x = 0.0;
            wall_point.y = 1.0;
            
            // for the min length nozzle case, every wavefront has a wall node
            // impacted by the last node in the wavefront 
            for (size_t i = 0; i < net.wavefronts.size(); i++) {
                auto& wavefront = net.wavefronts[i];
                const CharacteristicPoint& parent_point = wavefront.back();
                wall_point = solve_wall_point(parent_point, wall_point, static_cast<int>(i));
                net.wall_x.push_back(wall_point.x);
                net.wall_y.push_back(wall_point.y);
            }
            break;
        }
        case MocMode::DESIGN_RAO: {
            break;
        }
        case MocMode::ANALYSIS: {
            break;
        }
    }
}


CharacteristicPoint MocNozzle::solve_wall_point(
    const CharacteristicPoint& interior_parent,
    const CharacteristicPoint& previous_wall_point,
    int wall_point_index)
{
    switch (m_options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            // for minimum length design: theta is determined by theta schedule
            double theta_wall = m_options.theta_max - m_theta_schedule[wall_point_index];
            
            return solve_wall_point_design(interior_parent, previous_wall_point, theta_wall);
        }
        case MocMode::DESIGN_RAO: {
            throw NotImplementedError("Rao nozzle design is not implemented.");
        }
        case MocMode::ANALYSIS: {
            auto [x_wall, y_wall] = intersect_characteristic_with_wall(
                interior_parent, m_options.nozzle_profile);
            double theta_wall = m_options.nozzle_profile.theta_at(x_wall);
            return solve_wall_point_analysis(interior_parent, theta_wall, x_wall, y_wall);
        }
    }
    // unreachable
    throw std::runtime_error("Invalid MocMode in solve_wall_point");
}

CharacteristicPoint MocNozzle::solve_axis_point(const CharacteristicPoint& off_axis_parent) {
    CharacteristicPoint axis_point;
    axis_point.y = 0.0;
    axis_point.theta = 0.0; // symmetry condition

    if (m_options.chemistry == MocChemistry::PERFECT_GAS 
        && m_options.flow_type == MocFlowKind::PLANAR) 
    {
        double gamma = m_options.gamma;
        axis_point.gamma_s = gamma;
        axis_point.K_minus = off_axis_parent.K_minus;
        axis_point.nu = axis_point.K_minus - axis_point.theta;
        axis_point.K_plus = axis_point.theta - axis_point.nu;
        axis_point.mach = mach_from_prandtl_meyer(axis_point.nu, gamma, off_axis_parent.mach);
        axis_point.mu = asin(1.0/axis_point.mach);

        double off_axis_stag = stagnation_factor(off_axis_parent.mach, gamma);
        double axis_stag = stagnation_factor(axis_point.mach, gamma);
        double stagnation_ratio = off_axis_stag/axis_stag;
        axis_point.temperature = off_axis_parent.temperature * stagnation_ratio;
        axis_point.pressure = off_axis_parent.pressure * pow(stagnation_ratio, gamma/(gamma-1.0));

        double c_minus_angle = 0.5 * (off_axis_parent.theta + axis_point.theta) 
            - 0.5 * (off_axis_parent.mu + axis_point.mu);
        
        axis_point.x = off_axis_parent.x - off_axis_parent.y/tan(c_minus_angle);
        

    } else if (m_options.flow_type == MocFlowKind::PLANAR) {
        throw NotImplementedError("Flow solver for non-perfect gases.");
    } else {
        // estimate dtheta/dy from parent point (will be used when implemented)
        // double dtheta_dy = off_axis_parent.theta / off_axis_parent.y;
        throw NotImplementedError("Axis flow solver for axisymmetric flow.");
    }

    return axis_point;
}

std::pair<double, double> MocNozzle::intersect_characteristic_with_wall(
        const CharacteristicPoint& parent,
        const NozzleProfile& wall) {
    
    double char_angle = parent.theta + parent.mu;
    double char_slope = std::tan(char_angle);
    
    // we use a basic predictor-corrector scheme to find the wall intersection
    // predictor: straight-line intersection
    // linear scan is appropriate for tens to low hundreds of points. 
    // if npoints grows beyond that, switch to bisection algorithm
    double x_hit, y_hit;
    for (size_t k = 0; k < wall.x.size() - 1; k++) {
        double wall_slope = (wall.y[k+1] - wall.y[k])
            / (wall.x[k+1] - wall.x[k]);
        double wall_intercept = wall.y[k] - wall_slope * wall.x[k];
        double char_intercept = parent.y - char_slope * parent.x;

        double x_intercept = (wall_intercept - char_intercept) / (char_slope - wall_slope);

        // check if intersection is within the wall segment
        if (x_intercept >= wall.x[k] && x_intercept <= wall.x[k+1]) {
            x_hit = x_intercept;
            y_hit = parent.y + char_slope * (x_hit - parent.x);
            break;
        }
    }

    // corrector (for curved characteristics):
    // solve the wall point flow using computed theta
    // then recompute characteristic slope as average of parent and wall-point slopes.
    // one correction is usually sufficient.
    double wall_theta = wall.theta_at(x_hit);
    CharacteristicPoint wall_point = solve_wall_flow(parent, wall_theta);
    double wall_char_angle = wall_point.theta + wall_point.mu;
    double average_char_angle = 0.5*(wall_char_angle + char_angle);
    double corrected_char_slope = std::tan(average_char_angle);

    for (size_t k = 0; k < wall.x.size() - 1; k++) {
        double wall_slope = (wall.y[k+1] - wall.y[k])
            / (wall.x[k+1] - wall.x[k]);
        double wall_intercept = wall.y[k] - wall_slope * wall.x[k];
        double char_intercept = parent.y - corrected_char_slope * parent.x;

        double x_intercept = (wall_intercept - char_intercept) / (corrected_char_slope - wall_slope);

        // check if intersection is within the wall segment
        if (x_intercept >= wall.x[k] && x_intercept <= wall.x[k+1]) {
            x_hit = x_intercept;
            y_hit = parent.y + corrected_char_slope * (x_hit - parent.x);
            break;
        }
    }

    return {x_hit, y_hit};
}

CharacteristicPoint MocNozzle::solve_interior_point(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2)
{
    if (m_options.chemistry == MocChemistry::PERFECT_GAS 
        && m_options.flow_type == MocFlowKind::PLANAR) 
    {
        return solve_interior_point_algebraic(p1, p2);
    } 
    else {
        return solve_interior_point_iterative(p1, p2);
    }
} 

CharacteristicPoint MocNozzle::solve_interior_point_algebraic(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2) 
{

    double gamma = m_options.gamma;

    CharacteristicPoint p3;
    p3.K_minus = p1.K_minus;
    p3.K_plus = p2.K_plus;
    p3.theta = 0.5*(p1.K_minus + p2.K_plus);
    p3.nu = 0.5*(p1.K_minus - p2.K_plus);
    p3.gamma_s = gamma;

    double mach_guess = (p1.mach + p2.mach)/2.0;
    
    p3.mach = mach_from_prandtl_meyer(p3.nu, gamma, mach_guess);
    p3.mu = asin(1.0/p3.mach);

    // we arbitrarily use point 1 as the reference stagnation point.
    // this makes no difference due to Crocco's theorem
    double stagnation_factor_1 = stagnation_factor(p1.mach, gamma);
    double stagnation_factor_3 = stagnation_factor(p3.mach, gamma);
    
    p3.temperature = p1.temperature * stagnation_factor_1/stagnation_factor_3;
    p3.pressure = p1.pressure * pow(stagnation_factor_1/stagnation_factor_3, gamma/(gamma-1.0));

    // compute intersection point
    // assumimg characteristics are straight lines.
    double c_minus_angle = 0.5 * (p1.theta + p3.theta) - 0.5 * (p1.mu + p3.mu);
    double c_plus_angle = 0.5 * (p2.theta + p3.theta) + 0.5 * (p2.mu + p3.mu);

    auto [x, y] = characteristic_intersection_coordinates(p1, p2, c_minus_angle, c_plus_angle);
    p3.x = x;
    p3.y = y;

    // sanity check
    if (p3.x < p1.x || p3.x < p2.x) {
        p3.mach = -1.0; // invalid point
    }
    
    return p3;
}

CharacteristicPoint MocNozzle::solve_interior_point_iterative(
    const CharacteristicPoint& /*p1*/,
    const CharacteristicPoint& /*p2*/)
{
    throw NotImplementedError("Iterative interior points are not implemented.");

}

CharacteristicPoint MocNozzle::solve_wall_flow(
    const CharacteristicPoint& interior_parent,
    double theta_wall)
{
    double gamma = m_options.gamma;

    CharacteristicPoint wall_point;
    wall_point.theta = theta_wall;
    
    // Flow solve: compatibility equation along C+ from interior parent
    // K+ is preserved: theta - nu = const along C+
    wall_point.K_plus = interior_parent.K_plus;
    wall_point.nu = theta_wall - interior_parent.K_plus;
    wall_point.mach = mach_from_prandtl_meyer(
        wall_point.nu, gamma, interior_parent.mach);
    wall_point.mu = std::asin(1.0/ wall_point.mach);
    wall_point.K_minus = wall_point.theta + wall_point.nu;

    // Isentropic relations for pressure, temperature
    double stagnation_factor_parent = stagnation_factor(interior_parent.mach, gamma);
    double stagnation_factor_wall = stagnation_factor(wall_point.mach, gamma);
    double temperature_ratio = stagnation_factor_parent / stagnation_factor_wall;

    wall_point.temperature = interior_parent.temperature * temperature_ratio;
    wall_point.pressure = interior_parent.pressure * std::pow(temperature_ratio, gamma/(gamma-1.0));

    return wall_point;
}

CharacteristicPoint MocNozzle::solve_wall_point_design(
    const CharacteristicPoint& interior_parent,
    const CharacteristicPoint& previous_wall_point,
    double theta_wall)
{
    CharacteristicPoint wall_point = solve_wall_flow(interior_parent, theta_wall);

    // the angle between the interior point and the wall point 
    // is given by the C+ characteristic.
    double c_plus_angle = 0.5 * (interior_parent.theta + wall_point.theta) 
        + 0.5 * (interior_parent.mu + wall_point.mu);

    // the angle between the previous wall point and the current wall point
    // is given by wall_theta. 
    double prev_wall_angle = 0.5*(previous_wall_point.theta + wall_point.theta);

    auto [x,y] = characteristic_intersection_coordinates(
        interior_parent, previous_wall_point, 
        c_plus_angle, prev_wall_angle);

    wall_point.x = x; 
    wall_point.y = y;

    return wall_point;
}

CharacteristicPoint MocNozzle::solve_wall_point_analysis(
        const CharacteristicPoint& interior_parent,
        double theta_wall, double x_wall, double y_wall) 
{
    CharacteristicPoint wall_point = solve_wall_flow(interior_parent, theta_wall);
    
    wall_point.x = x_wall;
    wall_point.y = y_wall;
    return wall_point;
}

double MocNozzle::get_gamma_s() {
    switch (m_options.chemistry) {
        case MocChemistry::PERFECT_GAS: {
            return m_options.gamma;
        }
        case MocChemistry::FROZEN: {
            return m_gas->thermo()->cp_mass()/m_gas->thermo()->cv_mass();
        }
        case MocChemistry::EQUILIBRIUM: {
            return get_equilibrium_gamma(*m_gas->thermo());
        }
    }
    // unreachable
    return m_options.gamma;
}


} // namespace Goddard