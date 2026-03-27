#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <sstream>
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

NozzleProfile NozzleProfile::load_profile_csv(const std::string& filename) {
    NozzleProfile profile;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    std::string line;
    // skip header line
    std::getline(file, line);
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        double x_val, y_val;
        char comma;
        if (iss >> x_val >> comma >> y_val) {
            profile.x.push_back(x_val);
            profile.y.push_back(y_val);
        }
    }
    return profile;
}

void NozzleProfile::save_profile_csv(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    file << std::setprecision(15);
    file << "x,y\n";
    for (size_t i = 0; i < x.size(); i++) {
        file << x[i] << "," << y[i] << "\n";
    }
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
        m_S_ref = throat.S_inlet;

        double a_throat = gas_sonic_velocity(*m_gas->thermo(), nozzle->get_gamma_s(*m_gas->thermo()));
        
        pm_table.build_table(
            *m_gas->thermo(),
            m_options.chemistry == MocChemistry::EQUILIBRIUM,
            throat.S_inlet,
            throat.H_stagnation,
            a_throat);
                
        data_line = generate_initial_data_line(throat, m_options.geometry, m_options.num_characteristics);
    }
    net.num_c_minus = data_line.size();
    net.num_c_plus = data_line.size();
    net.wavefronts.push_back(data_line);
    net.wall_x.push_back(0.0);
    net.wall_y.push_back(1.0);

    solve_kernel_region(net);
    solve_wall_region(net);

    // Check for invalid points
    bool all_valid = true;
    for (const auto& wf : net.wavefronts) {
        for (const auto& pt : wf) {
            if (pt.mach < 0.0) {
                all_valid = false;
                m_messages.push_back("Invalid point detected: non-downstream intersection or subsonic pocket");
                break;
            }
        }
        if (!all_valid) break;
    }

    MocResult result;
    result.net = net;
    result.messages = m_messages;
    result.converged = all_valid;

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

    // Exit plane extraction:
    // For a min-length nozzle, the exit plane is the last wavefront + last wall point.
    // The last wavefront contains the axis point, and each preceding wavefront's
    // last point was absorbed into wall calculations.
    if (!net.wavefronts.empty()) {
        const auto& last_wf = net.wavefronts.back();
        for (const auto& pt : last_wf) {
            result.exit_plane.y.push_back(pt.y);
            result.exit_plane.mach.push_back(pt.mach);
            result.exit_plane.theta.push_back(pt.theta);
            result.exit_plane.pressure.push_back(pt.pressure);
            result.exit_plane.temperature.push_back(pt.temperature);
        }
        // Add the last wall point if available
        if (!net.wall_points.empty()) {
            const auto& wp = net.wall_points.back();
            result.exit_plane.y.push_back(wp.y);
            result.exit_plane.mach.push_back(wp.mach);
            result.exit_plane.theta.push_back(wp.theta);
            result.exit_plane.pressure.push_back(wp.pressure);
            result.exit_plane.temperature.push_back(wp.temperature);
        }
    }

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
    sonic_point.temperature = thermo->temperature()/m_T_ref; // should be 1.0 by definition
    sonic_point.pressure = thermo->pressure()/m_P_ref;
    sonic_point.cantera_state = throat.state;
    sonic_point.mach = 1.0; // by definition
    sonic_point.theta = 0.0; // by definition
    sonic_point.nu = 0.0;
    sonic_point.x = 0.0;
    sonic_point.y = 1.0; // dimensionless throat radius
    sonic_point.K_minus = 0.0;
    sonic_point.K_plus = 0.0;

    sonic_point.gamma_s = gamma_s_from_nu(sonic_point.nu);
    
    
    switch (m_options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            // the first theta should be small to minimize approximation error.
            double dtheta_initial = m_options.theta_max/(num_points*10);
            double dtheta = (m_options.theta_max-dtheta_initial)/(num_points - 1);
            CharacteristicPoint upstream_point;
            
            //for a minimum length nozzle, the sonic line is straight.
            //the initial data line is generated through a single Prandtl-Meyer expansion at the throat.
            for (size_t i = 0; i < num_points; i++){
                CharacteristicPoint expansion_point;
                expansion_point.x = sonic_point.x;
                expansion_point.y = sonic_point.y; // dimensionless throat radius
                expansion_point.theta = dtheta_initial + i*dtheta;
                // save theta for later use in wall solver
                m_theta_schedule.push_back(expansion_point.theta);
                // for a centered expansion fan, nu = theta by definition
                update_thermodynamic_state_from_nu(expansion_point, expansion_point.theta, 1.0);
                expansion_point.K_plus = expansion_point.theta - expansion_point.nu;
                expansion_point.K_minus = expansion_point.theta + expansion_point.nu;
                
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

    std::vector<CharacteristicPoint> data_line;
    data_line.reserve(num_points);

    CharacteristicPoint upstream_point;

    for (size_t i = 0; i < num_points; i++) {
        // Each C- characteristic from the centered expansion fan
        CharacteristicPoint expansion_point{};
        expansion_point.x = sonic_point.x;
        expansion_point.y = sonic_point.y;
        expansion_point.theta = m_theta_schedule[i];
        update_thermodynamic_state_from_nu(expansion_point, expansion_point.theta, 1.0); // centered fan: nu = theta
        expansion_point.K_plus = 0.0; // theta - nu = 0 for all fan rays
        expansion_point.K_minus = expansion_point.theta + expansion_point.nu; // theta + nu

        if (i == 0) {
            upstream_point = solve_initial_axis_point(expansion_point);            
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
            // For a minimum length nozzle, the expansion region is infinitely small
            // so we are guaranteed to not have Mach wave reflections off the upper wall.
            // This means the number of characteristics is known upfront.
            size_t num_characteristics = static_cast<size_t>(net.num_c_plus);
            // the last wavefront at num_characteristics-1 should be skipped.
            for (size_t i = 0; i < num_characteristics-1; i++) {
                CharacteristicNet::Wavefront next_wavefront;
                auto& curr_wavefront = net.wavefronts[i];
                // skip the first node in a wavefront -- it is the centerline node which doesn't impact
                // downstream nodes. The C- characteristic is irrelevant due to symmetry, 
                // and the C+ characteristic impacts nodes on the same wavefront.
                CharacteristicPoint prev_point;
                for (size_t j = 1; j < curr_wavefront.size(); j++) {
                    const CharacteristicPoint& parent_point = curr_wavefront[j];
                    // first non-centerline node impacts the downstream centerline node.
                    // all other nodes impact interior downstream nodes.
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
                net.wall_points.push_back(wall_point);
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

CharacteristicPoint MocNozzle::solve_initial_axis_point(const CharacteristicPoint& expansion_point) {
    CharacteristicPoint point;
    point.y = 0.0; //point always lies on axis
    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: {
            point.theta = expansion_point.theta;
            // nu = expansion point nu follows from geometric analysis
            update_thermodynamic_state_from_nu(point, expansion_point.nu, expansion_point.mach);
            point.K_plus = point.theta - point.nu;
            point.K_minus = expansion_point.K_minus; 
            
            double c_minus_angle = average_cminus_angle(expansion_point, point);
            point.x = expansion_point.x - expansion_point.y / tan(c_minus_angle);
            break;
        }
        case MocFlowKind::AXISYMMETRIC: {
            // Same approach as planar for the initial axis point:
            // assign theta from the expansion fan to bootstrap marching.
            // The axisymmetric source term is not applied here because the
            // initial expansion fan is modeled as a centered expansion at a point,
            // and the first axis point is a direct consequence of this geometric
            // construction. Source terms enter in subsequent solve_axis_point calls.
            point.theta = expansion_point.theta;
            update_thermodynamic_state_from_nu(point, expansion_point.nu, expansion_point.mach);
            point.K_plus = point.theta - point.nu;
            point.K_minus = expansion_point.K_minus;

            double c_minus_angle = average_cminus_angle(expansion_point, point);
            point.x = expansion_point.x - expansion_point.y / tan(c_minus_angle);
            break;
        }
        default:
            throw std::runtime_error("Invalid flow type specified.");
    }

    return point;
}

CharacteristicPoint MocNozzle::solve_axis_point(const CharacteristicPoint& off_axis_parent) {
    CharacteristicPoint axis_point;
    axis_point.y = 0.0;
    axis_point.theta = 0.0; // symmetry condition

    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: {

            axis_point.K_minus = off_axis_parent.K_minus;
            update_thermodynamic_state_from_nu(
                axis_point, 
                axis_point.K_minus - axis_point.theta,
                off_axis_parent.mach);

            axis_point.K_plus = axis_point.theta - axis_point.nu;

            double c_minus_angle = 0.5 * (off_axis_parent.theta + axis_point.theta) 
                - 0.5 * (off_axis_parent.mu + axis_point.mu);
            
            axis_point.x = off_axis_parent.x - off_axis_parent.y/tan(c_minus_angle);
            break;
        }
        case MocFlowKind::AXISYMMETRIC: {
            // The source term sin(theta)/y is 0/0 at y=0.
            // By L'Hopital's rule: lim_{y->0} sin(theta)/y = dtheta/dy,
            // estimated from the off-axis parent point.
            double dtheta_dy = off_axis_parent.theta / off_axis_parent.y;

            // Predictor: planar solution as initial guess
            axis_point.K_minus = off_axis_parent.K_minus;
            double nu_pred = axis_point.K_minus; // theta=0 => nu = K_minus
            update_thermodynamic_state_from_nu(axis_point, nu_pred, off_axis_parent.mach);
            axis_point.K_plus = -axis_point.nu;

            double c_minus_angle = average_cminus_angle(off_axis_parent, axis_point);
            axis_point.x = off_axis_parent.x - off_axis_parent.y / tan(c_minus_angle);

            // Corrector: apply limiting source term
            // The C- compatibility equation at the axis becomes:
            // dtheta + dnu = S_cminus * ds
            // where S_cminus_limit = dtheta_dy / (M * sin(mu))
            // and ds is the arc length along the C- characteristic.
            double ds = off_axis_parent.y / std::abs(sin(c_minus_angle));
            double mu_avg = 0.5 * (off_axis_parent.mu + axis_point.mu);
            double M_avg = 0.5 * (off_axis_parent.mach + axis_point.mach);
            double source_limit = dtheta_dy / (M_avg * sin(mu_avg)) * ds;

            // Corrected nu: K_minus from parent, minus source contribution
            double nu_corrected = off_axis_parent.K_minus + source_limit;
            update_thermodynamic_state_from_nu(axis_point, nu_corrected, axis_point.mach);
            axis_point.K_minus = axis_point.nu; // theta=0
            axis_point.K_plus = -axis_point.nu;

            // Recompute position with corrected slope
            c_minus_angle = average_cminus_angle(off_axis_parent, axis_point);
            axis_point.x = off_axis_parent.x - off_axis_parent.y / tan(c_minus_angle);
            break;
        }
        default:
            throw std::runtime_error("Invalid flow type specified.");
    }

    return axis_point;
}

std::pair<double, double> MocNozzle::intersect_characteristic_with_wall(
        const CharacteristicPoint& parent,
        const NozzleProfile& wall) {
    
    double char_angle = parent.theta + parent.mu;
    double char_slope = tan(char_angle);
    
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
    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: {
            return solve_interior_point_algebraic(p1, p2);
        }
        case MocFlowKind::AXISYMMETRIC: {
            return solve_interior_point_iterative(p1, p2);
        }
        default:
            throw std::runtime_error("Invalid MoC flow type specified.");
    }
} 

CharacteristicPoint MocNozzle::solve_interior_point_algebraic(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2)
{
    CharacteristicPoint p3;
    p3.theta = 0.5*(p1.K_minus + p2.K_plus);
    double mach_guess = 0.5*(p1.mach + p2.mach);
    
    update_thermodynamic_state_from_nu(p3, 0.5*(p1.K_minus - p2.K_plus), mach_guess);
    p3.K_minus = p3.theta + p3.nu;
    p3.K_plus = p3.theta - p3.nu;

    // compute intersection point
    // assumimg characteristics are straight lines.
    double c_minus_angle = average_cminus_angle(p1, p3);
    double c_plus_angle = average_cplus_angle(p2, p3);

    auto [x, y] = characteristic_intersection_coordinates(p1, p2, c_minus_angle, c_plus_angle);
    p3.x = x;
    p3.y = y;

    // validity checks
    if (p3.x < p1.x || p3.x < p2.x) {
        p3.mach = -1.0; // invalid point
    }
    if (p3.mach < 1.0 && p3.mach > 0.0) {
        p3.mach = -1.0; // subsonic point -- characteristics undefined
    }

    return p3;
}

CharacteristicPoint MocNozzle::solve_interior_point_iterative(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2)
{
    // General procedure:
    // 1. Receive theta and nu from the two incoming characteristics
    // 2. Estimate characteristic slopes at parent nodes
    // 3. Intersect straight line characteristics
    // 4. Compute source terms (if applicable)
    // 5. Solve nonlinear equation in V using upstream states only
    // 6. Get theta from either compatibility equation 
    // 7. Correct characteristic slope by averaging with endpoint 
    // 8. Correct V by averaging with endpoint 
    // 9. Recompute theta and source terms at corrected location 
    // 10. [optional] One additional pass if source term changed significantly 
    CharacteristicPoint p3;
    double c_minus_angle = p1.theta - p1.mu;
    double c_plus_angle = p2.theta + p2.mu;

    auto [x, y] = characteristic_intersection_coordinates(p1, p2, c_minus_angle, c_plus_angle);
    p3.x = x;
    p3.y = y;
    double S1 = cminus_source_term(p1, y);
    double S2 = cplus_source_term(p2, y);

    update_thermodynamic_state_from_mach(p3, find_node_mach(p1, p2, S1 - S2, 0.5 * (p1.mach + p2.mach)));

    // from compatibility equation.
    p3.theta = S1 + p1.theta - (p3.nu - p1.nu);
    p3.K_plus = p3.theta - p3.nu;
    p3.K_minus = p3.theta + p3.nu;

    // corrector step
    const int max_iters = 3;

    for (int i = 0; i < max_iters; i++) {
        c_minus_angle = average_cminus_angle(p1, p3);
        c_plus_angle = average_cplus_angle(p2, p3);
        auto [new_x, new_y] = characteristic_intersection_coordinates(p1, p2, c_minus_angle, c_plus_angle);
        p3.x = new_x;
        p3.y = new_y;
        double S1_new = cminus_source_term(p1, new_y);
        double S2_new = cplus_source_term(p2, new_y);
        update_thermodynamic_state_from_mach(p3, find_node_mach(p1, p2, S1_new - S2_new, p3.mach));

        p3.theta = S1_new + p1.theta - (p3.nu - p1.nu);
        p3.K_plus = p3.theta - p3.nu;
        p3.K_minus = p3.theta + p3.nu;
        double residual = std::max(std::abs(S1_new - S1), std::abs(S2_new - S2));
        if (residual < m_options.abstol) break;
        S1 = S1_new;
        S2 = S2_new;
    }
    return p3;
}

CharacteristicPoint MocNozzle::solve_wall_flow(
    const CharacteristicPoint& interior_parent,
    double theta_wall)
{
    
    CharacteristicPoint wall_point;
    // geometric constraint: the flow at the wall must follow the wall curvature
    wall_point.theta = theta_wall;

    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: {
            // Flow solve: compatibility equation along C+ from interior parent
            // K+ is preserved: theta - nu = const along C+
            // this only applies for planar flow.
            wall_point.K_plus = interior_parent.K_plus;
            update_thermodynamic_state_from_nu(
                wall_point, 
                wall_point.theta - wall_point.K_plus, 
                interior_parent.mach);
            wall_point.K_minus = wall_point.theta + wall_point.nu;
            break;
        }
        case MocFlowKind::AXISYMMETRIC: {
            // For axisymmetric flow, K+ is not preserved along C+.
            // Use predictor-corrector: start with planar approximation,
            // then correct using the source term.
            // Predictor: planar K+ preservation
            wall_point.K_plus = interior_parent.K_plus;
            update_thermodynamic_state_from_nu(
                wall_point,
                wall_point.theta - wall_point.K_plus,
                interior_parent.mach);
            wall_point.K_minus = wall_point.theta + wall_point.nu;

            // The source term correction requires y, which is computed later
            // in solve_wall_point_design/analysis. We store the predictor result
            // and the corrector is applied in those methods.
            // For now, this gives a first-order approximation.
            break;
        }
        default:
            throw std::runtime_error("Invalid flow type specified.");
    }

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

double MocNozzle::gamma_s_from_mach(double mach) const {
    switch (m_options.chemistry) {
        case MocChemistry::PERFECT_GAS: {
            return m_options.gamma;
        }
        case MocChemistry::FROZEN: 
        case MocChemistry::EQUILIBRIUM: {
            return pm_table.interpolate_gamma_s_from_mach(mach);
        }
        default:
            throw std::runtime_error("Invalid value of MocChemistry specified.");
    }
    // unreachable
    return m_options.gamma;
}

double MocNozzle::gamma_s_from_nu(double nu) const {
    switch (m_options.chemistry) {
        case MocChemistry::PERFECT_GAS: {
            return m_options.gamma;
        }
        case MocChemistry::FROZEN: 
        case MocChemistry::EQUILIBRIUM: {
            return pm_table.interpolate_gamma_s_from_nu(nu);
        }
        default:
            throw std::runtime_error("Invalid value of MocChemistry specified.");
    }
    // unreachable
    return m_options.gamma;
}

double MocNozzle::mach_from_nu(const CharacteristicPoint& point, double mach_guess) const {
    switch (m_options.chemistry) {
        case MocChemistry::PERFECT_GAS: {
            return mach_from_prandtl_meyer(point.nu, point.gamma_s, mach_guess);
        }
        case MocChemistry::FROZEN: 
        case MocChemistry::EQUILIBRIUM: {
            return pm_table.interpolate_mach(point.nu);
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of MocChemistry specified.");
    }
}

double MocNozzle::nu_from_mach(double mach) const {
    switch (m_options.chemistry) {
        case MocChemistry::PERFECT_GAS: {
            return prandtl_meyer(mach, m_options.gamma);
        }
        case MocChemistry::FROZEN:
        case MocChemistry::EQUILIBRIUM: {
            return pm_table.interpolate_nu_from_mach(mach);
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of MocChemistry specified.");
    }
}

double MocNozzle::find_node_mach(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2,
    double source_delta,
    double mach_guess) const
{

    // initial guess: mach of upstream
    double mach = (mach_guess > 1.0) ? mach_guess : p1.mach;

    double delta_theta = p1.theta - p2.theta;
    const int max_iter = 15;
    // Newton's method to find root
    // Residual function obtained from combining both compatibility equations:
    //
    //$$\epsilon (V_P) = (\theta_{A}- \theta_{B}) + \Delta \nu (V_{A} \rightarrow V_{P})
    // + \Delta \nu (V_{B} \rightarrow V_{P}) - (S_{A}- S_{B})$$
    // where A and B indicate upstream points.
    for (int i = 0; i < max_iter; i++) {
        double nu3, derivative;

        if (m_options.chemistry == MocChemistry::PERFECT_GAS) {
            nu3 = prandtl_meyer(mach, m_options.gamma);
            // dnu/dM for perfect gas
            derivative = 2.0 * prandtl_meyer_derivative(mach, m_options.gamma);
        } else {
            auto [idx, weight] = pm_table.find_mach_index_and_weight(mach);
            nu3 = pm_table.interpolate_at_index(idx, weight, pm_table.nus);
            double V = pm_table.interpolate_at_index(idx, weight, pm_table.velocities);
            derivative = 2.0 * sqrt(mach * mach - 1.0) / V;
        }

        double residual = delta_theta - source_delta + (nu3 - p1.nu) + (nu3 - p2.nu);
        if (std::abs(residual) < m_options.abstol) return mach;

        mach -= residual / derivative;
    }

    // rootfinding has failed
    return -1.0;
}

void MocNozzle::update_thermodynamic_state(CharacteristicPoint& point) {
    switch (m_options.chemistry) {
        case MocChemistry::PERFECT_GAS: {
            // the choice of upstream point can be arbitrary due to Crocco's theorem
            // stagnation factor at throat is = 1 by definition
            double stagnation_ratio = stagnation_factor(point.mach, point.gamma_s);
            point.temperature = m_T_ref / stagnation_ratio;
            point.pressure = m_P_ref / pow(stagnation_ratio, point.gamma_s / (point.gamma_s - 1.0));
            break;
        }
        case MocChemistry::FROZEN:
        case MocChemistry::EQUILIBRIUM: {
            m_gas->thermo()->restoreState(point.cantera_state);
            point.temperature = m_gas->thermo()->temperature() / m_T_ref;
            point.pressure = m_gas->thermo()->pressure() / m_P_ref;
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of MocChemistry specified.");
    }
}

void MocNozzle::update_thermodynamic_state_from_nu(
    CharacteristicPoint& point, 
    double nu, double mach_guess) {
    point.nu = nu;
    switch (m_options.chemistry) {
        case MocChemistry::PERFECT_GAS: {
            point.gamma_s = m_options.gamma;
            point.mach = mach_from_prandtl_meyer(nu, point.gamma_s, mach_guess);
            break;            
        };
        case MocChemistry::FROZEN:
        case MocChemistry::EQUILIBRIUM: {
            auto [idx, weight] = pm_table.find_nu_index_and_weight(nu);
            point.gamma_s = pm_table.interpolate_at_index(idx, weight, pm_table.gamma_s);
            point.mach = pm_table.interpolate_at_index(idx, weight, pm_table.machs);
            point.cantera_state = pm_table.interpolate_state_at_index(idx, weight);
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of MocChemistry specified.");
    }
    update_thermodynamic_state(point);
    point.mu = mach_to_mu(point.mach);
} 

void MocNozzle::update_thermodynamic_state_from_mach(CharacteristicPoint& point, double mach) {
    point.mach = mach;

    switch (m_options.chemistry) {
        case MocChemistry::PERFECT_GAS: {
            point.gamma_s = m_options.gamma;
            point.nu = prandtl_meyer(point.mach, point.gamma_s);
            break;
        }
        case MocChemistry::FROZEN:
        case MocChemistry::EQUILIBRIUM: {
            auto [idx, weight] = pm_table.find_mach_index_and_weight(mach);
            point.gamma_s = pm_table.interpolate_at_index(idx, weight, pm_table.gamma_s);
            point.nu = pm_table.interpolate_at_index(idx, weight, pm_table.nus);
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of MocChemistry specified.");
    }
    update_thermodynamic_state(point);
    point.mu = mach_to_mu(point.mach);
} 

double MocNozzle::cplus_source_term(
    const CharacteristicPoint& p, 
    double new_y) const 
{
    
    double y_avg = 0.5 * (p.y + new_y);
    double dy = new_y - p.y;
    return sin(p.theta)/(p.mach * sin(p.theta - p.mu)) * dy/y_avg;
}

double MocNozzle::cminus_source_term(
    const CharacteristicPoint& p, 
    double new_y) const 
{
    
    double y_avg = 0.5 * (p.y + new_y);
    double dy = new_y - p.y;
    return -sin(p.theta)/(p.mach * sin(p.theta + p.mu)) * dy/y_avg;
}


} // namespace Goddard