#include <cmath>
#include <algorithm>
#include <optional>
#include <limits>
#include <format>
#include "goddard/equilibrium.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/moc_initialization.hpp"
#include "goddard/prandtlmeyer.hpp"
#include "goddard/error.hpp"

namespace Goddard {

// -- Logging helpers --
void MocNozzle::log_warning(const std::string& msg) {
    m_messages.push_back("Warning: " + msg);
}

void MocNozzle::log_info(const std::string& msg) {
    m_messages.push_back("Info: " + msg);
}

bool MocNozzle::is_solved() const {
    return m_is_solved;
}

// -- MocNozzle --
MocResult MocNozzle::solve() {
    m_messages.clear();
    m_theta_schedule.clear();
    m_is_solved = false;

    CharacteristicNet net;
    std::vector<CharacteristicPoint> data_line;

    if (m_options.chemistry == GasChemistry::PERFECT_GAS) {
        // Pure algebraic path: no Cantera dependency
        m_L_ref = m_options.geometry.throat_radius;
        m_P_ref = 1.0; // dimensionless stagnation pressure
        m_T_ref = 1.0; // dimensionless stagnation temperature

        data_line = generate_initial_data_line_perfect_gas(
            m_options.num_characteristics);
    } 
    else {
        // Cantera-backed path for frozen/equilibrium chemistry
        NozzleOptions nozzle_opts;
        nozzle_opts.chemistry = m_options.chemistry;
        Nozzle nozzle(*m_gas, nozzle_opts);
        auto throat = nozzle.solve_throat_conditions();
        m_gas->restore_state(throat.state);

        m_L_ref = m_options.geometry.throat_radius;
        m_P_ref = throat.P_inlet;
        m_T_ref = m_gas->thermo()->temperature();
        m_S_ref = throat.S_inlet;

        double a_throat = m_gas->speed_of_sound();
        
        pm_table.build_table(
            *m_gas->thermo(),
            m_options.chemistry == GasChemistry::EQUILIBRIUM,
            throat.S_inlet,
            throat.H_stagnation,
            a_throat);
                
        data_line = generate_initial_data_line(throat, m_options.geometry, m_options.num_characteristics);
    }
    
    if (m_options.mode == MocMode::DESIGN_MIN_LENGTH) {
        // Seed the throat-lip wall point at (0, 1). It is the origin of the centered
        // expansion fan and the start of the wall contour; the first computed wall point
        // is built relative to it. Without it leading_wall_point() has nothing to anchor.
        CharacteristicPoint throat_lip{};
        throat_lip.x = 0.0;
        throat_lip.y = 1.0;
        throat_lip.theta = m_options.theta_max;
        net.seed_wall_point(throat_lip);

        net.add_initial_characteristic(data_line);
    }
    else {
        net.add_initial_data_line(data_line);
    }

    //TODO: needs to be reworked. In the general case, we cannot separately solve the wall and internal regions.
    //this only applies to the special case of minimum length nozzles, where reflected characteristics are
    //cancelled at the wall.
    solve_characteristic_kernel(net);
    

    // Check for invalid points (details already logged by individual solvers)
    bool all_valid = true;
    for (const auto& pt : net.points) {
        if (pt.mach < 0.0) {
            all_valid = false;
            break;
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
    if (!net.points.empty()) {
        result.exit_mach = net.points.back().mach;
    }

    // Build wall profile from computed wall coordinates
    result.profile.x = net.wall_x;
    result.profile.y = net.wall_y;


    // TODO: this needs to be updated to be generalized to all nozzles, not just min length ones.

    // Exit plane extraction:
    // For a min-length nozzle, the exit plane is the last wavefront + last wall point.
    // The last wavefront contains the axis point, and each preceding wavefront's
    // last point was absorbed into wall calculations.
    if (!net.points.empty()) {
        const auto& outflows = net.outflow_points();
        for (const auto& pt : outflows) {
            result.exit_plane.y.push_back(pt.y);
            result.exit_plane.mach.push_back(pt.mach);
            result.exit_plane.theta.push_back(pt.theta);
            result.exit_plane.pressure.push_back(pt.pressure);
            result.exit_plane.temperature.push_back(pt.temperature);
            result.exit_plane.gamma_s.push_back(pt.gamma_s);
            result.exit_plane.velocity.push_back(pt.V);
        }
        // Add the last wall point for min length nozzle
        if (m_options.mode == MocMode::DESIGN_MIN_LENGTH) {
            const auto& last_wall_pt = net.leading_wall_point();
            result.exit_plane.y.push_back(last_wall_pt.y);
            result.exit_plane.mach.push_back(last_wall_pt.mach);
            result.exit_plane.theta.push_back(last_wall_pt.theta);
            result.exit_plane.pressure.push_back(last_wall_pt.pressure);
            result.exit_plane.temperature.push_back(last_wall_pt.temperature);
            result.exit_plane.gamma_s.push_back(last_wall_pt.gamma_s);
            result.exit_plane.velocity.push_back(last_wall_pt.V);
        }
    }
    m_is_solved = true;

    return result;
}

std::vector<CharacteristicPoint> MocNozzle::generate_initial_data_line(
    const ThroatCondition& throat,
    const ThroatGeometry& geometry,
    size_t num_points)
{
    auto thermo = m_gas->thermo();
    thermo->restoreState(throat.state);

    std::vector<CharacteristicPoint> data_line;
    data_line.reserve(num_points);
    
    ThermodynamicContext context = build_thermo_context();
    MocInitialization initializer{geometry, context, m_options};
    
    switch (m_options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            // The minimum length nozzle is a special case as the sonic line is straight, and
            // a centered expansion fan is the "exact" solution.
            // Therefore all initialization methods simplify to straight line initialization.
    
            auto expansion_line = initializer.initialize_centered_expansion(throat);
            CharacteristicPoint upstream_point;
            //for a minimum length nozzle, the sonic line is straight.
            //the initial data line is generated through a single Prandtl-Meyer expansion at the throat.
            for (size_t i = 0; i < expansion_line.size(); i++){
                const auto& expansion_point = expansion_line[i];
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
        case MocMode::DESIGN_RAO:
        case MocMode::ANALYSIS: {
             // the first theta should be small to minimize approximation error.
            CharacteristicPoint upstream_point;
            auto expansion_line = initializer.initialize_kliegel_levine(throat);
            
            //for a minimum length nozzle, the sonic line is straight.
            //the initial data line is generated through a single Prandtl-Meyer expansion at the throat.
            for (size_t i = 0; i < expansion_line.size(); i++){
                const auto& expansion_point = expansion_line[i];
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
        case MocMode::DESIGN_CENTERLINE: {
            throw NotImplementedError("Centerline nozzle initialization not implemented.");
        }
    }
    return data_line;
}

std::vector<CharacteristicPoint> MocNozzle::generate_initial_data_line_perfect_gas(
    size_t num_points)
{
    double gamma = m_options.gamma;

    // For analysis mode, derive theta_max from wall contour
    double theta_max = m_options.theta_max;
    if (m_options.mode == MocMode::ANALYSIS) {
        const auto& wall = m_options.nozzle_profile;
        if (wall.size() < 2) {
            throw std::runtime_error("Analysis mode requires a wall profile with at least 2 points.");
        }
        theta_max = wall.max_theta();
        if (theta_max <= 0.0) {
            throw std::runtime_error("Wall profile must have a positive expansion angle for analysis mode.");
        }
    }

    // Build or use theta schedule
    if (!m_options.theta_schedule.empty()) {
        m_theta_schedule = m_options.theta_schedule;
        num_points = m_theta_schedule.size();
    } else {
        // Auto-generate: small first step, then uniform spacing
        double dtheta_initial = theta_max / (num_points * 10);
        double dtheta = (theta_max - dtheta_initial) / (num_points - 1);
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

void MocNozzle::solve_characteristic_kernel(CharacteristicNet& net) {
    using Family = ChainMetadata::Family;
    LeadingEdgeView plus_edges = leading_edges(net, Family::PLUS);
    LeadingEdgeView minus_edges = leading_edges(net, Family::MINUS);
    
    std::vector<std::pair<CharacteristicPoint, PointMembership>> intersections;
    std::vector<size_t> paired_cminus;
    std::vector<bool> cminus_is_intersected(net.c_chains.size(), false);

    intersections.reserve(net.c_chains.size());
    paired_cminus.reserve(net.c_chains.size());


    // A correctly converging 2D MoC kernel needs O(N) marching passes. This cap is a
    // safety bound so a non-converging/runaway net terminates promptly instead of
    // spinning indefinitely; reaching it indicates a kernel that has not converged.
    const int maxiter = 2000;
    int iters = 0;

    while (net.has_active_chains() && iters < maxiter) {
        intersections.clear();
        paired_cminus.clear();
        // Sized to the *current* minus-edge view: reflections add chains over time, so a
        // one-shot allocation sized to the initial chain count would be indexed out of range.
        cminus_is_intersected.assign(minus_edges.chain_indices.size(), false);
        
        /* We use a simple geometric approach to determine which points are intersecting.
        This requires no assumptions about net topology or precomputed traversal maps and is
        fairly efficient except for very large N (100,000+) which are unrealistic for 2D MoC methods.
        
        1. For each C+ leading point we find the corresponding C- leading point with the 
        smallest y (height) that is larger than the C+ y and pair them. These points will intersect
        on the next pass.
        2. After all C+'s are paired, the bottommost unpaired C- pairs with the axis if it exists
        3. Wall interaction the topmost C+: either reflection or absorption depending on circumstances.
        4. Track active-chain count, terminate when it reaches 0 or all leading edges are past the outflow

        Invariant: each C- index should appear at most only once per set of intersections
        */
        for (size_t i = 0; i < plus_edges.chain_indices.size(); i++) {
            double plus_y = plus_edges.y_values[i];
            double best_dy = std::numeric_limits<double>::max();
            std::optional<size_t> best_partner = std::nullopt;
            size_t best_partner_pt_idx = 0;
            size_t best_partner_edgevec_idx = 0;

            // Find the closest C- above 
            for (size_t j = 0; j < minus_edges.chain_indices.size(); j++) {
                if (cminus_is_intersected[j]) continue; //skip if already paired
                
                double dy = minus_edges.y_values[j] - plus_y;
                if (dy > 0 && dy < best_dy) {
                    best_dy = dy;
                    best_partner = minus_edges.chain_indices[j];
                    best_partner_pt_idx = minus_edges.leading_pt_indices[j];
                    best_partner_edgevec_idx = j;
                }
                // edge case: multiple points at literally the same y due to an expansion fan
                // the tiebreaker is determined by the characteristic angle theta-mu
                if (best_partner.has_value() && best_dy == dy) [[unlikely]] {
                    const CharacteristicPoint& old_pt = net.points[best_partner_pt_idx];
                    const CharacteristicPoint& new_pt = net.points[minus_edges.leading_pt_indices[j]];
                    double old_angle = old_pt.theta - old_pt.mu;
                    double new_angle = new_pt.theta - new_pt.mu;
                    if (new_angle < old_angle) {
                        best_partner = minus_edges.chain_indices[j];
                        best_partner_pt_idx = minus_edges.leading_pt_indices[j];
                        best_partner_edgevec_idx = j;
                    }
                }
            }
            // intersect C+ with closest C- above it.
            if (best_partner.has_value()) [[likely]] {
                intersections.push_back({
                    solve_interior_point(
                        net.points[best_partner_pt_idx],
                        net.points[plus_edges.leading_pt_indices[i]]
                    ),
                    PointMembership {
                        .c_plus_chain_idx = plus_edges.chain_indices[i], 
                        .c_minus_chain_idx = *best_partner
                    }
                });
                paired_cminus.push_back(*best_partner);
                cminus_is_intersected[best_partner_edgevec_idx] = true;
            }
            else [[unlikely]] { // C+ intersects with wall
                intersections.push_back({
                    solve_wall_point(
                        net.points[plus_edges.leading_pt_indices[i]],
                        net.leading_wall_point(),
                        net.wall_point_indices.size() - 1
                    ),
                    PointMembership {
                        .c_plus_chain_idx = plus_edges.chain_indices[i],
                        .c_minus_chain_idx = std::nullopt
                    }
                });
            }
        }

        // find C- characteristic reflecting off axis.
        double min_cminus_y = std::numeric_limits<double>::max();
        std::optional<size_t> min_y_cminus_idx = std::nullopt;
        for (size_t i = 0; i < minus_edges.chain_indices.size(); i++) {
            if (minus_edges.y_values[i] < min_cminus_y) {
                min_cminus_y = minus_edges.y_values[i];
                min_y_cminus_idx = minus_edges.chain_indices[i];
            }
        }
        // if the lowest C- characteristic is unpaired, reflect it off the axis
        if (min_y_cminus_idx.has_value() &&
            std::none_of(
                paired_cminus.cbegin(), 
                paired_cminus.cend(), 
                [min_y_cminus_idx](size_t x) {return x == *min_y_cminus_idx;})
            )
        {
            intersections.push_back({
                solve_axis_point(net.leading_point(*min_y_cminus_idx)),
                PointMembership {
                    .c_plus_chain_idx = std::nullopt,
                    .c_minus_chain_idx = min_y_cminus_idx
                }
            });
        }
        // insert all intersections into net
        for (auto&& [pt, mem] : intersections) {
            if (m_options.mode == MocMode::DESIGN_MIN_LENGTH) { //min length nozzle is a special case
                if (!mem.c_minus_chain_idx.has_value()) {
                    // C+ reaches the wall: absorb it (no reflected wave). Cancelling the
                    // reflection here is exactly what produces a minimum-length contour.
                    net.terminate_c_plus_at_wall(*mem.c_plus_chain_idx, pt);
                    net.wall_x.push_back(pt.x);
                    net.wall_y.push_back(pt.y);
                }
                else if (!mem.c_plus_chain_idx.has_value()) {
                    // C- reaches the axis: reflect it into a new upward-marching C+.
                    net.reflect_c_minus_off_axis(*mem.c_minus_chain_idx, pt);
                }
                else {
                    net.add_point(pt, mem);
                }
            }
            else if (maximum_nozzle_length() >= pt.x) { //standard path
                if (!mem.c_minus_chain_idx.has_value()) {
                    // wall intersection
                    net.reflect_c_plus_off_wall(*mem.c_plus_chain_idx, pt);
                    net.wall_x.push_back(pt.x);
                    net.wall_y.push_back(pt.y);
                }
                else if (!mem.c_plus_chain_idx.has_value()) {
                    // axis intersection
                    net.reflect_c_minus_off_axis(*mem.c_minus_chain_idx, pt);
                } else {
                    net.add_point(pt, mem);
                }
            }
            else { // intersection falls outside of nozzle bounds
                if (mem.c_minus_chain_idx.has_value())
                    net.terminate_chain(*mem.c_minus_chain_idx, ChainMetadata::TerminationType::OUTFLOW);
                if (mem.c_plus_chain_idx.has_value())
                    net.terminate_chain(*mem.c_plus_chain_idx, ChainMetadata::TerminationType::OUTFLOW);
            }
        }
        update_leading_edges(plus_edges, net, Family::PLUS);
        update_leading_edges(minus_edges, net, Family::MINUS);

        iters++;
    }
    
}

CharacteristicPoint MocNozzle::solve_wall_point(
    const CharacteristicPoint& interior_parent,
    const CharacteristicPoint& previous_wall_point,
    int wall_point_index)
{
    switch (m_options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            // for minimum length design: theta is determined by theta schedule.
            // Guard against an empty schedule (the Cantera min-length path does not
            // populate m_theta_schedule) and against indices beyond it (the kernel
            // currently produces more wall hits than scheduled characteristics).
            double theta_wall = m_options.theta_max;
            if (!m_theta_schedule.empty()) {
                size_t k = std::min(static_cast<size_t>(wall_point_index),
                                    m_theta_schedule.size() - 1);
                theta_wall = m_options.theta_max - m_theta_schedule[k];
            }

            return solve_wall_point_design(interior_parent, previous_wall_point, theta_wall);
        }
        case MocMode::DESIGN_RAO: 
        case MocMode::ANALYSIS: {
            
            return solve_wall_point_analysis(interior_parent);
        }
        case MocMode::DESIGN_CENTERLINE: {
            throw NotImplementedError("DESIGN_CENTERLINE mode is not implemented.");
        }
        default: {
            throw std::runtime_error("Invalid MocMode in solve_wall_point");
        }
    }
    // unreachable
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
            // dtheta + dnu = S_cminus * dy
            // where S_cminus_limit = dtheta_dy / (M * sin(mu))
            double dy = off_axis_parent.y;
            double theta_avg = 0.5 * off_axis_parent.theta;
            double mu_avg = 0.5 * (off_axis_parent.mu + axis_point.mu);
            double M_avg = 0.5 * (off_axis_parent.mach + axis_point.mach);
            double source_limit = dtheta_dy / (M_avg * sin(theta_avg - mu_avg)) * dy;

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
    
    // We use a basic predictor-corrector scheme to find the wall intersection
    // Predictor: straight-line intersection
    // Linear scan is appropriate for tens to low hundreds of points. 
    // if npoints grows beyond that, switch to bisection algorithm
    return find_wall_hit(parent, wall, char_angle);
}

CharacteristicPoint MocNozzle::solve_interior_point(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2)
{
    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: {
            return solve_interior_point_planar(p1, p2);
        }
        case MocFlowKind::AXISYMMETRIC: {
            return solve_interior_point_axisymmetric(p1, p2);
        }
        default:
            throw std::runtime_error("Invalid MoC flow type specified.");
    }
} 

CharacteristicPoint MocNozzle::solve_interior_point_planar(
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

    auto [x, y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
    p3.x = x;
    p3.y = y;

    // validity checks
    if (p3.x < p1.x || p3.x < p2.x) {
        log_warning("Non-downstream intersection at ({}, {}). Parents at x=({}, {}).",
            p3.x, p3.y, p1.x, p2.x);
        p3.mach = -1.0; // invalid point
    }
    if (p3.mach < 1.0 && p3.mach > 0.0) {
        log_warning("Subsonic Mach {} at ({}, {}). Possible shock formation.", p3.mach, p3.x, p3.y);
        p3.mach = -1.0; // subsonic point -- characteristics undefined
    }

    return p3;
}
CharacteristicPoint MocNozzle::solve_interior_point_axisymmetric(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2)
{
    CharacteristicPoint p3;
    
    // Exact solution for axisymmetric flow
    // Based on derivation by Guentert & Neumann (1959)

    // Predictor: assume straight characteristics
    double c_minus_angle = p1.theta - p1.mu;
    double c_plus_angle = p2.theta + p2.mu;

    auto [x, y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
    p3.x = x;
    p3.y = y;

    double cotmu1 = 1.0/tan(p1.mu);
    double cotmu2 = 1.0/tan(p2.mu);
    // C+ source term
    // NOTE: watch out for singularities when p2 is on the centerline. 
    // The averaging of y should prevent issues.
    double L = tan(p2.mu)*sin(p2.mu)*sin(p2.theta)/(0.5*(p2.y+p3.y) * cos(p2.theta + p2.mu));
    // C- source term
    double M = tan(p1.mu)*sin(p1.mu)*sin(p1.theta)/(0.5*(p1.y+p3.y) * cos(p1.theta - p1.mu));
    update_thermodynamic_state_from_V(
        p3,
        1.0/(cotmu2/p2.V + cotmu1/p1.V) 
        * (cotmu2*(1+L*(p3.x - p2.x)) 
            + cotmu1*(1+M*(p3.x - p1.x)) 
            + p1.theta - p2.theta)
    );
    
    p3.theta = p2.theta + cotmu2 * (p3.V-p2.V)/p2.V - L*(p3.x - p2.x);

    // corrector: adjust for the characteristic angles at p3
    // ultimately h^2 error
    double c_minus_angle_corrected = average_cminus_angle(p1, p3);
    double c_plus_angle_corrected = average_cplus_angle(p2, p3);
    auto [x_corrected, y_corrected] = characteristic_intersection_with_angle(
        p1, p2, 
        c_minus_angle_corrected, 
        c_plus_angle_corrected);
    p3.x = x_corrected;
    p3.y = y_corrected;
    double mu1_avg = 0.5*(p1.mu + p3.mu);
    double mu2_avg = 0.5*(p2.mu + p3.mu);
    double theta1_avg = 0.5*(p1.theta + p3.theta);
    double theta2_avg = 0.5*(p2.theta + p3.theta);
    double M_avg = tan(mu1_avg)*sin(mu1_avg)*sin(theta1_avg)/(0.5*(p1.y+p3.y) * cos(c_minus_angle_corrected));
    double L_avg = tan(mu2_avg)*sin(mu2_avg)*sin(theta2_avg)/(0.5*(p2.y+p3.y) * cos(c_plus_angle_corrected));

    double cotmu1_avg = 1.0/tan(0.5*(p3.mu + p1.mu));
    double cotmu2_avg = 1.0/tan(0.5*(p3.mu + p2.mu));
    double V1_avg = 0.5*(p1.V + p3.V);
    double V2_avg = 0.5*(p2.V + p3.V);

    update_thermodynamic_state_from_V(
        p3,
        1.0/(cotmu2_avg/V2_avg + cotmu1_avg/V1_avg)
        * (cotmu2_avg*(p2.V/V2_avg + L_avg * (p3.x - p2.x))
            + cotmu1_avg*(p1.V/V1_avg + M_avg * (p3.x - p1.x))
            + p1.theta - p2.theta)
    );
    p3.theta = p2.theta + cotmu2_avg * ((p3.V - p2.V)/V2_avg - L_avg * (p3.x  - p2.x));
    p3.K_minus = p3.theta + p3.nu;
    p3.K_plus = p3.theta - p3.nu;

    // validity checks
    if (p3.x < p1.x || p3.x < p2.x) {
        log_warning("Non-downstream intersection at ({}, {}). Parents at x=({}, {}).",
            p3.x, p3.y, p1.x, p2.x);
        p3.mach = -1.0; // invalid point
    } else if (p3.mach < 1.0) {
        log_warning("Subsonic Mach {} at ({}, {}). Possible shock formation.", p3.mach, p3.x, p3.y);
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

    auto [x, y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
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
        auto [new_x, new_y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
        p3.x = new_x;
        p3.y = new_y;
        double S1_new = cminus_source_term(p1, new_y);
        double S2_new = cplus_source_term(p2, new_y);
        update_thermodynamic_state_from_mach(p3, find_node_mach(p1, p2, S1_new - S2_new, p3.mach));

        p3.theta = S1_new + p1.theta - (p3.nu - p1.nu);
        p3.K_plus = p3.theta - p3.nu;
        p3.K_minus = p3.theta + p3.nu;
        double residual = std::max(std::abs(S1_new - S1), std::abs(S2_new - S2));
        if (residual < m_options.solver_options.abstol) break;
        S1 = S1_new;
        S2 = S2_new;

        if (i == max_iters - 1 && residual >= m_options.solver_options.abstol) {
            log_warning("Iterative interior solver did not converge."
                "Residual={} after {} iterations at ({}, {}).", residual, max_iters, p3.x, p3.y);
        }
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
    // Predictor: planar K+ preservation
    // Flow solve: compatibility equation along C+ from interior parent
    // K+ is preserved: theta - nu = const along C+
    // this only applies for planar flow.

    // For axisymmetric flow, K+ is not preserved along C+.
    // Use predictor-corrector: start with planar approximation,
    // then correct using the source term.
    // The source term correction requires y, which is computed later
    // in solve_wall_point_design/analysis. We store the predictor result
    // and the corrector is applied in those methods.
    // For now, this gives a first-order approximation.

    wall_point.K_plus = interior_parent.K_plus;
    update_thermodynamic_state_from_nu(
        wall_point, 
        wall_point.theta - wall_point.K_plus, 
        interior_parent.mach);
    wall_point.K_minus = wall_point.theta + wall_point.nu;

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
    double c_plus_angle = average_cplus_angle(interior_parent, wall_point); 

    // the angle between the previous wall point and the current wall point
    // is given by wall_theta. 
    double prev_wall_angle = 0.5*(previous_wall_point.theta + wall_point.theta);

    auto [x,y] = characteristic_intersection_with_angle(
        interior_parent, previous_wall_point, 
        c_plus_angle, prev_wall_angle);

    wall_point.x = x; 
    wall_point.y = y;

    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: { break; }
        case MocFlowKind::AXISYMMETRIC: {
            // Wall corrector step
            // Unlike the case of an interior point, we only have one characteristic line
            // and theta is already known (specified by wall).
            // The general compatibility equation for the C+ characteristic is:
            // $theta_P - theta_B - S_{BP} = nu_P - nu_B$
            // where S is the source term, B is the interior point and P is the wall point.

            // The solution procedure is:
            // 1. Calculate source term from estimated y.
            // 2. Calculate nu from axisymmetric compatibility equation.
            // 3. Obtain thermodynamic values for given value of nu.
            // 4. Recompute new position based on new C+ angle
            // 5. Repeat until source residual is satisfactory. 
            // In practice only 1-2 iterations should be needed.
            double old_S = 0.0;
            double residual = 0.0;
            const int max_iter = 4;
            for (int i = 0; i < max_iter; i++){

                double S = cplus_source_term(interior_parent, wall_point);
                residual = std::abs(S - old_S);
                old_S = S;
                if (residual < m_options.solver_options.abstol) break;

                update_thermodynamic_state_from_nu(
                    wall_point,
                    wall_point.theta-interior_parent.theta - S + interior_parent.nu,
                    interior_parent.mach);

                wall_point.K_plus = wall_point.theta - wall_point.nu;
                wall_point.K_minus = wall_point.theta + wall_point.nu;
                double corrected_cplus_angle = average_cplus_angle(interior_parent, wall_point);
                auto [new_x,new_y] = characteristic_intersection_with_angle(
                    interior_parent,
                    previous_wall_point,
                    corrected_cplus_angle,
                    prev_wall_angle);

                wall_point.x = new_x;
                wall_point.y = new_y;

                if (i == max_iter - 1 && residual >= m_options.solver_options.abstol) {
                    log_warning("Wall design source term iteration did not converge." 
                        "Residual={} at ({}, {}).", residual, wall_point.x, wall_point.y);
                }
            }
            break;
        }
        default:
            throw std::runtime_error("Invalid flow type specified.");
    }

    return wall_point;
}

CharacteristicPoint MocNozzle::solve_wall_point_analysis(
        const CharacteristicPoint& interior_parent)
{
    const NozzleProfile& wall = m_options.nozzle_profile;
    auto [x_wall, y_wall] = intersect_characteristic_with_wall(
                interior_parent, wall);
            

    double wall_theta = wall.theta_at(x_wall);
    double char_angle = interior_parent.theta + interior_parent.mu;
    CharacteristicPoint wall_point = solve_wall_flow(interior_parent, wall_theta);
    
    wall_point.x = x_wall;
    wall_point.y = y_wall;
    // Corrector (for curved characteristics):
    // Solve the wall point flow using computed theta
    // then recompute characteristic slope as average of parent and wall-point slopes.
    // One correction is usually sufficient.
    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: {

            //for planar flow, solve_wall_flow equations are exact.
            double wall_char_angle = wall_point.theta + wall_point.mu;
            double average_char_angle = 0.5*(wall_char_angle + char_angle);

            auto [x_corrected, y_corrected] = find_wall_hit(wall_point, wall, average_char_angle);
            wall_point.x = x_corrected;
            wall_point.y = y_corrected;
            wall_theta = wall.theta_at(wall_point.x);
            wall_point = solve_wall_flow(interior_parent, wall_theta);
            break;
        }
        case MocFlowKind::AXISYMMETRIC: {
            // For axisymmetric flow, we must account for the source term.
            double old_S = 0.0;
            const int max_iter = 4;
            for (int i = 0; i < max_iter; i++){

                double c_plus_angle = average_cplus_angle(interior_parent, wall_point);
                auto [x_corrected, y_corrected] = find_wall_hit(interior_parent, wall, c_plus_angle);
                wall_point.x = x_corrected;
                wall_point.y = y_corrected;
                wall_theta = wall.theta_at(wall_point.x);

                double S = cplus_source_term(interior_parent, wall_point);
                double residual = std::abs(S - old_S);
                if (residual < m_options.solver_options.abstol) break;
                old_S = S;

                update_thermodynamic_state_from_nu(
                    wall_point,
                    wall_point.theta-interior_parent.theta - S + interior_parent.nu,
                    interior_parent.mach);

                if (i == max_iter - 1 && residual >= m_options.solver_options.abstol) {
                    log_warning("Wall design source term iteration did not converge." 
                        "Residual={} at ({}, {}).", residual, wall_point.x, wall_point.y);
                }
            }
            wall_point.K_plus = wall_point.theta - wall_point.nu;
            wall_point.K_minus = wall_point.theta + wall_point.nu;
            break;
        }
        default:
            throw std::runtime_error("Invalid flow type specified.");
    }
    return wall_point;
}

double MocNozzle::maximum_nozzle_length() const {
    return m_options.nozzle_profile.length();
}

double MocNozzle::gamma_s_from_mach(double mach) const {
    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            return m_options.gamma;
        }
        case GasChemistry::FROZEN: 
        case GasChemistry::EQUILIBRIUM: {
            return pm_table.interpolate_gamma_s_from_mach(mach);
        }
        default:
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    // unreachable
    return m_options.gamma;
}

double MocNozzle::gamma_s_from_nu(double nu) const {
    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            return m_options.gamma;
        }
        case GasChemistry::FROZEN: 
        case GasChemistry::EQUILIBRIUM: {
            return pm_table.interpolate_gamma_s_from_nu(nu);
        }
        default:
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    // unreachable
    return m_options.gamma;
}

double MocNozzle::mach_from_nu(const CharacteristicPoint& point, double mach_guess) const {
    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            return mach_from_prandtl_meyer(point.nu, point.gamma_s, mach_guess);
        }
        case GasChemistry::FROZEN: 
        case GasChemistry::EQUILIBRIUM: {
            return pm_table.interpolate_mach(point.nu);
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
}

double MocNozzle::nu_from_mach(double mach) const {
    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            return prandtl_meyer(mach, m_options.gamma);
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            return pm_table.interpolate_nu_from_mach(mach);
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
}

double MocNozzle::find_node_mach(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2,
    double source_delta,
    double mach_guess)
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

        if (m_options.chemistry == GasChemistry::PERFECT_GAS) {
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
        if (std::abs(residual) < m_options.solver_options.abstol) return mach;

        mach -= residual / derivative;
    }

    // rootfinding has failed
    log_warning("find_node_mach did not converge after {} iterations. "
        "Last Mach={}, delta_theta={}.", max_iter, mach, delta_theta);
    return -1.0;
}

void MocNozzle::update_thermodynamic_state(CharacteristicPoint& point) {
    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            // the choice of upstream point can be arbitrary due to Crocco's theorem
            // stagnation factor at throat is = 1 by definition
            double stagnation_ratio = stagnation_factor(point.mach, point.gamma_s);
            point.temperature = m_T_ref / stagnation_ratio;
            point.pressure = m_P_ref / pow(stagnation_ratio, point.gamma_s / (point.gamma_s - 1.0));
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            // Guard against points without a stored Cantera state (e.g. the bare
            // throat-lip seed); restoring an empty state vector aborts inside Cantera.
            if (point.cantera_state.empty()) break;
            m_gas->thermo()->restoreState(point.cantera_state);
            point.temperature = m_gas->thermo()->temperature() / m_T_ref;
            point.pressure = m_gas->thermo()->pressure() / m_P_ref;
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
}

void MocNozzle::update_thermodynamic_state_from_nu(
    CharacteristicPoint& point,
    double nu, double mach_guess)
{
    if (nu < 0.0) {
        log_warning("Negative Prandtl-Meyer angle nu = {}. Subsonic flow or numerical error.", nu);
    }
    point.nu = nu;
    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            point.gamma_s = m_options.gamma;
            point.mach = mach_from_prandtl_meyer(nu, point.gamma_s, mach_guess);
            point.V = point.mach;
            break;            
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            auto [idx, weight] = pm_table.find_nu_index_and_weight(nu);
            point.V = pm_table.interpolate_at_index(idx, weight, pm_table.velocities);
            point.gamma_s = pm_table.interpolate_at_index(idx, weight, pm_table.gamma_s);
            point.mach = pm_table.interpolate_at_index(idx, weight, pm_table.machs);
            point.cantera_state = pm_table.interpolate_state_at_index(idx, weight);
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    update_thermodynamic_state(point);
    point.mu = mach_to_mu(point.mach);
} 

void MocNozzle::update_thermodynamic_state_from_mach(CharacteristicPoint& point, double mach) {
    point.mach = mach;

    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            point.gamma_s = m_options.gamma;
            point.nu = prandtl_meyer(point.mach, point.gamma_s);
            point.V = mach;
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            auto [idx, weight] = pm_table.find_mach_index_and_weight(mach);
            point.V = pm_table.interpolate_at_index(idx, weight, pm_table.velocities);
            point.gamma_s = pm_table.interpolate_at_index(idx, weight, pm_table.gamma_s);
            point.nu = pm_table.interpolate_at_index(idx, weight, pm_table.nus);
            point.cantera_state = pm_table.interpolate_state_at_index(idx, weight);
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    update_thermodynamic_state(point);
    point.mu = mach_to_mu(point.mach);
}

void MocNozzle::update_thermodynamic_state_from_V(CharacteristicPoint& point, double V) {
    point.V = V;

    switch (m_options.chemistry) {
        case GasChemistry::PERFECT_GAS: {
            point.gamma_s = m_options.gamma;
            point.mach = V; // for a perfect gas, the velocity is kept dimensionless
            point.nu = prandtl_meyer(point.mach, point.gamma_s);
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            auto [idx, weight] = pm_table.find_V_index_and_weight(V);
            point.mach = pm_table.interpolate_at_index(idx, weight, pm_table.machs);
            point.gamma_s = pm_table.interpolate_at_index(idx, weight, pm_table.gamma_s);
            point.nu = pm_table.interpolate_at_index(idx, weight, pm_table.nus);
            point.cantera_state = pm_table.interpolate_state_at_index(idx, weight);
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    update_thermodynamic_state(point);
    point.mu = mach_to_mu(point.mach);
}

ThermodynamicContext MocNozzle::build_thermo_context() {
    return ThermodynamicContext{
        .gas = m_gas, 
        .table = pm_table, 
        .T_ref = m_T_ref, 
        .P_ref = m_P_ref, 
        .gamma_s = m_options.gamma};
}

std::pair<double,double> MocNozzle::find_wall_hit(
    const CharacteristicPoint& p, 
    const NozzleProfile& wall, 
    double char_angle) const
{
    double c_plus_slope = tan(char_angle);
    double x,y;
    bool found = false;

    for (size_t k = 0; k < wall.x.size() - 1; k++) {
        double wall_slope = (wall.y[k+1] - wall.y[k])
            / (wall.x[k+1] - wall.x[k]);
        double wall_intercept = wall.y[k] - wall_slope * wall.x[k];
        double char_intercept = p.y - c_plus_slope * p.x;

        double x_intercept = (wall_intercept - char_intercept) / (c_plus_slope - wall_slope);

        // check if intersection is within the wall segment
        if (x_intercept >= wall.x[k] && x_intercept <= wall.x[k+1]) {
            x = x_intercept;
            y = p.y + c_plus_slope * (x - p.x);
            found = true;
            break;
        }
    }
    if (!found) { 
        // if the intersection cannot be found within the nozzle profile, 
        // extrapolate from the furthest point. 
        size_t last = wall.x.size() - 1;
        double wall_slope = (wall.y[last] - wall.y[last-1])/(wall.x[last] - wall.x[last - 1]);
        double wall_intercept = wall.y[last] - wall_slope * wall.x[last];
        double char_intercept = p.y - c_plus_slope * p.x;
        x = (wall_intercept - char_intercept) / (c_plus_slope - wall_slope);
        y = p.y + c_plus_slope * (x - p.x);
        // log_warning("Extrapolating wall intercept to x = {}, y = {}", x, y);
    }
    return {x,y};
}

double MocNozzle::cplus_source_term(
    const CharacteristicPoint& p, 
    double new_y) const 
{
    double y_avg = 0.5 * (p.y + new_y);
    double dy = new_y - p.y;
    return sin(p.theta)/(p.mach * sin(p.theta - p.mu)) * dy/y_avg;
}

double MocNozzle::cplus_source_term(
    const CharacteristicPoint& p1, 
    const CharacteristicPoint& p3) const 
{
    double y_avg = 0.5 * (p1.y + p3.y);
    double dy = p3.y - p1.y;
    double theta_avg = 0.5 * (p1.theta + p3.theta);
    double theta_plus_mu_avg = 0.5 * ((p1.theta + p1.mu) + (p3.theta + p3.mu));
    double mach_avg = 0.5 * (p1.mach + p3.mach);
    return sin(theta_avg)/(mach_avg * sin(theta_plus_mu_avg)) * dy/y_avg;
}

double MocNozzle::cminus_source_term(
    const CharacteristicPoint& p, 
    double new_y) const 
{
    
    double y_avg = 0.5 * (p.y + new_y);
    double dy = new_y - p.y;
    return -sin(p.theta)/(p.mach * sin(p.theta - p.mu)) * dy/y_avg;
}

MocNozzle::LeadingEdgeView MocNozzle::leading_edges(const CharacteristicNet& net, ChainMetadata::Family fam) const {
    LeadingEdgeView view;
    view.family = fam;
    for (size_t i = 0; i < net.chain_metadata.size(); i++) {
        const ChainMetadata& meta = net.chain_metadata[i];
        if (meta.active && meta.family == fam) {
            view.y_values.push_back(net.points[meta.latest_point_idx].y);
            view.chain_indices.push_back(i);
            view.leading_pt_indices.push_back(meta.latest_point_idx);
        }
    }
    return view;
}

void MocNozzle::update_leading_edges(LeadingEdgeView& view, const CharacteristicNet& net, ChainMetadata::Family family) const {
    view.y_values.clear();
    view.chain_indices.clear();
    view.leading_pt_indices.clear();
    for (size_t i = 0; i < net.chain_metadata.size(); i++) {
        const ChainMetadata& meta = net.chain_metadata[i];
        if (meta.active && meta.family == view.family) {
            view.y_values.push_back(net.points[meta.latest_point_idx].y);
            view.chain_indices.push_back(i);
            view.leading_pt_indices.push_back(meta.latest_point_idx);
        }
    }
}


} // namespace Goddard