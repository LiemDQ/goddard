#include <cmath>
#include <algorithm>
#include <numeric>
#include <optional>
#include <limits>
#include <format>
#include <iostream>
#include "goddard/equilibrium.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/moc_nozzle.hpp"
#include "goddard/moc_initialization.hpp"
#include "goddard/prandtlmeyer.hpp"
#include "goddard/error.hpp"

namespace Goddard {

namespace {
// Unwraps a PointResult produced while constructing the initial data line. Unlike
// the kernel (which communicates a unit-process failure back to solve() via a
// returned MocFailure), the initial data line runs before the kernel and
// communicates failure through solve()'s exception boundary: throwing here is
// caught by MocNozzle::solve() and converted to MocFailure{INITIALIZATION_FAILED}.
CharacteristicPoint unwrap_or_throw(const PointResult& result, std::string_view context) {
    if (result.error != MocErrorCode::NONE) {
        throw ConvergenceError(std::format(
            "{}: unit process failed at ({}, {}) with error {}.",
            context, result.point.x, result.point.y, to_string(result.error)));
    }
    return result.point;
}
} // namespace

// -- Logging helpers --
void MocNozzle::log_warning(const std::string& msg) {
    m_messages.push_back("Warning: " + msg);
}

void MocNozzle::log_info(const std::string& msg) {
    m_messages.push_back("Info: " + msg);
}

void MocNozzle::log_debug(const std::string& msg) {
    if (m_options.log_level != MocLogLevel::DEBUG) return;
    std::string full = "Debug: " + msg;
    m_messages.push_back(full);
    std::cerr << full << std::endl;
}

bool MocNozzle::is_solved() const {
    return m_is_solved;
}

// -- MocNozzle --
MocResult MocNozzle::solve() {
    m_messages.clear();
    m_theta_schedule.clear();
    m_initial_line_family.reset();
    m_is_solved = false;

    CharacteristicNet net;
    std::vector<CharacteristicPoint> data_line;
    m_options.nozzle_profile = setup_nozzle_profile(m_options.geometry);

    ThroatCondition throat{};

    if (m_options.chemistry == GasChemistry::PERFECT_GAS) {
        // Pure algebraic path: no Cantera dependency
        m_L_ref = m_options.geometry.throat_radius;
        m_P_ref = 1.0; // dimensionless stagnation pressure
        m_T_ref = 1.0; // dimensionless stagnation temperature

        // throat populated with dummy values; they are not used in the perfect gas case.
        throat = ThroatCondition{
            .converged = true,
            .speed_of_sound = 1.0,
            .H_stagnation = 1.0,
            .P_inlet = m_P_ref,
            .S_inlet = 1.0,
            .gamma_s = m_options.gamma,
            .dlV_dlP_T = -1.0,
            .dlV_dlT_P = 1.0,
            .state = {}
        };
    }
    else {
        // Cantera-backed path for frozen/equilibrium chemistry
        NozzleOptions nozzle_opts;
        nozzle_opts.chemistry = m_options.chemistry;
        Nozzle nozzle(*m_gas, nozzle_opts);
        throat = nozzle.solve_throat_conditions();
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
    }

    MocResult result;

    // The initial data line (a Newton solve for the Kliegel-Levine transonic line, a
    // centered-fan Prandtl-Meyer inversion, or a PrandtlMeyerTable lookup) can fail
    // numerically before the kernel ever starts marching. Converting those failures
    // here, rather than letting them propagate as exceptions, is what makes solve()
    // never throw for a numerical/convergence failure.
    try {
        data_line = generate_initial_data_line(throat, m_options.geometry, m_options.num_characteristics);

        // The initial-data-line generators declared whether the line lies along a single
        // characteristic (m_initial_line_family). A centered expansion fan does, and needs a
        // wall anchor at the throat lip (0, 1): it seeds leading_wall_point() for the first wall
        // solve and gives the area ratio its throat reference. A transonic start line spans
        // axis-to-wall and already carries its own wall point, so no separate anchor is seeded.
        if (m_initial_line_family.has_value()) {
            CharacteristicPoint throat_lip{};
            throat_lip.x = 0.0;
            throat_lip.y = 1.0;
            // The lip angle anchors the first wall solve in DESIGN_MIN_LENGTH, where it equals
            // theta_max. In analysis the wall angles come from the profile and the lip is purely
            // an anchor, so the topmost fan angle is a harmless stand-in.
            throat_lip.theta = m_options.mode == MocMode::DESIGN_MIN_LENGTH
                ? m_options.theta_max
                : (data_line.empty() ? 0.0 : data_line.back().theta);
            net.seed_wall_point(throat_lip);
        }
        net.add_initial_data_line(data_line, m_initial_line_family);
    }
    catch (const NotImplementedError&) {
        // Programmer error (unimplemented mode): not a numerical failure, keep throwing.
        throw;
    }
    catch (const std::out_of_range& e) {
        result.failure = MocFailure{
            MocErrorCode::TABLE_RANGE_EXCEEDED,
            std::string("Initial data line construction failed: ") + e.what(),
            0.0, 0.0, -1
        };
    }
    catch (const std::exception& e) {
        result.failure = MocFailure{
            MocErrorCode::INITIALIZATION_FAILED,
            std::string("Initial data line construction failed: ") + e.what(),
            0.0, 0.0, -1
        };
    }

    if (result.failure.code != MocErrorCode::NONE) {
        result.converged = false;
        result.messages = m_messages;
        m_is_solved = true;
        return result;
    }

    std::optional<MocFailure> kernel_failure = solve_characteristic_kernel(net);

    result.net = net;
    result.messages = m_messages;
    if (kernel_failure.has_value()) {
        result.failure = *kernel_failure;
        result.converged = false;
    } else {
        result.converged = true;
    }

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

    // Exit Mach:
    // For a min-length nozzle the exit flow is uniform, so the last computed point
    // (on the final characteristic) is representative.
    // In analysis mode the outflow boundary is a ragged staircase of terminated
    // chains and the last computed point is arbitrary; report the centerline exit
    // Mach (most downstream axis point) instead.
    if (!net.points.empty()) {
        if (m_options.mode != MocMode::DESIGN_MIN_LENGTH && !net.axis_point_indices.empty()) {
            result.exit_mach = net.leading_axis_point().mach;
        } else {
            result.exit_mach = net.points.back().mach;
        }
    }

    // Build wall profile from computed wall coordinates
    result.profile.x = net.wall_x;
    result.profile.y = net.wall_y;



    // Exit plane extraction:
    // For a min-length nozzle, the exit plane is the last wavefront + last wall point.
    // The last wavefront contains the axis point, and each preceding wavefront's
    // last point was absorbed into wall calculations.
    if (!net.points.empty()) {
        const auto& outflows = net.outflow_points();
        // Add the last characteristic for min length nozzle. Guarded by
        // emptiness: a kernel that now fails fast can abort before any axis/wall
        // point has been recorded, where leading_axis_point()/leading_wall_point()
        // would otherwise index an empty vector.
        if (m_options.mode == MocMode::DESIGN_MIN_LENGTH) {
            if (!net.axis_point_indices.empty()) {
                const auto& last_axis_pt = net.leading_axis_point();
                result.exit_plane.y.push_back(last_axis_pt.y);
                result.exit_plane.mach.push_back(last_axis_pt.mach);
                result.exit_plane.theta.push_back(last_axis_pt.theta);
                result.exit_plane.pressure.push_back(last_axis_pt.pressure);
                result.exit_plane.temperature.push_back(last_axis_pt.temperature);
                result.exit_plane.gamma_s.push_back(last_axis_pt.gamma_s);
                result.exit_plane.velocity.push_back(last_axis_pt.V);
            }

            if (!net.wall_point_indices.empty()) {
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
        else {
            for (const auto& pt : outflows) {
                result.exit_plane.y.push_back(pt.y);
                result.exit_plane.mach.push_back(pt.mach);
                result.exit_plane.theta.push_back(pt.theta);
                result.exit_plane.pressure.push_back(pt.pressure);
                result.exit_plane.temperature.push_back(pt.temperature);
                result.exit_plane.gamma_s.push_back(pt.gamma_s);
                result.exit_plane.velocity.push_back(pt.V);
            }
        }
    }
    m_is_solved = true;

    return result;
}


auto MocNozzle::setup_nozzle_profile(const NozzleGeometry& geometry) -> NozzleProfile {
    switch (m_options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            return m_options.nozzle_profile;
        }
        case MocMode::DESIGN_RAO: {
            return NozzleProfile::generate_Rao_TOP_nozzle(geometry.expansion_ratio, 1.0, geometry.length_fraction);
        }
        case MocMode::ANALYSIS: {
            return m_options.nozzle_profile;
        }
        case MocMode::DESIGN_CENTERLINE: {
            throw NotImplementedError("DESIGN_CENTERLINE is not implemented.");
        }
        default: {
            //unreachable
            throw std::runtime_error("Invalid value for MocMode specified.");
        }
    }
}


std::vector<CharacteristicPoint> MocNozzle::generate_initial_data_line(
    const ThroatCondition& throat,
    const NozzleGeometry& geometry,
    size_t num_points)
{
    std::vector<CharacteristicPoint> data_line;
    data_line.reserve(num_points);
    if (!m_options.theta_schedule.empty()) {
        m_theta_schedule = m_options.theta_schedule;
    }
    
    if (m_gas.has_value()) {
        m_gas->thermo()->restoreState(throat.state);
    }

    ThermodynamicContext context = build_thermo_context();

    // The Kliegel-Levine series implemented here is the axisymmetric transonic
    // solution; it is selected by a positive downstream wall curvature radius but
    // is not valid for planar flow, which falls back to a centered fan.
    bool use_centered_fan =
        m_options.geometry.downstream_wall_curvature_radius <= 0.0 ||
        m_options.flow_type == MocFlowKind::PLANAR;
    if (use_centered_fan &&
        m_options.flow_type == MocFlowKind::PLANAR &&
        m_options.geometry.downstream_wall_curvature_radius > 0.0 &&
        m_options.mode != MocMode::DESIGN_MIN_LENGTH)
    {
        log_info("Kliegel-Levine initialization is only valid for axisymmetric flow; "
                 "using centered-fan initialization for planar flow.");
    }

    // In analysis mode the expansion is set by the wall contour, not by a design
    // theta_max (which is meaningless there and typically left unset). When the
    // centered-fan initializer will be used, derive theta_max from the profile.
    MocOptions init_options = m_options;
    if (m_options.mode == MocMode::ANALYSIS && use_centered_fan)
    {
        const NozzleProfile& wall = m_options.nozzle_profile;
        if (wall.size() < 2) {
            throw std::runtime_error("Analysis mode requires a wall profile with at least 2 points.");
        }
        init_options.theta_max = wall.max_theta();
        if (init_options.theta_max <= 0.0) {
            throw std::runtime_error("Wall profile must have a positive expansion angle for analysis mode.");
        }
    }
    MocInitialization initializer{geometry, context, init_options};

    switch (m_options.mode) {
        case MocMode::DESIGN_MIN_LENGTH: {
            // The minimum length nozzle is a special case as the sonic line is straight, and
            // a centered expansion fan is the "exact" solution.
            // Therefore all initialization methods simplify to straight line initialization.

            // A centered expansion fan is collinear along a single C+ characteristic.
            m_initial_line_family = ChainMetadata::Family::PLUS;

            auto expansion_line = initializer.initialize_centered_expansion(throat);
            // Mirror the fan angles (user-supplied or auto-generated inside the
            // initializer) into m_theta_schedule: the min-length wall solve reads
            // theta_wall = theta_max - m_theta_schedule[k] for each wall point.
            m_theta_schedule.resize(expansion_line.size());
            for (size_t i = 0; i < expansion_line.size(); i++) {
                m_theta_schedule[i] = expansion_line[i].theta;
            }
            CharacteristicPoint upstream_point;
            //for a minimum length nozzle, the sonic line is straight.
            //the initial data line is generated through a single Prandtl-Meyer expansion at the throat.
            for (size_t i = 0; i < expansion_line.size(); i++){
                const auto& expansion_point = expansion_line[i];
                // special case: the first point is on the centerline
                if (i == 0) {
                    upstream_point = unwrap_or_throw(
                        solve_initial_axis_point_centered_exp(expansion_point),
                        "Initial data line (min-length axis point)");
                }
                else {
                    upstream_point = unwrap_or_throw(
                        solve_interior_point(expansion_point, upstream_point),
                        "Initial data line (min-length interior point)");
                }
                data_line.push_back(upstream_point);
            }
            break;
        }
        case MocMode::DESIGN_RAO:
        case MocMode::ANALYSIS: {
            // The Kliegel-Levine transonic start line crosses many characteristics, so it is
            // seeded as a generic data line (m_initial_line_family stays empty).
             // the first theta should be small to minimize approximation error.
            CharacteristicPoint upstream_point;
            if (use_centered_fan) {

                m_initial_line_family = ChainMetadata::Family::PLUS;
                auto expansion_line = initializer.initialize_centered_expansion(throat);
                m_theta_schedule.resize(expansion_line.size());
                for (size_t i = 0; i < expansion_line.size(); i++) {
                    m_theta_schedule[i] = expansion_line[i].theta;
                }
                //for a minimum length nozzle, the sonic line is straight.
                //the initial data line is generated through a single Prandtl-Meyer expansion at the throat.
                for (size_t i = 0; i < expansion_line.size(); i++){
                    const auto& expansion_point = expansion_line[i];
                    // special case: the first point is on the centerline
                    if (i == 0) {
                        upstream_point = unwrap_or_throw(
                            solve_initial_axis_point_centered_exp(expansion_point),
                            "Initial data line (centered-fan axis point)");
                        data_line.push_back(upstream_point);
                    }
                    else {
                        upstream_point = unwrap_or_throw(
                            solve_interior_point(expansion_point, upstream_point),
                            "Initial data line (centered-fan interior point)");
                        data_line.push_back(upstream_point);
                    }
                }
            }
            else {
                // The Kliegel-Levine transonic line already spans axis-to-wall with full
                // thermodynamic state and K+/K- set: it IS the initial data line. Unlike
                // the centered fan (whose rays all emanate from the throat lip and must be
                // marched into the flow field), re-marching these points pairwise would
                // intersect characteristics *behind* their parents. See
                // CharacteristicNet::add_initial_data_line's nullopt branch for how this
                // non-collinear line is actually seeded into the net without that re-marching.
                data_line = initializer.initialize_kliegel_levine(throat);
            }
            break;
        }
        case MocMode::DESIGN_CENTERLINE: {
            throw NotImplementedError("Centerline nozzle initialization not implemented.");
        }
    }

    for (size_t i = 0; i < data_line.size(); i++) {
        const auto& pt = data_line[i];
        log_debug("INIT[{}] x={:.6f} y={:.6f} theta={:.6f} nu={:.6f} mach={:.6f} mu={:.6f}",
            i, pt.x, pt.y, pt.theta, pt.nu, pt.mach, pt.mu);
        // The Kliegel-Levine and Sauer transonic-line initializers compute Mach
        // directly (not via Prandtl-Meyer inversion) and are not funneled through a
        // unit process that would otherwise validate them; check them here so a
        // breakdown of the small-perturbation series (e.g. a subsonic or non-finite
        // result near the wall) is reported as MocFailure::INITIALIZATION_FAILED
        // instead of an invalid point silently entering the net.
        MocErrorCode validity = check_point_validity(pt, m_options.solver_options.abstol);
        if (validity != MocErrorCode::NONE) {
            throw ConvergenceError(std::format(
                "Initial data line point {} at (x={}, y={}) failed validity check: {}.",
                i, pt.x, pt.y, to_string(validity)));
        }
    }
    return data_line;
}

std::optional<MocFailure> MocNozzle::solve_characteristic_kernel(CharacteristicNet& net) {
    using Family = ChainMetadata::Family;
    LeadingEdgeView plus_edges = leading_edges(net, Family::PLUS);
    sort_plus_edges_by_proximity(plus_edges);
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

    // Populated the instant any unit process reports a numerical failure; the march
    // aborts immediately rather than continuing to build on top of an invalid point
    // (which is what let corrupted points silently propagate through the net before).
    std::optional<MocFailure> failure;

    while (net.has_active_chains() && iters < maxiter && !failure.has_value()) {
        log_debug("--- kernel pass {}: {} active C+, {} active C- ---",
            iters, plus_edges.chain_indices.size(), minus_edges.chain_indices.size());
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
        for (size_t i = 0; i < plus_edges.chain_indices.size() && !failure.has_value(); i++) {
            double plus_y = plus_edges.y_values[i];
            double best_dy = std::numeric_limits<double>::max();
            std::optional<size_t> best_partner = std::nullopt;
            size_t best_partner_pt_idx = 0;
            size_t best_partner_edgevec_idx = 0;
            bool any_cminus_above = false;

            // Find the closest C- above
            for (size_t j = 0; j < minus_edges.chain_indices.size(); j++) {
                if (minus_edges.y_values[j] - plus_y > 0) any_cminus_above = true;
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
                const CharacteristicPoint& minus_pt = net.points[best_partner_pt_idx];
                const CharacteristicPoint& plus_pt = net.points[plus_edges.leading_pt_indices[i]];
                PointResult result = solve_interior_point(minus_pt, plus_pt);
                if (result.error != MocErrorCode::NONE) {
                    failure = MocFailure{
                        result.error,
                        std::format("Interior point solve failed ({}) pairing "
                            "minus(x={:.6f},y={:.6f}) with plus(x={:.6f},y={:.6f}).",
                            to_string(result.error), minus_pt.x, minus_pt.y, plus_pt.x, plus_pt.y),
                        result.point.x, result.point.y,
                        iters
                    };
                    break;
                }
                log_debug("PAIR minus(x={:.6f},y={:.6f},th={:.6f},mu={:.6f}) "
                    "plus(x={:.6f},y={:.6f},th={:.6f},mu={:.6f}) -> (x={:.6f},y={:.6f},mach={:.6f})",
                    minus_pt.x, minus_pt.y, minus_pt.theta, minus_pt.mu,
                    plus_pt.x, plus_pt.y, plus_pt.theta, plus_pt.mu,
                    result.point.x, result.point.y, result.point.mach);
                intersections.push_back({
                    result.point,
                    PointMembership {
                        .c_plus_chain_idx = plus_edges.chain_indices[i],
                        .c_minus_chain_idx = *best_partner
                    }
                });
                paired_cminus.push_back(*best_partner);
                cminus_is_intersected[best_partner_edgevec_idx] = true;
            }
            else if (any_cminus_above) {
                // A C- does exist above this C+, but a closer competitor already claimed it
                // this pass (e.g. many individual C+ chains from a Kliegel-Levine transonic
                // line, competing for a single C- freshly born from a wall reflection). Leave
                // this chain active and retry once that C- has advanced on the next pass,
                // rather than wrongly treating it as a wall hit.
                const CharacteristicPoint& plus_pt = net.points[plus_edges.leading_pt_indices[i]];
                log_debug("SKIP plus(x={:.6f},y={:.6f}) waiting for scarce C- partner "
                    "(already claimed this pass)", plus_pt.x, plus_pt.y);
            }
            else [[unlikely]] { // C+ intersects with wall

                const CharacteristicPoint& plus_pt = net.points[plus_edges.leading_pt_indices[i]];
                std::optional<PointResult> maybe_result = solve_wall_point(
                    plus_pt,
                    net.leading_wall_point(),
                    static_cast<int>(net.wall_point_indices.size()) - 1
                );

                if (maybe_result.has_value()) {
                    if (maybe_result->error != MocErrorCode::NONE) {
                        failure = MocFailure{
                            maybe_result->error,
                            std::format("Wall point solve failed ({}) for plus(x={:.6f},y={:.6f}).",
                                to_string(maybe_result->error), plus_pt.x, plus_pt.y),
                            maybe_result->point.x, maybe_result->point.y,
                            iters
                        };
                        break;
                    }
                    log_debug("WALL plus(x={:.6f},y={:.6f}) -> wall hit at (x={:.6f},y={:.6f})",
                        plus_pt.x, plus_pt.y, maybe_result->point.x, maybe_result->point.y);
                    intersections.push_back({
                        maybe_result->point,
                        PointMembership {
                            .c_plus_chain_idx = plus_edges.chain_indices[i],
                            .c_minus_chain_idx = std::nullopt
                        }
                    });
                }
                else {
                    // if no intersection found, this point intersects beyond the profile boundary
                    log_debug("OUTFLOW plus(x={:.6f},y={:.6f}) terminates: "
                        "no wall hit found within profile bounds", plus_pt.x, plus_pt.y);
                    net.terminate_chain(
                        plus_edges.chain_indices[i],
                        ChainMetadata::TerminationType::OUTFLOW
                    );
                }
            }
        }

        if (failure.has_value()) break;

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
            const CharacteristicPoint& minus_pt = net.leading_point(*min_y_cminus_idx);
            PointResult axis_result = solve_axis_point(minus_pt);
            if (axis_result.error != MocErrorCode::NONE) {
                failure = MocFailure{
                    axis_result.error,
                    std::format("Axis point solve failed ({}) for minus(x={:.6f},y={:.6f}).",
                        to_string(axis_result.error), minus_pt.x, minus_pt.y),
                    axis_result.point.x, axis_result.point.y,
                    iters
                };
                break;
            }
            log_debug("AXIS minus(x={:.6f},y={:.6f}) reflects off axis", minus_pt.x, minus_pt.y);
            intersections.push_back({
                axis_result.point,
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
                log_debug("OUTFLOW pt(x={:.6f},y={:.6f}) terminates: beyond "
                    "maximum_nozzle_length() ({:.6f})", pt.x, pt.y, maximum_nozzle_length());
                if (mem.c_minus_chain_idx.has_value())
                    net.terminate_chain(*mem.c_minus_chain_idx, ChainMetadata::TerminationType::OUTFLOW);
                if (mem.c_plus_chain_idx.has_value())
                    net.terminate_chain(*mem.c_plus_chain_idx, ChainMetadata::TerminationType::OUTFLOW);
            }
        }
        update_leading_edges(plus_edges, net, Family::PLUS);
        sort_plus_edges_by_proximity(plus_edges);
        update_leading_edges(minus_edges, net, Family::MINUS);

        iters++;
    }

    if (failure.has_value()) {
        return failure;
    }

    if (iters >= maxiter) {
        return MocFailure{
            MocErrorCode::MAX_ITERATIONS_REACHED,
            std::format("Kernel reached the iteration safety cap ({}) without all "
                "characteristics terminating.", maxiter),
            0.0, 0.0, iters
        };
    }

    return std::nullopt;
}

std::optional<PointResult> MocNozzle::solve_wall_point(
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

PointResult MocNozzle::solve_initial_axis_point_centered_exp(const CharacteristicPoint& expansion_point) {
    CharacteristicPoint point{};
    point.y = 0.0; //point always lies on axis
    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: {
            point.theta = expansion_point.theta;
            // nu = expansion point nu follows from geometric analysis
            MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
                point, expansion_point.nu, expansion_point.mach);
            if (thermo_error != MocErrorCode::NONE) {
                point.x = expansion_point.x;
                return {point, thermo_error};
            }
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
            MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
                point, expansion_point.nu, expansion_point.mach);
            if (thermo_error != MocErrorCode::NONE) {
                point.x = expansion_point.x;
                return {point, thermo_error};
            }
            point.K_plus = point.theta - point.nu;
            point.K_minus = expansion_point.K_minus;

            double c_minus_angle = average_cminus_angle(expansion_point, point);
            point.x = expansion_point.x - expansion_point.y / tan(c_minus_angle);
            break;
        }
        default:
            throw std::runtime_error("Invalid flow type specified.");
    }

    MocErrorCode validity = check_point_validity(point, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        return {point, validity};
    }
    return {point, MocErrorCode::NONE};
}

PointResult MocNozzle::solve_axis_point(const CharacteristicPoint& off_axis_parent) {
    CharacteristicPoint axis_point{};
    axis_point.y = 0.0;
    axis_point.theta = 0.0; // symmetry condition

    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: {

            axis_point.K_minus = off_axis_parent.K_minus;
            MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
                axis_point,
                axis_point.K_minus - axis_point.theta,
                off_axis_parent.mach);
            if (thermo_error != MocErrorCode::NONE) {
                axis_point.x = off_axis_parent.x;
                return {axis_point, thermo_error};
            }

            axis_point.K_plus = axis_point.theta - axis_point.nu;

            double c_minus_angle = 0.5 * (off_axis_parent.theta + axis_point.theta)
                - 0.5 * (off_axis_parent.mu + axis_point.mu);

            axis_point.x = off_axis_parent.x - off_axis_parent.y/tan(c_minus_angle);
            break;
        }
        case MocFlowKind::AXISYMMETRIC: {

            const double theta = off_axis_parent.theta;
            const double mu = off_axis_parent.mu;
            // Predictor method recommended by Zucrow & Hoffman, ch. 17
            // Use parent point as initial guess for source term.
            // The C- compatibility per dy is d(theta+nu) = [sin(theta)/(M sin(theta-mu))] (dy/y).
            // The descent from the parent to the axis has dy = -y_parent (signed), which makes
            // the source contribution positive (sin(theta-mu) < 0).
            axis_point.K_minus = off_axis_parent.K_minus;

            const double dy = -off_axis_parent.y; // y_axis - y_parent
            double source_pred = sin(theta) / (off_axis_parent.mach * sin(theta - mu))
                * dy / off_axis_parent.y;

            // dtheta + dnu = S => nu_2 - nu_1 = S + theta_1 = S + K_minus_1
            double nu_pred = axis_point.K_minus + source_pred;

            MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
                axis_point,
                nu_pred,
                off_axis_parent.mach);
            if (thermo_error != MocErrorCode::NONE) {
                axis_point.x = off_axis_parent.x;
                return {axis_point, thermo_error};
            }

            axis_point.K_plus = axis_point.theta - axis_point.nu;

            double c_minus_angle = average_cminus_angle(off_axis_parent, axis_point);
            axis_point.x = off_axis_parent.x - off_axis_parent.y / tan(c_minus_angle);

            // Corrector: apply the averaged source term with the same signed dy:
            // dy/y_avg = (0 - y_parent)/(y_parent/2) = -2, and theta_avg = theta_parent/2,
            // which together reproduce the finite sin(theta)/y limit at the axis.
            double dy_over_y_avg = -2.0;
            double theta_avg = 0.5 * off_axis_parent.theta;
            double mu_avg = 0.5 * (off_axis_parent.mu + axis_point.mu);
            double M_avg = 0.5 * (off_axis_parent.mach + axis_point.mach);
            double source_corrected = sin(theta_avg)/ (M_avg * sin(theta_avg - mu_avg)) * dy_over_y_avg;

            // Corrected nu: K_minus from parent, + source contribution
            double nu_corrected = off_axis_parent.K_minus + source_corrected;
            thermo_error = update_thermodynamic_state_from_nu(axis_point, nu_corrected, axis_point.mach);
            if (thermo_error != MocErrorCode::NONE) {
                return {axis_point, thermo_error};
            }
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

    // solve_axis_point has no non-downstream guard today (a documented gap): the
    // resulting axis point's x must not precede its parent's.
    if (axis_point.x < off_axis_parent.x) {
        log_warning("Non-downstream axis intersection at ({}, {}). Parent at x={}.",
            axis_point.x, axis_point.y, off_axis_parent.x);
        return {axis_point, MocErrorCode::NON_DOWNSTREAM_POINT};
    }

    MocErrorCode validity = check_point_validity(axis_point, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        log_warning("Invalid axis point at ({}, {}): error code {}.",
            axis_point.x, axis_point.y, to_string(validity));
        return {axis_point, validity};
    }

    return {axis_point, MocErrorCode::NONE};
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

PointResult MocNozzle::solve_interior_point(
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

PointResult MocNozzle::solve_interior_point_planar(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2)
{
    CharacteristicPoint p3{};
    p3.theta = 0.5*(p1.K_minus + p2.K_plus);
    double mach_guess = 0.5*(p1.mach + p2.mach);

    // Position (x, y) is computed after the thermo update below, so a thermo
    // failure here leaves p3.x/p3.y unset; fall back to a parent's coordinates so
    // the returned point still carries a meaningful location for diagnostics.
    MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
        p3, 0.5*(p1.K_minus - p2.K_plus), mach_guess);
    if (thermo_error != MocErrorCode::NONE) {
        p3.x = p1.x;
        p3.y = p1.y;
        return {p3, thermo_error};
    }
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
        return {p3, MocErrorCode::NON_DOWNSTREAM_POINT};
    }

    MocErrorCode validity = check_point_validity(p3, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        log_warning("Invalid interior point at ({}, {}): error code {}.",
            p3.x, p3.y, to_string(validity));
        return {p3, validity};
    }

    return {p3, MocErrorCode::NONE};
}
PointResult MocNozzle::solve_interior_point_axisymmetric(
    const CharacteristicPoint& p1,
    const CharacteristicPoint& p2)
{
    CharacteristicPoint p3{};

    // Axisymmetric interior point in the (theta, nu) formulation, consistent with the planar,
    // wall, and axis solvers. The Riemann invariants K+ = theta - nu (constant along C+) and
    // K- = theta + nu (constant along C-) pick up source terms in axisymmetric flow:
    //   along C+:  d(theta - nu) = -L dx,   L = sin(mu) sin(theta) / (y cos(theta + mu))
    //   along C-:  d(theta + nu) = +M dx,   M = sin(mu) sin(theta) / (y cos(theta - mu))
    // With L = M = 0 this reduces exactly to the planar solver. Working in nu (rather than the
    // velocity form) keeps it exact for perfect gas, where nu = nu(M); the velocity form would
    // require the *true* velocity for cot(mu) dV/V = dnu to hold, but V is stored as M here.

    // Predictor: straight characteristics, source terms evaluated at the parents.
    double c_minus_angle = p1.theta - p1.mu;
    double c_plus_angle = p2.theta + p2.mu;

    auto [x, y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
    p3.x = x;
    p3.y = y;

    // NOTE: watch out for singularities when a parent is on the centerline.
    // The averaging of y should prevent issues.
    double L = sin(p2.mu)*sin(p2.theta)/(0.5*(p2.y+p3.y) * cos(p2.theta + p2.mu));
    double M = sin(p1.mu)*sin(p1.theta)/(0.5*(p1.y+p3.y) * cos(p1.theta - p1.mu));

    double K_plus = (p2.theta - p2.nu) - L*(p3.x - p2.x);
    double K_minus = (p1.theta + p1.nu) + M*(p3.x - p1.x);
    p3.theta = 0.5*(K_minus + K_plus);
    MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
        p3, 0.5*(K_minus - K_plus), 0.5*(p1.mach + p2.mach));
    if (thermo_error != MocErrorCode::NONE) {
        return {p3, thermo_error};
    }

    // Corrector: average the characteristic angles and source terms across the step.
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
    double M_avg = sin(mu1_avg)*sin(theta1_avg)/(0.5*(p1.y+p3.y) * cos(c_minus_angle_corrected));
    double L_avg = sin(mu2_avg)*sin(theta2_avg)/(0.5*(p2.y+p3.y) * cos(c_plus_angle_corrected));

    K_plus = (p2.theta - p2.nu) - L_avg*(p3.x - p2.x);
    K_minus = (p1.theta + p1.nu) + M_avg*(p3.x - p1.x);
    p3.theta = 0.5*(K_minus + K_plus);
    thermo_error = update_thermodynamic_state_from_nu(p3, 0.5*(K_minus - K_plus), p3.mach);
    if (thermo_error != MocErrorCode::NONE) {
        return {p3, thermo_error};
    }
    p3.K_minus = p3.theta + p3.nu;
    p3.K_plus = p3.theta - p3.nu;

    // validity checks
    if (p3.x < p1.x || p3.x < p2.x) {
        log_warning("Non-downstream intersection at ({}, {}). Parents at x=({}, {}).",
            p3.x, p3.y, p1.x, p2.x);
        return {p3, MocErrorCode::NON_DOWNSTREAM_POINT};
    }

    MocErrorCode validity = check_point_validity(p3, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        log_warning("Invalid interior point at ({}, {}): error code {}.",
            p3.x, p3.y, to_string(validity));
        return {p3, validity};
    }

    return {p3, MocErrorCode::NONE};
}

PointResult MocNozzle::solve_interior_point_iterative(
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
    CharacteristicPoint p3{};
    double c_minus_angle = p1.theta - p1.mu;
    double c_plus_angle = p2.theta + p2.mu;

    auto [x, y] = characteristic_intersection_with_angle(p1, p2, c_minus_angle, c_plus_angle);
    p3.x = x;
    p3.y = y;
    double S1 = cminus_source_term(p1, y);
    double S2 = cplus_source_term(p2, y);

    double mach = find_node_mach(p1, p2, S1 - S2, 0.5 * (p1.mach + p2.mach));
    if (mach < 1.0) {
        return {p3, MocErrorCode::PM_INVERSION_FAILED};
    }
    MocErrorCode thermo_error = update_thermodynamic_state_from_mach(p3, mach);
    if (thermo_error != MocErrorCode::NONE) {
        return {p3, thermo_error};
    }

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
        double mach_new = find_node_mach(p1, p2, S1_new - S2_new, p3.mach);
        if (mach_new < 1.0) {
            return {p3, MocErrorCode::PM_INVERSION_FAILED};
        }
        thermo_error = update_thermodynamic_state_from_mach(p3, mach_new);
        if (thermo_error != MocErrorCode::NONE) {
            return {p3, thermo_error};
        }

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

    MocErrorCode validity = check_point_validity(p3, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        return {p3, validity};
    }
    return {p3, MocErrorCode::NONE};
}

PointResult MocNozzle::solve_wall_flow(
    const CharacteristicPoint& interior_parent,
    double theta_wall)
{
    CharacteristicPoint wall_point{};
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
    MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
        wall_point,
        wall_point.theta - wall_point.K_plus,
        interior_parent.mach);
    if (thermo_error != MocErrorCode::NONE) {
        return {wall_point, thermo_error};
    }
    wall_point.K_minus = wall_point.theta + wall_point.nu;

    return {wall_point, MocErrorCode::NONE};
}

PointResult MocNozzle::solve_wall_point_design(
    const CharacteristicPoint& interior_parent,
    const CharacteristicPoint& previous_wall_point,
    double theta_wall)
{
    PointResult flow_result = solve_wall_flow(interior_parent, theta_wall);
    if (flow_result.error != MocErrorCode::NONE) {
        flow_result.point.x = interior_parent.x;
        flow_result.point.y = interior_parent.y;
        return flow_result;
    }
    CharacteristicPoint wall_point = flow_result.point;

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

                // C+ compatibility nu_P = theta_P - theta_B + nu_B + S, matching the interior
                // solver's K+ = theta - nu decreasing by the source along C+ (d(theta-nu) = -S).
                // The source enters with a +S sign here, NOT -S.
                MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
                    wall_point,
                    wall_point.theta-interior_parent.theta + S + interior_parent.nu,
                    interior_parent.mach);
                if (thermo_error != MocErrorCode::NONE) {
                    return {wall_point, thermo_error};
                }

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

    MocErrorCode validity = check_point_validity(wall_point, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        log_warning("Invalid wall point at ({}, {}): error code {}.",
            wall_point.x, wall_point.y, to_string(validity));
        return {wall_point, validity};
    }

    return {wall_point, MocErrorCode::NONE};
}

std::optional<PointResult> MocNozzle::solve_wall_point_analysis(
        const CharacteristicPoint& interior_parent)
{
    const NozzleProfile& wall = m_options.nozzle_profile;
    auto [x_wall, y_wall] = intersect_characteristic_with_wall(
                interior_parent, wall);

    // no wall hit found
    if (x_wall < 0.0 && y_wall < 0.0) {
        return std::nullopt;
    }

    // theta_at() throws (std::runtime_error, "Queried x ... larger than nozzle
    // profile") when the predicted C+ hits beyond the last profile station -- a
    // known consequence of the analysis-outflow gap (characteristics leaving the
    // contour are not yet terminated as OUTFLOW before the wall solve runs; see
    // instructions/moc_algorithm.md Sec. 9.1, a Phase 2 fix). Converting that here
    // instead of letting it escape is exactly what keeps solve() from throwing.
    double wall_theta;
    try {
        wall_theta = wall.theta_at(x_wall);
    } catch (const std::exception& e) {
        log_warning("Wall profile query out of bounds at x={}: {}", x_wall, e.what());
        CharacteristicPoint failed_point{};
        failed_point.x = interior_parent.x;
        failed_point.y = interior_parent.y;
        return PointResult{failed_point, MocErrorCode::WALL_QUERY_OUT_OF_BOUNDS};
    }

    double char_angle = interior_parent.theta + interior_parent.mu;
    PointResult flow_result = solve_wall_flow(interior_parent, wall_theta);
    if (flow_result.error != MocErrorCode::NONE) {
        flow_result.point.x = interior_parent.x;
        flow_result.point.y = interior_parent.y;
        return flow_result;
    }
    CharacteristicPoint wall_point = flow_result.point;

    wall_point.x = x_wall;
    wall_point.y = y_wall;
    // Corrector (for curved characteristics):
    // Solve the wall point flow using computed theta
    // then recompute characteristic slope as average of parent and wall-point slopes.
    // One correction is usually sufficient.
    switch (m_options.flow_type) {
        case MocFlowKind::PLANAR: {

            //for planar flow, solve_wall_flow equations are exact.
            // Corrector: re-intersect the C+ ray from the interior parent using the average
            // of the parent and predicted wall-point characteristic angles.
            double wall_char_angle = wall_point.theta + wall_point.mu;
            double average_char_angle = 0.5*(wall_char_angle + char_angle);

            auto [x_corrected, y_corrected] = find_wall_hit(interior_parent, wall, average_char_angle);
            try {
                wall_theta = wall.theta_at(x_corrected);
            } catch (const std::exception& e) {
                log_warning("Wall profile query out of bounds at x={}: {}", x_corrected, e.what());
                wall_point.x = interior_parent.x;
                wall_point.y = interior_parent.y;
                return PointResult{wall_point, MocErrorCode::WALL_QUERY_OUT_OF_BOUNDS};
            }
            // Recompute the flow at the corrected wall angle, then restore the position:
            // solve_wall_flow returns flow properties only and leaves x,y at their defaults.
            flow_result = solve_wall_flow(interior_parent, wall_theta);
            if (flow_result.error != MocErrorCode::NONE) {
                flow_result.point.x = interior_parent.x;
                flow_result.point.y = interior_parent.y;
                return flow_result;
            }
            wall_point = flow_result.point;
            wall_point.x = x_corrected;
            wall_point.y = y_corrected;
            break;
        }
        case MocFlowKind::AXISYMMETRIC: {
            // For axisymmetric flow, we must account for the source term.
            double old_S = 0.0;
            const int max_iter = 4;
            for (int i = 0; i < max_iter; i++){

                double c_plus_angle = average_cplus_angle(interior_parent, wall_point);
                auto [x_corrected, y_corrected] = find_wall_hit(interior_parent, wall, c_plus_angle);
                // The corrected C+ angle can miss the contour even when the predictor hit
                // it (e.g. near the exit lip). Terminate the characteristic as outflow rather
                // than marching on with a (-1, -1) miss sentinel.
                if (x_corrected < 0.0 && y_corrected < 0.0) {
                    return std::nullopt;
                }
                wall_point.x = x_corrected;
                wall_point.y = y_corrected;
                // In analysis mode the wall point's flow angle is set by the wall, not
                // prescribed as in design mode. Update it to the wall angle at the corrected
                // position before applying the C+ compatibility relation; otherwise the
                // reflected wave is computed from a stale predictor angle.
                try {
                    wall_theta = wall.theta_at(wall_point.x);
                } catch (const std::exception& e) {
                    log_warning("Wall profile query out of bounds at x={}: {}", wall_point.x, e.what());
                    return PointResult{wall_point, MocErrorCode::WALL_QUERY_OUT_OF_BOUNDS};
                }
                wall_point.theta = wall_theta;

                double S = cplus_source_term(interior_parent, wall_point);
                double residual = std::abs(S - old_S);
                if (residual < m_options.solver_options.abstol) break;
                old_S = S;

                // C+ compatibility nu_P = theta_P - theta_B + nu_B + S (see solve_wall_point_design).
                MocErrorCode thermo_error = update_thermodynamic_state_from_nu(
                    wall_point,
                    wall_point.theta-interior_parent.theta + S + interior_parent.nu,
                    interior_parent.mach);
                if (thermo_error != MocErrorCode::NONE) {
                    return PointResult{wall_point, thermo_error};
                }

                if (i == max_iter - 1 && residual >= m_options.solver_options.abstol) {
                    log_warning("Wall analysis source term iteration did not converge."
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

    MocErrorCode validity = check_point_validity(wall_point, m_options.solver_options.abstol);
    if (validity != MocErrorCode::NONE) {
        log_warning("Invalid analysis wall point at ({}, {}): error code {}.",
            wall_point.x, wall_point.y, to_string(validity));
        return PointResult{wall_point, validity};
    }

    return PointResult{wall_point, MocErrorCode::NONE};
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

MocErrorCode MocNozzle::update_thermodynamic_state_from_nu(
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
            // mach_from_prandtl_meyer's documented failure sentinel: the Newton
            // solve for the inverse Prandtl-Meyer function did not converge. Report
            // it as data instead of letting -1.0 flow into mach_to_mu() and
            // corrector-step Mach averages downstream.
            if (point.mach < 0.0) {
                log_warning("Prandtl-Meyer inversion failed to converge for nu = {}.", nu);
                return MocErrorCode::PM_INVERSION_FAILED;
            }
            point.V = point.mach;
            break;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::EQUILIBRIUM: {
            try {
                auto [idx, weight] = pm_table.find_nu_index_and_weight(nu);
                point.V = pm_table.interpolate_at_index(idx, weight, pm_table.velocities);
                point.gamma_s = pm_table.interpolate_at_index(idx, weight, pm_table.gamma_s);
                point.mach = pm_table.interpolate_at_index(idx, weight, pm_table.machs);
                point.cantera_state = pm_table.interpolate_state_at_index(idx, weight);
            } catch (const std::out_of_range& e) {
                log_warning("Prandtl-Meyer table lookup by nu = {} out of range: {}", nu, e.what());
                return MocErrorCode::TABLE_RANGE_EXCEEDED;
            }
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    update_thermodynamic_state(point);
    point.mu = mach_to_mu(point.mach);
    return MocErrorCode::NONE;
}

MocErrorCode MocNozzle::update_thermodynamic_state_from_mach(CharacteristicPoint& point, double mach) {
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
            try {
                auto [idx, weight] = pm_table.find_mach_index_and_weight(mach);
                point.V = pm_table.interpolate_at_index(idx, weight, pm_table.velocities);
                point.gamma_s = pm_table.interpolate_at_index(idx, weight, pm_table.gamma_s);
                point.nu = pm_table.interpolate_at_index(idx, weight, pm_table.nus);
                point.cantera_state = pm_table.interpolate_state_at_index(idx, weight);
            } catch (const std::out_of_range& e) {
                log_warning("Prandtl-Meyer table lookup by mach = {} out of range: {}", mach, e.what());
                return MocErrorCode::TABLE_RANGE_EXCEEDED;
            }
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    update_thermodynamic_state(point);
    point.mu = mach_to_mu(point.mach);
    return MocErrorCode::NONE;
}

MocErrorCode MocNozzle::update_thermodynamic_state_from_V(CharacteristicPoint& point, double V) {
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
            try {
                auto [idx, weight] = pm_table.find_V_index_and_weight(V);
                point.mach = pm_table.interpolate_at_index(idx, weight, pm_table.machs);
                point.gamma_s = pm_table.interpolate_at_index(idx, weight, pm_table.gamma_s);
                point.nu = pm_table.interpolate_at_index(idx, weight, pm_table.nus);
                point.cantera_state = pm_table.interpolate_state_at_index(idx, weight);
            } catch (const std::out_of_range& e) {
                log_warning("Prandtl-Meyer table lookup by V = {} out of range: {}", V, e.what());
                return MocErrorCode::TABLE_RANGE_EXCEEDED;
            }
            break;
        }
        default: //unreachable
            throw std::runtime_error("Invalid value of GasChemistry specified.");
    }
    update_thermodynamic_state(point);
    point.mu = mach_to_mu(point.mach);
    return MocErrorCode::NONE;
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

    // A wall hit needs at least one segment; a degenerate profile would underflow the
    // loop bound and the extrapolation indices below.
    if (wall.x.size() < 2) {
        throw std::runtime_error("find_wall_hit: nozzle profile has fewer than two points.");
    }

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
        // trun negative values
        return {-1.0, -1.0};
        // size_t last = wall.x.size() - 1;
        // double wall_slope = (wall.y[last] - wall.y[last-1])/(wall.x[last] - wall.x[last - 1]);
        // double wall_intercept = wall.y[last] - wall_slope * wall.x[last];
        // double char_intercept = p.y - c_plus_slope * p.x;
        // x = (wall_intercept - char_intercept) / (c_plus_slope - wall_slope);
        // y = p.y + c_plus_slope * (x - p.x);
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
        if (meta.active && meta.family == family) {
            view.y_values.push_back(net.points[meta.latest_point_idx].y);
            view.chain_indices.push_back(i);
            view.leading_pt_indices.push_back(meta.latest_point_idx);
        }
    }
}

void MocNozzle::sort_plus_edges_by_proximity(LeadingEdgeView& view) const {
    std::vector<size_t> order(view.chain_indices.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return view.y_values[a] > view.y_values[b];
    });

    std::vector<double> y_sorted(order.size());
    std::vector<size_t> chain_sorted(order.size());
    std::vector<size_t> pt_sorted(order.size());
    for (size_t k = 0; k < order.size(); k++) {
        y_sorted[k] = view.y_values[order[k]];
        chain_sorted[k] = view.chain_indices[order[k]];
        pt_sorted[k] = view.leading_pt_indices[order[k]];
    }
    view.y_values = std::move(y_sorted);
    view.chain_indices = std::move(chain_sorted);
    view.leading_pt_indices = std::move(pt_sorted);
}


} // namespace Goddard