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

bool MocNozzle::is_solved() const {
    return m_is_solved;
}

// -- MocNozzle --
MocResult MocNozzle::solve() {
    // Re-checked here, not just in the constructor: m_options is public and callers do
    // adjust it between construction and solve().
    validate_moc_options(m_options);

    m_log = MocLog(m_options.log_level);
    m_is_solved = false;
    m_pass_diagnostics.clear();
    m_reference_spacing = 0.0;
    m_throat_radius = 1.0;

    CharacteristicNet net;
    StartLine line;
    m_options.nozzle_profile = setup_nozzle_profile(m_options.geometry);

    ThroatCondition throat{};

    if (m_options.chemistry == GasChemistry::PERFECT_GAS) {
        // Pure algebraic path: no Cantera dependency
        // throat populated with dummy values; they are not used in the perfect gas case.
        throat = ThroatCondition{
            .converged = true,
            .speed_of_sound = 1.0,
            .H_stagnation = 1.0,
            .P_inlet = 1.0, // dimensionless stagnation pressure
            .S_inlet = 1.0,
            .gamma_s = m_options.gamma,
            .dlV_dlP_T = -1.0,
            .dlV_dlT_P = 1.0,
            .state = {}
        };
        m_thermo = MocThermo::perfect_gas(m_options.gamma);
    }
    else {
        // Cantera-backed path for frozen/equilibrium chemistry
        NozzleOptions nozzle_opts;
        nozzle_opts.chemistry = m_options.chemistry;
        Nozzle nozzle(*m_gas, nozzle_opts);
        throat = nozzle.solve_throat_conditions();
        m_thermo = MocThermo::tabulated(*m_gas, m_options.chemistry, throat);
    }

    // Everything a kernel or unit process needs for this solve, besides the points it works
    // on: m_options.nozzle_profile is still where the resolved contour lives at this step
    // (see setup_nozzle_profile above); a later step moves it out of MocOptions.
    m_context.emplace(MocSolveContext{m_options, m_options.nozzle_profile, *m_thermo, m_log});

    MocResult result;

    // The initial data line (a Newton solve for the Kliegel-Levine transonic line, a
    // centered-fan Prandtl-Meyer inversion, or a PrandtlMeyerTable lookup) can fail
    // numerically before the kernel ever starts marching. Converting those failures
    // here, rather than letting them propagate as exceptions, is what makes solve()
    // never throw for a numerical/convergence failure.
    try {
        if (m_inverse_front_override.has_value()) {
            // Test-only path (see m_inverse_front_override): skip the normal throat/KL/fan
            // construction entirely -- it need not be consistent with a hand-built front --
            // and seed the inverse kernel directly from the override. Treated as a generic
            // (Kliegel-Levine-like) line for build_inverse_initial_front and
            // measure_start_line: it is not collinear along a single characteristic.
            line.points = *m_inverse_front_override;
            line.family.reset();
            line.used = MocStartLine::KLIEGEL_LEVINE;
        } else {
            line = build_start_line(*m_context, throat);
        }

        if (m_options.mode == MocMode::DESIGN_MIN_LENGTH) {
            // The start line declares whether it lies along a single characteristic
            // (line.family). A centered expansion fan does, and needs a wall anchor at the
            // throat lip (0, 1): it seeds leading_wall_point() for the first wall solve and
            // gives the area ratio its throat reference. A transonic start line spans
            // axis-to-wall and already carries its own wall point, so no separate anchor is seeded.
            if (line.family.has_value()) {
                CharacteristicPoint throat_lip{};
                throat_lip.x = 0.0;
                throat_lip.y = 1.0;
                // The lip angle anchors the first wall solve in DESIGN_MIN_LENGTH, where it equals
                // theta_max. In analysis the wall angles come from the profile and the lip is purely
                // an anchor, so the topmost fan angle is a harmless stand-in.
                throat_lip.theta = m_options.mode == MocMode::DESIGN_MIN_LENGTH
                    ? m_options.theta_max
                    : (line.points.empty() ? 0.0 : line.points.back().theta);
                net.seed_wall_point(throat_lip);
            }
            net.add_initial_data_line(line.points, line.family);

            // Anchor mesh control to the geometry, once, from the line actually seeded. Both
            // init paths put a wall point at the throat lip first, so wall_y.front() is the
            // throat radius in whatever units the net is carrying; taking it from the net rather
            // than from geometry.throat_radius keeps the two consistent even when the initial
            // line is built in normalized coordinates.
            if (!net.wall_y.empty() && net.wall_y.front() > 0.0 &&
                m_options.num_characteristics > 1) {
                m_throat_radius = net.wall_y.front();
                m_reference_spacing =
                    m_throat_radius / static_cast<double>(m_options.num_characteristics - 1);
            }
        }
        else {
            // INVERSE: no throat-lip anchor is seeded (see build_inverse_initial_front /
            // CharacteristicNet::add_front); wall_x.front() is F_0's own wall point, not
            // (0, 1), so the throat radius for area_ratio/mesh purposes comes from the
            // geometry directly.
            std::vector<CharacteristicPoint> front0 = build_inverse_initial_front(line);
            if (front0.size() < 2) {
                throw ConvergenceError(
                    "Inverse march: initial front has fewer than two points.");
            }
            net.add_front(front0);
            m_throat_radius = (m_options.geometry.throat_radius > 0.0)
                ? m_options.geometry.throat_radius : 1.0;
        }
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
        result.messages = m_log.messages;
        m_is_solved = true;
        return result;
    }

    // Measured before the march so it is available even when the kernel fails partway --
    // an initialization defect is exactly the case where the march does not finish.
    result.init_diagnostics = measure_start_line(line, *m_context, m_reference_spacing, m_throat_radius);

    std::optional<MocFailure> kernel_failure = (m_options.mode == MocMode::DESIGN_MIN_LENGTH)
        ? solve_characteristic_kernel(net, line)
        : solve_inverse_characteristic_kernel(net);

    // Minimum flow angle and its location, for every scheme: see MocResult::min_theta.
    if (!net.points.empty()) {
        const CharacteristicPoint* lowest = &net.points.front();
        for (const CharacteristicPoint& pt : net.points) {
            if (pt.theta < lowest->theta) lowest = &pt;
        }
        result.min_theta = lowest->theta;
        result.min_theta_x = lowest->x;
        result.min_theta_y = lowest->y;
        if (result.min_theta < -0.5 * M_PI / 180.0) {
            m_log.warning("Flow angle reaches {:.2f} deg at (x={:.3f}, y={:.3f}): a compression is "
                "converging on the axis (a forming shock); the isentropic solution downstream "
                "of it is approximate.", result.min_theta * 180.0 / M_PI,
                result.min_theta_x, result.min_theta_y);
        }
    }

    result.net = net;
    result.messages = m_log.messages;
    if (kernel_failure.has_value()) {
        result.failure = *kernel_failure;
        result.converged = false;
    } else {
        result.converged = true;
    }

    // Populate performance fields
    result.crossings = find_like_characteristic_crossings(net);
    result.pass_diagnostics = m_pass_diagnostics;

    // How much of the requested exit radius the outflow boundary actually reached. The
    // boundary is a ragged staircase of independently terminated chains, so even a healthy
    // march ends up to about one characteristic spacing short; the allowance below is two
    // spacings. This is reported, not folded into `converged`: measured coverage does not
    // separate a coarse-but-healthy solve from a truncated one (planar N=8 and the
    // axisymmetric AR=4 N=8 truncation both sit at 0.92), so the meaningful test is whether
    // the shortfall shrinks under refinement, which only a grid sweep can see.
    if (m_options.mode == MocMode::DESIGN_MIN_LENGTH) {
        result.exit_coverage = 1.0;
        result.reached_exit_plane = true;
    }
    else {
        // The inverse march always lands its last front exactly on the exit plane (see
        // solve_inverse_characteristic_kernel's MocStepLimiter::EXIT step) -- there is no
        // ragged outflow staircase to fall short of. A march that failed partway did not
        // reach it; report how far it got instead.
        if (!kernel_failure.has_value()) {
            result.exit_coverage = 1.0;
            result.reached_exit_plane = true;
        }
        else if (!m_options.nozzle_profile.y.empty() && m_options.nozzle_profile.y.back() > 0.0
                 && !net.wall_y.empty()) {
            result.exit_coverage = net.wall_y.back() / m_options.nozzle_profile.y.back();
            result.reached_exit_plane = false;
        }
    }
    if (!net.wall_x.empty()) {
        result.nozzle_length = net.wall_x.back();
    }
    if (m_options.mode != MocMode::DESIGN_MIN_LENGTH) {
        // net.wall_x.front()/wall_y.front() are F_0's own wall point, not the throat lip
        // (0, 1) DIRECT seeds -- see the INVERSE seeding branch above -- so the throat
        // radius reference comes from the geometry directly instead.
        const double r_throat = (m_options.geometry.throat_radius > 0.0)
            ? m_options.geometry.throat_radius : 1.0;
        if (!net.wall_y.empty()) {
            if (m_options.flow_type == MocFlowKind::PLANAR) {
                result.area_ratio = net.wall_y.back() / r_throat;
            } else {
                const double y_ratio = net.wall_y.back() / r_throat;
                result.area_ratio = y_ratio * y_ratio;
            }
        }
    }
    else if (!net.wall_y.empty() && net.wall_y.front() > 0.0) {
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
        if (m_options.mode != MocMode::DESIGN_MIN_LENGTH) {
            // The exit plane is exactly the last front the kernel built (or, on a partial
            // failure, whatever front it last completed) -- axis to wall, as stored.
            if (!net.fronts.empty()) {
                for (size_t idx : net.fronts.back()) {
                    const CharacteristicPoint& pt = net.points[idx];
                    result.exit_plane.y.push_back(pt.y);
                    result.exit_plane.mach.push_back(pt.mach);
                    result.exit_plane.theta.push_back(pt.theta);
                    result.exit_plane.pressure.push_back(pt.pressure);
                    result.exit_plane.temperature.push_back(pt.temperature);
                    result.exit_plane.gamma_s.push_back(pt.gamma_s);
                    result.exit_plane.velocity.push_back(pt.V);
                }
            }
            m_is_solved = true;
            return result;
        }

        // Add the last characteristic for min length nozzle. Guarded by
        // emptiness: a kernel that now fails fast can abort before any axis/wall
        // point has been recorded, where leading_axis_point()/leading_wall_point()
        // would otherwise index an empty vector.
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


std::optional<MocFailure> MocNozzle::solve_characteristic_kernel(CharacteristicNet& net, const StartLine& line) {
    using Family = CharacteristicFamily;
    LeadingEdgeView plus_edges = leading_edges(net, Family::PLUS);
    sort_plus_edges_by_proximity(plus_edges, net);
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
        m_log.debug("--- kernel pass {}: {} active C+, {} active C- ---",
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
            PairSearch match = find_pair_partner(
                plus_edges.y_values[i], minus_edges, cminus_is_intersected, net);

            // intersect C+ with closest C- above it.
            if (match.chain_idx.has_value()) [[likely]] {
                const CharacteristicPoint& minus_pt = net.points[match.pt_idx];
                const CharacteristicPoint& plus_pt = net.points[plus_edges.leading_pt_indices[i]];

                PointResult result = solve_interior_point(*m_context, minus_pt, plus_pt);
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
                m_log.debug("PAIR minus(x={:.6f},y={:.6f},th={:.6f},mu={:.6f}) "
                    "plus(x={:.6f},y={:.6f},th={:.6f},mu={:.6f}) -> (x={:.6f},y={:.6f},mach={:.6f})",
                    minus_pt.x, minus_pt.y, minus_pt.theta, minus_pt.mu,
                    plus_pt.x, plus_pt.y, plus_pt.theta, plus_pt.mu,
                    result.point.x, result.point.y, result.point.mach);
                intersections.push_back({
                    result.point,
                    PointMembership {
                        .c_plus_chain_idx = plus_edges.chain_indices[i],
                        .c_minus_chain_idx = *match.chain_idx
                    }
                });
                paired_cminus.push_back(*match.chain_idx);
                cminus_is_intersected[match.edgevec_idx] = true;
            }
            else if (match.any_cminus_above) {
                // A C- does exist above this C+, but a closer competitor already claimed it
                // this pass (e.g. many individual C+ chains from a Kliegel-Levine transonic
                // line, competing for a single C- freshly born from a wall reflection). Leave
                // this chain active and retry once that C- has advanced on the next pass,
                // rather than wrongly treating it as a wall hit.
                const CharacteristicPoint& plus_pt = net.points[plus_edges.leading_pt_indices[i]];
                m_log.debug("SKIP plus(x={:.6f},y={:.6f}) waiting for scarce C- partner "
                    "(already claimed this pass)", plus_pt.x, plus_pt.y);
            }
            else [[unlikely]] { // C+ intersects with wall

                const CharacteristicPoint& plus_pt = net.points[plus_edges.leading_pt_indices[i]];
                PointResult result = solve_wall_point(
                    plus_pt,
                    net.leading_wall_point(),
                    static_cast<int>(net.wall_point_indices.size()) - 1,
                    line
                );

                if (result.error != MocErrorCode::NONE) {
                    failure = MocFailure{
                        result.error,
                        std::format("Wall point solve failed ({}) for plus(x={:.6f},y={:.6f}).",
                            to_string(result.error), plus_pt.x, plus_pt.y),
                        result.point.x, result.point.y,
                        iters
                    };
                    break;
                }
                m_log.debug("WALL plus(x={:.6f},y={:.6f}) -> wall hit at (x={:.6f},y={:.6f})",
                    plus_pt.x, plus_pt.y, result.point.x, result.point.y);
                intersections.push_back({
                    result.point,
                    PointMembership {
                        .c_plus_chain_idx = plus_edges.chain_indices[i],
                        .c_minus_chain_idx = std::nullopt
                    }
                });
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
            PointResult axis_result = solve_axis_point(*m_context, minus_pt);
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
            m_log.debug("AXIS minus(x={:.6f},y={:.6f}) reflects off axis", minus_pt.x, minus_pt.y);
            intersections.push_back({
                axis_result.point,
                PointMembership {
                    .c_plus_chain_idx = std::nullopt,
                    .c_minus_chain_idx = min_y_cminus_idx
                }
            });
        }
        // insert all intersections into net. Minimum-length design only: the wall
        // absorbs each C+ (no reflected wave), which is exactly what produces a
        // minimum-length contour.
        for (auto&& [pt, mem] : intersections) {
            if (!mem.c_minus_chain_idx.has_value()) {
                // C+ reaches the wall: absorb it (no reflected wave). terminate_c_plus_at_wall
                // pushes wall_x/wall_y itself.
                net.terminate_c_plus_at_wall(*mem.c_plus_chain_idx, pt);
            }
            else if (!mem.c_plus_chain_idx.has_value()) {
                // C- reaches the axis: reflect it into a new upward-marching C+.
                net.reflect_c_minus_off_axis(*mem.c_minus_chain_idx, pt);
            }
            else {
                net.add_point(pt, mem);
            }
        }
        update_leading_edges(plus_edges, net, Family::PLUS);
        sort_plus_edges_by_proximity(plus_edges, net);
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

PointResult MocNozzle::solve_wall_point(
    const CharacteristicPoint& interior_parent,
    const CharacteristicPoint& previous_wall_point,
    int wall_point_index,
    const StartLine& line)
{
    // Minimum length design only: theta is determined by theta schedule.
    // Guard against an empty schedule (the Cantera min-length path does not
    // populate line.theta_schedule) and against indices beyond it (the kernel
    // currently produces more wall hits than scheduled characteristics).
    double theta_wall = m_options.theta_max;
    if (!line.theta_schedule.empty()) {
        size_t k = std::min(static_cast<size_t>(wall_point_index),
                            line.theta_schedule.size() - 1);
        theta_wall = m_options.theta_max - line.theta_schedule[k];
    }

    return solve_wall_point_design(*m_context, interior_parent, previous_wall_point, theta_wall);
}

MocNozzle::LeadingEdgeView MocNozzle::leading_edges(const CharacteristicNet& net, CharacteristicFamily fam) const {
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

void MocNozzle::update_leading_edges(LeadingEdgeView& view, const CharacteristicNet& net, CharacteristicFamily family) const {
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

void MocNozzle::sort_plus_edges_by_proximity(LeadingEdgeView& view, const CharacteristicNet& net) const {
    std::vector<size_t> order(view.chain_indices.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        if (view.y_values[a] != view.y_values[b]) return view.y_values[a] > view.y_values[b];
        return net.points[view.leading_pt_indices[a]].x < net.points[view.leading_pt_indices[b]].x;
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

MocNozzle::PairSearch MocNozzle::find_pair_partner(
    double plus_y,
    const LeadingEdgeView& minus_edges,
    const std::vector<bool>& claimed,
    const CharacteristicNet& net) const
{
    PairSearch search;
    double best_dy = std::numeric_limits<double>::max();

    // Nearest C- above regardless of claim status -- tracked separately from the best
    // *available* partner so a C+ whose true nearest partner was already claimed by a
    // closer competitor is told to wait rather than settling for a farther one.
    double nearest_dy = std::numeric_limits<double>::max();
    bool nearest_is_claimed = false;

    for (size_t j = 0; j < minus_edges.chain_indices.size(); j++) {
        double dy = minus_edges.y_values[j] - plus_y;
        if (dy <= 0) continue;
        search.any_cminus_above = true;

        if (dy < nearest_dy) {
            nearest_dy = dy;
            nearest_is_claimed = claimed[j];
        }
        if (claimed[j]) continue; //skip if already paired

        if (dy < best_dy) {
            best_dy = dy;
            search.chain_idx = minus_edges.chain_indices[j];
            search.pt_idx = minus_edges.leading_pt_indices[j];
            search.edgevec_idx = j;
        }
        // edge case: multiple points at literally the same y due to an expansion fan
        // the tiebreaker is determined by the characteristic angle theta-mu
        if (search.chain_idx.has_value() && best_dy == dy) [[unlikely]] {
            const CharacteristicPoint& old_pt = net.points[search.pt_idx];
            const CharacteristicPoint& new_pt = net.points[minus_edges.leading_pt_indices[j]];
            if (new_pt.theta - new_pt.mu < old_pt.theta - old_pt.mu) {
                search.chain_idx = minus_edges.chain_indices[j];
                search.pt_idx = minus_edges.leading_pt_indices[j];
                search.edgevec_idx = j;
            }
        }
    }

    if (nearest_is_claimed) search.chain_idx = std::nullopt;
    search.dy = best_dy;
    return search;
}

} // namespace Goddard
