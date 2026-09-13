// Chain-pairing kernel for minimum-length design (see moc_direct_march.hpp).
#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <numeric>
#include "goddard/moc_direct_march.hpp"

namespace Goddard {

void DirectMarch::seed() {
    // The start line declares whether it lies along a single characteristic
    // (m_start.family). A centered expansion fan does, and needs a wall anchor at the
    // throat lip (0, 1): it seeds leading_wall_point() for the first wall solve and gives
    // the area ratio its throat reference. A transonic start line spans axis-to-wall and
    // already carries its own wall point, so no separate anchor is seeded.
    if (m_start.family.has_value()) {
        CharacteristicPoint throat_lip{};
        throat_lip.x = 0.0;
        throat_lip.y = 1.0;
        // The lip angle anchors the first wall solve in DESIGN_MIN_LENGTH, where it equals
        // theta_max. In analysis the wall angles come from the profile and the lip is purely
        // an anchor, so the topmost fan angle is a harmless stand-in.
        throat_lip.theta = m_ctx.options.mode == MocMode::DESIGN_MIN_LENGTH
            ? m_ctx.options.theta_max
            : (m_start.points.empty() ? 0.0 : m_start.points.back().theta);
        m_net.seed_wall_point(throat_lip);
    }
    m_net.add_initial_data_line(m_start.points, m_start.family);
}

std::optional<MocFailure> DirectMarch::run() {
    using Family = CharacteristicFamily;
    LeadingEdgeView plus_edges = leading_edges(Family::PLUS);
    sort_plus_edges_by_proximity(plus_edges);
    LeadingEdgeView minus_edges = leading_edges(Family::MINUS);

    std::vector<std::pair<CharacteristicPoint, PointMembership>> intersections;
    std::vector<size_t> paired_cminus;
    std::vector<bool> cminus_is_intersected(m_net.c_chains.size(), false);

    intersections.reserve(m_net.c_chains.size());
    paired_cminus.reserve(m_net.c_chains.size());


    // A correctly converging 2D MoC kernel needs O(N) marching passes. This cap is a
    // safety bound so a non-converging/runaway net terminates promptly instead of
    // spinning indefinitely; reaching it indicates a kernel that has not converged.
    const int maxiter = 2000;
    int iters = 0;

    // Populated the instant any unit process reports a numerical failure; the march
    // aborts immediately rather than continuing to build on top of an invalid point, which
    // would let a corrupted point silently propagate through the rest of the net.
    std::optional<MocFailure> failure;

    while (m_net.has_active_chains() && iters < maxiter && !failure.has_value()) {
        m_ctx.log.debug("--- kernel pass {}: {} active C+, {} active C- ---",
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
                plus_edges.y_values[i], minus_edges, cminus_is_intersected);

            // intersect C+ with closest C- above it.
            if (match.chain_idx.has_value()) [[likely]] {
                const CharacteristicPoint& minus_pt = m_net.points[match.pt_idx];
                const CharacteristicPoint& plus_pt = m_net.points[plus_edges.leading_pt_indices[i]];

                PointResult result = solve_interior_point(m_ctx, minus_pt, plus_pt);
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
                m_ctx.log.debug("PAIR minus(x={:.6f},y={:.6f},th={:.6f},mu={:.6f}) "
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
                const CharacteristicPoint& plus_pt = m_net.points[plus_edges.leading_pt_indices[i]];
                m_ctx.log.debug("SKIP plus(x={:.6f},y={:.6f}) waiting for scarce C- partner "
                    "(already claimed this pass)", plus_pt.x, plus_pt.y);
            }
            else [[unlikely]] { // C+ intersects with wall

                const CharacteristicPoint& plus_pt = m_net.points[plus_edges.leading_pt_indices[i]];
                PointResult result = solve_wall_point(
                    plus_pt,
                    m_net.leading_wall_point(),
                    static_cast<int>(m_net.wall_point_indices.size()) - 1
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
                m_ctx.log.debug("WALL plus(x={:.6f},y={:.6f}) -> wall hit at (x={:.6f},y={:.6f})",
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
            const CharacteristicPoint& minus_pt = m_net.leading_point(*min_y_cminus_idx);
            PointResult axis_result = solve_axis_point(m_ctx, minus_pt);
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
            m_ctx.log.debug("AXIS minus(x={:.6f},y={:.6f}) reflects off axis", minus_pt.x, minus_pt.y);
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
                m_net.terminate_c_plus_at_wall(*mem.c_plus_chain_idx, pt);
            }
            else if (!mem.c_plus_chain_idx.has_value()) {
                // C- reaches the axis: reflect it into a new upward-marching C+.
                m_net.reflect_c_minus_off_axis(*mem.c_minus_chain_idx, pt);
            }
            else {
                m_net.add_point(pt, mem);
            }
        }
        plus_edges = leading_edges(Family::PLUS);
        sort_plus_edges_by_proximity(plus_edges);
        minus_edges = leading_edges(Family::MINUS);

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

PointResult DirectMarch::solve_wall_point(
    const CharacteristicPoint& interior_parent,
    const CharacteristicPoint& previous_wall_point,
    int wall_point_index)
{
    // Minimum length design only: theta is determined by theta schedule.
    // Guard against an empty schedule (the Cantera min-length path does not
    // populate m_start.theta_schedule) and against indices beyond it (the kernel
    // currently produces more wall hits than scheduled characteristics).
    double theta_wall = m_ctx.options.theta_max;
    if (!m_start.theta_schedule.empty()) {
        size_t k = std::min(static_cast<size_t>(wall_point_index),
                            m_start.theta_schedule.size() - 1);
        theta_wall = m_ctx.options.theta_max - m_start.theta_schedule[k];
    }

    return solve_wall_point_design(m_ctx, interior_parent, previous_wall_point, theta_wall);
}

DirectMarch::LeadingEdgeView DirectMarch::leading_edges(CharacteristicFamily family) const {
    LeadingEdgeView view;
    view.family = family;
    for (size_t i = 0; i < m_net.chain_metadata.size(); i++) {
        const ChainMetadata& meta = m_net.chain_metadata[i];
        if (meta.active && meta.family == family) {
            view.y_values.push_back(m_net.points[meta.latest_point_idx].y);
            view.chain_indices.push_back(i);
            view.leading_pt_indices.push_back(meta.latest_point_idx);
        }
    }
    return view;
}

void DirectMarch::sort_plus_edges_by_proximity(LeadingEdgeView& view) const {
    std::vector<size_t> order(view.chain_indices.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        if (view.y_values[a] != view.y_values[b]) return view.y_values[a] > view.y_values[b];
        return m_net.points[view.leading_pt_indices[a]].x < m_net.points[view.leading_pt_indices[b]].x;
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

DirectMarch::PairSearch DirectMarch::find_pair_partner(
    double plus_y,
    const LeadingEdgeView& minus_edges,
    const std::vector<bool>& claimed) const
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
            const CharacteristicPoint& old_pt = m_net.points[search.pt_idx];
            const CharacteristicPoint& new_pt = m_net.points[minus_edges.leading_pt_indices[j]];
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
