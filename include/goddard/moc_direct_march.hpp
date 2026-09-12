#pragma once
#include <optional>
#include <vector>
#include "goddard/characteristic_net.hpp"
#include "goddard/moc.hpp"
#include "goddard/moc_context.hpp"
#include "goddard/moc_initialization.hpp"
#include "goddard/moc_unit_processes.hpp"

namespace Goddard {

/**
 * Chain-pairing kernel for minimum-length design: advances a ladder of C+/C- chain leading
 * edges, pairing the nearest unclaimed opposite-family edge each pass (or absorbing a C+ at
 * the wall, or reflecting a C- off the axis), until every chain has terminated.
 */
class DirectMarch {
public:
    DirectMarch(const MocSolveContext& ctx, const StartLine& start, CharacteristicNet& net)
        : m_ctx(ctx), m_start(start), m_net(net) {}

    /**
     * Seed the net: the throat-lip anchor (0, 1) when the start line is a centered fan (a
     * transonic start line spans axis-to-wall and already carries its own wall point, so no
     * separate anchor is seeded), then the start line itself via
     * CharacteristicNet::add_initial_data_line.
     */
    void seed();

    /**
     * Advance the ladder to completion.
     *
     * @return the failure that aborted the march (a unit-process error, or the iteration
     *         safety cap being reached), or std::nullopt if every chain terminated cleanly.
     */
    std::optional<MocFailure> run();

private:
    /** Struct-of-array of leading edge values for a specified family. */
    struct LeadingEdgeView {
        CharacteristicFamily family;
        std::vector<double> y_values;
        std::vector<size_t> chain_indices;
        std::vector<size_t> leading_pt_indices;
    };

    /** Outcome of searching for a C+ leading edge's C- pairing partner. */
    struct PairSearch {
        bool any_cminus_above = false; ///< A C- exists above this C+ at all. False means the C+ is wall-bound.
        std::optional<size_t> chain_idx = std::nullopt; ///< Partner C- chain; empty means SKIP (or wall-bound, when any_cminus_above is false).
        size_t pt_idx = 0;      ///< Partner's leading point index.
        size_t edgevec_idx = 0; ///< Partner's index within the C- leading-edge view (the claim key).
        double dy = 0.0;        ///< Height from the C+ to the partner.
    };

    /** Struct-of-array of leading edge values for `family`, built fresh from the net's current
     *  chain metadata (used both to seed a view and to refresh one after a pass). */
    LeadingEdgeView leading_edges(CharacteristicFamily family) const;

    /**
     * Reorder a C+ leading-edge view by descending y (closest to the wall first), tie-broken
     * by ascending x.
     *
     * When several individual C+ chains compete for the same partner in a single kernel
     * pass (e.g. many chains seeded from a Kliegel-Levine transonic line, all wanting the
     * one C- freshly born from a wall reflection), the per-pass search must resolve the
     * geometrically closest competitor first. Iterating in chain-creation order instead
     * lets a far-away C+ (e.g. the transonic line's axis point) claim a partner meant for
     * a much closer chain, producing a physically invalid, oversized jump.
     *
     * Multiple C+ chains can also sit at exactly y=0 simultaneously (e.g. several axis
     * reflections in flight at once, which dual-family KL-init seeding makes common): the
     * y-only comparator leaves their relative order unspecified under `std::sort`
     * (unstable). Ties are broken by ascending x -- the chain whose leading point is
     * further upstream reaches its next partner first, so it must claim before a chain
     * that is already further downstream.
     */
    void sort_plus_edges_by_proximity(LeadingEdgeView& view) const;

    /**
     * Find the C- leading edge that a C+ at `plus_y` pairs with this pass: the nearest one
     * above it that no closer competitor has already claimed.
     *
     * When the truly-nearest C- above is already claimed, the result is empty rather than
     * the next-nearest unclaimed one -- settling for a farther partner would violate lattice
     * adjacency. A claim implies someone else progressed this pass, so waiting cannot
     * deadlock.
     */
    PairSearch find_pair_partner(
        double plus_y,
        const LeadingEdgeView& minus_edges,
        const std::vector<bool>& claimed) const;

    /**
     * Minimum-length design only. Check PointResult::error for a numerical failure.
     * Resolves the wall angle from m_start.theta_schedule, then delegates to
     * solve_wall_point_design (moc_unit_processes.hpp).
     */
    PointResult solve_wall_point(
        const CharacteristicPoint& interior_parent,
        const CharacteristicPoint& previous_wall_point,
        int wall_point_index);

    const MocSolveContext& m_ctx;
    const StartLine& m_start;
    CharacteristicNet& m_net;
};

} // namespace Goddard
