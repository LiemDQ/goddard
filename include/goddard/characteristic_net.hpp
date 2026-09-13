#pragma once
#include <vector>
#include <optional>
#include <utility>
#include <ranges>
#include "goddard/characteristics.hpp"

namespace Goddard {

/** How a characteristic chain stopped being marched. */
enum class ChainTermination {
    NOT_TERMINATED, WALL, AXIS
};

/**
 * Bookkeeping for one characteristic chain in a CharacteristicNet.
 *
 * A chain is a single characteristic line, stored as an ordered list of point indices. It is
 * born at the initial data line or at a wall/axis reflection, and stops when it terminates.
 */
struct ChainMetadata {
    /// False once the chain has terminated; only active chains are still marched.
    bool active = true;
    /// How this chain stopped being marched.
    ChainTermination termination = ChainTermination::NOT_TERMINATED;
    /// Which characteristic family this chain belongs to.
    CharacteristicFamily family;
    /// Index into CharacteristicNet::points of the chain's first point.
    size_t origin_point_idx;
    /// Index into CharacteristicNet::points of the chain's leading (most downstream) point.
    size_t latest_point_idx;
};

/**
 * The chains one point belongs to.
 *
 * An interior point lies on exactly one characteristic of each family, so it leads both. A wall,
 * axis, or initial-line point may belong to only one.
 */
struct PointMembership {
    /// Index into CharacteristicNet::c_chains of the C+ chain through this point, if any.
    std::optional<size_t> c_plus_chain_idx;
    /// Index into CharacteristicNet::c_chains of the C- chain through this point, if any.
    std::optional<size_t> c_minus_chain_idx;
};

/**
 * The method-of-characteristics flow field mesh: every solved point, the characteristic
 * chains and marching fronts that connect them, and the wall/axis boundary bookkeeping.
 * Built up incrementally by the two kernels (the chain-pairing ladder for minimum-length
 * design, the reference-plane march for analysis and Rao design) via the mutation methods
 * below; queried afterward via the columnar point data and the boundary/chain accessors.
 */
class CharacteristicNet {

    public:
    /**
     * Every point in the net, in creation order. The index is opaque; the flow-field topology
     * lives in `c_chains`, and the boundaries in `wall_point_indices`/`axis_point_indices`.
     */
    std::vector<CharacteristicPoint> points;

    
    /**
     * Each entry is the chain of point indices along one characteristic.
     * Note that wall and axis points terminate a chain, and also start the next chain.
     */
    std::vector<std::vector<size_t>> c_chains;

    /**
     * Given a point index, which chains does it belong to?
     * `membership[i]` describes `points[i]`.
     */
    std::vector<PointMembership> membership;

    /** Metadata for each chain in `c_chains`, in the same order. */
    std::vector<ChainMetadata> chain_metadata;
    
    /** Indices into `points` of the points lying on the nozzle wall, in march order. */
    std::vector<size_t> wall_point_indices;
    /** Indices into `points` of the points lying on the centerline, in march order. */
    std::vector<size_t> axis_point_indices;

    /**
     * Axial coordinates of the wall points, in length units. The net owns the invariant
     * that `wall_x`/`wall_y` stay parallel to `wall_point_indices` (same length, same
     * order); every method that appends to `wall_point_indices` appends to these too.
     */
    std::vector<double> wall_x;
    /** Radial coordinates of the wall points, in length units. See `wall_x`. */
    std::vector<double> wall_y;

    /**
     * Point indices of every marching front, axis to wall, in order. Populated by the
     * front-based nets (analysis and Rao design): `fronts.front()` is the initial front F_0
     * and `fronts.back()` is the exit plane once the march has completed. Empty for the
     * chain-pairing ladder (minimum-length design), which has no synchronized front -- its
     * topology lives entirely in `c_chains`.
     */
    std::vector<std::vector<size_t>> fronts;

    /** True when the net holds no points at all. */
    bool empty() const;

    /**
     * Get the leading point of the chain at the given index.
     */
    CharacteristicPoint& leading_point(size_t chain_idx);
    const CharacteristicPoint& leading_point(size_t chain_idx) const;

    /**
     * Get the leading wall point.
     */
    CharacteristicPoint& leading_wall_point();
    const CharacteristicPoint& leading_wall_point() const;

    /**
     * Get the leading axis point.
     */
    CharacteristicPoint& leading_axis_point();
    const CharacteristicPoint& leading_axis_point() const;

    /** Add a point and its membership to list of points, and
     * appends it to its member characteristic chains. 
     * @return index of the added point
     */
    size_t add_point(CharacteristicPoint pt, PointMembership m);

    /** Add a point along an initial data line. 
     * @return index of the added point and corresponding characteristics 
     */
    size_t add_initialization_point(CharacteristicPoint pt, bool cminus, bool cplus);
    
    /** Seed the net from an initial data line.
     *
     * Handles both data-line topologies. When `on_characteristic` is empty the points are
     * assumed to lie on distinct characteristics (a non-collinear transonic start line, e.g.
     * Kliegel-Levine): every interior point seeds both a C+ and a C- chain, mirroring the
     * fan-init topology, and the last point (already at the wall) immediately reflects into
     * the first C- chain, exactly as a wall reflection would during marching. The first
     * point (the axis bootstrap) seeds only a C+: giving it a C- would immediately
     * re-reflect it off the axis as a degenerate point. When `on_characteristic` is set the
     * points are collinear along a single characteristic of that family (e.g. a centered
     * expansion fan), and are seeded as one shared chain via add_initial_characteristic.
     *
     * The caller must supply the topology; it cannot be inferred from the points, since a
     * Riemann invariant is constant along a characteristic only for planar flow.
     */
    void add_initial_data_line(
        const std::vector<CharacteristicPoint>& points,
        std::optional<CharacteristicFamily> on_characteristic = std::nullopt);

    /** Add points from an initial data line that happens to be along a characteristic.
     * This is primarily used when the initial data line originates from a centered expansion.
    */
    void add_initial_characteristic(const std::vector<CharacteristicPoint>& points,
                                     CharacteristicFamily family = CharacteristicFamily::PLUS);

    /** Add a marching front's points (axis to wall) to the net: no chain membership (the
     * front-based kernels do not use c_chains), but registers the front's axis and wall
     * points and records it in `fronts`.
     * @return Indices into `points` of the added front, axis to wall.
     */
    std::vector<size_t> add_front(const std::vector<CharacteristicPoint>& front_points);

    /** Create a new characteristic with `pt_idx` as the starting point,
     * and add it to the tracking lists.
     *
     * @return Index of the chain
     */
    size_t create_chain(size_t pt_idx, CharacteristicFamily family);
    
    /** Add a characteristic chain and its metadata to be tracked. 
     * @return Index of the chain
     */
    size_t push_chain(std::vector<size_t>&& chain, ChainMetadata metadata);

    /** Terminate a C+ characteristic off a wall and generate a new reflected C- characteristic.
     * 
     * @return Index of wall point, index of new C- chain
     */
    std::pair<size_t,size_t> reflect_c_plus_off_wall(size_t chain_idx, const CharacteristicPoint& pt);

    /** Terminate a C- characteristic off the central axis and generate a new reflected C+ characteristic.
     * 
     * @return Index of axis point, index of new C+ chain
     */
    std::pair<size_t,size_t> reflect_c_minus_off_axis(size_t chain_idx, const CharacteristicPoint& pt);

    /** Add a point where the C+ characteristic hits the wall without emitting a reflected C- characteristic.
     * This is used when designing minimum length nozzles.
     */
    size_t terminate_c_plus_at_wall(size_t chain_idx, const CharacteristicPoint& pt);

    /** Seed an initial wall point (e.g. the throat lip) that anchors the wall march.
     * The point owns no characteristic chain; it only bootstraps leading_wall_point()
     * and the wall coordinate lists.
     * @return index of the seeded point
     */
    size_t seed_wall_point(const CharacteristicPoint& pt);

    /** Mark a chain as inactive. */
    void terminate_chain(size_t chain_idx, ChainTermination termtype);

    /** Mark a chain as inactive while updating the last point. */
    void update_and_terminate_chain(size_t chain_idx, size_t last_pt_idx, ChainTermination termtype);

    /** True while at least one chain is still being marched. */
    bool has_active_chains() const;

    /** Every point at `axis_point_indices`, in march order. Kernel-independent: populated
     *  for both the chain-pairing ladder and the front-based (analysis/Rao) kernels. */
    std::vector<CharacteristicPoint> axis_points() const;
    /** Every point at `wall_point_indices`, in march order. Kernel-independent: populated
     *  for both the chain-pairing ladder and the front-based (analysis/Rao) kernels. */
    std::vector<CharacteristicPoint> wall_points() const;

    /** View of the metadata of active chains that belong to a given family (or all active
     *  chains, if `family` is UNSPECIFIED). */
    auto active_chains(CharacteristicFamily family = CharacteristicFamily::UNSPECIFIED) {
        auto is_active = [family](const ChainMetadata& meta) {
            if (family != CharacteristicFamily::UNSPECIFIED) {
                return meta.active && (meta.family == family);
            }
            else return meta.active;
        };

        return chain_metadata | std::views::filter(is_active);
    }
};

} //namespace Goddard