#pragma once
#include <vector>
#include <optional>
#include <utility>
#include <ranges>
#include "goddard/characteristics.hpp"

namespace Goddard {

/**
 * Bookkeeping for one characteristic chain in a CharacteristicNet.
 *
 * A chain is a single characteristic line, stored as an ordered list of point indices. It is
 * born at the initial data line or at a wall/axis reflection, and stops when it terminates.
 */
struct ChainMetadata {
    /// False once the chain has terminated; only active chains are still marched.
    bool active = true;
    /// How a chain stopped being marched.
    enum class TerminationType {
        NOT_TERMINATED, WALL, AXIS, OUTFLOW, CORNER_FAN_ORIGIN,
        // Retired by mesh control because the front crowded around it (see
        // MocNozzle::control_front_spacing). Deliberately distinct from OUTFLOW: only
        // OUTFLOW-terminated chains contribute to outflow_points(), so a merged chain
        // must not be mistaken for one that reached the exit plane.
        MERGED
    } termination = TerminationType::NOT_TERMINATED;
    /// Which characteristic family the chain belongs to.
    enum class Family {
        UNSPECIFIED, PLUS, MINUS
    } family;
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

class CharacteristicNet {

    public:
    using Family = ChainMetadata::Family;
    using TerminationType = ChainMetadata::TerminationType;
    
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

    /** Axial coordinates of the wall points, in length units. */
    std::vector<double> wall_x;
    /** Radial coordinates of the wall points, in length units. */
    std::vector<double> wall_y;

    /**
     * Inverse march only (MocMarchScheme::INVERSE): point indices of every marching front,
     * axis to wall, in order. `fronts.front()` is the initial front F_0 and `fronts.back()`
     * is the exit plane once the march has completed. Empty for the DIRECT kernel, which
     * has no synchronized front -- its topology lives entirely in `c_chains`.
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
     * the first C- chain, exactly as a wall reflection would during marching. This requires
     * the line to be lifted off the raw sonic (v=0) locus first (see
     * MocInitialization::initialize_kliegel_levine and
     * MocOptions::initial_line_axial_shift) -- pairing two points still on the raw sonic
     * locus directly can land behind both parents, since mu -> 90 deg there near the axis.
     * The first point (the axis bootstrap) seeds only a C+: giving it a C- would immediately
     * re-reflect it off the axis as a degenerate point. When `on_characteristic` is set the
     * points are collinear along a single characteristic of that family (e.g. a centered
     * expansion fan), and are seeded as one shared chain via add_initial_characteristic.
     *
     * The caller must supply the topology; it cannot be inferred from the points, since a
     * Riemann invariant is constant along a characteristic only for planar flow.
     */
    void add_initial_data_line(
        const std::vector<CharacteristicPoint>& points,
        std::optional<Family> on_characteristic = std::nullopt);

    /** Add points from an initial data line that happens to be along a characteristic.
     * This is primarily used when the initial data line originates from a centered expansion.
    */
    void add_initial_characteristic(const std::vector<CharacteristicPoint>& points, Family family = Family::PLUS);

    /** Create a new characteristic with `pt_idx` as the starting point, 
     * and add it to the tracking lists.
     * 
     * @return Index of the chain
     */
    size_t create_chain(size_t pt_idx, Family family);
    
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

    /** Point index and the two chain indices created by insert_rung. */
    struct InsertedRung { size_t point_idx; size_t c_plus_chain_idx; size_t c_minus_chain_idx; };

    /** Insert a new mesh point into the marching front, owning a fresh chain of each family.
     *
     * Every interior point of the net is simultaneously the leading edge of one C+ chain and
     * one C- chain, and it is that pairing the kernel marches. Refining the front therefore
     * means adding a point that owns both, so the oversized step it splits becomes two
     * ordinary unit processes on the next pass. See MocNozzle::refine_front for when and why
     * the front needs refining.
     *
     * @return index of the new point and of the C+ and C- chains it originates
     */
    InsertedRung insert_rung(const CharacteristicPoint& pt);

    /** Terminate the C+ and C- chains led by `point_idx`, removing it from the marching front.
     *
     * The inverse of insert_rung: where insert_rung splits an over-stretched front segment,
     * this retires a rung the front has crowded around, so a compression region cannot drive
     * adjacent front points together without bound. The point itself stays in `points` --
     * every index in the net is permanent, and the retired point remains a valid interior
     * node of the chains that already passed through it. Only its two *leading* chains stop.
     *
     * @return the retired C+ and C- chain indices, or nullopt when `point_idx` does not
     *         currently lead an active chain of each family (in which case nothing changes).
     */
    std::optional<std::pair<size_t, size_t>> retire_rung(size_t point_idx);

    /** Seed an initial wall point (e.g. the throat lip) that anchors the wall march.
     * The point owns no characteristic chain; it only bootstraps leading_wall_point()
     * and the wall coordinate lists.
     * @return index of the seeded point
     */
    size_t seed_wall_point(const CharacteristicPoint& pt);

    // Mark a chain as inactive.
    void terminate_chain(size_t chain_idx, TerminationType termtype);
    
    // Mark a chain as inactive while updating the last point. 
    void update_and_terminate_chain(size_t chain_idx, size_t last_pt_idx, TerminationType termtype);

    /** True while at least one chain is still being marched. */
    bool has_active_chains() const;

    /**
     * The leading points of every chain that terminated by flowing out of the domain.
     *
     * Only OUTFLOW-terminated chains contribute, so chains retired by mesh control
     * (ChainMetadata::TerminationType::MERGED) are excluded -- they never reached the exit
     * plane and must not be mistaken for points that did.
     */
    std::vector<CharacteristicPoint> outflow_points() const;
    /** Every point lying on the centerline, in march order. */
    std::vector<CharacteristicPoint> axis_points() const;
    /** Every point lying on the nozzle wall, in march order. */
    std::vector<CharacteristicPoint> wall_points() const;

    // Generate view of metadata of active chains that belong to a given family.
    auto active_chains(Family family = Family::UNSPECIFIED) {
        auto is_active = [family](const ChainMetadata& meta) { 
            if (family != Family::UNSPECIFIED) {
                return meta.active && (meta.family == family);
            }
            else return meta.active;
        };

        return chain_metadata | std::views::filter(is_active);
    }
};

} //namespace Goddard