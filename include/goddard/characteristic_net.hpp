#pragma once
#include <vector>
#include <optional>
#include <utility>
#include <ranges>
#include "goddard/characteristics.hpp"

namespace Goddard {

struct ChainMetadata {
    bool active = true;
    enum class TerminationType {
        NOT_TERMINATED, WALL, AXIS, OUTFLOW, CORNER_FAN_ORIGIN
    } termination = TerminationType::NOT_TERMINATED;
    enum class Family {
        UNSPECIFIED, PLUS, MINUS
    } family;
    size_t origin_point_idx;
    size_t latest_point_idx; //indicates leading edge
};

struct PointMembership {std::optional<size_t> c_plus_chain_idx; std::optional<size_t> c_minus_chain_idx; };

class CharacteristicNet {

    public:
    using Family = ChainMetadata::Family;
    using TerminationType = ChainMetadata::TerminationType;
    
    // Index is opaque; access via family chains
    std::vector<CharacteristicPoint> points;

    
    // Each entry is the chain of point indices along one characteristic
    // Note that wall and axis points terminate a chain, and also start the next chain
    std::vector<std::vector<size_t>> c_chains;

    // Given a point index, which chains does it belong to?
    // membership at index i describes point i
    std::vector<PointMembership> membership;

    std::vector<ChainMetadata> chain_metadata;
    
    // Wall and axis lists for boundary tracking
    std::vector<size_t> wall_point_indices;
    std::vector<size_t> axis_point_indices;

    // wall coordinates
    std::vector<double> wall_x;
    std::vector<double> wall_y;

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
     * Get the predicted x-coordinate of the intersection between a C+ characteristic
     * and a C- characteristic.
     * @return (x,y) coordinates of intersection point
     */
    std::pair<double, double> get_predicted_chain_intersection(size_t plus_idx, size_t minus_idx) const;

    /** Add a point and its membership to list of points, and 
     * appends it to its member characteristic chains. 
     * @return index of the added point
     */
    size_t add_point(CharacteristicPoint pt, PointMembership m);

    /** Add a point along an initial data line. 
     * @return index of the added point and corresponding characteristics 
     */
    size_t add_initialization_point(CharacteristicPoint pt, bool cminus, bool cplus);
    
    /** Add points from an initial data line, where the points are not 
     * coincident along a single characteristic.
     */
    void add_initial_data_line(const std::vector<CharacteristicPoint>& points);

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

    // Mark a chain as inactive.
    void terminate_chain(size_t chain_idx, TerminationType termtype);
    
    // Mark a chain as inactive while updating the last point. 
    void update_and_terminate_chain(size_t chain_idx, size_t last_pt_idx, TerminationType termtype);

    bool has_active_chains() const;

    std::vector<CharacteristicPoint> outflow_points() const;
    std::vector<CharacteristicPoint> axis_points() const;
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