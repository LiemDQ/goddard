#include <utility>
#include "goddard/characteristic_net.hpp"

namespace Goddard {


bool CharacteristicNet::empty() const {
    return points.empty();
}

CharacteristicPoint& CharacteristicNet::leading_point(size_t chain_idx) {
    return points[chain_metadata[chain_idx].latest_point_idx];
}

const CharacteristicPoint& CharacteristicNet::leading_point(size_t chain_idx) const {
    return points[chain_metadata[chain_idx].latest_point_idx];
}

CharacteristicPoint& CharacteristicNet::leading_wall_point() {
    return points[wall_point_indices.back()];
}

const CharacteristicPoint& CharacteristicNet::leading_wall_point() const {
    return points[wall_point_indices.back()];
}

CharacteristicPoint& CharacteristicNet::leading_axis_point() {
    return points[axis_point_indices.back()];
}

const CharacteristicPoint& CharacteristicNet::leading_axis_point() const {
    return points[axis_point_indices.back()];
}

size_t CharacteristicNet::add_point(CharacteristicPoint pt, PointMembership m)
{
    points.push_back(std::move(pt));
    size_t idx = points.size() - 1;
    if (m.c_plus_chain_idx.has_value()) {
        c_chains[*m.c_plus_chain_idx].push_back(idx);
        chain_metadata[*m.c_plus_chain_idx].latest_point_idx = idx;
    }
    if (m.c_minus_chain_idx.has_value()) {
        c_chains[*m.c_minus_chain_idx].push_back(idx);
        chain_metadata[*m.c_minus_chain_idx].latest_point_idx = idx;
    }
    membership.push_back(std::move(m));
    return idx;
}

size_t CharacteristicNet::add_initialization_point(
    CharacteristicPoint pt, bool cminus, bool cplus)
{
    points.push_back(std::move(pt));
    size_t idx = points.size() - 1;
    PointMembership m;
    if (cminus) {
        size_t minus_idx = create_chain(idx, CharacteristicFamily::MINUS);
        m.c_minus_chain_idx = minus_idx;
    }
    if (cplus) {
        size_t plus_idx = create_chain(idx, CharacteristicFamily::PLUS);
        m.c_plus_chain_idx = plus_idx;
    }
    membership.push_back(m);
    return idx;
}

void CharacteristicNet::add_initial_data_line(
    const std::vector<CharacteristicPoint>& init_pts,
    std::optional<CharacteristicFamily> on_characteristic)
{
    // A data line collinear along a single characteristic (e.g. a centered expansion fan)
    // has a fundamentally different chain topology: there is only one characteristic of the
    // given family, not one per point. Treating it as a generic data line would pair adjacent
    // collinear points as opposing-family parents and produce degenerate intersections.
    if (on_characteristic.has_value()) {
        add_initial_characteristic(init_pts, *on_characteristic);
        return;
    }

    // A non-collinear transonic start line (e.g. Kliegel-Levine) crosses many
    // characteristics, so every interior point seeds both a C+ and a C- chain, exactly like
    // a fan-init data line -- this fills the near-axis void that a C+-only line would leave.
    // The axis bootstrap point (i=0) still seeds only a C+, mirroring
    // add_initial_characteristic's own axis special case: giving it a C- would immediately
    // re-reflect it off the axis as a degenerate point. The last point IS the wall point and
    // immediately reflects into the first C- chain, exactly like a wall reflection
    // encountered during marching.
    size_t last = init_pts.size() - 1;
    for (size_t i = 0; i < last; i++) {
        bool cminus = (i != 0);
        size_t idx = add_initialization_point(init_pts[i], /*cminus=*/cminus, /*cplus=*/true);
        if (i == 0) {
            axis_point_indices.push_back(idx);
        }
    }

    const CharacteristicPoint& wall_pt = init_pts[last];
    points.push_back(wall_pt);
    size_t wall_idx = points.size() - 1;
    membership.push_back(PointMembership{});
    wall_point_indices.push_back(wall_idx);
    wall_x.push_back(wall_pt.x);
    wall_y.push_back(wall_pt.y);

    ChainMetadata minus_metadata;
    minus_metadata.family = CharacteristicFamily::MINUS;
    minus_metadata.origin_point_idx = wall_idx;
    minus_metadata.latest_point_idx = wall_idx;
    size_t minus_chain_idx = push_chain({wall_idx}, minus_metadata);
    membership[wall_idx].c_minus_chain_idx = minus_chain_idx;
}

void CharacteristicNet::add_initial_characteristic(
    const std::vector<CharacteristicPoint>& init_pts, CharacteristicFamily fam)
{
    bool cplus = fam != CharacteristicFamily::PLUS;
    bool cminus = fam != CharacteristicFamily::MINUS;
    ChainMetadata init_metadata{
        .family = fam,
        .origin_point_idx = 0,
        .latest_point_idx = 0
    };
    size_t init_chain_index = push_chain({}, init_metadata);

    for (size_t i = 0; i < init_pts.size(); i++) {
        // The first point lies on the axis (the bootstrap reflection of the first ray).
        // It must not start a minor-family chain, or it would immediately re-reflect off
        // the axis as a degenerate point. It belongs only to the shared initial chain.
        bool make_minor = (i != 0);
        size_t idx = add_initialization_point(init_pts[i], cminus && make_minor, cplus && make_minor);
        c_chains[init_chain_index].push_back(idx);
        chain_metadata[init_chain_index].latest_point_idx = idx;
        PointMembership& mem = membership[idx];
        if (!cplus) {
            mem.c_plus_chain_idx = init_chain_index;
        }
        if (!cminus) {
            mem.c_minus_chain_idx = init_chain_index;
        }
    }
}

std::pair<size_t,size_t> CharacteristicNet::reflect_c_plus_off_wall(
    size_t c_plus_chain_idx, const CharacteristicPoint& pt)
{
    // Append the wall point to the incoming C+ chain only. The reflected C- chain
    // is created afterward, once it has a valid index, then recorded in membership.
    size_t idx = add_point(pt, {.c_plus_chain_idx = c_plus_chain_idx, .c_minus_chain_idx = std::nullopt});

    // cap off the incoming C+ chain at the wall
    update_and_terminate_chain(c_plus_chain_idx, idx, ChainTermination::WALL);

    // start a new reflected C- chain originating at the wall point
    ChainMetadata minus_metadata;
    minus_metadata.family = CharacteristicFamily::MINUS;
    minus_metadata.origin_point_idx = idx;
    minus_metadata.latest_point_idx = idx;
    size_t minus_chain_idx = push_chain({idx}, minus_metadata);
    membership[idx].c_minus_chain_idx = minus_chain_idx;

    wall_point_indices.push_back(idx);
    wall_x.push_back(pt.x);
    wall_y.push_back(pt.y);

    return {idx, minus_chain_idx};
}

std::pair<size_t,size_t> CharacteristicNet::reflect_c_minus_off_axis(size_t c_minus_chain_idx, const CharacteristicPoint& pt) {
    // Append the axis point to the incoming C- chain only. The reflected C+ chain
    // is created afterward, once it has a valid index, then recorded in membership.
    size_t idx = add_point(pt, {.c_plus_chain_idx = std::nullopt, .c_minus_chain_idx = c_minus_chain_idx});

    // cap off the incoming C- chain at the axis
    update_and_terminate_chain(c_minus_chain_idx, idx, ChainTermination::AXIS);

    // start a new reflected C+ chain originating at the axis point
    ChainMetadata plus_metadata;
    plus_metadata.family = CharacteristicFamily::PLUS;
    plus_metadata.origin_point_idx = idx;
    plus_metadata.latest_point_idx = idx;
    size_t plus_chain_idx = push_chain({idx}, plus_metadata);
    membership[idx].c_plus_chain_idx = plus_chain_idx;

    axis_point_indices.push_back(idx);

    return {idx, plus_chain_idx};
}

size_t CharacteristicNet::terminate_c_plus_at_wall(size_t c_plus_chain_idx, const CharacteristicPoint& pt) {
    // add_point already appends idx to the C+ chain and updates its leading edge.
    size_t idx = add_point(pt, {.c_plus_chain_idx = c_plus_chain_idx, .c_minus_chain_idx = std::nullopt});

    update_and_terminate_chain(c_plus_chain_idx, idx, ChainTermination::WALL);
    // advance the leading wall point so the next wall solve uses this point as its predecessor
    wall_point_indices.push_back(idx);
    wall_x.push_back(pt.x);
    wall_y.push_back(pt.y);
    return idx;
}

size_t CharacteristicNet::seed_wall_point(const CharacteristicPoint& pt) {
    // Seed an initial wall point (e.g. the throat lip) that bootstraps the wall march.
    // It owns no characteristic chain; it only anchors leading_wall_point() and the
    // wall coordinate lists.
    points.push_back(pt);
    size_t idx = points.size() - 1;
    membership.push_back(PointMembership{});
    wall_point_indices.push_back(idx);
    wall_x.push_back(pt.x);
    wall_y.push_back(pt.y);
    return idx;
}

size_t CharacteristicNet::create_chain(size_t pt_idx, CharacteristicFamily family) {
    ChainMetadata metadata;
    metadata.family = family;
    metadata.origin_point_idx = pt_idx;
    metadata.latest_point_idx = pt_idx;

    return push_chain({pt_idx}, metadata);
}

size_t CharacteristicNet::push_chain(std::vector<size_t>&& chain, ChainMetadata metadata) {
    c_chains.push_back(chain);
    chain_metadata.push_back(metadata);

    return c_chains.size() - 1;
}

void CharacteristicNet::terminate_chain(size_t chain_idx, ChainTermination termtype) {
    ChainMetadata& meta = chain_metadata[chain_idx];
    meta.active = false;
    meta.termination = termtype;
}

void CharacteristicNet::update_and_terminate_chain(size_t chain_idx, size_t last_pt_idx, ChainTermination termtype) {
    ChainMetadata& meta = chain_metadata[chain_idx];
    meta.active = false;
    meta.termination = termtype;
    meta.latest_point_idx = last_pt_idx;
}

std::vector<size_t> CharacteristicNet::add_front(const std::vector<CharacteristicPoint>& front_points) {
    std::vector<size_t> indices;
    indices.reserve(front_points.size());
    for (const CharacteristicPoint& pt : front_points) {
        points.push_back(pt);
        membership.push_back(PointMembership{});
        indices.push_back(points.size() - 1);
    }
    axis_point_indices.push_back(indices.front());
    wall_point_indices.push_back(indices.back());
    wall_x.push_back(front_points.back().x);
    wall_y.push_back(front_points.back().y);
    fronts.push_back(indices);
    return indices;
}

std::vector<CharacteristicPoint> CharacteristicNet::axis_points() const {
    std::vector<CharacteristicPoint> axis_pts;
    axis_pts.reserve(axis_point_indices.size());
    for (size_t idx : axis_point_indices) {
        axis_pts.push_back(points[idx]);
    }
    return axis_pts;
}

std::vector<CharacteristicPoint> CharacteristicNet::wall_points() const {
    std::vector<CharacteristicPoint> wall_pts;
    wall_pts.reserve(wall_point_indices.size());
    for (size_t idx : wall_point_indices) {
        wall_pts.push_back(points[idx]);
    }
    return wall_pts;
}

bool CharacteristicNet::has_active_chains() const {
    bool active = false;
    for (const auto& meta : chain_metadata) {
        active |= meta.active;
        if (active) return true;
    }
    return false;
}
    
} //namespace Goddard