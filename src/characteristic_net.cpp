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
        size_t minus_idx = create_chain(idx, Family::MINUS);
        m.c_minus_chain_idx = minus_idx;
    }
    if (cplus) {
        size_t plus_idx = create_chain(idx, Family::PLUS);
        m.c_plus_chain_idx = plus_idx;
    }
    membership.push_back(m);
    return idx;
}

void CharacteristicNet::add_initial_data_line(const std::vector<CharacteristicPoint>& init_pts) {
    // The data line is ordered from the axis (i == 0) up to the wall (i == num_pts).
    // The kernel pairs each C+ leading point with the nearest C- leading point above it
    // (i.e. the lower point contributes the C+, the upper point the C-). Therefore the
    // bottommost (axis) point owns only a C+ chain and the topmost (wall) point owns only
    // a C- chain; all interior points own both.
    size_t num_pts = init_pts.size() - 1;
    for (size_t i = 0; i <= num_pts; i++) {
        bool cplus = i != num_pts;
        bool cminus = i != 0;
        size_t idx = add_initialization_point(init_pts[i], cminus, cplus);
        if (i == 0) {
            axis_point_indices.push_back(idx);
        }
        if (i == num_pts) {
            wall_point_indices.push_back(idx);
        }
    }
}

void CharacteristicNet::add_initial_characteristic(const std::vector<CharacteristicPoint>& init_pts, Family fam) {
    bool cplus = fam != Family::PLUS;
    bool cminus = fam != Family::MINUS;
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
    update_and_terminate_chain(c_plus_chain_idx, idx, TerminationType::WALL);

    // start a new reflected C- chain originating at the wall point
    ChainMetadata minus_metadata;
    minus_metadata.family = Family::MINUS;
    minus_metadata.origin_point_idx = idx;
    minus_metadata.latest_point_idx = idx;
    size_t minus_chain_idx = push_chain({idx}, minus_metadata);
    membership[idx].c_minus_chain_idx = minus_chain_idx;

    wall_point_indices.push_back(idx);

    return {idx, minus_chain_idx};
}

std::pair<size_t,size_t> CharacteristicNet::reflect_c_minus_off_axis(size_t c_minus_chain_idx, const CharacteristicPoint& pt) {
    // Append the axis point to the incoming C- chain only. The reflected C+ chain
    // is created afterward, once it has a valid index, then recorded in membership.
    size_t idx = add_point(pt, {.c_plus_chain_idx = std::nullopt, .c_minus_chain_idx = c_minus_chain_idx});

    // cap off the incoming C- chain at the axis
    update_and_terminate_chain(c_minus_chain_idx, idx, TerminationType::AXIS);

    // start a new reflected C+ chain originating at the axis point
    ChainMetadata plus_metadata;
    plus_metadata.family = Family::PLUS;
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

    update_and_terminate_chain(c_plus_chain_idx, idx, ChainMetadata::TerminationType::WALL);
    // advance the leading wall point so the next wall solve uses this point as its predecessor
    wall_point_indices.push_back(idx);
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

size_t CharacteristicNet::create_chain(size_t pt_idx, Family family) {
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

void CharacteristicNet::terminate_chain(size_t chain_idx, TerminationType termtype) {
    ChainMetadata& meta = chain_metadata[chain_idx];
    meta.active = false;
    meta.termination = termtype;
}

void CharacteristicNet::update_and_terminate_chain(size_t chain_idx, size_t last_pt_idx, TerminationType termtype) {
    ChainMetadata& meta = chain_metadata[chain_idx];
    meta.active = false;
    meta.termination = termtype;
    meta.latest_point_idx = last_pt_idx;
}

std::vector<CharacteristicPoint> CharacteristicNet::outflow_points() const {
    std::vector<CharacteristicPoint> outflow_pts;
    for (auto&& meta: chain_metadata) {
        if (!meta.active && meta.termination == TerminationType::OUTFLOW) {
            outflow_pts.push_back(points[meta.latest_point_idx]);
        }
    }
    return outflow_pts;
}

std::vector<CharacteristicPoint> CharacteristicNet::axis_points() const {
    std::vector<CharacteristicPoint> axis_pts;
    for (auto&& meta: chain_metadata) {
        if (!meta.active && meta.termination == TerminationType::AXIS) {
            axis_pts.push_back(points[meta.latest_point_idx]);
        }
    }
    return axis_pts;
}

std::vector<CharacteristicPoint> CharacteristicNet::wall_points() const {
    std::vector<CharacteristicPoint> wall_pts;
    for (auto&& meta: chain_metadata) {
        if (!meta.active && meta.termination == TerminationType::WALL) {
            wall_pts.push_back(points[meta.latest_point_idx]);
        }
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