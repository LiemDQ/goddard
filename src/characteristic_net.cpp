#include <utility>
#include "goddard/characteristic_net.hpp"

namespace Goddard {

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

std::pair<double, double> CharacteristicNet::get_predicted_chain_intersection(size_t plus_idx, size_t minus_idx) const {
    const CharacteristicPoint& plus_pt = leading_point(plus_idx);
    const CharacteristicPoint& minus_pt = leading_point(minus_idx);
    
    return characteristic_intersection(minus_pt, plus_pt);
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
        chain_metadata[*m.c_plus_chain_idx].latest_point_idx = idx;
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
        *m.c_minus_chain_idx = minus_idx;
    }
    if (cplus) {
        size_t plus_idx = create_chain(idx, Family::PLUS);
        *m.c_plus_chain_idx = plus_idx;
    }
    membership.push_back(m);
    return idx;
}

void CharacteristicNet::add_initial_data_line(const std::vector<CharacteristicPoint>& init_pts) {
    bool cminus = false;
    bool cplus = false;
    size_t num_pts = init_pts.size() - 1;
    for (size_t i = 0; i <= num_pts; i++) {
        cplus = i != 0;
        cminus = i != num_pts;
        size_t idx = add_initialization_point(init_pts[i], cminus, cplus);
        if (i == 0) {
            wall_point_indices.push_back(idx);
        }
        if (i == num_pts) {
            axis_point_indices.push_back(idx);
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

    for (const auto& pt : init_pts) {
        size_t idx = add_initialization_point(pt, cminus, cplus);
        c_chains[init_chain_index].push_back(idx);
        chain_metadata[init_chain_index].latest_point_idx = idx;
        PointMembership& mem = membership[idx];
        if (!cplus) {
            *mem.c_plus_chain_idx = init_chain_index;
        }
        if (!cminus) {
            *mem.c_minus_chain_idx = init_chain_index;
        }
    }
}

std::pair<size_t,size_t> CharacteristicNet::reflect_c_plus_off_wall(
    size_t c_plus_chain_idx, const CharacteristicPoint& pt) 
{
    const size_t new_c_minus_chain_idx = c_chains.size();
    size_t idx = add_point(pt, {.c_plus_chain_idx = c_plus_chain_idx, .c_minus_chain_idx = new_c_minus_chain_idx});

    // cap off the chain and create a new one
    c_chains[c_plus_chain_idx].push_back(idx);
    terminate_chain(c_plus_chain_idx, idx, TerminationType::WALL);

    ChainMetadata minus_metadata;
    minus_metadata.family = Family::MINUS;
    minus_metadata.origin_point_idx = idx;
    minus_metadata.latest_point_idx = idx;
    size_t minus_chain_idx = push_chain({idx}, minus_metadata);

    wall_point_indices.push_back(idx);
    
    return {idx, minus_chain_idx};
}

std::pair<size_t,size_t> CharacteristicNet::reflect_c_minus_off_axis(size_t c_minus_chain_idx, const CharacteristicPoint& pt) {
    const size_t new_c_plus_chain_idx = c_chains.size();
    size_t idx = add_point(pt, {.c_plus_chain_idx = new_c_plus_chain_idx, .c_minus_chain_idx = c_minus_chain_idx});

    // cap off the chain and create a new one
    
    c_chains[c_minus_chain_idx].push_back(idx);
    terminate_chain(c_minus_chain_idx, idx,TerminationType::AXIS);

    ChainMetadata plus_metadata;
    plus_metadata.family = Family::PLUS;
    plus_metadata.origin_point_idx = idx;
    plus_metadata.latest_point_idx = idx;

    size_t plus_chain_idx = push_chain({idx}, plus_metadata);

    axis_point_indices.push_back(idx);
    
    return {idx, plus_chain_idx};
}

size_t CharacteristicNet::terminate_c_plus_at_wall(size_t c_plus_chain_idx, const CharacteristicPoint& pt) {
    size_t idx = add_point(pt, {.c_plus_chain_idx = c_plus_chain_idx, .c_minus_chain_idx = std::nullopt});

    c_chains[c_plus_chain_idx].push_back(idx);
    terminate_chain(c_plus_chain_idx, idx, ChainMetadata::TerminationType::WALL);
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

void CharacteristicNet::terminate_chain(size_t chain_idx, size_t last_pt_idx, TerminationType termtype) {
    ChainMetadata& meta = chain_metadata[chain_idx];
    meta.active = false;
    meta.termination = termtype;
    meta.latest_point_idx = last_pt_idx;
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