#pragma once
#include <vector>
#include <utility>
#include <cmath>
#include <optional>

#include "goddard/prandtlmeyer.hpp"
#include "goddard/gas.hpp"
namespace Goddard {

struct ThermodynamicContext {
    std::optional<Gas> gas;
    const PrandtlMeyerTable& table;
    double T_ref = 273.15;
    double P_ref = 101325.0;
    double gamma_s = 1.4;
};
class CharacteristicPoint {
    public:
    double theta; // flow angle
    double nu; // Prandtl-Meyer angle (or generalized PM function)
    double pressure;
    double temperature;
    double gamma_s; // local isentropic gamma
    double mach;
    double V; //velocity
    // Characteristic slopes
    double mu; // Mach angle = asin(1/M)
    // Riemann invariants
    double K_plus; 
    double K_minus;
    // geometry
    double x;
    double y;
    // for chemistry
    std::vector<double> cantera_state; 

    void update_thermodynamic_state_from_nu(ThermodynamicContext& ctxt, double nu, double mach_guess = 1.0);
    void update_thermodynamic_state_from_mach(ThermodynamicContext& ctxt, double mach);
    void update_thermodynamic_state_from_V(ThermodynamicContext& ctxt, double V);
    void update_Ks();

    private:
    void update_thermodynamic_state(ThermodynamicContext& ctxt);
};

constexpr double average_angle(double angle1, double angle2) {
    return 0.5*(angle1+angle2);
}

constexpr double average_cminus_angle(
    const CharacteristicPoint& p1, const CharacteristicPoint& p2)
{
    return average_angle(p1.theta-p1.mu, p2.theta-p2.mu);
}

constexpr double average_cplus_angle(
    const CharacteristicPoint& p1, const CharacteristicPoint& p2) 
{
    return average_angle(p1.theta + p1.mu, p2.theta + p2.mu);
}

/**
 * Set temperature and pressure based on isentropic relations and thermodynamic
 * state of upstream characteristic node.
 */
void characteristic_isentropic_PT_from_parent(
    CharacteristicPoint& point, const CharacteristicPoint& parent);

/**
 * Get the coordinates of a downstream characteristic, 
 * based on the intersection of the characteristics of two upstream parent points.
 */
std::pair<double, double> characteristic_intersection_coordinates(
    const CharacteristicPoint& p1, 
    const CharacteristicPoint& p2,
    double angle1,
    double angle2);


struct ChainMetadata {
    bool active = true;
    enum class TerminationType {
        NOT_TERMINATED, WALL, AXIS, OUTFLOW, CORNER_FAN_ORIGIN
    } termination = TerminationType::NOT_TERMINATED;
    size_t origin_point_idx;
    size_t latest_point_idx;
    size_t left_neighbor_chain_idx;
    size_t right_neighbor_chain_idx;
};
class CharacteristicNet {
    
    public:
    
    // Index is opaque; access via family chains
    std::vector<CharacteristicPoint> points;
    
    using Wavefront = std::vector<CharacteristicPoint>;
    // access point at (i_plus, j_minus)
    std::vector<Wavefront> wavefronts;

    

    // Each entry is the chain of point indices along one characteristic
    // Note that wall and axis points terminate a chain, and also start the next chain
    std::vector<std::vector<size_t>> c_plus_chains;
    std::vector<std::vector<size_t>> c_minus_chains;

    // Given a point index, which chains does it belong to?
    struct PointMembership {size_t c_plus_chain_idx; size_t c_minus_chain_idx; };
    // membership at index i describes point i
    std::vector<PointMembership> membership;

    std::vector<ChainMetadata> c_plus_metadata;
    std::vector<ChainMetadata> c_minus_metadata;
    
    // Wall contour (profile output)
    std::vector<double> wall_x;
    std::vector<double> wall_y;
    
    // Wall and axis lists for boundary tracking
    std::vector<size_t> wall_point_indices;
    std::vector<size_t> axis_point_indices;

    size_t add_point(const CharacteristicPoint& pt, PointMembership m);
    size_t add_starting_point(const CharacteristicPoint& pt);

    size_t add_plus_chain(std::vector<size_t> chain, ChainMetadata metadata);
    
    size_t reflect_c_plus_off_wall(size_t chain_idx, const CharacteristicPoint& pt);
    size_t reflect_c_minus_off_axis(size_t chain_idx, const CharacteristicPoint& pt);

    // Add a point where the C+ characteristic hits the wall without emitting a reflected C- characteristic.
    size_t terminate_c_plus_at_wall(size_t chain_idx, const CharacteristicPoint& pt);
};

} // namespace Goddard