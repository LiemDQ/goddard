#pragma once
#include <vector>
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

/** Which family a characteristic belongs to: C+ (PLUS) or C- (MINUS). */
enum class CharacteristicFamily {
    UNSPECIFIED, PLUS, MINUS
};

/**
 * One node of the characteristic mesh: the flow state at a point, plus the two Riemann
 * invariants carried through it.
 *
 * For perfect-gas solves the thermodynamic quantities are normalized by their stagnation
 * values and `V` holds the Mach number; for frozen and equilibrium chemistry they are
 * dimensional SI.
 */
class CharacteristicPoint {
    public:
    double theta;       ///< Flow angle, in radians.
    double nu;          ///< Prandtl-Meyer angle (or generalized PM function), in radians.
    double pressure;    ///< Static pressure; Pa, or normalized by the stagnation pressure for perfect gas.
    double temperature; ///< Static temperature; K, or normalized by the stagnation temperature for perfect gas.
    double gamma_s;     ///< Local isentropic exponent.
    double mach;        ///< Local Mach number.
    double V;           ///< Velocity in m/s; equal to the Mach number for perfect gas.
    /// Mach angle asin(1/M), in radians. Sets the characteristic slopes theta +/- mu.
    double mu;
    double K_plus;      ///< Riemann invariant carried along the C+ characteristic.
    double K_minus;     ///< Riemann invariant carried along the C- characteristic.
    double x;           ///< Axial position, in length units.
    double y;           ///< Radial (or transverse) position, in length units.
    /// Serialized Cantera state; empty for perfect-gas solves. Internal detail.
    std::vector<double> cantera_state; 

    void update_thermodynamic_state_from_nu(ThermodynamicContext& ctxt, double nu, double mach_guess = 1.0);
    void update_thermodynamic_state_from_mach(ThermodynamicContext& ctxt, double mach);
    void update_thermodynamic_state_from_V(ThermodynamicContext& ctxt, double V);
    void update_Ks();

    private:
    void update_thermodynamic_state(ThermodynamicContext& ctxt);

    // FROZEN/EQUILIBRIUM counterpart to update_thermodynamic_state(): sets
    // temperature and pressure directly from the PrandtlMeyerTable's
    // temperature/pressure columns at the (idx, weight) pair already found by
    // the caller for the other interpolated columns (V, gamma_s, mach/nu,
    // cantera_state), instead of restoring a Cantera state per point.
    void update_thermodynamic_state_from_table(ThermodynamicContext& ctxt, size_t idx, double weight);
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
 * based on the intersection of the characteristics of two upstream parent points
 * and specified characteristic angles.
 * 
 * @param p1 Point along C- characteristic
 * @param p2 Point along C+ characteristic
 * @param angle1 C- characteristic angle
 * @param angle2 C+ characteristic angle
 */
std::pair<double, double> characteristic_intersection_with_angle(
    const CharacteristicPoint& p1, 
    const CharacteristicPoint& p2,
    double angle1,
    double angle2);


} // namespace Goddard