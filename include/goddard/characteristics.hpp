#pragma once
#include <vector>
#include <optional>
#include <string_view>

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
};

/**
 * Classification of the numerical/physical failure modes that can occur while
 * resolving a single unit process (interior/wall/axis point solve) or marching
 * the method-of-characteristics net as a whole.
 */
enum class MocErrorCode {
    NONE,                      ///< No error; the point/solve is valid.
    NEGATIVE_NU,               ///< Prandtl-Meyer angle (or generalized PM function) nu is negative beyond tolerance.
    NEGATIVE_THETA,            ///< Flow angle theta is negative beyond tolerance.
    SUBSONIC_MACH,             ///< Mach number is below 1.0; the supersonic compatibility relations no longer apply.
    NONFINITE_VALUE,           ///< One of the point's numeric fields (x, y, theta, nu, mach, mu) is NaN or infinite.
    PM_INVERSION_FAILED,       ///< The Prandtl-Meyer inversion (nu -> Mach), or an equivalent Mach rootfind, failed to converge.
    TABLE_RANGE_EXCEEDED,      ///< A PrandtlMeyerTable lookup (by nu, Mach, or velocity) fell outside the tabulated range.
    NON_DOWNSTREAM_POINT,      ///< The computed intersection lies at or behind (upstream of) one of its parent points.
    WALL_QUERY_OUT_OF_BOUNDS,  ///< A wall-profile query (e.g. theta_at) fell outside the profile's domain.
    INITIALIZATION_FAILED,     ///< Construction of the initial data line (transonic start line) failed to converge.
    MAX_ITERATIONS_REACHED,    ///< The characteristic kernel reached its iteration safety cap before all chains terminated.
    INCOMPLETE_MARCH           ///< Every chain terminated without the net reaching the exit plane; the field covers only part of the nozzle.
};

/**
 * Human-readable name for a MocErrorCode, for log/diagnostic messages.
 */
std::string_view to_string(MocErrorCode code);

/**
 * A CharacteristicPoint together with the error code (if any) produced while
 * computing it. Unit-process solvers return this instead of a bare
 * CharacteristicPoint so a numerical failure is reported as data -- instead of
 * being laundered into the point's fields (e.g. a sentinel Mach number) or
 * thrown as an exception.
 */
struct PointResult {
    CharacteristicPoint point;
    MocErrorCode error = MocErrorCode::NONE;
};

/**
 * Check a computed characteristic point for physical/numerical validity.
 *
 * @param pt  Point to validate.
 * @param tol Absolute tolerance for the nu/theta non-negativity checks (axis points
 *            legitimately pin theta to exactly 0.0, and interior points can carry
 *            small negative roundoff).
 * @return MocErrorCode::NONE if the point is valid, else the first violated
 *         condition in priority order: nu, theta, mach, finiteness.
 */
/**
 * Validity of a computed point: finite fields, supersonic, non-negative Prandtl-Meyer angle,
 * and (by default) non-negative flow angle.
 *
 * The flow-angle requirement is what the DIRECT kernel uses to catch a march that has gone
 * wrong. It is optional because a negative flow angle is physically admissible in the
 * analysis of an arbitrary contour (any locally converging wall) and because, wherever the
 * exact angle is zero -- the uniform exit region of a minimum-length nozzle, the axis
 * neighbourhood -- a computed angle is discretization noise of either sign, which a
 * tolerance of the solver's abstol cannot absorb.
 */
MocErrorCode check_point_validity(const CharacteristicPoint& pt, double tol,
                                  bool require_nonnegative_theta = true);

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