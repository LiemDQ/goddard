#pragma once
#include <string>

#include "goddard/chemistry.hpp"
#include "goddard/characteristics.hpp"
#include "goddard/moc.hpp"   // MocErrorCode
#include "goddard/prandtlmeyer.hpp"

namespace Goddard {

class Gas;
struct ThroatCondition;

/**
 * The isentrope every point of a MoC solve lies on: state as a function of nu, Mach or
 * velocity, for either a perfect gas or a tabulated frozen/equilibrium real gas.
 *
 * Points are nondimensional: pressure P/P0, temperature T/T0 and velocity V/a0, all relative
 * to the stagnation state.
 *
 * The single thermodynamic dispatch shared by initialization, both kernels, and every unit
 * process. Failure is reported as a MocErrorCode return rather than by throwing, so a caller
 * can attach point/pass context before surfacing a MocFailure. Immutable once built.
 */
class MocThermo {
public:
    /** Perfect gas with constant gamma. */
    static MocThermo perfect_gas(double gamma);
    /**
     * Frozen or equilibrium real gas: builds the Prandtl-Meyer table from the throat state.
     * `gas` must be at the stagnation state on entry, which sets `T_ref`, `P_ref` and `V_ref`,
     * and is left at it on return.
     *
     * The table ends before the first state far outside the temperature range of the species
     * data (`Gas::check_temperature_range()` fails there; see `table_end()`), so that a point
     * beyond it fails as a lookup outside the table.
     */
    static MocThermo tabulated(Gas& gas, GasChemistry chemistry, const ThroatCondition& throat);

    /** Which dispatch branch set_state_from_*() takes. */
    GasChemistry chemistry() const;
    double gamma() const;                     ///< Perfect gas only: the constant gamma.
    const PrandtlMeyerTable& table() const;    ///< Real gas only.
    /**
     * Why the table ends before the pressure ratio it was built to: the failure of the range
     * check at the first state dropped. Empty when the table is complete.
     */
    const std::string& table_end() const;
    /** Stagnation temperature T0 [K] that point temperatures are divided by; 1 for a perfect gas. */
    double T_ref() const;
    /** Stagnation pressure P0 [Pa] that point pressures are divided by; 1 for a perfect gas. */
    double P_ref() const;
    /**
     * Stagnation speed of sound a0 [m/s], in the chemistry's composition model, that point
     * velocities are divided by; 1 for a perfect gas.
     */
    double V_ref() const;

    // Convenience accessors; not called by the kernels themselves, which go through the
    // set_state_from_* chokepoints below.
    /** Prandtl-Meyer angle (or generalized PM function), in radians, at a given Mach number. */
    double nu_from_mach(double mach) const;
    /** Mach number at a given Prandtl-Meyer angle (radians); `mach_guess` seeds perfect-gas Newton iteration. */
    double mach_from_nu(double nu, double mach_guess = 0.0) const;
    /** Local isentropic exponent at a given Mach number. */
    double gamma_s_from_mach(double mach) const;
    /** Local isentropic exponent at a given Prandtl-Meyer angle (radians). */
    double gamma_s_from_nu(double nu) const;

    // The chokepoints. On failure the point's fields are unspecified and the code says why.
    /** Set `point`'s full state from a Prandtl-Meyer angle nu (radians); `mach_guess` seeds perfect-gas Newton iteration. */
    MocErrorCode set_state_from_nu(CharacteristicPoint& point, double nu, double mach_guess = 0.0) const;
    /** Set `point`'s full state from a Mach number. */
    MocErrorCode set_state_from_mach(CharacteristicPoint& point, double mach) const;
    /** Set `point`'s full state from a velocity V/a0 [-]. */
    MocErrorCode set_state_from_V(CharacteristicPoint& point, double V) const;

private:
    GasChemistry m_chemistry = GasChemistry::PERFECT_GAS;
    double m_gamma = 1.4;
    PrandtlMeyerTable m_table;
    std::string m_table_end;
    double m_T_ref = 1.0;
    double m_P_ref = 1.0;
    double m_V_ref = 1.0;

    // Sets point.temperature/point.pressure from the isentropic stagnation relations
    // (perfect gas only; FROZEN/EQUILIBRIUM points go through set_pressure_temperature_from_table).
    void set_pressure_temperature(CharacteristicPoint& point) const;

    // FROZEN/EQUILIBRIUM counterpart: sets temperature and pressure directly from the
    // PrandtlMeyerTable's temperature/pressure columns at the (idx, weight) pair already
    // found by the caller for the other interpolated columns (V, gamma_s, mach/nu,
    // cantera_state), instead of restoring a Cantera state per point.
    void set_pressure_temperature_from_table(CharacteristicPoint& point, size_t idx, double weight) const;
};

} // namespace Goddard
