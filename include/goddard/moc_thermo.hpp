#pragma once
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
 * The single thermodynamic dispatch shared by initialization, both kernels, and every unit
 * process. Failure is reported as a MocErrorCode return rather than by throwing, so a caller
 * can attach point/pass context before surfacing a MocFailure. Immutable once built.
 */
class MocThermo {
public:
    /** Perfect gas with constant gamma; T and P normalized by stagnation values (T_ref = P_ref = 1). */
    static MocThermo perfect_gas(double gamma);
    /**
     * Frozen or equilibrium real gas: builds the Prandtl-Meyer table from the throat state.
     * `gas` is left at the throat state on return.
     */
    static MocThermo tabulated(Gas& gas, GasChemistry chemistry, const ThroatCondition& throat);

    /** Which dispatch branch set_state_from_*() takes. */
    GasChemistry chemistry() const;
    double gamma() const;                     ///< Perfect gas only: the constant gamma.
    const PrandtlMeyerTable& table() const;    ///< Real gas only.
    /** Real gas only: throat temperature (K) that table T is normalized against by set_pressure_temperature_from_table. */
    double T_ref() const;
    /** Real gas only: throat pressure (Pa) that table P is normalized against by set_pressure_temperature_from_table. */
    double P_ref() const;

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
    /** Set `point`'s full state from a velocity (m/s for real gas; Mach number for perfect gas). */
    MocErrorCode set_state_from_V(CharacteristicPoint& point, double V) const;

private:
    GasChemistry m_chemistry = GasChemistry::PERFECT_GAS;
    double m_gamma = 1.4;
    PrandtlMeyerTable m_table;
    double m_T_ref = 1.0;
    double m_P_ref = 1.0;

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
