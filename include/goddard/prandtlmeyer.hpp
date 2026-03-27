#pragma once
#include "cantera/core.h"
#include <vector>
#include <memory>
#include <utility>

namespace Goddard {

// Prandtl-Meyer angle function
double prandtl_meyer(double mach, double gamma); 

//derivative of Prandtl-Meyer function w.r.t. mach number
double prandtl_meyer_derivative(double mach, double gamma); 

// inverse function: obtain Mach number from Prandtl-meter angle using Newton's method
double mach_from_prandtl_meyer(double nu, double gamma, 
                                double mach_guess = 0.0,
                                double tol = 1e-10,
                                int max_iter = 20);

/**
 * Generalized Prandtl-Meyer angle for a non-calorically perfect gas.
 */
double frozen_prandtl_meyer(Cantera::ThermoPhase& thermo, double mach);

/**
 * Table of thermodynamic data for Prandtl-Meyer expansion fans
 */
class PrandtlMeyerTable {
public:
    
    std::vector<double> velocities;
    std::vector<double> nus;
    std::vector<double> machs;
    std::vector<double> enthalpies;
    std::vector<double> dnu_dV;
    std::vector<double> gamma_s;
    std::vector<std::vector<double>> states;

    /**
     * Build table of (velocity, nu) pairs along an isentropic expansion.
     * Steps in pressure from throat conditions to a specified pressure ratio.
     *
     * @param thermo     Cantera ThermoPhase, set to throat conditions on entry
     * @param equilibrium  If true, equilibrate at each point; otherwise frozen composition
     * @param s0         Stagnation entropy (J/kg/K)
     * @param h0         Stagnation enthalpy (J/kg)
     * @param a_throat   Sound speed at the throat (m/s)
     * @param pressure_ratio  Ratio P_exit/P_throat to expand to (default: 1e-4)
     * @param num_points Number of table points
     */
    void build_table(
        Cantera::ThermoPhase& thermo,
        bool equilibrium,
        double s0,
        double h0,
        double a_throat,
        double pressure_ratio = 1e-4,
        size_t num_points = 500);
    
    bool is_built() const { return built; }

    bool is_equilibrium() const { return equil; }

    // Get velocity from nu.
    double interpolate_V(double nu) const;
    
    // Get velocity from mach number.
    double interpolate_V_from_mach(double mach) const;

    // Get mach number from nu.
    double interpolate_mach(double nu) const;

    // Get static enthalpy from nu.
    double interpolate_h_from_nu(double nu) const;
    
    // Get static enthalpy from mach number.
    double interpolate_h_from_mach(double mach) const;
    
    // Get static enthalpy from velocity.
    double interpolate_h(double V) const;

    // Get adiabatic index from nu.
    double interpolate_gamma_s_from_nu(double nu) const;

    // Get adiabatic index from mach number.
    double interpolate_gamma_s_from_mach(double mach) const;
        
    // Get nu from velocity.
    double interpolate_nu(double V) const;

    // Get nu form mach number.
    double interpolate_nu_from_mach(double mach) const;

    // Get dnu_dV from velocity.
    double interpolate_dnu(double V) const;

    std::vector<double> interpolate_state_from_mach(double mach) const;
    
    using IdxWeight = std::pair<size_t, double>;

    // Get index containing value closest to input Mach number.
    IdxWeight find_mach_index_and_weight(double mach) const;
    
    // Get index containing value closest to input nu.
    IdxWeight find_nu_index_and_weight(double nu) const;

    /**
     * Linear interpolation over vals with a specified weight.
     */
    double interpolate_at_index(size_t idx, 
        double weight, 
        const std::vector<double>& vals) const;
    
    /**
     * Linear interpolation over vals with a specified weight.
     */
    std::vector<double> interpolate_state_at_index(size_t idx, double weight) const;
    
    /**
     * Generic linear interpolation with binary search.
     */
    double interp(double query, 
                    const std::vector<double>& keys,
                    const std::vector<double>& vals) const;
    
    IdxWeight index_and_weight(double query, const std::vector<double>& keys) const;
private:
    bool built = false;
    bool equil = false;
    
    /**
     * Find the index closest to a queried value using binary search.
     * This works for nu, Mach number, and velocity 
     * as they monotonically increase with expansion.
     */
    size_t find_closest_nMv_index(double query, const std::vector<double>& keys) const;
                    
    std::vector<double> interp_state_vector(double query, const std::vector<double>& keys) const;

};



} // namespace Goddard