#pragma once

#include <vector>
#include <string>

#include "cea_loader.hpp"
#include "goddard/problem.hpp"

/**
 * CEA Test Utilities - helper functions for comparing Goddard vs CEA results
 */
class CEATestUtils {
public:
    struct ComparisonTolerance {
        double temperature_rel = 1e-2;      // 1% relative error for temperature
        double pressure_rel = 1e-2;         // 0.1% relative error for pressure
        double density_rel = 1e-2;          // 1% relative error for density
        double enthalpy_rel = 1e-2;         // 1% relative error for enthalpy
        double entropy_rel = 1e-2;          // 1% relative error for entropy
        double molecular_weight_rel = 1e-3; // 0.1% relative error for molecular weight
        double gamma_rel = 1e-2;            // 1% relative error for gamma
        double sound_speed_rel = 1e-2;      // 1% relative error for sound speed
        double thermo_deriv_rel = 5e-2;     // 5% relative error for dlV/dlP, dlV/dlT
        double mass_fraction_abs = 1e-4;    // Absolute error for mass fractions
    };
    
    struct ComparisonResult {
        bool passed;
        std::string parameter;
        double goddard_value;
        double cea_value;
        double relative_error;
        double tolerance_used;
        std::string message;
    };
    
    /**
     * Compare Goddard combustion result against CEA chamber state
     * @param goddard_state Goddard combustion result
     * @param cea_state CEA chamber state
     * @param tolerance Comparison tolerances
     * @return Vector of comparison results for each parameter
     */
    static std::vector<ComparisonResult> compare_chamber_states(
        const Goddard::ThermodynamicState& goddard_state,
        const CEAState& cea_state,
        const ComparisonTolerance& tolerance
    );

    // Overload with default tolerance
    static std::vector<ComparisonResult> compare_chamber_states(
        const Goddard::ThermodynamicState& goddard_state,
        const CEAState& cea_state
    );
    
    /**
     * Check if relative error is within tolerance
     * @param actual Actual value
     * @param expected Expected value
     * @param relative_tolerance Relative tolerance (e.g., 0.01 for 1%)
     * @return true if within tolerance
     */
    static bool within_relative_tolerance(double actual, double expected, double relative_tolerance);
    
    /**
     * Check if absolute error is within tolerance
     * @param actual Actual value
     * @param expected Expected value  
     * @param absolute_tolerance Absolute tolerance
     * @return true if within tolerance
     */
    static bool within_absolute_tolerance(double actual, double expected, double absolute_tolerance);
    
    /**
     * Calculate relative error between two values
     * @param actual Actual value
     * @param expected Expected value
     * @return Relative error (|actual - expected| / |expected|)
     */
    static double relative_error(double actual, double expected);
    
    /**
     * Print comparison results summary
     * @param results Vector of comparison results
     * @param test_name Name of test case
     */
    static void print_comparison_summary(
        const std::vector<ComparisonResult>& results,
        const std::string& test_name
    );
};
