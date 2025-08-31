#include "cea_comparer.hpp"
#include <iostream>

// CEATestUtils implementation
std::vector<CEATestUtils::ComparisonResult> CEATestUtils::compare_chamber_states(
    const Goddard::StateInfo& goddard_state,
    const CEAState& cea_state,
    const ComparisonTolerance& tolerance) {
    
    std::vector<ComparisonResult> results;
    
    // Temperature comparison
    if (cea_state.temperature_k > 0) {
        ComparisonResult temp_result;
        temp_result.parameter = "temperature_k";
        temp_result.goddard_value = goddard_state.temperature;  // Assuming this field exists
        temp_result.cea_value = cea_state.temperature_k;
        temp_result.relative_error = relative_error(temp_result.goddard_value, temp_result.cea_value);
        temp_result.tolerance_used = tolerance.temperature_rel;
        temp_result.passed = within_relative_tolerance(
            temp_result.goddard_value, temp_result.cea_value, tolerance.temperature_rel);
        
        std::ostringstream msg;
        msg << "Temperature: Goddard=" << temp_result.goddard_value 
            << "K, CEA=" << temp_result.cea_value 
            << "K, RelErr=" << temp_result.relative_error;
        temp_result.message = msg.str();
        
        results.push_back(temp_result);
    }
    
    // Pressure comparison
    if (cea_state.pressure_bar > 0) {
        ComparisonResult press_result;
        press_result.parameter = "pressure_bar";
        press_result.goddard_value = goddard_state.pressure / 1e5;  // Convert Pa to bar
        press_result.cea_value = cea_state.pressure_bar;
        press_result.relative_error = relative_error(press_result.goddard_value, press_result.cea_value);
        press_result.tolerance_used = tolerance.pressure_rel;
        press_result.passed = within_relative_tolerance(
            press_result.goddard_value, press_result.cea_value, tolerance.pressure_rel);
        
        std::ostringstream msg;
        msg << "Pressure: Goddard=" << press_result.goddard_value 
            << "bar, CEA=" << press_result.cea_value 
            << "bar, RelErr=" << press_result.relative_error;
        press_result.message = msg.str();
        
        results.push_back(press_result);
    }
    
    // Molecular weight comparison
    if (cea_state.molecular_weight > 0) {
        ComparisonResult mw_result;
        mw_result.parameter = "molecular_weight";
        mw_result.goddard_value = goddard_state.molecular_weight;
        mw_result.cea_value = cea_state.molecular_weight;
        mw_result.relative_error = relative_error(mw_result.goddard_value, mw_result.cea_value);
        mw_result.tolerance_used = tolerance.molecular_weight_rel;
        mw_result.passed = within_relative_tolerance(
            mw_result.goddard_value, mw_result.cea_value, tolerance.molecular_weight_rel);
        
        std::ostringstream msg;
        msg << "Molecular Weight: Goddard=" << mw_result.goddard_value 
            << "g/mol, CEA=" << mw_result.cea_value 
            << "g/mol, RelErr=" << mw_result.relative_error;
        mw_result.message = msg.str();
        
        results.push_back(mw_result);
    }
    
    // TODO: Add more parameter comparisons as needed
    
    return results;
}

bool CEATestUtils::within_relative_tolerance(double actual, double expected, double relative_tolerance) {
    if (std::abs(expected) < 1e-12) {
        return std::abs(actual - expected) < relative_tolerance;
    }
    return std::abs(actual - expected) / std::abs(expected) <= relative_tolerance;
}

bool CEATestUtils::within_absolute_tolerance(double actual, double expected, double absolute_tolerance) {
    return std::abs(actual - expected) <= absolute_tolerance;
}

double CEATestUtils::relative_error(double actual, double expected) {
    if (std::abs(expected) < 1e-12) {
        return std::abs(actual - expected);
    }
    return std::abs(actual - expected) / std::abs(expected);
}

void CEATestUtils::print_comparison_summary(
    const std::vector<ComparisonResult>& results,
    const std::string& test_name) {
    
    std::cout << "\n=== CEA Comparison Results: " << test_name << " ===" << std::endl;
    
    int passed = 0;
    int total = results.size();
    
    for (const auto& result : results) {
        std::string status = result.passed ? "✓" : "✗";
        std::cout << status << " " << result.message << std::endl;
        if (result.passed) passed++;
    }
    
    std::cout << "Summary: " << passed << "/" << total << " tests passed" << std::endl;
}
