#include "cea_comparer.hpp"
#include <iostream>

// CEATestUtils implementation - overload with default tolerance
std::vector<CEATestUtils::ComparisonResult> CEATestUtils::compare_chamber_states(
    const Goddard::ThermodynamicState& goddard_state,
    const CEAState& cea_state) {
    return compare_chamber_states(goddard_state, cea_state, ComparisonTolerance{});
}

std::vector<CEATestUtils::ComparisonResult> CEATestUtils::compare_chamber_states(
    const Goddard::ThermodynamicState& goddard_state,
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

    // Density comparison (both in kg/m3)
    if (cea_state.density_kg_m3 > 0) {
        ComparisonResult density_result;
        density_result.parameter = "density";
        density_result.goddard_value = goddard_state.density;
        density_result.cea_value = cea_state.density_kg_m3;
        density_result.relative_error = relative_error(density_result.goddard_value, density_result.cea_value);
        density_result.tolerance_used = tolerance.density_rel;
        density_result.passed = within_relative_tolerance(
            density_result.goddard_value, density_result.cea_value, tolerance.density_rel);

        std::ostringstream msg;
        msg << "Density: Goddard=" << density_result.goddard_value
            << "kg/m3, CEA=" << density_result.cea_value
            << "kg/m3, RelErr=" << density_result.relative_error;
        density_result.message = msg.str();

        results.push_back(density_result);
    }

    // Enthalpy comparison (Goddard J/kg, CEA kJ/kg - convert)
    if (std::abs(cea_state.enthalpy_kj_kg) > 0) {
        ComparisonResult h_result;
        h_result.parameter = "enthalpy";
        h_result.goddard_value = goddard_state.enthalpy / 1000.0;  // J/kg -> kJ/kg
        h_result.cea_value = cea_state.enthalpy_kj_kg;
        h_result.relative_error = relative_error(h_result.goddard_value, h_result.cea_value);
        h_result.tolerance_used = tolerance.enthalpy_rel;
        h_result.passed = within_relative_tolerance(
            h_result.goddard_value, h_result.cea_value, tolerance.enthalpy_rel);

        std::ostringstream msg;
        msg << "Enthalpy: Goddard=" << h_result.goddard_value
            << "kJ/kg, CEA=" << h_result.cea_value
            << "kJ/kg, RelErr=" << h_result.relative_error;
        h_result.message = msg.str();

        results.push_back(h_result);
    }

    // Entropy comparison (Goddard J/(kg*K), CEA kJ/(kg*K) - convert)
    if (cea_state.entropy_kj_kg_k > 0) {
        ComparisonResult s_result;
        s_result.parameter = "entropy";
        s_result.goddard_value = goddard_state.entropy / 1000.0;  // J/(kg*K) -> kJ/(kg*K)
        s_result.cea_value = cea_state.entropy_kj_kg_k;
        s_result.relative_error = relative_error(s_result.goddard_value, s_result.cea_value);
        s_result.tolerance_used = tolerance.entropy_rel;
        s_result.passed = within_relative_tolerance(
            s_result.goddard_value, s_result.cea_value, tolerance.entropy_rel);

        std::ostringstream msg;
        msg << "Entropy: Goddard=" << s_result.goddard_value
            << "kJ/(kg*K), CEA=" << s_result.cea_value
            << "kJ/(kg*K), RelErr=" << s_result.relative_error;
        s_result.message = msg.str();

        results.push_back(s_result);
    }

    // Gamma (specific heat ratio) comparison
    if (cea_state.gamma > 0) {
        ComparisonResult gamma_result;
        gamma_result.parameter = "gamma";
        gamma_result.goddard_value = goddard_state.gamma_s;
        gamma_result.cea_value = cea_state.gamma;
        gamma_result.relative_error = relative_error(gamma_result.goddard_value, gamma_result.cea_value);
        gamma_result.tolerance_used = tolerance.gamma_rel;
        gamma_result.passed = within_relative_tolerance(
            gamma_result.goddard_value, gamma_result.cea_value, tolerance.gamma_rel);

        std::ostringstream msg;
        msg << "Gamma: Goddard=" << gamma_result.goddard_value
            << ", CEA=" << gamma_result.cea_value
            << ", RelErr=" << gamma_result.relative_error;
        gamma_result.message = msg.str();

        results.push_back(gamma_result);
    }

    // Speed of sound comparison (both in m/s)
    if (cea_state.sound_speed_m_s > 0) {
        ComparisonResult sound_result;
        sound_result.parameter = "sound_speed";
        sound_result.goddard_value = goddard_state.speed_of_sound;
        sound_result.cea_value = cea_state.sound_speed_m_s;
        sound_result.relative_error = relative_error(sound_result.goddard_value, sound_result.cea_value);
        sound_result.tolerance_used = tolerance.sound_speed_rel;
        sound_result.passed = within_relative_tolerance(
            sound_result.goddard_value, sound_result.cea_value, tolerance.sound_speed_rel);

        std::ostringstream msg;
        msg << "Sound Speed: Goddard=" << sound_result.goddard_value
            << "m/s, CEA=" << sound_result.cea_value
            << "m/s, RelErr=" << sound_result.relative_error;
        sound_result.message = msg.str();

        results.push_back(sound_result);
    }

    // dlV/dlP_T comparison (thermodynamic derivative)
    if (std::abs(cea_state.dlv_dlp_t) > 1e-10) {
        ComparisonResult dlvp_result;
        dlvp_result.parameter = "dlV_dlP_T";
        dlvp_result.goddard_value = goddard_state.dlV_dlP_T;
        dlvp_result.cea_value = cea_state.dlv_dlp_t;
        dlvp_result.relative_error = relative_error(dlvp_result.goddard_value, dlvp_result.cea_value);
        dlvp_result.tolerance_used = tolerance.thermo_deriv_rel;
        dlvp_result.passed = within_relative_tolerance(
            dlvp_result.goddard_value, dlvp_result.cea_value, tolerance.thermo_deriv_rel);

        std::ostringstream msg;
        msg << "dlV/dlP_T: Goddard=" << dlvp_result.goddard_value
            << ", CEA=" << dlvp_result.cea_value
            << ", RelErr=" << dlvp_result.relative_error;
        dlvp_result.message = msg.str();

        results.push_back(dlvp_result);
    }

    // dlV/dlT_P comparison (thermodynamic derivative)
    if (std::abs(cea_state.dlv_dlt_p) > 1e-10) {
        ComparisonResult dlvt_result;
        dlvt_result.parameter = "dlV_dlT_P";
        dlvt_result.goddard_value = goddard_state.dlV_dlT_P;
        dlvt_result.cea_value = cea_state.dlv_dlt_p;
        dlvt_result.relative_error = relative_error(dlvt_result.goddard_value, dlvt_result.cea_value);
        dlvt_result.tolerance_used = tolerance.thermo_deriv_rel;
        dlvt_result.passed = within_relative_tolerance(
            dlvt_result.goddard_value, dlvt_result.cea_value, tolerance.thermo_deriv_rel);

        std::ostringstream msg;
        msg << "dlV/dlT_P: Goddard=" << dlvt_result.goddard_value
            << ", CEA=" << dlvt_result.cea_value
            << ", RelErr=" << dlvt_result.relative_error;
        dlvt_result.message = msg.str();

        results.push_back(dlvt_result);
    }

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

    size_t passed = 0;
    size_t total = results.size();

    for (const auto& result : results) {
        std::string status = result.passed ? "[PASS]" : "[FAIL]";
        std::cout << status << " " << result.message << std::endl;
        if (result.passed) passed++;
    }

    std::cout << "Summary: " << passed << "/" << total << " tests passed" << std::endl;
}
