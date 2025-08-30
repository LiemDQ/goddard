#include "cea_loader.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <sstream>


// CEAResult convenience methods
const CEAState* CEAResult::find_chamber_state() const {
    for (const auto& state : equilibrium_states) {
        if (state.location == "CHAMBER") {
            return &state;
        }
    }
    return nullptr;
}

const CEAState* CEAResult::find_throat_state() const {
    for (const auto& state : equilibrium_states) {
        if (state.location == "THROAT") {
            return &state;
        }
    }
    return nullptr;
}

std::vector<const CEAState*> CEAResult::find_exit_states() const {
    std::vector<const CEAState*> exits;
    for (const auto& state : equilibrium_states) {
        if (state.location.find("EXIT") != std::string::npos) {
            exits.push_back(&state);
        }
    }
    return exits;
}

// CEADataLoader implementation
std::unique_ptr<CEAResult> CEADataLoader::load_from_file(const std::string& json_path) const {
    try {
        std::ifstream file(json_path);
        if (!file.is_open()) {
            std::cerr << "Failed to open CEA JSON file: " << json_path << std::endl;
            return nullptr;
        }
        
        nlohmann::json json;
        file >> json;
        
        auto result = std::make_unique<CEAResult>();
        
        // Parse conditions
        if (json.contains("conditions")) {
            result->conditions = parse_conditions(json["conditions"]);
        }
        
        // Parse equilibrium states
        if (json.contains("equilibrium_states")) {
            for (const auto& state_json : json["equilibrium_states"]) {
                result->equilibrium_states.push_back(parse_state(state_json));
            }
        }
        
        // Parse frozen states
        if (json.contains("frozen_states")) {
            for (const auto& state_json : json["frozen_states"]) {
                result->frozen_states.push_back(parse_state(state_json));
            }
        }
        
        // Parse performance (if needed)
        if (json.contains("equilibrium_performance")) {
            for (const auto& perf_json : json["equilibrium_performance"]) {
                result->equilibrium_performance.push_back(parse_performance(perf_json));
            }
        }
        
        if (json.contains("frozen_performance")) {
            for (const auto& perf_json : json["frozen_performance"]) {
                result->frozen_performance.push_back(parse_performance(perf_json));
            }
        }
        
        return result;
        
    } catch (const std::exception& e) {
        std::cerr << "Error parsing CEA JSON file " << json_path << ": " << e.what() << std::endl;
        return nullptr;
    }
}

std::unordered_map<std::string, std::unique_ptr<CEAResult>> CEADataLoader::load_from_directory(
    const std::string& directory_path) const {
    
    std::unordered_map<std::string, std::unique_ptr<CEAResult>> results;
    
    try {
        for (const auto& entry : std::filesystem::directory_iterator(directory_path)) {
            if (entry.path().extension() == ".json") {
                std::string filename = entry.path().stem();
                auto cea_result = load_from_file(entry.path().string());
                if (cea_result) {
                    results[filename] = std::move(cea_result);
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error loading CEA files from " << directory_path << ": " << e.what() << std::endl;
    }
    
    return results;
}

CEAConditions CEADataLoader::parse_conditions(const nlohmann::json& json) const {
    CEAConditions conditions;
    
    conditions.pressure_psia = json.value("pressure_psia", 0.0);
    conditions.of_ratio = json.value("of_ratio", 0.0);
    conditions.fuel = json.value("fuel", "");
    conditions.oxidizer = json.value("oxidizer", "");
    conditions.fuel_temp = json.value("fuel_temp", 0.0);
    conditions.fuel_energy = json.value("fuel_energy", 0.0);
    conditions.oxidizer_temp = json.value("oxidizer_temp", 0.0);
    conditions.oxidizer_energy = json.value("oxidizer_energy", 0.0);
    conditions.case_name = json.value("case_name", "");
    conditions.problem_type = json.value("problem_type", "");
    
    
    return conditions;
}

CEAState CEADataLoader::parse_state(const nlohmann::json& json) const {
    CEAState state;
    
    state.location = json.value("location", "");
    state.pressure_ratio = json.value("pressure_ratio", 0.0);
    state.pressure_bar = json.value("pressure_bar", 0.0);
    state.temperature_k = json.value("temperature_k", 0.0);
    state.density_kg_m3 = json.value("density_kg_m3", 0.0);
    state.enthalpy_kj_kg = json.value("enthalpy_kj_kg", 0.0);
    state.internal_energy_kj_kg = json.value("internal_energy_kj_kg", 0.0);
    state.gibbs_kj_kg = json.value("gibbs_kj_kg", 0.0);
    state.entropy_kj_kg_k = json.value("entropy_kj_kg_k", 0.0);
    state.molecular_weight = json.value("molecular_weight", 0.0);
    state.dlv_dlp_t = json.value("dlv_dlp_t", 0.0);
    state.dlv_dlt_p = json.value("dlv_dlt_p", 0.0);
    state.cp_kj_kg_k = json.value("cp_kj_kg_k", 0.0);
    state.gamma = json.value("gamma", 0.0);
    state.sound_speed_m_s = json.value("sound_speed_m_s", 0.0);
    state.mach_number = json.value("mach_number", 0.0);
    
    if (json.contains("mass_fractions")) {
        state.mass_fractions = parse_mass_fractions(json["mass_fractions"]);
    }
    
    return state;
}

CEAPerformance CEADataLoader::parse_performance(const nlohmann::json& json) const {
    CEAPerformance perf;
    
    perf.cstar_m_s = json.value("cstar_m_s", 0.0);
    perf.cf = json.value("cf", 0.0);
    perf.isp_m_s = json.value("isp_m_s", 0.0);
    perf.ivac_m_s = json.value("ivac_m_s", 0.0);
    
    return perf;
}

CEAComposition CEADataLoader::parse_mass_fractions(const nlohmann::json& json) const {
    std::unordered_map<std::string, double> mass_fractions;
    
    for (auto it = json.begin(); it != json.end(); ++it) {
        mass_fractions[it.key()] = it.value().get<double>();
    }
    
    return mass_fractions;
}
