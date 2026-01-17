#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <memory>
#include <optional>
#include <nlohmann/json.hpp>

#include "goddard/problem.hpp"

/**
 * CEA test data structures for validation against NASA CEA results
 */

/**
 * @brief Input conditions for a CEA problem.
 */
struct CEAConditions {
    double pressure_psia;
    double of_ratio;
    std::string fuel;
    std::string oxidizer;
    double fuel_temp;
    double fuel_energy;
    double oxidizer_temp;
    double oxidizer_energy;
    std::vector<double> area_ratios;
    std::string case_name;
    std::string problem_type;
};


using CEAComposition = std::unordered_map<std::string, double>;

/**
 * @brief Thermodynamic state at a specific point
 */
struct CEAState {
    std::string location; //chamber, throat, exit
    double pressure_ratio;
    double pressure_bar;
    double temperature_k;
    double density_kg_m3;
    double enthalpy_kj_kg;
    double internal_energy_kj_kg;
    double gibbs_kj_kg;
    double entropy_kj_kg_k;
    double molecular_weight;
    double dlv_dlp_t;
    double dlv_dlt_p;
    double cp_kj_kg_k;
    double gamma;
    double sound_speed_m_s;
    double mach_number;
    CEAComposition mass_fractions;
};

/**
 * @brief Rocket performance parameters
 */
struct CEAPerformance {
    double area_ratio;
    double cstar_m_s;
    double cf;
    double isp_m_s;
    double ivac_m_s;
};

struct CEAResult {
    CEAConditions conditions;
    std::vector<CEAState> equilibrium_states;
    std::vector<CEAState> frozen_states;
    std::vector<CEAPerformance> equilibrium_performance;
    std::vector<CEAPerformance> frozen_performance;
    
    // Convenience methods to find specific states
    const CEAState* find_chamber_state() const;
    const CEAState* find_throat_state() const;
    std::vector<const CEAState*> find_exit_states() const;
};


/**
 * 
 * CEA Data Loader - reads JSON files produced by cea_parser.py
 * 
 * Usage:
 * 
 *   CEADataLoader loader;
 *   auto cea_result = loader.load_from_file("data/cea_results/h2.json");
 *   if (cea_result) {
 *       // Use CEA data for testing
 *       auto chamber = cea_result->find_chamber_state();
 *       if (chamber) {
 *           // Compare against Goddard results
 *       }
 *   }
 */
class CEADataLoader {
public:
    /**
     * Load CEA data from JSON file
     * @param json_path Path to JSON file created by cea_parser.py
     * @return CEAResult or nullptr if failed to load
     */
    std::unique_ptr<CEAResult> load_from_file(const std::string& json_path) const;
    
    /**
     * Load all CEA files from a directory
     * @param directory_path Path to directory containing *.json files
     * @return Map of filename -> CEAResult
     */
    std::unordered_map<std::string, std::unique_ptr<CEAResult>> load_from_directory(
        const std::string& directory_path) const;

private:
    CEAConditions parse_conditions(const nlohmann::json& json) const;
    CEAState parse_state(const nlohmann::json& json) const;
    CEAPerformance parse_performance(const nlohmann::json& json) const;
    CEAComposition parse_mass_fractions(const nlohmann::json& json) const;
};

