#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <memory>
#include <optional>
#include "cantera/core.h"
#include "goddard/combustor.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/thermo.hpp"
#include "goddard/rocket_results.hpp"

namespace Goddard {

struct RocketCaseParameters {
    std::string name;
    std::string problem_type;
    CombustorOptions combustor_options;
    NozzleOptions nozzle_options;
};

struct Speciation { //TODO: implement speciation functionality
    std::unordered_set<std::string> species;
    std::unordered_set<std::string> elements;
    int max_species = 10;
};

struct ChemicalParameters {
    std::string thermo_file;
    std::unordered_set<std::string> species;
    ThermodynamicState cantera_fuel_state;
    ThermodynamicState cantera_oxidizer_state;
    std::vector<double> OF_ratios;
    std::vector<double> phi_ratios;
    std::vector<double> fuel_weight_percentages;
    std::vector<double> valance_equivalences;
};

struct RocketState {
    std::string name;
    std::vector<double> cantera_state;
    double pressure_ratio;
    double area_ratio;
    double dlv_dlp_t;
    double dlv_dlt_p;
    double gamma_s;
    double speed_of_sound;
};



// Forward declaration
class RocketProblemResults;

// TODO: implement this function
// RocketPerformance performance_from_state();


class RocketProblem {

    public:
    
    RocketProblem(const ChemicalParameters& chem_params,
        const std::vector<RocketCaseParameters>& cases, 
        const std::string& name = "", 
        bool transport = false,
        bool ionized_species = false,
        double trace = 1e-6);
    
    RocketProblemResults solve();
    
    inline std::shared_ptr<Cantera::Solution> solution() { return m_sln; }
    inline std::shared_ptr<Cantera::ThermoPhase> thermo() { return m_sln->thermo(); }
    bool include_transport = false; //TODO: implement transport functionality
    bool include_ionized_species = false; //TODO: implement ionized species
    double trace_cutoff = 1e-6;
    std::vector<RocketCaseParameters> problem_cases;
    ChemicalParameters chemical_params;
    
    private:
    std::shared_ptr<Cantera::Solution> m_sln;

};



} //namespace Goddard