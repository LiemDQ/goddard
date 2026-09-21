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

struct ChemicalParameters {
    std::string thermo_file;
    std::unordered_set<std::string> species;
    /**
     * Fuel stream state. `composition` is always a **mole-fraction** map. Its keys are product
     * species names when `reactant_file` is empty, and species of `reactant_file` otherwise
     * (e.g. `H2(L)`, `RP-1`). CEA-style weight percentages must be converted to mole fractions
     * by the caller.
     */
    PhaseSpecification cantera_fuel_state;
    /** Oxidizer stream state; see `cantera_fuel_state` for how `composition` is interpreted. */
    PhaseSpecification cantera_oxidizer_state;
    MixtureRatioType mixture_type;
    std::vector<double> mixtures;
    std::vector<double> OF_ratios;
    std::vector<double> phi_ratios;
    std::vector<double> fuel_weight_percentages;
    /**
     * Cantera YAML file holding the reactant species, e.g. `data/nasa9_reactants.yaml`.
     *
     * When set, the fuel and oxidizer are built as separate `Gas` streams from this file and
     * their element amounts and enthalpies are transferred to the products, so reactants need not
     * be product species. When empty, the reactant compositions refer to product species and the
     * combustor blends them in the product phase, as before.
     */
    std::string reactant_file;
    /** Cantera YAML file holding candidate condensed product species, e.g. `data/nasa9_condensed.yaml`. */
    std::string condensed_file;
    /** Condensed species of `condensed_file` to offer as candidates. Ignored if `all_condensed_species`. */
    std::unordered_set<std::string> condensed_species;
    /** Offer every species of `condensed_file` whose elements the product phase has. */
    bool all_condensed_species = false;
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
    /**
     * Product `Gas` for one solver stage. Every returned `Gas` references the same `Solution` and
     * the same candidate condensed species set, so attaching candidates once in the constructor is
     * enough for the combustor and the nozzle to see them. The chemistry mode is left at its
     * default: the combustor does not read it, and the nozzle and the results set their own.
     */
    Gas product_gas() const;

    /** Reactant stream `Gas` built from `ChemicalParameters::reactant_file` and set to `state`. */
    Gas reactant_gas(const PhaseSpecification& state) const;

    std::shared_ptr<Cantera::Solution> m_sln;
    /**
     * Product gas carrying the candidate condensed species, built once so that every `Gas` handed
     * to a solver shares the same set. Empty when no `condensed_file` was given.
     */
    std::optional<Gas> m_condensed_prototype;

};



} //namespace Goddard