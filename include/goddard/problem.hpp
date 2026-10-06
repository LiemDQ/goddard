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
    PhaseSpecification fuel_state;
    /** Oxidizer stream state; see `fuel_state` for how `composition` is interpreted. */
    PhaseSpecification oxidizer_state;
    /**
     * Mixture ratios of the problem, interpreted through `CombustorOptions::mixture_type` of each
     * case (O/F ratios by default).
     */
    std::vector<double> OF_ratios;
    /** Not implemented: `RocketProblem::solve` raises `NotImplementedError` if non-empty. */
    std::vector<double> mixtures;
    /** Not implemented: `RocketProblem::solve` raises `NotImplementedError` if non-empty. */
    std::vector<double> phi_ratios;
    /** Not implemented: `RocketProblem::solve` raises `NotImplementedError` if non-empty. */
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

// Forward declaration
class RocketProblemResults;

// TODO: implement this function
// RocketPerformance performance_from_state();


class RocketProblem {

    public:
    
    /**
     * @param chem_params Thermodynamic data, reactants and mixture ratios.
     * @param cases Combustor and nozzle settings of each case.
     * @param phase_name Phase of `chem_params.thermo_file` to build the products from; empty
     * selects the first phase in the file.
     * @param transport Not implemented; `solve` raises NotImplementedError if true.
     * @param ionized_species Not implemented; `solve` raises NotImplementedError if true.
     * @param trace Not implemented; `solve` raises NotImplementedError if non-zero.
     */
    RocketProblem(const ChemicalParameters& chem_params,
        const std::vector<RocketCaseParameters>& cases, 
        const std::string& phase_name = "", 
        bool transport = false,
        bool ionized_species = false,
        double trace = 0.0);
    
    RocketProblemResults solve();
    
    inline std::shared_ptr<Cantera::Solution> solution() { return m_sln; }
    inline std::shared_ptr<Cantera::ThermoPhase> thermo() { return m_sln->thermo(); }
    /** Not implemented: `solve` raises `NotImplementedError` if true. */
    bool include_transport = false;
    /** Not implemented: `solve` raises `NotImplementedError` if true. */
    bool include_ionized_species = false;
    /** Not implemented: `solve` raises `NotImplementedError` if non-zero. */
    double trace_cutoff = 0.0;
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