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

/** Combustor and nozzle settings of one case of a `RocketProblem`. */
struct RocketCaseParameters {
    /// Case name, the key of the case in `RocketProblemResults`.
    std::string name;
    /// Chamber model, pressures [Pa] and combustion process.
    CombustorOptions combustor_options;
    /// Nozzle chemistry, freezing station and exit stations.
    NozzleOptions nozzle_options;
};

/** Thermodynamic data, reactants and mixture ratios of a `RocketProblem`. */
struct ChemicalParameters {
    /// Cantera YAML file holding the gas-phase product species, e.g. `nasa9_gas.yaml`.
    std::string thermo_file;
    /**
     * Gas-phase product species of `thermo_file` to include, names matched exactly. Every product
     * species must be listed: species not named are dropped from the phase.
     *
     * In Python, assign a whole set: adding to the attribute acts on a copy.
     */
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
     *
     * In Python, assign a whole list: item assignment and append act on a copy.
     */
    std::vector<double> OF_ratios;
    /**
     * Not implemented: `RocketProblem::solve` raises `NotImplementedError` if non-empty.
     *
     * In Python, assign a whole list: item assignment and append act on a copy.
     */
    std::vector<double> mixtures;
    /**
     * Not implemented: `RocketProblem::solve` raises `NotImplementedError` if non-empty.
     *
     * In Python, assign a whole list: item assignment and append act on a copy.
     */
    std::vector<double> phi_ratios;
    /**
     * Not implemented: `RocketProblem::solve` raises `NotImplementedError` if non-empty.
     *
     * In Python, assign a whole list: item assignment and append act on a copy.
     */
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
    /**
     * Condensed species of `condensed_file` to offer as candidates. Ignored if
     * `all_condensed_species`.
     *
     * In Python, assign a whole set: adding to the attribute acts on a copy.
     */
    std::unordered_set<std::string> condensed_species;
    /** Offer every species of `condensed_file` whose elements the product phase has. */
    bool all_condensed_species = false;
};

// Forward declaration
class RocketProblemResults;

// TODO: implement this function
// RocketPerformance performance_from_state();


/**
 * CEA-style rocket problem: combustion followed by an isentropic nozzle expansion, for every
 * case, mixture ratio and chamber pressure.
 */
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
    
    /**
     * Solve every case.
     *
     * @return Stations and performance of every case, mixture ratio and pressure.
     * @throws NotImplementedError if an unimplemented option is set (`include_transport`,
     * `include_ionized_species`, `trace_cutoff`, `mixtures`, `phi_ratios`,
     * `fuel_weight_percentages`).
     */
    RocketProblemResults solve();
    
    inline std::shared_ptr<Cantera::Solution> solution() { return m_sln; }
    inline std::shared_ptr<Cantera::ThermoPhase> thermo() { return m_sln->thermo(); }
    /** Not implemented: `solve` raises `NotImplementedError` if true. */
    bool include_transport = false;
    /** Not implemented: `solve` raises `NotImplementedError` if true. */
    bool include_ionized_species = false;
    /** Not implemented: `solve` raises `NotImplementedError` if non-zero. */
    double trace_cutoff = 0.0;
    /**
     * Cases to solve.
     *
     * In Python, assign a whole list: item assignment and append act on a copy.
     */
    std::vector<RocketCaseParameters> problem_cases;
    /**
     * Thermodynamic data, reactants and mixture ratios. The product phase is built from it by
     * the constructor; later changes affect only the reactants and mixture ratios.
     */
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