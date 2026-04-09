#include "goddard/problem.hpp"
#include "goddard/combustor.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/gas.hpp"
#include "goddard/gas_dynamics.hpp"
#include "goddard/error.hpp"
#include "goddard/thermo.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/thermoarray.hpp"
#include "goddard/utils.hpp"
#include "goddard/speciate.hpp"
#include "goddard/format.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <exception>
#include <memory>
#include <set>
#include <sstream>
#include <utility>
#include <vector>
#include <unordered_set>
#include "eigen3/Eigen/Dense"
#include "cantera/core.h"

namespace Goddard {

std::shared_ptr<RocketProblem> create(const ChemicalParameters& chem_params,
        std::vector<RocketCaseParameters>& cases,
        const std::string& name,
        bool include_transport,
        bool include_ionized_species,
        double trace_cutoff) {
    
    
    return std::make_shared<RocketProblem>(
        chem_params,
        cases,
        name,
        include_transport,
        include_ionized_species,
        trace_cutoff
    );
} 

RocketProblem::RocketProblem(const ChemicalParameters& chem_params,
        const std::vector<RocketCaseParameters>& cases,
        const std::string& name,
        bool transport,
        bool ionized_species,
        double trace):
    include_transport(transport),
    include_ionized_species(ionized_species),
    trace_cutoff(trace),
    problem_cases(cases),
    chemical_params(chem_params) {

    auto root_node = select_species(chemical_params.thermo_file, chemical_params.species);
    const Cantera::AnyMap& phase_node = root_node.at("phases").getMapWhere("name", name);
    m_sln = Cantera::newSolution(phase_node, root_node);
    m_sln->setSource(chemical_params.thermo_file);
}

RocketProblemResults RocketProblem::solve() {
    std::unordered_map<std::string, RocketProblemCaseResult> case_results;
    
    Eigen::ArrayXd OFs = vector_to_eigenarray(chemical_params.OF_ratios);
    std::shared_ptr<Cantera::ThermoPhase> thermo = m_sln->thermo();
    std::vector<double> state(m_sln->thermo()->stateSize());

    const auto& fuel_input = chemical_params.cantera_fuel_state;
    const auto& ox_input = chemical_params.cantera_oxidizer_state;

    // Parse fuel and oxidizer compositions into Cantera Composition maps
    thermo->setState_TPX(fuel_input.T, fuel_input.P, fuel_input.composition);
    Cantera::Composition fuel_comp = thermo->getMoleFractionsByName();
    double fuel_T = fuel_input.T;

    thermo->setState_TPX(ox_input.T, ox_input.P, ox_input.composition);
    Cantera::Composition ox_comp = thermo->getMoleFractionsByName();
    double ox_T = ox_input.T;

    for (RocketCaseParameters& params : problem_cases) {
        Eigen::ArrayXd pressures = vector_to_eigenarray(params.combustor_options.pressures);

        Gas combustor_gas(m_sln);
        Combustor combustor(combustor_gas, fuel_comp, ox_comp);

        // TODO: workaround for bug in Cantera SolutionArray
        // prevents incorrect values from being written in RocketProblemResults.
        // this will have to remain in place until Cantera merges a fix.
        // see https://github.com/Cantera/cantera/issues/2067
        if (pressures.size() == 1 && OFs.size() == 1) {
            pressures.conservativeResize(2);
            pressures(1) = pressures(0);
        }
        ThermoArray combustion_states = combustor.solve(fuel_T, ox_T, pressures, OFs, params.combustor_options);

        std::vector<NozzleResults> expansion_results;
        expansion_results.reserve(static_cast<std::size_t>(combustion_states.size()));

        Gas gas(m_sln, params.nozzle_options.chemistry);
        Nozzle nozzle(gas, params.nozzle_options.chemistry);

        for (int i = 0; i < combustion_states.size(); i++) {
            state = combustion_states.get_state(i);
            nozzle.set_inlet_state(state);
            NozzleResults expansions = nozzle.solve(params.nozzle_options.expansion_type, params.nozzle_options.expansion_ratios);
            expansion_results.push_back(expansions);
        }

        case_results.emplace(params.name, RocketProblemCaseResult{
            params.problem_type,
            std::move(combustion_states),
            params.nozzle_options.chemistry,
            expansion_results,
            chemical_params.OF_ratios,
            params.combustor_options.pressures,
            params.nozzle_options.expansion_ratios,
            params.nozzle_options.expansion_type
        });
    }
    
    return {std::move(case_results), m_sln};
}


} //namespace Goddard