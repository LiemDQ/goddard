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

RocketProblem::RocketProblem(const ChemicalParameters& chem_params,
        const std::vector<RocketCaseParameters>& cases,
        const std::string& phase_name,
        bool transport,
        bool ionized_species,
        double trace):
    include_transport(transport),
    include_ionized_species(ionized_species),
    trace_cutoff(trace),
    problem_cases(cases),
    chemical_params(chem_params) {

    auto root_node = select_species(chemical_params.thermo_file, chemical_params.species);
    const Cantera::AnyMap& phase_node = root_node.at("phases").getMapWhere("name", phase_name);
    m_sln = Cantera::newSolution(phase_node, root_node);
    m_sln->setSource(chemical_params.thermo_file);

    if (!chemical_params.condensed_file.empty()) {
        Gas prototype(m_sln);
        if (chemical_params.all_condensed_species) {
            prototype.add_all_condensed_species(chemical_params.condensed_file);
        } else {
            std::vector<std::string> names(chemical_params.condensed_species.begin(),
                chemical_params.condensed_species.end());
            std::sort(names.begin(), names.end());
            prototype.add_condensed_species(chemical_params.condensed_file, names);
        }
        m_condensed_prototype = std::move(prototype);
    }
}

Gas RocketProblem::product_gas() const {
    if (m_condensed_prototype) {
        // A copy shares the `Solution` and the candidate condensed species set of the prototype.
        return *m_condensed_prototype;
    }
    return Gas(m_sln);
}

Gas RocketProblem::reactant_gas(const PhaseSpecification& state) const {
    std::vector<std::string> species;
    species.reserve(state.composition.size());
    for (const auto& entry : state.composition) {
        species.push_back(entry.first);
    }
    Gas stream = Gas::create_from_species(chemical_params.reactant_file, "reactants", species);
    stream.set_state_TPX(state.T, state.P, state.composition);
    return stream;
}

RocketProblemResults RocketProblem::solve() {
    if (include_transport) {
        throw NotImplementedError("RocketProblem: transport properties are not implemented.");
    }
    if (include_ionized_species) {
        throw NotImplementedError("RocketProblem: ionized species are not implemented.");
    }
    if (trace_cutoff != 0.0) {
        throw NotImplementedError("RocketProblem: a trace cutoff is not implemented.");
    }
    if (!chemical_params.mixtures.empty() || !chemical_params.phi_ratios.empty()
            || !chemical_params.fuel_weight_percentages.empty()) {
        throw NotImplementedError(
            "RocketProblem: mixtures, phi_ratios and fuel_weight_percentages are not implemented; "
            "set OF_ratios and CombustorOptions.mixture_type instead.");
    }

    std::unordered_map<std::string, RocketProblemCaseResult> case_results;
    
    Eigen::ArrayXd OFs = vector_to_eigenarray(chemical_params.OF_ratios);
    std::shared_ptr<Cantera::ThermoPhase> thermo = m_sln->thermo();
    std::vector<double> state(m_sln->thermo()->stateSize());

    const auto& fuel_input = chemical_params.fuel_state;
    const auto& ox_input = chemical_params.oxidizer_state;
    const bool reactant_streams = !chemical_params.reactant_file.empty();

    Cantera::Composition fuel_comp;
    Cantera::Composition ox_comp;
    if (!reactant_streams) {
        // Parse fuel and oxidizer compositions into Cantera Composition maps
        thermo->setState_TPX(fuel_input.T, fuel_input.P, fuel_input.composition);
        fuel_comp = thermo->getMoleFractionsByName();

        thermo->setState_TPX(ox_input.T, ox_input.P, ox_input.composition);
        ox_comp = thermo->getMoleFractionsByName();
    }
    double fuel_T = fuel_input.T;
    double ox_T = ox_input.T;

    for (RocketCaseParameters& params : problem_cases) {
        Eigen::ArrayXd pressures = vector_to_eigenarray(params.combustor_options.pressures);

        Gas combustor_gas = product_gas();

        ThermoArray combustion_states = reactant_streams
            ? Combustor(combustor_gas, reactant_gas(fuel_input), reactant_gas(ox_input))
                  .solve(pressures, OFs, params.combustor_options)
            : Combustor(combustor_gas, fuel_comp, ox_comp)
                  .solve(fuel_T, ox_T, pressures, OFs, params.combustor_options);

        std::vector<NozzleResults> expansion_results;
        expansion_results.reserve(static_cast<std::size_t>(combustion_states.size()));
        std::vector<FiniteAreaChamber> finite_area_chambers;

        Gas gas = product_gas();
        Nozzle nozzle(gas, params.nozzle_options);

        for (int i = 0; i < combustion_states.size(); i++) {
            state = combustion_states.get_state(i);

            switch (params.combustor_options.type) {
                case CombustorType::INFINITE_AREA: {
                    nozzle.set_inlet_state(state);
                    NozzleResults expansions = nozzle.solve(
                        params.nozzle_options.expansion_type, params.nozzle_options.expansion_ratios);
                    expansion_results.push_back(std::move(expansions));
                    break;
                }
                case CombustorType::FINITE_MASS_FLUX:
                case CombustorType::FINITE_CONTRACTION_RATIO: {
                    const double value =
                        params.combustor_options.type == CombustorType::FINITE_CONTRACTION_RATIO
                            ? params.combustor_options.contraction_ratio
                            : params.combustor_options.mass_flux;
                    FiniteAreaChamber fac = nozzle.solve_finite_area_chamber(
                        state, params.combustor_options.type, value);
                    std::vector<NozzleStation> stations = nozzle.solve_stations(fac,
                        params.nozzle_options.expansion_type, params.nozzle_options.expansion_ratios);

                    expansion_results.push_back(
                        NozzleResults{fac.injector, fac.throat, std::move(stations)});
                    finite_area_chambers.push_back(std::move(fac));
                    break;
                }
            }
        }

        case_results.emplace(params.name, RocketProblemCaseResult{
            std::move(combustion_states),
            params.nozzle_options.chemistry,
            expansion_results,
            chemical_params.OF_ratios,
            params.combustor_options.pressures,
            params.nozzle_options.expansion_ratios,
            params.nozzle_options.expansion_type,
            params.combustor_options.process,
            params.combustor_options.type,
            std::move(finite_area_chambers),
            params.nozzle_options.frozen_NFZ,
            params.combustor_options.mixture_type
        });
    }
    
    return {std::move(case_results), product_gas()};
}


} //namespace Goddard