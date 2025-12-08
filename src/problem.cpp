#include "goddard/problem.hpp"
#include "goddard/combustor.hpp"
#include "goddard/nozzle.hpp"
#include "goddard/nozzle_utils.hpp"
#include "goddard/mixture_ratio.hpp"
#include "goddard/error.hpp"
#include "goddard/thermo.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/thermoarray.hpp"
#include "goddard/utils.hpp"
#include "goddard/speciate.hpp"
#include <cmath>
#include <iostream>
#include <exception>
#include <memory>
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

    for (RocketCaseParameters& params : problem_cases) {
        Eigen::ArrayXd pressures = vector_to_eigenarray(params.combustor_options.pressures);
        
        double M_fuel = molar_mass_from_composition(*thermo, chemical_params.cantera_fuel_state);
        double M_oxidizer = molar_mass_from_composition(*thermo, chemical_params.cantera_oxidizer_state);
        MixtureRatios MRs(OFs, M_fuel, M_oxidizer);
        Combustor combustor(m_sln, chemical_params.cantera_fuel_state, chemical_params.cantera_oxidizer_state);
        ThermoArray combustion_states = combustor.solve(pressures, MRs, params.combustor_options);

        std::vector<NozzleResults> expansion_results;
        expansion_results.reserve(static_cast<std::size_t>(combustion_states.size()));

        std::unique_ptr<NozzleBase> nozzle;

        switch (params.nozzle_options.chemistry) {
            case NozzleChemistryType::FROZEN: {
                nozzle = std::make_unique<FrozenNozzle>(*m_sln);
                break;
            }
            case NozzleChemistryType::EQUILIBRIUM: {
                nozzle = std::make_unique<EquilibriumNozzle>(*m_sln);
                break;
            }
            default: throw NotImplementedError("Nozzle type is not implemented.");
            //TODO: separate execution path for "NONE" nozzle chemistry
        }

        for (int i = 0; i < combustion_states.size(); i++) {
            state = combustion_states.get_state(i);
            nozzle->set_inlet_state(state);
            NozzleResults expansions = nozzle->solve(params.nozzle_options.expansion_type, params.nozzle_options.expansion_ratios);
            expansion_results.push_back(expansions);
        }

        case_results.emplace(params.name, RocketProblemCaseResult{
            params.problem_type,
            std::move(combustion_states),
            params.nozzle_options.chemistry,
            expansion_results
        });
    }
    
    return {std::move(case_results), m_sln};
}

RocketProblemResults::RocketProblemResults(
    std::unordered_map<std::string, RocketProblemCaseResult>&& case_results, 
    std::shared_ptr<Cantera::Solution> sln) 
    : cases(case_results), m_sln(std::move(sln)) {
}

std::vector<ThermoStateInfo> RocketProblemResults::extract_thermo_info(const std::string& case_name, std::size_t index){
    RocketProblemCaseResult& case_result = cases.at(case_name);
    auto tmo = thermo();

    if (index >= static_cast<std::size_t>(case_result.inlet_states.size())) {
        throw std::runtime_error("Provided index exceeds length of ThermoArray.");
    }

    std::vector<double> state = case_result.inlet_states.get_state(static_cast<int>(index));
    NozzleResults nozzle_results = case_result.nozzle_states[index];
    

    std::unordered_map<std::string, double> inlet_compositions; 
    std::vector<std::string> names = tmo->speciesNames();
    

    for (std::size_t i = 0; i < tmo->nSpecies(); i++) {
        inlet_compositions[names[i]] = state[i+2];
    }

    tmo->restoreState(state);
    double inlet_gamma = 0;
    double inlet_dlogV_dlogP_T = 0;
    double inlet_dlogV_dlogT_P = 0;

    switch (case_result.chemistry) {
        case NozzleChemistryType::EQUILIBRIUM: {
            auto props = get_thermo_equilibrium_properties(*tmo);
            inlet_gamma = props.gamma_s;
            inlet_dlogV_dlogP_T = props.dlogV_dlogP_T;
            inlet_dlogV_dlogT_P = props.dlogV_dlogT_P;
            break;
        }
        case NozzleChemistryType::FROZEN: {
            inlet_gamma = tmo->cp_mass()/tmo->cv_mass();
            inlet_dlogV_dlogP_T = -1.0; //by definition -- see NASA reference publication 1311 eq 6.38
            inlet_dlogV_dlogT_P = 1.0; //by definition -- see NASA reference publication 1311 eq 6.37
            break;
        }
        default: throw NotImplementedError("This chemistry type has not been implemented.");
    }
        
    ThermoStateInfo inlet_state{
        tmo->pressure(),
        tmo->temperature(),
        tmo->density(),
        tmo->enthalpy_mass(),
        tmo->intEnergy_mass(),
        tmo->gibbs_mass(),
        tmo->entropy_mass(),
        tmo->meanMolecularWeight(),
        inlet_gamma,
        inlet_dlogV_dlogP_T,
        inlet_dlogV_dlogT_P,
        gas_sonic_velocity(*tmo, inlet_gamma),
        inlet_compositions
    };

    //throat
    ThroatCondition throat = nozzle_results.throat;
    std::unordered_map<std::string, double> throat_compositions; 
    tmo->restoreState(throat.state);

    for (std::size_t i = 0; i < tmo->nSpecies(); i++) {
        throat_compositions[names[i]] = state[i+2];
    }

    
    ThermoStateInfo throat_state {
        tmo->pressure(),
        tmo->temperature(),
        tmo->density(),
        tmo->enthalpy_mass(),
        tmo->intEnergy_mass(),
        tmo->gibbs_mass(),
        tmo->entropy_mass(),
        tmo->meanMolecularWeight(),
        throat.gamma_s,
        throat.dlV_dlP_T,
        throat.dlV_dlT_P,
        gas_sonic_velocity(*tmo, throat.gamma_s),
        throat_compositions
    };

    std::vector<ThermoStateInfo> state_info {inlet_state, throat_state}; 

    for (auto&& exp: nozzle_results.expansions) {
        tmo->restoreState(exp.state);
        std::unordered_map<std::string, double> compositions;
        for (std::size_t i = 0; i < tmo->nSpecies(); i++) {
           compositions[names[i]] = state[i+2];
        }

        ThermoStateInfo exp_state {
            tmo->pressure(),
            tmo->temperature(),
            tmo->density(),
            tmo->enthalpy_mass(),
            tmo->intEnergy_mass(),
            tmo->gibbs_mass(),
            tmo->entropy_mass(),
            tmo->meanMolecularWeight(),
            exp.gamma_s,
            exp.dlV_dlP_T,
            exp.dlV_dlT_P,
            gas_sonic_velocity(*tmo, throat.gamma_s),
            compositions
        };

        state_info.push_back(exp_state);

    }

    return state_info;
}

std::optional<ThermoStateInfo> RocketProblemResults::get_chamber_state(
    const std::string& case_name, std::size_t of_index) {

    auto it = cases.find(case_name);
    if (it == cases.end()) {
        return std::nullopt;
    }

    RocketProblemCaseResult& case_result = it->second;
    if (of_index >= static_cast<std::size_t>(case_result.inlet_states.size())) {
        return std::nullopt;
    }

    auto all_states = extract_thermo_info(case_name, of_index);
    if (all_states.empty()) {
        return std::nullopt;
    }

    return all_states[0]; // Chamber is first state
}

std::optional<ThermoStateInfo> RocketProblemResults::get_throat_state(
    const std::string& case_name, std::size_t of_index) {

    auto it = cases.find(case_name);
    if (it == cases.end()) {
        return std::nullopt;
    }

    RocketProblemCaseResult& case_result = it->second;
    if (of_index >= static_cast<std::size_t>(case_result.inlet_states.size())) {
        return std::nullopt;
    }

    auto all_states = extract_thermo_info(case_name, of_index);
    if (all_states.size() < 2) {
        return std::nullopt;
    }

    return all_states[1]; // Throat is second state
}

std::vector<ThermoStateInfo> RocketProblemResults::get_exit_states(
    const std::string& case_name, std::size_t of_index) {

    auto it = cases.find(case_name);
    if (it == cases.end()) {
        return {};
    }

    RocketProblemCaseResult& case_result = it->second;
    if (of_index >= static_cast<std::size_t>(case_result.inlet_states.size())) {
        return {};
    }

    auto all_states = extract_thermo_info(case_name, of_index);
    if (all_states.size() <= 2) {
        return {}; // No exit states
    }

    // Exit states start at index 2
    return std::vector<ThermoStateInfo>(all_states.begin() + 2, all_states.end());
}

RocketPerformance RocketProblemResults::calculate_performance(
    const ThermoStateInfo& chamber,
    const ThermoStateInfo& throat,
    const ThermoStateInfo& exit) {

    // Pressure ratios
    double pressure_ratio = chamber.pressure / exit.pressure;
    double throat_pressure_ratio = chamber.pressure / throat.pressure;

    // Area ratio from isentropic flow relations
    // A/A* = (1/M) * [(2/(gamma+1)) * (1 + (gamma-1)/2 * M^2)]^((gamma+1)/(2*(gamma-1)))
    // For now, estimate from density ratio (approximate)
    double area_ratio = (throat.density * throat.speed_of_sound) /
                        (exit.density * exit.speed_of_sound) *
                        (throat_pressure_ratio / pressure_ratio);

    // Characteristic velocity (c*)
    // c* = P_c * A_t / m_dot = sqrt(gamma * R * T_c) / gamma * sqrt((2/(gamma+1))^((gamma+1)/(gamma-1)))
    // Simplified: c* = throat.speed_of_sound / sqrt(throat.gamma_s) * factor
    double gamma = throat.gamma_s;
    double cstar = throat.speed_of_sound * std::sqrt(
        std::pow(2.0 / (gamma + 1.0), (gamma + 1.0) / (gamma - 1.0)) / gamma
    );

    // Exit velocity from enthalpy difference
    // v_e = sqrt(2 * (h_chamber - h_exit))
    double exit_velocity = std::sqrt(2.0 * (chamber.enthalpy - exit.enthalpy));

    // Thrust coefficient CF = v_e / c* + (P_e - P_amb) * A_e / (P_c * A_t)
    // For vacuum: CF_vac = v_e / c* + P_e * A_e / (P_c * A_t)
    // Simplified (assuming matched nozzle): CF ≈ v_e / c*
    double CF = exit_velocity / cstar;

    // Specific impulse Isp = v_e / g0
    constexpr double g0 = 9.80665; // m/s^2
    double isp = exit_velocity / g0;

    // Vacuum specific impulse (includes pressure thrust term)
    // Ivac = Isp + P_e * A_e / (m_dot * g0)
    // Approximation: Ivac ≈ Isp * (1 + P_e/(P_c) * area_ratio * some_factor)
    double ivac = isp + (exit.pressure / chamber.pressure) * area_ratio * cstar / g0;

    // Mach number at exit (approximate from speed of sound)
    double mach_number = exit_velocity / exit.speed_of_sound;

    return RocketPerformance{
        pressure_ratio,
        area_ratio,
        mach_number,
        cstar,
        CF,
        isp,
        ivac
    };
}

} //namespace Goddard