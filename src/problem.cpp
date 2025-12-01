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
    m_problem_cases(cases), m_chemical_params(chem_params), 
    include_transport(transport), include_ionized_species(ionized_species), 
    m_trace_cutoff(trace) {

    auto root_node = select_species(m_chemical_params.thermo_file, m_chemical_params.species);
    const Cantera::AnyMap& phase_node = root_node.at("phases").getMapWhere("name", name);
    m_solution = Cantera::newSolution(phase_node, root_node);
    m_solution->setSource(m_chemical_params.thermo_file);
}

RocketProblemResults RocketProblem::solve() {
    std::unordered_map<std::string, RocketProblemCaseResult> case_results;
    
    Eigen::ArrayXd OFs = vector_to_eigenarray(m_chemical_params.OF_ratios);
    std::shared_ptr<Cantera::ThermoPhase> thermo = m_solution->thermo();
    std::vector<double> state(m_solution->thermo()->stateSize());

    for (RocketCaseParameters& params : m_problem_cases) {
        Eigen::ArrayXd pressures = vector_to_eigenarray(params.combustor_options.pressures);
        
        double M_fuel = molar_mass_from_composition(*thermo, m_chemical_params.cantera_oxidizer_state);
        double M_oxidizer = molar_mass_from_composition(*thermo, m_chemical_params.cantera_oxidizer_state);
        MixtureRatios MRs(OFs, M_fuel, M_oxidizer);
        Combustor combustor(m_solution, m_chemical_params.cantera_fuel_state, m_chemical_params.cantera_oxidizer_state);
        ThermoArray combustion_states = combustor.solve(pressures, MRs, params.combustor_options);

        std::vector<NozzleResults> expansion_results(static_cast<std::size_t>(combustion_states.size()));

        NozzleBase* nozzle;

        switch (params.nozzle_options.chemistry) {
            case NozzleChemistryType::FROZEN: {
                nozzle = new FrozenNozzle(*m_solution);
                break;
            }
            case NozzleChemistryType::EQUILIBRIUM: {
                nozzle = new EquilibriumNozzle(*m_solution);
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

        case_results[params.name] = {
            params.problem_type, 
            combustion_states, 
            params.nozzle_options.chemistry, 
            expansion_results
        };
    }
    
    return {shared_from_this(), std::move(case_results)};
}



std::vector<ThermoStateInfo> RocketProblemResults::extract_thermo_info(const std::string& case_name, std::size_t index){
    RocketProblemCaseResult& case_result = m_cases[case_name];
    auto tmo = thermo();
    std::vector<double> state(tmo->stateSize());
    
    
    if (index > case_result.inlet_states.size()) {
        throw std::runtime_error("Provided index exceeds length of ThermoArray.");
    }
    
    std::vector<double> state = case_result.inlet_states.get_state(index);
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

} //namespace Goddard