#include "goddard/problem.hpp"
#include "goddard/error.hpp"
#include <exception>

namespace Goddard {

    StateInfo RocketProblemResults::get_state_info(const std::string& case_name, const std::string& state_name){
        RocketProblemCaseResult& caseresult = cases[case_name];
        auto thermo = params.solution->thermo();
        RocketState state;

        bool found = false;
        if (caseresult.inlet_state.name == state_name){
            state = caseresult.inlet_state;
            found = true;
        }
        else {
            for (RocketState& s : caseresult.nozzle_states) {
                if (s.name == state_name) {
                    state = s;
                    found = true;
                    break;
                }
            }
        }
        if (!found) {
            throw std::runtime_error("State name not found in case.");
        }

        std::unordered_map<std::string, double> compositions; 
        std::vector<std::string> names = thermo->speciesNames();

        for (std::size_t i = 0; i < thermo->nSpecies(); i++) {
            compositions[names[i]] = state.ct_state[i+2];
        }

        return StateInfo{
            thermo->pressure(),
            thermo->temperature(),
            thermo->density(),
            thermo->enthalpy_mass(),
            thermo->intEnergy_mass(),
            thermo->gibbs_mass(),
            thermo->entropy_mass(),
            thermo->meanMolecularWeight(),
            state.dlv_dlp_t,
            state.dlv_dlt_p,
            state.gamma_s,
            state.speed_of_sound,
            compositions
        };
    }

} //namespace Goddard