#include "goddard/thermo.hpp"

namespace Goddard {

// double molar_mass_from_composition(Cantera::ThermoPhase& thermo, const std::vector<double>& state) {
//     std::vector<double> old_state(thermo.stateSize());

//     thermo.saveState(old_state);
//     thermo.restoreState(state);
//     double molar_mass = thermo.meanMolecularWeight();
//     thermo.restoreState(old_state);

//     return molar_mass;
// }

std::vector<double> ThermodynamicState::to_vector(Cantera::ThermoPhase& thermo){
    std::vector<double> out(thermo.stateSize());

    thermo.setState_TPX(T, P, composition);
    thermo.saveState(out);\

    return out;
}

std::vector<double> ThermodynamicState::to_mass_vector(Cantera::ThermoPhase& thermo){
    std::vector<double> out(thermo.stateSize());

    thermo.setState_TPY(T, P, composition);
    thermo.saveState(out);
    
    return out;
}

} //namespace goddard