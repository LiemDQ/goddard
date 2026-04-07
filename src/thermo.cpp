#include "goddard/thermo.hpp"

namespace Goddard {

std::vector<double> ThermodynamicState::to_vector(Cantera::ThermoPhase& thermo){
    std::vector<double> out(thermo.stateSize());

    thermo.setState_TPX(T, P, composition);
    thermo.saveState(out);

    return out;
}

std::vector<double> ThermodynamicState::to_mass_vector(Cantera::ThermoPhase& thermo){
    std::vector<double> out(thermo.stateSize());

    thermo.setState_TPY(T, P, composition);
    thermo.saveState(out);
    
    return out;
}

} //namespace goddard