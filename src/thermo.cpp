#include "goddard/thermo.hpp"
#include "goddard/error.hpp"
namespace Goddard {

size_t ThermodynamicState::state_size() const {
    return 2 + composition.size();
}

std::vector<double> ThermodynamicState::to_mole_vector(Cantera::ThermoPhase& thermo) const {
    std::vector<double> out(state_size());

    thermo.setState_TPX(temperature, pressure, composition);
    thermo.saveState(out);

    return out;
}

std::vector<double> ThermodynamicState::to_vector() const {
    std::vector<double> out(state_size());

    out[0] = temperature;
    out[1] = density;
    size_t i = 2;
    for (const auto& species: composition) {
        if (i >= out.size())
            throw std::runtime_error("Size of output vector was misspecified.");
        out[i] = species.second;
        i++;
    }
    return out;
}

std::vector<double> InputState::to_vector(Cantera::ThermoPhase& thermo) const {
    std::vector<double> out(thermo.stateSize());

    thermo.setState_TPX(T, P, composition);
    thermo.saveState(out);

    return out;
}

std::vector<double> InputState::to_mass_vector(Cantera::ThermoPhase& thermo) const {
    std::vector<double> out(thermo.stateSize());

    thermo.setState_TPY(T, P, composition);
    thermo.saveState(out);

    return out;
}


} //namespace goddard