#include "goddard/thermo.hpp"

double molar_mass_from_composition(Cantera::ThermoPhase& thermo, const std::vector<double>& state) {
    std::vector<double> old_state(thermo.stateSize());

    thermo.saveState(old_state);
    thermo.restoreState(state);
    double molar_mass = thermo.meanMolecularWeight();
    thermo.restoreState(old_state);

    return molar_mass;
}