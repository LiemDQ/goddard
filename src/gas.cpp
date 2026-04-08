#include "goddard/gas.hpp"
#include "goddard/gas_dynamics.hpp"

namespace Goddard {

Gas::Gas(std::shared_ptr<Cantera::Solution> gas, GasChemistry chem)
    : chemistry(chem), m_sol(gas->clone()) {
    m_H_stagnation = m_sol->thermo()->enthalpy_mass();
    m_S0 = m_sol->thermo()->entropy_mass();
}

Gas::Gas(Cantera::Solution& gas, GasChemistry chem)
    : Gas(gas.clone(), chem) {}

Gas::Gas(Cantera::Solution&& gas, GasChemistry chem)
    : chemistry(chem), m_sol(gas.shared_from_this()) {}

Gas::Gas(const Gas& gas) 
    : chemistry(gas.chemistry), m_sol(gas.m_sol->clone()), 
      m_H_stagnation(gas.m_H_stagnation), m_S0(gas.m_S0) {}

Gas Gas::operator=(const Gas& gas) {
    return Gas(gas);
}
// State setters

void Gas::set_state_TP(double T, double P) {
    m_sol->thermo()->setState_TP(T, P);
}

void Gas::set_state_TPX(double T, double P, const std::string& composition) {
    m_sol->thermo()->setState_TPX(T, P, composition);
}

void Gas::set_state_HP(double H, double P) {
    m_sol->thermo()->setState_HP(H, P);
}

void Gas::set_state_SP(double S, double P) {
    m_sol->thermo()->setState_SP(S, P);
}

void Gas::save_state(std::vector<double>& state) const {
    state.resize(m_sol->thermo()->stateSize());
    m_sol->thermo()->saveState(state);
}

void Gas::restore_state(const std::vector<double>& state) {
    m_sol->thermo()->restoreState(state);
}

// Basic thermodynamic properties

double Gas::temperature() const { return m_sol->thermo()->temperature(); }
double Gas::pressure() const { return m_sol->thermo()->pressure(); }
double Gas::density() const { return m_sol->thermo()->density(); }
double Gas::enthalpy_mass() const { return m_sol->thermo()->enthalpy_mass(); }
double Gas::entropy_mass() const { return m_sol->thermo()->entropy_mass(); }
double Gas::cp_mass() const { return m_sol->thermo()->cp_mass(); }
double Gas::cv_mass() const { return m_sol->thermo()->cv_mass(); }
double Gas::molecular_weight() const { return m_sol->thermo()->meanMolecularWeight(); }

// Chemistry-aware derived properties

double Gas::gamma_s() const {
    switch (chemistry) {
        case GasChemistry::PERFECT_GAS:
        case GasChemistry::FROZEN:
        case GasChemistry::KINETIC:
            return m_sol->thermo()->cp_mass() / m_sol->thermo()->cv_mass();
        case GasChemistry::EQUILIBRIUM:
            return get_thermo_equilibrium_properties(*m_sol->thermo()).gamma_s;
    }
    return m_sol->thermo()->cp_mass() / m_sol->thermo()->cv_mass(); // unreachable
}

double Gas::speed_of_sound() const {
    return gas_sonic_velocity(*m_sol->thermo(), gamma_s());
}

double Gas::stagnation_enthalpy(double velocity) const {
    return gas_stagnation_enthalpy(*m_sol->thermo(), velocity);
}

double Gas::stagnation_pressure(double velocity) const {
    return gas_stagnation_pressure(*m_sol->thermo(), velocity);
}

double Gas::isenthalpic_velocity(double H_stagnation) const {
    return gas_isenthalpic_velocity(*m_sol->thermo(), H_stagnation);
}

double Gas::isenthalpic_velocity() const {
    return gas_isenthalpic_velocity(*m_sol->thermo(), m_H_stagnation);
}

double Gas::mach(double velocity) const {
    return velocity / speed_of_sound();
}

// Expansion properties

ExpansionProperties Gas::expansion_properties() const {
    switch (chemistry) {
        case GasChemistry::EQUILIBRIUM:
            return get_thermo_equilibrium_properties(*m_sol->thermo());
        case GasChemistry::PERFECT_GAS:
        case GasChemistry::FROZEN:
        case GasChemistry::KINETIC: {
            double cp = m_sol->thermo()->cp_mass();
            double g = cp / m_sol->thermo()->cv_mass();
            return {1.0, -1.0, cp, g};
        }
    }
    return {1.0, -1.0, m_sol->thermo()->cp_mass(),
            m_sol->thermo()->cp_mass() / m_sol->thermo()->cv_mass()}; // unreachable
}

void Gas::equilibrate(const std::string& XY) {
    m_sol->thermo()->equilibrate(XY, "gibbs");
}

// Snapshot

ThermoStateInfo Gas::snapshot() const {
    auto t = m_sol->thermo();

    ThermoStateInfo info;
    info.pressure = t->pressure();
    info.temperature = t->temperature();
    info.density = t->density();
    info.enthalpy = t->enthalpy_mass();
    info.internal_energy = t->intEnergy_mass();
    info.gibbs = t->gibbs_mass();
    info.entropy = t->entropy_mass();
    info.molecular_weight = t->meanMolecularWeight();
    info.cp = t->cp_mass();
    info.gamma_s = gamma_s();
    info.stagnation_enthalpy = m_H_stagnation;

    auto props = expansion_properties();
    info.dlV_dlP_T = props.dlogV_dlogP_T;
    info.dlV_dlT_P = props.dlogV_dlogT_P;
    info.speed_of_sound = gas_sonic_velocity(*t, info.gamma_s);

    // Composition
    std::vector<double> mole_fractions(t->nSpecies());
    t->getMoleFractions(mole_fractions.data());
    for (size_t i = 0; i < t->nSpecies(); ++i) {
        if (mole_fractions[i] > 0.0) {
            info.composition[t->speciesName(i)] = mole_fractions[i];
        }
    }

    return info;
}

// Reference state

void Gas::set_stagnation_enthalpy(double H) { m_H_stagnation = H; }
double Gas::get_stagnation_enthalpy() const { return m_H_stagnation; }
void Gas::set_reference_entropy(double S) { m_S0 = S; }
double Gas::get_reference_entropy() const { return m_S0; }

// Accessors

std::shared_ptr<Cantera::Solution> Gas::solution() const { return m_sol; }
std::shared_ptr<Cantera::ThermoPhase> Gas::thermo() const { return m_sol->thermo(); }
std::shared_ptr<Cantera::Kinetics> Gas::kinetics() const { return m_sol->kinetics(); }
std::shared_ptr<Cantera::Transport> Gas::transport() const { return m_sol->transport(); }

} // namespace Goddard
