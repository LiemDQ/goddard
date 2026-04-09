#include "goddard/gas.hpp"
#include "goddard/gas_dynamics.hpp"

namespace Goddard {

Gas::Gas(std::shared_ptr<Cantera::Solution> gas, GasChemistry chem)
    : chemistry(chem), m_sol(gas) 
{
    set_current_state_as_reference();
}

Gas::Gas(const Cantera::Solution& gas, GasChemistry chem)
    : Gas(gas.clone(), chem) 
{
    set_current_state_as_reference();
}

Gas::Gas(Cantera::Solution&& gas, GasChemistry chem)
    : chemistry(chem), m_sol(gas.shared_from_this()) {}

Gas::Gas(const Gas& gas) 
    : chemistry(gas.chemistry), m_sol(gas.m_sol->clone()), 
      m_H_stagnation(gas.m_H_stagnation), m_S0(gas.m_S0) {}

Gas Gas::operator=(const Gas& gas) {
    return Gas(gas);
}

Gas Gas::create(const std::string& filename, GasChemistry chemistry, const std::string& phase_name) {
    return Gas(Cantera::newSolution(filename, phase_name), chemistry);
}
// State setters

void Gas::set_state_TD(double T, double D) {
    m_sol->thermo()->setState_TD(T,D);
}

void Gas::set_state_TP(double T, double P) {
    m_sol->thermo()->setState_TP(T, P);
}

void Gas::set_state_TPX(double T, double P, const std::string& composition) {
    m_sol->thermo()->setState_TPX(T, P, composition);
}

void Gas::set_state_TPX(double T, double P, const Composition& composition) {
    m_sol->thermo()->setState_TPX(T, P, composition);
}

void Gas::set_state_TPY(double T, double P, const std::string& composition) {
    m_sol->thermo()->setState_TPY(T, P, composition);
}

void Gas::set_state_TPY(double T, double P, const Composition& composition) {
    m_sol->thermo()->setState_TPY(T, P, composition);
}

void Gas::set_state_HP(double H, double P) {
    m_sol->thermo()->setState_HP(H, P);
}

void Gas::set_state_SP(double S, double P) {
    m_sol->thermo()->setState_SP(S, P);
}

void Gas::set_state_UV(double U, double V) {
    m_sol->thermo()->setState_UV(U,V);
}

// 

std::vector<double> Gas::save_state() const {
    std::vector<double> out(thermo()->stateSize());
    m_sol->thermo()->saveState(out);
    return out;
}

void Gas::copy_state(std::vector<double>& state) const {
    state.resize(m_sol->thermo()->stateSize());
    m_sol->thermo()->saveState(state);
}

void Gas::restore_state(const std::vector<double>& state) {
    m_sol->thermo()->restoreState(state);
}

ThermodynamicState Gas::snapshot() const {
    auto t = m_sol->thermo();

    ThermodynamicState info;
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

// Basic thermodynamic properties

double Gas::temperature() const { return m_sol->thermo()->temperature(); }
double Gas::pressure() const { return m_sol->thermo()->pressure(); }
double Gas::density() const { return m_sol->thermo()->density(); }
double Gas::enthalpy_mass() const { return m_sol->thermo()->enthalpy_mass(); }
double Gas::entropy_mass() const { return m_sol->thermo()->entropy_mass(); }
double Gas::cp_mass() const { return m_sol->thermo()->cp_mass(); }
double Gas::cv_mass() const { return m_sol->thermo()->cv_mass(); }
double Gas::molecular_weight() const { return m_sol->thermo()->meanMolecularWeight(); }
size_t Gas::num_species() const { return m_sol->thermo()->nSpecies(); }
std::vector<std::string> Gas::species_names() const { return m_sol->thermo()->speciesNames(); }
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

// Mixture properties
double Gas::fuel_fraction(
    const std::string& fuel_comp, 
    const std::string& ox_comp,
    Cantera::ThermoBasis basis) const 
{
    return thermo()->mixtureFraction(fuel_comp, ox_comp, basis);
}
double Gas::fuel_fraction(
    const Composition& fuel_comp, 
    const Composition& ox_comp,
    Cantera::ThermoBasis basis) const 
{
    return thermo()->mixtureFraction(fuel_comp, ox_comp, basis);
}

double Gas::equivalence_ratio() const {
    return thermo()->equivalenceRatio();
}

double Gas::equivalence_ratio(
    const std::string& fuel_comp, 
    const std::string& ox_comp,
    Cantera::ThermoBasis basis) const 
{
    return thermo()->equivalenceRatio(fuel_comp, ox_comp, basis);
}
double Gas::equivalence_ratio(
    const Composition& fuel_comp, 
    const Composition& ox_comp,
    Cantera::ThermoBasis basis) const 
{
    return thermo()->equivalenceRatio(fuel_comp, ox_comp, basis);
}

double Gas::stoich_OF_ratio(
    const std::string& fuel_comp, 
    const std::string& ox_comp,
    Cantera::ThermoBasis basis) const 
{
    return thermo()->stoichAirFuelRatio(fuel_comp, ox_comp, basis);
}
double Gas::stoich_OF_ratio(
    const Composition& fuel_comp, 
    const Composition& ox_comp,
    Cantera::ThermoBasis basis) const 
{
    return thermo()->stoichAirFuelRatio(fuel_comp, ox_comp, basis);
}

void Gas::set_fuel_fraction(
    double fuel_frac,
    const std::string& fuel_comp, 
    const std::string& ox_comp,
    Cantera::ThermoBasis basis)
{
    thermo()->setMixtureFraction(fuel_frac, fuel_comp, ox_comp, basis);
    set_current_state_as_reference();
}
void Gas::set_fuel_fraction(
    double fuel_frac,
    const Composition& fuel_comp, 
    const Composition& ox_comp,
    Cantera::ThermoBasis basis)
{
    thermo()->setMixtureFraction(fuel_frac, fuel_comp, ox_comp, basis);
    set_current_state_as_reference();
}

void Gas::set_equivalence_ratio(
    double phi,
    const std::string& fuel_comp, 
    const std::string& ox_comp,
    Cantera::ThermoBasis basis)
{
    thermo()->setEquivalenceRatio(phi, fuel_comp, ox_comp, basis);
    set_current_state_as_reference();
}
void Gas::set_equivalence_ratio(
    double phi,
    const Composition& fuel_comp, 
    const Composition& ox_comp,
    Cantera::ThermoBasis basis)
{
    thermo()->setEquivalenceRatio(phi, fuel_comp, ox_comp, basis);
    set_current_state_as_reference();
}

void Gas::set_OF_ratio(
    double OF,
    const std::string& fuel_comp, 
    const std::string& ox_comp,
    Cantera::ThermoBasis basis)
{
    double frac = 1.0/(OF+1.0);
    thermo()->setMixtureFraction(frac, fuel_comp, ox_comp, basis);
    set_current_state_as_reference();
}
void Gas::set_OF_ratio(
    double OF,
    const Composition& fuel_comp, 
    const Composition& ox_comp,
    Cantera::ThermoBasis basis)
{
    double frac = 1.0/(OF+1.0);
    thermo()->setMixtureFraction(frac, fuel_comp, ox_comp, basis);
    set_current_state_as_reference();
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


// Reference state

void Gas::set_stagnation_enthalpy(double H) { m_H_stagnation = H; }
double Gas::get_stagnation_enthalpy() const { return m_H_stagnation; }
void Gas::set_reference_entropy(double S) { m_S0 = S; }
double Gas::get_reference_entropy() const { return m_S0; }

void Gas::set_current_state_as_reference() {
    m_H_stagnation = m_sol->thermo()->enthalpy_mass();
    m_S0 = m_sol->thermo()->entropy_mass();
}

// Accessors

std::shared_ptr<Cantera::Solution> Gas::solution() const { return m_sol; }
std::shared_ptr<Cantera::ThermoPhase> Gas::thermo() const { return m_sol->thermo(); }
std::shared_ptr<Cantera::Kinetics> Gas::kinetics() const { return m_sol->kinetics(); }
std::shared_ptr<Cantera::Transport> Gas::transport() const { return m_sol->transport(); }

std::string Gas::report(bool show_thermo, double threshold) const {
    return m_sol->thermo()->report(show_thermo, threshold);
}

} // namespace Goddard
