#include "goddard/utils.hpp"
#include "goddard/speciate.hpp"
#include "goddard/condensed.hpp"
#include "goddard/error.hpp"
#include "goddard/gas.hpp"
#include "goddard/gas_dynamics.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace Goddard {

Gas::Gas(std::shared_ptr<Cantera::Solution> gas, GasChemistry chem)
    : chemistry(chem), m_sol(gas)
{
    set_current_state_as_reference();
}

Gas::Gas(std::shared_ptr<Cantera::Solution> gas,
         std::shared_ptr<CondensedPhaseSet> condensed,
         GasChemistry chem)
    : chemistry(chem), m_sol(std::move(gas)), m_condensed(std::move(condensed))
{
    set_current_state_as_reference();
}

Gas::Gas(const Cantera::Solution& gas, GasChemistry chem)
    : Gas(gas.clone(), chem) 
{
    set_current_state_as_reference();
}

Gas::Gas(Cantera::Solution&& gas, GasChemistry chem)
    : chemistry(chem), m_sol(gas.shared_from_this()) 
{
    set_current_state_as_reference();
}

Gas::Gas(const std::string& infile, 
        const std::string& phase_name,
        GasChemistry chem) 
    : chemistry(chem), m_sol(Cantera::newSolution(infile, phase_name)) 
{
    set_current_state_as_reference();
}

Gas::Gas(double gamma) : m_gamma(gamma)
{
    chemistry = GasChemistry::PERFECT_GAS;
    set_current_state_as_reference();
}

Gas Gas::clone() const {
    Gas copy = *this;
    if (m_condensed) {
        copy.m_condensed = m_condensed->clone();
    }
    if (m_sol) {
        copy.m_sol = m_sol->clone();
        copy.restore_state(save_state());
    }
    return copy;
}

Gas Gas::create(const std::string& filename, const std::string& phase_name, GasChemistry chemistry) {
    return Gas(filename, phase_name, chemistry);
}

namespace {

/**
 * Drop the kinetics and transport entries of a generated phase node when the data file cannot
 * support them: thermodynamic databases such as the NASA9 files have no `reactions` section and
 * no species transport data, and Cantera refuses to build the phase otherwise.
 */
void adapt_phase_node_to_data(Cantera::AnyMap& phase_node, const Cantera::AnyMap& root_node) {
    if (!root_node.hasKey("reactions")) {
        phase_node.erase("kinetics");
        phase_node.erase("reactions");
    }
    if (root_node.hasKey("species")) {
        for (const auto& species_node : root_node.at("species").asVector<Cantera::AnyMap>()) {
            if (!species_node.hasKey("transport")) {
                phase_node.erase("transport");
                break;
            }
        }
    }
}

} // namespace

Gas Gas::create_from_elements(
    const std::string& infile,
    const std::string& name,
    const std::vector<std::string>& elements,
    GasChemistry chemistry)
{
    Cantera::AnyMap root_node = load_root_node(infile);
    Cantera::AnyMap phase_node = create_speciated_phase_node(name, elements);
    adapt_phase_node_to_data(phase_node, root_node);
    return Gas(Cantera::newSolution(phase_node, root_node), chemistry);
}

Gas Gas::create_from_species(
    const std::string& infile,
    const std::string& name,
    const std::vector<std::string>& species,
    GasChemistry chemistry)
{
    Cantera::AnyMap root_node = load_root_node(infile);
    Cantera::AnyMap phase_node = create_phase_node(name, species);
    adapt_phase_node_to_data(phase_node, root_node);
    return Gas(Cantera::newSolution(phase_node, root_node), chemistry);
}

// State setters
void Gas::set_state_TD(double T, double D) {
    thermo()->setState_TD(T,D);
}

void Gas::set_state_TP(double T, double P) {
    thermo()->setState_TP(T, P);
}

void Gas::set_state_TPX(double T, double P, const std::string& composition) {
    thermo()->setState_TPX(T, P, composition);
}

void Gas::set_state_TPX(double T, double P, const Composition& composition) {
    thermo()->setState_TPX(T, P, composition);
}

void Gas::set_state_TPX(double T, double P, const double* composition) {
    thermo()->setState_TPX(T, P, composition);
}

void Gas::set_state_TPY(double T, double P, const std::string& composition) {
    thermo()->setState_TPY(T, P, composition);
}

void Gas::set_state_TPY(double T, double P, const Composition& composition) {
    thermo()->setState_TPY(T, P, composition);
}

void Gas::set_state_TPY(double T, double P, const double* composition) {
    thermo()->setState_TPY(T, P, composition);
}

void Gas::set_state_HP(double H, double P) {
    thermo()->setState_HP(H, P);
}

void Gas::set_state_SP(double S, double P) {
    thermo()->setState_SP(S, P);
}

void Gas::set_state_UV(double U, double V) {
    thermo()->setState_UV(U,V);
}

// 

std::vector<double> Gas::save_state() const {
    std::vector<double> out;
    copy_state(out);
    return out;
}

void Gas::copy_state(std::vector<double>& state) const {
    const size_t cantera_size = thermo()->stateSize();
    const size_t n_condensed = m_condensed ? m_condensed->size() : 0;

    state.resize(cantera_size + n_condensed);
    thermo()->saveState(cantera_size, state.data());
    for (size_t k = 0; k < n_condensed; k++) {
        state[cantera_size + k] = m_condensed->moles[k];
    }
}

void Gas::restore_state(const std::vector<double>& state) {
    const size_t cantera_size = thermo()->stateSize();
    const size_t n_condensed = m_condensed ? m_condensed->size() : 0;

    if (state.size() == cantera_size) {
        thermo()->restoreState(state);
        if (m_condensed) {
            std::fill(m_condensed->moles.begin(), m_condensed->moles.end(), 0.0);
            m_condensed->pinned_group = -1;
        }
        return;
    }

    if (state.size() != cantera_size + n_condensed) {
        throw std::invalid_argument(
            "Gas::restore_state: state vector has length " + std::to_string(state.size())
            + ", expected " + std::to_string(cantera_size) + " (gas only) or "
            + std::to_string(cantera_size + n_condensed) + " (gas plus condensed species).");
    }

    thermo()->restoreState(cantera_size, state.data());
    for (size_t k = 0; k < n_condensed; k++) {
        m_condensed->moles[k] = state[cantera_size + k];
    }
}

ThermodynamicState Gas::snapshot() const {
    auto t = thermo();

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

std::string Gas::name() const {
    if (has_cantera_sln())  return m_sol->name();
    else return "";
}

void Gas::set_name(const std::string& name) {
    if (has_cantera_sln()) m_sol->setName(name);
}

// Basic thermodynamic properties

double Gas::temperature() const { return thermo()->temperature(); }
double Gas::pressure() const { return thermo()->pressure(); }
double Gas::density() const { return thermo()->density(); }
double Gas::enthalpy_mass() const { return thermo()->enthalpy_mass(); }
double Gas::entropy_mass() const { return thermo()->entropy_mass(); }
double Gas::cp_mass() const { return thermo()->cp_mass(); }
double Gas::cv_mass() const { return thermo()->cv_mass(); }
double Gas::molecular_weight() const { return thermo()->meanMolecularWeight(); }
std::vector<double> Gas::mole_fractions() const { 
    std::vector<double> fracs(num_species());
    thermo()->getMoleFractions(fracs.data());
    return fracs; 
}
std::vector<double> Gas::mass_fractions() const { 
    std::vector<double> fracs(num_species());
    thermo()->getMassFractions(fracs.data());
    return fracs; 
}

size_t Gas::num_species() const { return thermo()->nSpecies(); }
std::vector<std::string> Gas::species_names() const { return thermo()->speciesNames(); }
// Chemistry-aware derived properties

double Gas::gamma_s() const {
    switch (chemistry) {
        case GasChemistry::PERFECT_GAS: {
            if (has_cantera_sln()) {
                return thermo()->cp_mass() / thermo()->cv_mass();
            }
            else return m_gamma;
        }
        case GasChemistry::FROZEN:
        case GasChemistry::KINETIC:
            return thermo()->cp_mass() / thermo()->cv_mass();
        case GasChemistry::EQUILIBRIUM:
            return get_thermo_equilibrium_properties(*thermo()).gamma_s;
    }
    return thermo()->cp_mass() / thermo()->cv_mass(); // unreachable
}

double Gas::speed_of_sound() const {
    return gas_sonic_velocity(*thermo(), gamma_s());
}

double Gas::stagnation_enthalpy(double velocity) const {
    return gas_stagnation_enthalpy(*thermo(), velocity);
}

double Gas::stagnation_pressure(double velocity) const {
    return gas_stagnation_pressure(*thermo(), velocity);
}

double Gas::isenthalpic_velocity(double H_stagnation) const {
    return gas_isenthalpic_velocity(*thermo(), H_stagnation);
}

double Gas::isenthalpic_velocity() const {
    return gas_isenthalpic_velocity(*thermo(), m_H_stagnation);
}

double Gas::area_per_mdot(double velocity) const {
    return Goddard::area_per_mdot(*thermo(), velocity);
}

double Gas::mach(double velocity) const {
    return velocity / speed_of_sound();
}

double Gas::cstar() const {
    return Goddard::cstar(gamma_s(), temperature(), molecular_weight());
}

double Gas::isp() const {
    return isenthalpic_velocity();
}

double Gas::ivac() const {
    double exit_velocity = isp();
    return exit_velocity + temperature()*Cantera::GasConstant/(exit_velocity * molecular_weight());
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
            return get_thermo_equilibrium_properties(*thermo());
        case GasChemistry::PERFECT_GAS: {
            if (!has_cantera_sln()) return {1.0, -1.0, 0.0, m_gamma};
            else {
                double cp = thermo()->cp_mass();
                double g = cp / thermo()->cv_mass();
                return {1.0, -1.0, cp, g};
            }
        }
        case GasChemistry::FROZEN:
        case GasChemistry::KINETIC: {
            double cp = thermo()->cp_mass();
            double g = cp / thermo()->cv_mass();
            return {1.0, -1.0, cp, g};
        }
    }
    return {1.0, -1.0, thermo()->cp_mass(),
            thermo()->cp_mass() / thermo()->cv_mass()}; // unreachable
}

void Gas::equilibrate(const std::string& XY, const std::string& solver) {
    thermo()->equilibrate(XY, solver);
}

namespace {

/** The multiphase solver is work package B; until it lands the condensed paths are stubs. */
void reject_condensed_candidates(const Gas& gas, const char* function) {
    if (gas.has_condensed_candidates()) {
        throw NotImplementedError(std::string(function)
            + ": equilibrium with condensed species is not implemented yet (WP-B).");
    }
}

/**
 * Index of the species representing element `m` on its own in a basis composition: the
 * homonuclear diatomic if there is one, else the monatomic, else any single-element species.
 */
size_t basis_species_for_element(const Cantera::ThermoPhase& gas, size_t m) {
    size_t diatomic = Cantera::npos;
    size_t monatomic = Cantera::npos;
    size_t any = Cantera::npos;

    for (size_t j = 0; j < gas.nSpecies(); j++) {
        const double atoms = gas.nAtoms(j, m);
        if (atoms <= 0.0) {
            continue;
        }
        bool single_element = true;
        for (size_t i = 0; i < gas.nElements(); i++) {
            if (i != m && gas.nAtoms(j, i) != 0.0) {
                single_element = false;
                break;
            }
        }
        if (!single_element) {
            continue;
        }
        if (atoms == 2.0 && diatomic == Cantera::npos) {
            diatomic = j;
        } else if (atoms == 1.0 && monatomic == Cantera::npos) {
            monatomic = j;
        }
        if (any == Cantera::npos) {
            any = j;
        }
    }

    if (diatomic != Cantera::npos) return diatomic;
    if (monatomic != Cantera::npos) return monatomic;
    return any;
}

} // namespace

void Gas::equilibrate_TP(double T, double P) {
    reject_condensed_candidates(*this, "Gas::equilibrate_TP");
    thermo()->setState_TP(T, P);
    thermo()->equilibrate("TP", "gibbs");
}

void Gas::equilibrate_HP(double H, double P) {
    reject_condensed_candidates(*this, "Gas::equilibrate_HP");
    thermo()->setState_HP(H, P);
    thermo()->equilibrate("HP", "gibbs");
}

void Gas::equilibrate_SP(double S, double P) {
    reject_condensed_candidates(*this, "Gas::equilibrate_SP");
    thermo()->setState_SP(S, P);
    thermo()->equilibrate("SP", "gibbs");
}

void Gas::set_element_moles(const Eigen::ArrayXd& element_moles, double T, double P) {
    reject_condensed_candidates(*this, "Gas::set_element_moles");

    auto gas = thermo();
    const size_t n_elements = gas->nElements();
    if (static_cast<size_t>(element_moles.size()) != n_elements) {
        throw std::invalid_argument(
            "Gas::set_element_moles: expected " + std::to_string(n_elements)
            + " element amounts, got " + std::to_string(element_moles.size()) + ".");
    }

    Eigen::ArrayXd moles = Eigen::ArrayXd::Zero(static_cast<long>(gas->nSpecies()));
    for (size_t m = 0; m < n_elements; m++) {
        if (element_moles(static_cast<long>(m)) <= 0.0) {
            continue;
        }
        const size_t j = basis_species_for_element(*gas, m);
        if (j == Cantera::npos) {
            throw FmtError(
                "Gas::set_element_moles: the phase contains no species made up of element '{}' "
                "alone, so no basis composition can be formed.",
                gas->elementName(m));
        }
        moles(static_cast<long>(j)) += element_moles(static_cast<long>(m)) / gas->nAtoms(j, m);
    }

    // setState_TPX normalizes, so only the composition ratios of `moles` survive.
    gas->setState_TPX(T, P, moles.data());

    if (m_condensed) {
        std::fill(m_condensed->moles.begin(), m_condensed->moles.end(), 0.0);
        m_condensed->pinned_group = -1;
    }
}

// Condensed species

void Gas::add_condensed_species(const std::string& infile, const std::vector<std::string>& names) {
    Cantera::AnyMap root_node = load_root_node(infile);
    if (m_condensed) {
        m_condensed->add_species(root_node, names, *thermo());
    } else {
        m_condensed = std::make_shared<CondensedPhaseSet>(root_node, names, *thermo());
    }
}

void Gas::add_all_condensed_species(const std::string& infile) {
    Cantera::AnyMap root_node = load_root_node(infile);
    std::vector<std::string> names = CondensedPhaseSet::compatible_species(root_node, *thermo());
    if (m_condensed) {
        m_condensed->add_species(root_node, names, *thermo());
    } else {
        m_condensed = std::make_shared<CondensedPhaseSet>(root_node, names, *thermo());
    }
}

bool Gas::has_condensed_candidates() const {
    return m_condensed != nullptr && m_condensed->size() > 0;
}

bool Gas::has_condensed_phases() const {
    if (!m_condensed) return false;
    for (double n : m_condensed->moles) {
        if (n > 0.0) return true;
    }
    return false;
}

std::vector<std::string> Gas::condensed_species_names() const {
    if (!m_condensed) return {};
    return m_condensed->names();
}

std::vector<double> Gas::condensed_moles() const {
    if (!m_condensed) return {};
    return m_condensed->moles;
}

void Gas::set_condensed_moles(const std::vector<double>& moles) {
    const size_t expected = m_condensed ? m_condensed->size() : 0;
    if (moles.size() != expected) {
        throw std::invalid_argument(
            "Gas::set_condensed_moles: expected " + std::to_string(expected)
            + " values, got " + std::to_string(moles.size()) + ".");
    }
    if (m_condensed) {
        m_condensed->moles = moles;
    }
}

double Gas::gas_mass_fraction() const {
    if (!m_condensed) return 1.0;
    double condensed_mass = 0.0;
    for (size_t k = 0; k < m_condensed->size(); k++) {
        condensed_mass += m_condensed->moles[k] * m_condensed->species[k].molar_mass;
    }
    return 1.0 - condensed_mass;
}

double Gas::mixture_molecular_weight() const {
    const double gas_moles = gas_mass_fraction() / thermo()->meanMolecularWeight();
    double total_moles = gas_moles;
    if (m_condensed) {
        for (double n : m_condensed->moles) {
            total_moles += n;
        }
    }
    return 1.0 / total_moles;
}

bool Gas::at_phase_transition() const {
    return m_condensed != nullptr && m_condensed->pinned_group >= 0;
}

std::pair<long, long> Gas::pinned_polymorphs() const {
    if (!at_phase_transition()) {
        return {-1L, -1L};
    }

    // Map candidate indices onto the compacted list of species that are actually present.
    std::vector<long> present_index(m_condensed->size(), -1);
    long n_present = 0;
    for (size_t k = 0; k < m_condensed->size(); k++) {
        if (m_condensed->moles[k] > 0.0) {
            present_index[k] = n_present++;
        }
    }

    const PolymorphGroup& group = m_condensed->groups.at(static_cast<size_t>(m_condensed->pinned_group));
    std::pair<long, long> pair{-1L, -1L};
    for (size_t k : group.members) { // members are ordered by increasing T_min
        if (present_index[k] < 0) continue;
        if (pair.first < 0) {
            pair.first = present_index[k];
        } else if (pair.second < 0) {
            pair.second = present_index[k];
        }
    }
    return pair;
}

std::vector<double> Gas::mixture_mass_fractions() const {
    const double w_gas = gas_mass_fraction();
    std::vector<double> fractions = mass_fractions();
    for (double& y : fractions) {
        y *= w_gas;
    }
    if (m_condensed) {
        for (size_t k = 0; k < m_condensed->size(); k++) {
            fractions.push_back(m_condensed->moles[k] * m_condensed->species[k].molar_mass);
        }
    }
    return fractions;
}

namespace {

/** Candidate indices of the condensed species that are currently present. */
std::vector<size_t> present_condensed(const CondensedPhaseSet& set) {
    std::vector<size_t> present;
    for (size_t k = 0; k < set.size(); k++) {
        if (set.moles[k] > 0.0) present.push_back(k);
    }
    return present;
}

} // namespace

Eigen::ArrayXd Gas::condensed_enthalpy_RT() const {
    if (!m_condensed) return {};
    const std::vector<size_t> present = present_condensed(*m_condensed);
    const double T = thermo()->temperature();

    Eigen::ArrayXd values(static_cast<long>(present.size()));
    for (size_t i = 0; i < present.size(); i++) {
        values(static_cast<long>(i)) = m_condensed->enthalpy_RT(present[i], T);
    }
    return values;
}

Eigen::ArrayXd Gas::condensed_cp_R() const {
    if (!m_condensed) return {};
    const std::vector<size_t> present = present_condensed(*m_condensed);
    const double T = thermo()->temperature();

    Eigen::ArrayXd values(static_cast<long>(present.size()));
    for (size_t i = 0; i < present.size(); i++) {
        values(static_cast<long>(i)) = m_condensed->cp_R(present[i], T);
    }
    return values;
}

Eigen::ArrayXd Gas::condensed_molar_masses() const {
    if (!m_condensed) return {};
    const std::vector<size_t> present = present_condensed(*m_condensed);

    Eigen::ArrayXd values(static_cast<long>(present.size()));
    for (size_t i = 0; i < present.size(); i++) {
        values(static_cast<long>(i)) = m_condensed->species[present[i]].molar_mass;
    }
    return values;
}

Eigen::ArrayXXd Gas::condensed_stoich_coeffs() const {
    const long n_elements = static_cast<long>(thermo()->nElements());
    if (!m_condensed) return Eigen::ArrayXXd(0, n_elements);

    const std::vector<size_t> present = present_condensed(*m_condensed);
    Eigen::ArrayXXd coeffs(static_cast<long>(present.size()), n_elements);
    for (size_t i = 0; i < present.size(); i++) {
        coeffs.row(static_cast<long>(i)) = m_condensed->species[present[i]].element_atoms.transpose();
    }
    return coeffs;
}

std::vector<std::string> Gas::element_names() const {
    return thermo()->elementNames();
}

Eigen::ArrayXd Gas::element_moles() const {
    auto gas = thermo();
    const long n_elements = static_cast<long>(gas->nElements());
    const size_t n_species = gas->nSpecies();

    std::vector<double> mole_fracs(n_species);
    gas->getMoleFractions(mole_fracs.data());
    // Moles of gas per kg of mixture: w_g / M_gas.
    const double gas_moles = gas_mass_fraction() / gas->meanMolecularWeight();

    Eigen::ArrayXd amounts = Eigen::ArrayXd::Zero(n_elements);
    for (long m = 0; m < n_elements; m++) {
        double total = 0.0;
        for (size_t j = 0; j < n_species; j++) {
            total += mole_fracs[j] * gas->nAtoms(j, static_cast<size_t>(m));
        }
        amounts(m) = gas_moles * total;
    }

    if (m_condensed) {
        for (size_t k = 0; k < m_condensed->size(); k++) {
            if (m_condensed->moles[k] == 0.0) continue;
            amounts += m_condensed->moles[k] * m_condensed->species[k].element_atoms;
        }
    }
    return amounts;
}


// Reference state

void Gas::set_stagnation_enthalpy(double H) { m_H_stagnation = H; }
double Gas::get_stagnation_enthalpy() const { return m_H_stagnation; }
void Gas::set_reference_entropy(double S) { m_S0 = S; }
double Gas::get_reference_entropy() const { return m_S0; }

void Gas::set_current_state_as_reference() {
    if (m_sol) {
        m_H_stagnation = thermo()->enthalpy_mass();
        m_S0 = thermo()->entropy_mass();
    }
    else {
        m_H_stagnation = 1.0;
        m_S0 = 1.0;
    }
}

// Accessors

std::shared_ptr<Cantera::Solution> Gas::solution() const { 
    check_for_valid_cantera();
    return m_sol; 
}
std::shared_ptr<Cantera::ThermoPhase> Gas::thermo() const {
    check_for_valid_cantera();
    return m_sol->thermo();
}

std::shared_ptr<Cantera::Kinetics> Gas::kinetics() const { 
    check_for_valid_cantera();
    return m_sol->kinetics(); 
}
std::shared_ptr<Cantera::Transport> Gas::transport() const { 
    check_for_valid_cantera();
    return m_sol->transport(); 
}

std::string Gas::report(bool show_thermo, double threshold) const {
    return thermo()->report(show_thermo, threshold);
}

} // namespace Goddard
