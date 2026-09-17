#include "goddard/condensed.hpp"
#include "goddard/equilibrium.hpp"
#include "goddard/error.hpp"
#include "goddard/gas.hpp"

#include "cantera/base/AnyMap.h"
#include "cantera/core.h"
#include "cantera/equil/MultiPhase.h"
#include "cantera/thermo/ThermoFactory.h"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

using Cantera::AnyMap;
using Cantera::ThermoPhase;

namespace Goddard {

namespace {

/** Temperature window [K] within which two data ranges are considered contiguous. */
constexpr double TRANSITION_TOLERANCE = 1e-6;

/** Pressure [Pa] the candidate phases are held at while reference-state data is read. */
constexpr double REFERENCE_PRESSURE = 1.0e5;

/**
 * Largest mass imbalance [kg per kg of mixture] accepted from a converged multiphase solve.
 * `MultiPhaseEquil` conserves the element amounts to a few times its own relative tolerance, and
 * the mixture mass is a weighted sum of those, so this cannot be pushed to machine precision.
 */
constexpr double MASS_CLOSURE_TOLERANCE = 1e-8;

/** Largest relative change of any element amount accepted from a converged multiphase solve. */
constexpr double ELEMENT_CLOSURE_TOLERANCE = 1e-7;

/** Map from species name to its node, for every species in a data file root node. */
std::map<std::string, AnyMap> index_species_nodes(const AnyMap& root_node) {
    std::map<std::string, AnyMap> nodes;
    if (!root_node.hasKey("species")) {
        return nodes;
    }
    const auto& species_list = root_node.at("species").asVector<AnyMap>();
    for (const AnyMap& node : species_list) {
        nodes[node["name"].asString()] = node;
    }
    return nodes;
}

} // namespace

AnyMap make_condensed_phase_node(const std::string& species_name) {
    AnyMap phase;
    phase["name"] = species_name;
    phase["thermo"] = "fixed-stoichiometry";
    phase["species"] = std::vector<std::string>{species_name};
    return phase;
}

std::shared_ptr<ThermoPhase> build_condensed_phase(AnyMap species_node, double density) {
    const std::string name = species_node["name"].asString();

    // NASA-database condensed species carry no equation of state; StoichSubstance needs one.
    if (!species_node.hasKey("equation-of-state")) {
        AnyMap eos;
        eos["model"] = "constant-volume";
        eos["density"] = density;
        species_node["equation-of-state"] = eos;
    }

    AnyMap root;
    root["species"] = std::vector<AnyMap>{species_node};

    AnyMap phase_node = make_condensed_phase_node(name);
    return Cantera::newThermo(phase_node, root);
}

// ---- CondensedPhaseSet ----

CondensedPhaseSet::CondensedPhaseSet(const AnyMap& root_node,
                                     const std::vector<std::string>& names,
                                     const ThermoPhase& gas)
{
    m_element_names = gas.elementNames();
    add_species(root_node, names, gas);
}

std::vector<std::string> CondensedPhaseSet::compatible_species(const AnyMap& root_node,
                                                               const ThermoPhase& gas)
{
    std::unordered_set<std::string> gas_elements;
    for (const std::string& element : gas.elementNames()) {
        gas_elements.insert(element);
    }

    std::vector<std::string> compatible;
    if (!root_node.hasKey("species")) {
        return compatible;
    }
    const auto& species_list = root_node.at("species").asVector<AnyMap>();
    for (const AnyMap& node : species_list) {
        const auto composition = node["composition"].asMap<double>();
        bool usable = true;
        for (const auto& entry : composition) {
            if (gas_elements.count(entry.first) == 0) {
                usable = false;
                break;
            }
        }
        if (usable) {
            compatible.push_back(node["name"].asString());
        }
    }
    return compatible;
}

void CondensedPhaseSet::add_species(const AnyMap& root_node,
                                    const std::vector<std::string>& names,
                                    const ThermoPhase& gas)
{
    if (m_element_names.empty()) {
        m_element_names = gas.elementNames();
    }

    const std::map<std::string, AnyMap> nodes = index_species_nodes(root_node);
    const size_t n_elements = gas.nElements();

    for (const std::string& name : names) {
        if (species_index(name) != Cantera::npos) {
            continue; // already a candidate
        }

        auto found = nodes.find(name);
        if (found == nodes.end()) {
            throw FmtError("Condensed species '{}' was not found in the data file.", name);
        }

        CondensedSpecies entry;
        entry.name = name;
        entry.phase = build_condensed_phase(found->second);
        entry.molar_mass = entry.phase->molecularWeight(0);
        entry.T_min = entry.phase->minTemp();
        entry.T_max = entry.phase->maxTemp();

        entry.element_atoms = Eigen::ArrayXd::Zero(static_cast<long>(n_elements));
        const auto composition = found->second["composition"].asMap<double>();
        for (const auto& element : composition) {
            size_t m = gas.elementIndex(element.first);
            if (m == Cantera::npos) {
                throw FmtError(
                    "Condensed species '{}' contains element '{}', which is not present in the "
                    "gas phase.",
                    name, element.first);
            }
            entry.element_atoms(static_cast<long>(m)) = element.second;
        }

        species.push_back(std::move(entry));
        m_species_nodes.push_back(found->second);
        moles.push_back(0.0);
    }

    build_groups();
}

void CondensedPhaseSet::build_groups() {
    groups.clear();

    // Species with identical elemental composition are polymorphs of one another.
    std::map<std::vector<std::pair<long, double>>, size_t> group_of_composition;
    for (size_t k = 0; k < species.size(); k++) {
        std::vector<std::pair<long, double>> key;
        for (long m = 0; m < species[k].element_atoms.size(); m++) {
            if (species[k].element_atoms(m) != 0.0) {
                key.emplace_back(m, species[k].element_atoms(m));
            }
        }
        auto found = group_of_composition.find(key);
        if (found == group_of_composition.end()) {
            group_of_composition[key] = groups.size();
            species[k].group = static_cast<int>(groups.size());
            groups.push_back(PolymorphGroup{{k}, {}});
        } else {
            species[k].group = static_cast<int>(found->second);
            groups[found->second].members.push_back(k);
        }
    }

    for (PolymorphGroup& group : groups) {
        std::sort(group.members.begin(), group.members.end(),
                  [this](size_t a, size_t b) { return species[a].T_min < species[b].T_min; });

        group.transition_temperatures.clear();
        for (size_t i = 0; i + 1 < group.members.size(); i++) {
            const double upper = species[group.members[i]].T_max;
            const double lower = species[group.members[i + 1]].T_min;
            if (std::abs(upper - lower) < TRANSITION_TOLERANCE) {
                group.transition_temperatures.push_back(upper);
            }
        }
    }
}

std::shared_ptr<CondensedPhaseSet> CondensedPhaseSet::clone() const {
    // The copy constructor would share the phase objects, whose state is mutated when
    // reference-state data is read, so the phases are rebuilt from the stored species nodes.
    std::shared_ptr<CondensedPhaseSet> copy(new CondensedPhaseSet(*this));
    for (size_t k = 0; k < copy->species.size(); k++) {
        copy->species[k].phase = build_condensed_phase(copy->m_species_nodes[k]);
    }
    return copy;
}

size_t CondensedPhaseSet::species_index(const std::string& name) const {
    for (size_t k = 0; k < species.size(); k++) {
        if (species[k].name == name) {
            return k;
        }
    }
    return Cantera::npos;
}

std::vector<std::string> CondensedPhaseSet::names() const {
    std::vector<std::string> out;
    out.reserve(species.size());
    for (const CondensedSpecies& entry : species) {
        out.push_back(entry.name);
    }
    return out;
}

bool CondensedPhaseSet::in_range(size_t k, double T) const {
    return T >= species.at(k).T_min && T <= species.at(k).T_max;
}

size_t CondensedPhaseSet::offered_member(size_t group_index, double T) const {
    // Members are ordered by increasing T_min, so the first in-range member is the lower
    // polymorph of a coexisting pair at a transition temperature.
    for (size_t k : groups.at(group_index).members) {
        if (in_range(k, T)) {
            return k;
        }
    }
    return Cantera::npos;
}

std::vector<size_t> CondensedPhaseSet::offered_at(double T) const {
    std::vector<size_t> offered;
    for (size_t g = 0; g < groups.size(); g++) {
        const size_t k = offered_member(g, T);
        if (k != Cantera::npos) {
            offered.push_back(k);
        }
    }
    std::sort(offered.begin(), offered.end());
    return offered;
}

void CondensedPhaseSet::update_pinned_group() {
    pinned_group = -1;
    for (size_t g = 0; g < groups.size(); g++) {
        int present = 0;
        for (size_t k : groups[g].members) {
            if (moles[k] > 0.0) present++;
        }
        if (present > 1) {
            pinned_group = static_cast<int>(g);
            return;
        }
    }
}

double CondensedPhaseSet::enthalpy_RT(size_t k, double T) const {
    const auto& phase = species.at(k).phase;
    phase->setState_TP(T, REFERENCE_PRESSURE);
    double value = 0.0;
    phase->getEnthalpy_RT_ref(&value);
    return value;
}

double CondensedPhaseSet::entropy_R(size_t k, double T) const {
    const auto& phase = species.at(k).phase;
    phase->setState_TP(T, REFERENCE_PRESSURE);
    double value = 0.0;
    phase->getEntropy_R_ref(&value);
    return value;
}

double CondensedPhaseSet::cp_R(size_t k, double T) const {
    const auto& phase = species.at(k).phase;
    phase->setState_TP(T, REFERENCE_PRESSURE);
    double value = 0.0;
    phase->getCp_R_ref(&value);
    return value;
}


// ---- Multiphase equilibrium (the Gas members implemented on top of CondensedPhaseSet) ----

namespace {

/** Mixture property of `gas` selected by `property`, per kg of mixture. */
double mixture_value(const Gas& gas, EquilibriumProperty property) {
    return property == EquilibriumProperty::ENTHALPY ? gas.enthalpy_mass() : gas.entropy_mass();
}

/**
 * Reference-state molar value of candidate `k` at `T` selected by `property`
 * [J/kmol] or [J/(kmol.K)].
 */
double condensed_molar_value(const CondensedPhaseSet& set,
                             EquilibriumProperty property,
                             size_t k,
                             double T)
{
    if (property == EquilibriumProperty::ENTHALPY) {
        return set.enthalpy_RT(k, T) * Cantera::GasConstant * T;
    }
    return set.entropy_R(k, T) * Cantera::GasConstant;
}

/** One phase transition of one polymorph group, with the two polymorphs that meet there. */
struct GroupTransition {
    double temperature = 0.0;   //!< Transition temperature [K].
    size_t group = 0;           //!< Index into `CondensedPhaseSet::groups`.
    size_t low = 0;             //!< Candidate index of the lower-temperature polymorph.
    size_t high = 0;            //!< Candidate index of the higher-temperature polymorph.
};

/**
 * Set the gas composition to the monatomic species of each element, holding `element_moles`.
 *
 * `Gas::set_element_moles` prefers the homonuclear diatomic of an element, which is the better
 * physical guess but makes `MultiPhaseEquil` stall on some carbon-bearing mixtures. Starting from
 * free atoms is the guess the Python study found robust, so it is one rung of the restart ladder.
 */
void set_monatomic_basis(Cantera::ThermoPhase& gas, const Eigen::ArrayXd& element_moles,
                         double T, double P)
{
    Eigen::ArrayXd moles = Eigen::ArrayXd::Zero(static_cast<long>(gas.nSpecies()));
    for (size_t m = 0; m < gas.nElements(); m++) {
        const double amount = element_moles(static_cast<long>(m));
        if (amount <= 0.0) continue;

        size_t chosen = Cantera::npos;
        double chosen_atoms = 0.0;
        for (size_t j = 0; j < gas.nSpecies(); j++) {
            const double atoms = gas.nAtoms(j, m);
            if (atoms <= 0.0) continue;
            bool single_element = true;
            for (size_t i = 0; i < gas.nElements(); i++) {
                if (i != m && gas.nAtoms(j, i) != 0.0) { single_element = false; break; }
            }
            if (!single_element) continue;
            if (chosen == Cantera::npos || atoms < chosen_atoms) {
                chosen = j;
                chosen_atoms = atoms;
            }
        }
        if (chosen == Cantera::npos) {
            throw FmtError("no species made of element '{}' alone is available as a starting guess",
                           gas.elementName(m));
        }
        moles(static_cast<long>(chosen)) += amount / chosen_atoms;
    }
    gas.setState_TPX(T, P, moles.data());
}

/** All phase transitions of a set, sorted by increasing temperature. */
std::vector<GroupTransition> all_transitions(const CondensedPhaseSet& set) {
    std::vector<GroupTransition> transitions;
    for (size_t g = 0; g < set.groups.size(); g++) {
        const PolymorphGroup& group = set.groups[g];
        for (size_t i = 0; i < group.transition_temperatures.size(); i++) {
            // transition_temperatures[i] separates members[i] from members[i + 1].
            transitions.push_back(GroupTransition{group.transition_temperatures[i], g,
                                                  group.members[i], group.members[i + 1]});
        }
    }
    std::sort(transitions.begin(), transitions.end(),
              [](const GroupTransition& a, const GroupTransition& b) {
                  return a.temperature < b.temperature;
              });
    return transitions;
}

} // namespace

void Gas::solve_multiphase_TP(double T, double P, const std::vector<size_t>& offered_override) {
    CondensedPhaseSet& set = *m_condensed;

    const std::vector<size_t> offered =
        offered_override.empty() ? set.offered_at(T) : offered_override;

    std::vector<size_t> offered_of_group(set.groups.size(), Cantera::npos);
    for (size_t k : offered) {
        offered_of_group[static_cast<size_t>(set.species[k].group)] = k;
    }

    // Gordon & McBride replacement step: a group's whole amount moves onto the one polymorph that
    // is offered to the solver. Away from a transition that is the only in-range member; exactly at
    // one it is the lower polymorph, whose chemical potential equals the upper one's, so the gas
    // composition of the solve does not depend on the choice.
    for (size_t g = 0; g < set.groups.size(); g++) {
        double total = 0.0;
        for (size_t k : set.groups[g].members) {
            total += set.moles[k];
        }
        if (total <= 0.0) {
            continue;
        }
        const size_t target = offered_of_group[g];
        if (target == Cantera::npos) {
            const size_t present = set.groups[g].members.front();
            throw FmtError(
                "Gas::equilibrate_TP: condensed species '{}' is present ({:.6g} kmol/kg) but no "
                "polymorph of its group has thermodynamic data at T = {:.4f} K.",
                set.species[present].name, total, T);
        }
        for (size_t k : set.groups[g].members) {
            set.moles[k] = 0.0;
        }
        set.moles[target] = total;
    }
    set.pinned_group = -1;

    const Eigen::ArrayXd elements_before = element_moles();
    std::vector<double> gas_state_before(thermo()->stateSize());
    thermo()->saveState(gas_state_before);
    const std::vector<double> moles_before = set.moles;

    // Restart ladder of starting compositions, in increasing order of desperation. There is no
    // rung varying Cantera's `estimate_equil`: MultiPhaseEquil ignores it (checked against
    // Cantera 3.2), so such a rung would repeat the previous solve step for step.
    enum class StartingGuess { INCOMING, FREE_ATOMS, GAS_EQUILIBRIUM };
    const StartingGuess restarts[] = {
        StartingGuess::INCOMING,         // warm-started from the previous solve
        StartingGuess::FREE_ATOMS,       // no condensed phase
        StartingGuess::GAS_EQUILIBRIUM,  // gas-only equilibrium holding the same elements
    };

    std::string last_failure;
    for (const StartingGuess restart : restarts) {
        thermo()->restoreState(gas_state_before);
        set.moles = moles_before;
        set.pinned_group = -1;

        try {
            if (restart != StartingGuess::INCOMING) {
                std::fill(set.moles.begin(), set.moles.end(), 0.0);
                if (restart == StartingGuess::FREE_ATOMS) {
                    set_monatomic_basis(*thermo(), elements_before, T, P);
                } else {
                    set_element_moles(elements_before, T, P);
                    thermo()->equilibrate("TP", "gibbs");
                }
            }

            double condensed_mass = 0.0;
            for (size_t k = 0; k < set.size(); k++) {
                condensed_mass += set.moles[k] * set.species[k].molar_mass;
            }
            if (condensed_mass >= 1.0) {
                throw FmtError(
                    "Gas::equilibrate_TP: the condensed amounts weigh {:.6g} kg per kg of mixture, "
                    "leaving no gas phase.", condensed_mass);
            }

            thermo()->setState_TP(T, P);
            const double gas_moles = (1.0 - condensed_mass) / thermo()->meanMolecularWeight();

            Cantera::MultiPhase mixture;
            mixture.addPhase(thermo(), gas_moles);
            for (size_t k : offered) {
                // Phases offered with zero moles still grow if the minimization calls for them.
                mixture.addPhase(set.species[k].phase, set.moles[k]);
            }
            mixture.init();
            mixture.setState_TP(T, P);
            // "vcs" returns wrong answers with condensed phases, so the Gibbs solver is the only
            // one used here.
            mixture.equilibrate("TP", "gibbs", equilibrium_options.rtol,
                                equilibrium_options.max_steps, 100, 0, 0);

            std::fill(set.moles.begin(), set.moles.end(), 0.0);
            for (size_t i = 0; i < offered.size(); i++) {
                const double amount = mixture.phaseMoles(i + 1);
                set.moles[offered[i]] = amount > 0.0 ? amount : 0.0;
            }
            m_equilibrium_solve_count++;
            thermo()->setState_TP(T, P);

            double mass = mixture.phaseMoles(0) * thermo()->meanMolecularWeight();
            for (size_t k = 0; k < set.size(); k++) {
                mass += set.moles[k] * set.species[k].molar_mass;
            }
            if (std::abs(mass - 1.0) > MASS_CLOSURE_TOLERANCE) {
                throw FmtError("mass closure is off by {:.3e} kg per kg of mixture", mass - 1.0);
            }

            const Eigen::ArrayXd elements_after = element_moles();
            for (long m = 0; m < elements_before.size(); m++) {
                const double scale = std::abs(elements_before(m)) + 1e-30;
                if (std::abs(elements_after(m) - elements_before(m)) > ELEMENT_CLOSURE_TOLERANCE * scale) {
                    throw FmtError("element '{}' changed by {:.3e} relative",
                                   thermo()->elementName(static_cast<size_t>(m)),
                                   (elements_after(m) - elements_before(m)) / scale);
                }
            }
            return;
        } catch (const std::exception& error) {
            last_failure = error.what();
        }
    }

    thermo()->restoreState(gas_state_before);
    set.moles = moles_before;
    throw ConvergenceError(
        std::format("Gas::equilibrate_TP: multiphase equilibrium at T = {:.4f} K, P = {:.6g} Pa "
                    "failed from every starting guess. Last failure: {}", T, P, last_failure),
        static_cast<int>(std::size(restarts)), equilibrium_options.rtol, -1.0);
}

int Gas::solve_multiphase_XP(EquilibriumProperty property, double target, double P) {
    CondensedPhaseSet& set = *m_condensed;
    const EquilibriumOptions& options = equilibrium_options;
    m_equilibrium_solve_count = 0;

    // The residual is monotonically increasing in temperature (cp > 0) and jumps by the latent heat
    // of a transition wherever the group undergoing it is present.
    auto residual_at = [&](double T, const std::vector<size_t>& offer) {
        solve_multiphase_TP(T, P, offer);
        return mixture_value(*this, property) - target;
    };

    // Newton step in ln T from the state just solved, using the equilibrium heat capacity:
    // (dH/d ln T)_P = cp T and (dS/d ln T)_P = cp. NaN when that heat capacity is unusable.
    auto newton_log_T_step = [&](double residual) {
        double heat_capacity = 0.0;
        try {
            heat_capacity = get_thermo_equilibrium_properties(*this).spec_heat_p;
        } catch (const std::exception&) {
            // A singular derivative system leaves only the bracketing steps.
            return std::numeric_limits<double>::quiet_NaN();
        }
        if (!(std::isfinite(heat_capacity) && heat_capacity > 0.0)) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return property == EquilibriumProperty::ENTHALPY
            ? -residual / (heat_capacity * temperature())
            : -residual / heat_capacity;
    };

    const double log_T_min = std::log(options.T_min);
    const double log_T_max = std::log(options.T_max);

    double T_start = temperature();
    if (!(T_start >= options.T_min && T_start <= options.T_max)) {
        T_start = options.T_default;
    }
    const double log_T_start = std::clamp(std::log(T_start), log_T_min, log_T_max);

    double log_T_a = log_T_start;
    double residual_a = residual_at(std::exp(log_T_a), {});
    double log_T_b = log_T_a;
    double residual_b = residual_a;

    if (residual_a != 0.0) {
        const double direction = residual_a < 0.0 ? 1.0 : -1.0;
        // Each probe marches on from the previous one by an overshot Newton step, so a probe that
        // falls short of the root tightens the bracket for the next. The smallest step doubles
        // with every probe, so a Newton step that keeps falling short cannot stall the march;
        // without a usable heat capacity the step is that doubling one.
        double min_step = 1e-3;
        bool bracketed = false;
        for (int i = 0; i < options.max_bracket_steps; i++) {
            const double newton_step = newton_log_T_step(residual_b);
            double step = std::isfinite(newton_step)
                ? std::clamp(1.5 * std::abs(newton_step), min_step, std::max(min_step, 0.5))
                : 50.0 * min_step;
            min_step *= 2.0;

            log_T_a = log_T_b;
            residual_a = residual_b;
            log_T_b = std::clamp(log_T_a + direction * step, log_T_min, log_T_max);
            residual_b = residual_at(std::exp(log_T_b), {});
            if (residual_a * residual_b <= 0.0) {
                bracketed = true;
                break;
            }
            if (log_T_b <= log_T_min || log_T_b >= log_T_max) {
                break;
            }
        }
        if (!bracketed) {
            const double T_end = std::exp(log_T_b);
            throw ConvergenceError(
                std::format("Gas::equilibrate_{}P: could not bracket the temperature root between "
                            "{:.2f} K and {:.2f} K (residuals {:.6g} and {:.6g}).",
                            property == EquilibriumProperty::ENTHALPY ? "H" : "S",
                            T_start, T_end, residual_a, residual_b),
                m_equilibrium_solve_count, options.T_rel_tol, residual_b);
        }
    }

    double log_T_lo = log_T_a;
    double residual_lo = residual_a;
    double log_T_hi = log_T_b;
    double residual_hi = residual_b;
    if (residual_lo > 0.0) {
        std::swap(log_T_lo, log_T_hi);
        std::swap(residual_lo, residual_hi);
    }

    // Each transition inside the bracket is a jump in the residual. Evaluating just below and just
    // above it either places the target inside the jump (a pinned state) or says which side of the
    // transition the root is on.
    for (const GroupTransition& transition : all_transitions(set)) {
        // Endpoints count: a warm start from the previous station often sits exactly on the
        // transition, and the root of the next station may be the very same pinned temperature.
        const double log_T_tr = std::log(transition.temperature);
        if (log_T_tr < std::min(log_T_lo, log_T_hi) || log_T_tr > std::max(log_T_lo, log_T_hi)) {
            continue;
        }

        const double residual_below = residual_at(transition.temperature, {});
        const double group_moles = set.moles[transition.low];

        if (group_moles <= 0.0) {
            // The polymorphs have equal chemical potentials at the transition, so a group absent
            // just below it is absent just above it too: the residual does not jump there.
            if (residual_below < 0.0) {
                log_T_lo = log_T_tr;
                residual_lo = residual_below;
            } else {
                log_T_hi = log_T_tr;
                residual_hi = residual_below;
            }
            continue;
        }
        const std::vector<double> state_below = save_state();

        std::vector<size_t> offer_above = set.offered_at(transition.temperature);
        for (size_t& k : offer_above) {
            if (set.species[k].group == set.species[transition.high].group) {
                k = transition.high;
            }
        }
        const double residual_above = residual_at(transition.temperature, offer_above);

        if (group_moles > 0.0 && residual_below <= 0.0 && residual_above >= 0.0) {
            // Pinned: both polymorphs coexist at the transition. The gas composition is the one of
            // the lower-polymorph solve (their chemical potentials are equal there), and the split
            // follows from the linear enthalpy (entropy) balance across the latent heat.
            // Keep the lower-polymorph state itself and re-read its property: the split has to
            // reproduce the target property exactly.
            restore_state(state_below);
            const double value_low_only = mixture_value(*this, property);
            const double value_low = condensed_molar_value(set, property, transition.low,
                                                           transition.temperature);
            const double value_high = condensed_molar_value(set, property, transition.high,
                                                            transition.temperature);
            const double moles = set.moles[transition.low];
            double upper_fraction =
                (target - value_low_only) / (moles * (value_high - value_low));
            upper_fraction = std::clamp(upper_fraction, 0.0, 1.0);

            set.moles[transition.low] = (1.0 - upper_fraction) * moles;
            set.moles[transition.high] = upper_fraction * moles;
            // Callers (and the derivative code) treat a state as pinned only when both polymorphs
            // are really present; at the ends of the latent-heat band it is an ordinary state.
            if (set.moles[transition.low] > 0.0 && set.moles[transition.high] > 0.0) {
                set.pinned_group = static_cast<int>(transition.group);
            }
            return m_equilibrium_solve_count;
        }

        if (residual_above < 0.0) {
            log_T_lo = log_T_tr;
            residual_lo = residual_above;
        } else {
            log_T_hi = log_T_tr;
            residual_hi = residual_below;
        }
    }

    // Newton on ln T, safeguarded by the bracket: a Newton step that leaves the bracket, or follows
    // one that did not halve the residual, is replaced by an Illinois-modified regula falsi step.
    // The latter happens once the residual reaches the noise of the inner solves, where Newton steps
    // no longer shrink the bracket. The first Newton step starts from whatever state was solved
    // last, whose residual needs no further solve.
    double previous_residual = mixture_value(*this, property) - target;
    double predicted_log_T = std::log(temperature()) + newton_log_T_step(previous_residual);
    bool newton_productive = true;
    const int max_iterations = 200;
    int side = 0;
    for (int iteration = 0; iteration < max_iterations; iteration++) {
        double log_T = (log_T_lo * residual_hi - log_T_hi * residual_lo)
            / (residual_hi - residual_lo);
        if (!(log_T > std::min(log_T_lo, log_T_hi) && log_T < std::max(log_T_lo, log_T_hi))) {
            log_T = 0.5 * (log_T_lo + log_T_hi);
        }
        const bool newton = newton_productive
            && predicted_log_T > std::min(log_T_lo, log_T_hi)
            && predicted_log_T < std::max(log_T_lo, log_T_hi);
        if (newton) {
            log_T = predicted_log_T;
        }

        const double T = std::exp(log_T);
        const double residual = residual_at(T, {});
        const double scale =
            property == EquilibriumProperty::ENTHALPY ? cp_mass() * T : cp_mass();
        if (std::abs(residual) <= 1e-10 * scale) {
            return m_equilibrium_solve_count;
        }

        predicted_log_T = log_T + newton_log_T_step(residual);
        newton_productive = !newton || std::abs(residual) <= 0.5 * std::abs(previous_residual);
        previous_residual = residual;

        if (residual < 0.0) {
            log_T_lo = log_T;
            residual_lo = residual;
            if (!newton && side == -1) residual_hi *= 0.5;
            side = -1;
        } else {
            log_T_hi = log_T;
            residual_hi = residual;
            if (!newton && side == 1) residual_lo *= 0.5;
            side = 1;
        }

        if (std::abs(log_T_hi - log_T_lo) <= options.T_rel_tol) {
            return m_equilibrium_solve_count;
        }
    }

    throw ConvergenceError(
        std::format("Gas::equilibrate_{}P: the temperature root did not converge; the bracket is "
                    "still [{:.6f}, {:.6f}] K.",
                    property == EquilibriumProperty::ENTHALPY ? "H" : "S",
                    std::exp(std::min(log_T_lo, log_T_hi)),
                    std::exp(std::max(log_T_lo, log_T_hi))),
        max_iterations, options.T_rel_tol, residual_hi);
}

void Gas::solve_frozen_XP(EquilibriumProperty property, double target, double P) {
    CondensedPhaseSet& set = *m_condensed;

    // A frozen expansion keeps the condensed amounts, so the temperature cannot leave the data
    // range of any species that is present. CEA stops with an error there; so does Goddard.
    double T_lowest = 0.0;
    double T_highest = std::numeric_limits<double>::max();
    size_t low_limiter = Cantera::npos;
    size_t high_limiter = Cantera::npos;
    for (size_t k = 0; k < set.size(); k++) {
        if (set.moles[k] <= 0.0) continue;
        if (set.species[k].T_min > T_lowest) {
            T_lowest = set.species[k].T_min;
            low_limiter = k;
        }
        if (set.species[k].T_max < T_highest) {
            T_highest = set.species[k].T_max;
            high_limiter = k;
        }
    }

    double T = temperature();
    if (!(T > 0.0)) {
        T = equilibrium_options.T_default;
    }
    T = std::clamp(T, T_lowest, T_highest);

    const int max_iterations = 100;
    double residual = 0.0;
    for (int iteration = 0; iteration < max_iterations; iteration++) {
        thermo()->setState_TP(T, P);
        const double value = mixture_value(*this, property);
        const double heat_capacity = cp_mass();
        residual = value - target;

        const double scale = property == EquilibriumProperty::ENTHALPY
            ? std::abs(target) + heat_capacity * T
            : std::abs(target) + heat_capacity;
        if (std::abs(residual) <= 1e-12 * scale) {
            return;
        }

        const double step = property == EquilibriumProperty::ENTHALPY
            ? -residual / heat_capacity
            : -residual * T / heat_capacity;
        // Damp the step so a bad guess cannot leave the valid temperature window in one go.
        const double T_next = std::clamp(std::clamp(T + step, 0.5 * T, 2.0 * T),
                                         T_lowest, T_highest);

        if (T_next == T) {
            const size_t limiter = step < 0.0 ? low_limiter : high_limiter;
            if (limiter != Cantera::npos) {
                throw FmtError(
                    "Gas::set_state_{}P: the frozen state leaves the temperature range of "
                    "condensed species '{}' ({:.2f} - {:.2f} K) at T = {:.2f} K.",
                    property == EquilibriumProperty::ENTHALPY ? "H" : "S",
                    set.species[limiter].name, set.species[limiter].T_min,
                    set.species[limiter].T_max, T);
            }
            return;
        }
        T = T_next;
    }

    throw ConvergenceError("Gas::set_state_HP/SP: the frozen temperature iteration did not "
                           "converge with condensed species present.",
                           max_iterations, 1e-12, residual);
}

} // namespace Goddard
