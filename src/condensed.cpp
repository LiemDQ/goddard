#include "goddard/condensed.hpp"
#include "goddard/error.hpp"

#include "cantera/base/AnyMap.h"
#include "cantera/core.h"
#include "cantera/thermo/ThermoFactory.h"

#include <algorithm>
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

std::vector<size_t> CondensedPhaseSet::offered_at(double T) const {
    std::vector<size_t> offered;
    for (const PolymorphGroup& group : groups) {
        // Members are ordered by increasing T_min, so the first in-range member is the lower
        // polymorph of a coexisting pair at a transition temperature.
        for (size_t k : group.members) {
            if (in_range(k, T)) {
                offered.push_back(k);
                break;
            }
        }
    }
    std::sort(offered.begin(), offered.end());
    return offered;
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

} // namespace Goddard
