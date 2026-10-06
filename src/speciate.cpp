#include "goddard/speciate.hpp"
#include "goddard/utils.hpp"

#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <cassert>
#include <stdexcept>

using Cantera::AnyMap;
using Cantera::AnyValue;
using Cantera::Solution;

using std::shared_ptr;
using std::string;

namespace Goddard {

namespace {

/**
 * Throw if any requested species name is missing from the names found in an input file.
 *
 * @param requested Species names asked for.
 * @param found Species names present in the input file.
 * @throws std::invalid_argument listing the missing names, sorted.
 */
void check_all_species_found(const std::unordered_set<string>& requested,
                             const std::unordered_set<string>& found) {
    std::vector<string> missing;
    for (const string& name : requested) {
        if (!found.count(name)) {
            missing.push_back(name);
        }
    }
    if (missing.empty()) {
        return;
    }
    std::sort(missing.begin(), missing.end());
    string message = "Species not found in the input file:";
    for (const string& name : missing) {
        message += " " + name;
    }
    throw std::invalid_argument(message);
}

} // namespace

AnyMap speciate(const string& infile, const std::unordered_set<string>& elements) {
    auto root_node = load_root_node(infile);
    return speciate(root_node, elements);
}


AnyMap speciate(AnyMap& root_node, const std::unordered_set<string>& elements) {
    auto& species_list = root_node["species"].asVector<AnyMap>();
    //for the species to be included, it must contain one of the elements
    //provided AND it must not contain any elements that were not provided.
    bool contains_element;
    bool contains_excluded_element;
    
    for (auto it = species_list.begin(); it != species_list.end();) {
        contains_element = false;
        contains_excluded_element = false;
        const auto& composition = (*it)["composition"].asMap<double>();
        for (const auto& element : elements) {
            if (composition.count(element)) {
                contains_element = true;
                break;
            }
        }
        for (const auto& comp: composition) {
            if (!elements.count(comp.first)){
                contains_excluded_element = true;
            }
        }
        //erase invalidates the iterator and returns the next one
        if (!contains_element || contains_excluded_element){
            it = species_list.erase(it);
        }
        else {
            ++it;
        }
    }
    return root_node;
}

AnyMap speciate_from_yaml_string(const string& yaml_str, const std::unordered_set<string>& elements) {
    auto root_node = AnyMap::fromYamlString(yaml_str);
    return speciate(root_node, elements);
}

AnyMap select_species(const string& infile, const std::unordered_set<string>& species) {
    auto root_node = load_root_node(infile);
    return select_species(root_node, species);
}

AnyMap select_species(AnyMap& root_node, const std::unordered_set<string>& species) {
    auto& species_list = root_node["species"].asVector<AnyMap>();
    std::unordered_set<string> found;
    for (auto it = species_list.begin(); it != species_list.end();) {
        const string name = (*it)["name"].asString();
        if (!species.count(name)) {
            //erase invalidates the iterator and returns the next one
            it = species_list.erase(it);
        }
        else {
            found.insert(name);
            ++it;
        }
    }
    check_all_species_found(species, found);
    return root_node;
}

AnyMap select_species_from_yaml_string(const string& yaml_str, const std::unordered_set<string>& species) {
    auto root_node = AnyMap::fromYamlString(yaml_str);
    return select_species(root_node, species);
}

Cantera::AnyMap create_phase_node(
    const std::string& name, 
    const std::vector<std::string>& species, 
    double T, double P)
{
    AnyMap phase;
    phase["name"] = name;
    phase["thermo"] = "ideal-gas";
    phase["species"] = species;
    phase["transport"] = "mixture-averaged";
    phase["skip-undeclared-third-bodies"] = true;
    phase["kinetics"] = "gas";
    phase["reactions"] = "declared-species";
    phase["state"]["T"] = T;
    phase["state"]["P"] = P;

    return phase;
}

AnyMap create_speciated_phase_node(
    const std::string& name, 
    const std::vector<std::string>& elements, 
    double T, double P) 
{
    AnyMap phase;
    phase["name"] = name;
    phase["elements"] = elements;
    phase["thermo"] = "ideal-gas";
    phase["species"] = "all";
    phase["skip-undeclared-elements"] = true;
    phase["skip-undeclared-third-bodies"] = true;
    phase["transport"] = "mixture-averaged";
    phase["kinetics"] = "gas";
    phase["reactions"] = "declared-species";
    phase["state"]["T"] = T;
    phase["state"]["P"] = P;

    return phase;
}

} //namespace Goddard