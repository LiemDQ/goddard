#pragma once
#include "cantera/core.h"
#include "cantera/base/AnyMap.h"

#include <memory>
#include <unordered_set>
#include <string>

namespace Goddard {
/**
 * @brief Generate a solution object with all species from the input file containing only the provided elements.
 * 
 */
Cantera::AnyMap speciate(const std::string& infile, const std::unordered_set<std::string>& elements);

Cantera::AnyMap speciate(Cantera::AnyMap& root_node, const std::unordered_set<std::string>& elements);

Cantera::AnyMap speciate_from_yaml_string(const std::string& yaml_str, const std::unordered_set<std::string>& elements);

/**
 * @brief Select only the given species from the input file. 
 * Note that the species name must be an exact match.
 */
Cantera::AnyMap select_species(const std::string& infile, const std::unordered_set<std::string>& species);

Cantera::AnyMap select_species(Cantera::AnyMap& root_node, const std::unordered_set<std::string>& species);

Cantera::AnyMap select_species_from_yaml_string(const std::string& yaml_str, const std::unordered_set<std::string>& species);

} //namespace Goddard