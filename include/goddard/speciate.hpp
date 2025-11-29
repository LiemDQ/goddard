#pragma once
#include "cantera/core.h"
#include "cantera/base/AnyMap.h"

#include <memory>
#include <unordered_set>
#include <string>

namespace Goddard {
/**
 * @brief Generate a solution object containing all species containing only the provided elements.
 * 
 */
Cantera::AnyMap speciate(const std::string& infile, const std::unordered_set<std::string>& elements);

Cantera::AnyMap speciate(Cantera::AnyMap& root_node, const std::unordered_set<std::string>& elements);

Cantera::AnyMap speciate_from_yaml_string(const std::string& yaml_str, const std::unordered_set<std::string>& elements);



} //namespace Goddard