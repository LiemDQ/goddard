#pragma once
#include "cantera/core.h"
#include "cantera/base/AnyMap.h"

#include <unordered_set>
#include <string>
#include <vector>

namespace Goddard {
/**
 * Generate an AnyMap object with all species from the input file containing only the provided elements.
 * @warning For databases with large numbers of species (e.g. nasa_gas), this can lead to poor performance
 * due to large numbers of trace species.
 * 
 * @return root node with selected species.
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


/**
 * Generate a Phase node for use with `newSolution` with 
 * the provided species.
 * 
 * @return Cantera `AnyMap` object suitable for use as a phase node in Cantera factory functions. 
 */
Cantera::AnyMap create_phase_node(
    const std::string& name, 
    const std::vector<std::string>& species, 
    double T = 298.15, double P = 101325);

/**
 * Generate a Phase node for use with `newSolution` with all species from the input file containing only the provided
 * elements.
 * 
 * @return Cantera `AnyMap` object suitable for use as a phase node in Cantera factory functions. 
 */
Cantera::AnyMap create_speciated_phase_node(
    const std::string& name, 
    const std::vector<std::string>& elements, 
    double T = 298.15, double P = 101325);

struct Speciation { //TODO: implement speciation functionality
    std::unordered_set<std::string> elements;
    std::unordered_set<std::string> species;
    int max_species = 10;
};

} //namespace Goddard