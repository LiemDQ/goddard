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
 * Select only the given species from the input file.
 *
 * The species names must match exactly. The selected species keep their order in the file.
 *
 * @throws Cantera::CanteraError for a CTI or XML input file.
 * @throws std::invalid_argument if a requested species is not in the file, listing the missing
 * names.
 */
Cantera::AnyMap select_species(const std::string& infile, const std::unordered_set<std::string>& species);

/**
 * Remove every species not in `species` from the `species` section of a root node, in place.
 *
 * @return the modified root node.
 * @throws std::invalid_argument if a requested species is not in the root node, listing the
 * missing names.
 */
Cantera::AnyMap select_species(Cantera::AnyMap& root_node, const std::unordered_set<std::string>& species);

/**
 * Order species names as they appear in the `species` section of a root node.
 *
 * Building a phase from the result gives the same species indices whatever the order of
 * `species`, e.g. when it comes from an unordered set. If the root node has no `species` section
 * the names are sorted instead.
 *
 * @param root_node Root node of a Cantera YAML input file.
 * @param species Species names; a repeated name is kept once.
 * @return the names in file order.
 * @throws std::invalid_argument if a name is not in the `species` section, listing the missing
 * names.
 */
std::vector<std::string> species_in_file_order(
    const Cantera::AnyMap& root_node, const std::vector<std::string>& species);

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