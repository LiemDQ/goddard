#include "goddard/speciate.hpp"

#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"

#include <vector>
#include <string>
#include <iostream>
#include <cassert>

using Cantera::AnyMap;
using Cantera::AnyValue;
using Cantera::Solution;

using std::shared_ptr;
using std::string;

namespace Goddard {

AnyMap speciate(const string& infile, const std::unordered_set<string>& elements) {
    
    size_t dot = infile.find_last_of('c');
    string extension;
    if (dot != Cantera::npos) {
        extension = Cantera::toLowerCopy(infile.substr(dot+1));
    }

    if (extension == "cti" || extension == "xml") {
        throw Cantera::CanteraError("newSolution",
                           "The CTI and XML formats are no longer supported.");
    }
    
    auto root_node = AnyMap::fromYamlFile(infile);
    
    return speciate(root_node, elements);
}


AnyMap speciate(AnyMap& root_node, const std::unordered_set<string>& elements) {
    
    // const auto& species_node = root_node["species"].asMap("name");
    assert(root_node["species"].isVector<AnyMap>() && "Must be a vector of Maps");
    auto& species_list = root_node["species"].asVector<AnyMap>();
    std::cerr << "GOT SPECIES LIST\n";
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
        //iterator is invalidated when erase is called, 
        //so only increment when no erasure happens
        if (!contains_element || contains_excluded_element){
            species_list.erase(it);
        }
        else {
            ++it;
        }
    }
    
    return root_node;
}

AnyMap speciate_from_yaml_string(const std::string& yaml_str, const std::unordered_set<std::string>& elements) {
    
    auto root_node = AnyMap::fromYamlString(yaml_str);
    
    return speciate(root_node, elements);
}

} //namespace Goddard