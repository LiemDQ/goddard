#include "goddard/speciate.hpp"
#include "goddard/global.hpp"
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "cantera/base/AnyMap.h"
#include <vector>
#include <unordered_set>
#include <string>

using Cantera::AnyValue;
using Cantera::AnyMap;
using std::unordered_set;
using std::string;
using std::vector;

using ::testing::UnorderedElementsAreArray;

TEST(Speciation, speciesFilteredByElement) {
    Goddard::setup_defaults();
    unordered_set<string> elements{"H", "O"};

    AnyMap node = Goddard::speciate("nasa_reactants.yaml", elements);

    vector<AnyMap> species_maps = node["species"].asVector<AnyMap>();
    vector<string> species;
    for (const AnyMap& map : species_maps) {
        species.push_back(map["name"].asString());
    }
    
    vector<string> expected_species = {"H2(L)", "H2O2(L)", "O2(L)", "O3(L)"};
    
    EXPECT_EQ(species.size(), expected_species.size());
    
    EXPECT_THAT(species, UnorderedElementsAreArray(expected_species));
}


TEST(Speciation, speciesFilteredByElementH2O2) {
    Goddard::setup_defaults();
    unordered_set<string> elements{"H", "O"};

    AnyMap node = Goddard::speciate("h2o2.yaml", elements);

    vector<AnyMap> species_maps = node["species"].asVector<AnyMap>();
    vector<string> species;
    for (const AnyMap& map : species_maps) {
        species.push_back(map["name"].asString());
    }
    
    vector<string> expected_species = {"H2","H","O", "H2O2", "O2", "OH", "H2O", "HO2"};
    
    EXPECT_EQ(species.size(), expected_species.size());
    
    EXPECT_THAT(species, UnorderedElementsAreArray(expected_species));
}

TEST(Speciation, noSpeciesFoundH2O2) {
    Goddard::setup_defaults();
    unordered_set<string> elements{"F", "B"};

    AnyMap node = Goddard::speciate("h2o2.yaml", elements);

    vector<AnyMap> species_maps = node["species"].asVector<AnyMap>();
    vector<string> species;
    for (const AnyMap& map : species_maps) {
        species.push_back(map["name"].asString());
    }
    
    vector<string> expected_species = {};

    EXPECT_EQ(species.size(), 0);
    
    // EXPECT_THAT(species, UnorderedElementsAreArray(expected_species));
}