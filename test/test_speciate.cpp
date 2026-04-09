#include "goddard/speciate.hpp"
#include "goddard/global.hpp"
#include "goddard/utils.hpp"
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
using ::testing::NotNull;

TEST(Speciate, speciesFilteredByElement) {
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


TEST(Speciate, speciesFilteredByElementH2O2) {
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

TEST(Speciate, noSpeciesFoundH2O2) {
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
    
    EXPECT_THAT(species, UnorderedElementsAreArray(expected_species));
}

TEST(Speciate, canConstructCanteraSolutionFromFilePhase) {
    Goddard::setup_defaults();
    unordered_set<string> elements{"H", "O", "Ar", "N"};
    std::vector<std::string> expected_species{"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "AR", "N2"};

    AnyMap node = Goddard::speciate("h2o2.yaml", elements);

    auto phase_nodes = node["phases"].asVector<AnyMap>();
    auto phase_node = phase_nodes[0];

    auto sln = Cantera::newSolution(phase_node, node);
    ASSERT_THAT(sln, NotNull()) << "Solution object should be initialized.";
    EXPECT_STREQ(sln->name().c_str(), "ohmech");

    EXPECT_THAT(sln->thermo()->speciesNames(), UnorderedElementsAreArray(expected_species));

}

TEST(PhaseNodeCreation, canConstructCanteraSolutionFromSpecies) {
    Goddard::setup_defaults();
    std::vector<std::string> species{"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2"};

    AnyMap root_node = Goddard::load_root_node("h2o2.yaml");
    AnyMap phase_node = Goddard::create_phase_node("Phase0", species);

    auto sln = Cantera::newSolution(phase_node, root_node);
    ASSERT_THAT(sln, NotNull()) << "Solution object should be initialized.";
    EXPECT_STREQ(sln->name().c_str(), "Phase0");

    EXPECT_THAT(sln->thermo()->speciesNames(), UnorderedElementsAreArray(species));
}

TEST(PhaseNodeCreation, canConstructCanteraSolutionFromElements) {
    Goddard::setup_defaults();
    vector<string> elements{"H", "O"};

    AnyMap root_node = Goddard::load_root_node("h2o2.yaml");
    AnyMap phase_node = Goddard::create_speciated_phase_node("Phase0", elements);

    auto sln = Cantera::newSolution(phase_node, root_node);
    ASSERT_THAT(sln, NotNull()) << "Solution object should be initialized.";
    EXPECT_STREQ(sln->name().c_str(), "Phase0");

    vector<string> expected_species = {"H2","H","O", "H2O2", "O2", "OH", "H2O", "HO2"};

    EXPECT_THAT(sln->thermo()->speciesNames(), UnorderedElementsAreArray(expected_species));
}