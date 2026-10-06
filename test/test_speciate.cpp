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

namespace {

/** Species names of a root node's `species` section, in file order. */
vector<string> species_names_of(const AnyMap& root_node) {
    vector<string> names;
    for (const AnyMap& species : root_node.at("species").asVector<AnyMap>()) {
        names.push_back(species.at("name").asString());
    }
    return names;
}

} // namespace

TEST(SelectSpecies, keepsExactlyTheRequestedSpeciesInFileOrder) {
    // h2o2.yaml lists H2, H, O, O2, OH, H2O, HO2, H2O2, AR, N2. Selecting a non-contiguous subset
    // removes neighbouring entries, which exercises the erase loop.
    Goddard::setup_defaults();
    unordered_set<string> requested{"H2O", "O2", "H2", "OH", "N2"};

    AnyMap node = Goddard::select_species("h2o2.yaml", requested);

    vector<string> expected;
    for (const string& name : species_names_of(Goddard::load_root_node("h2o2.yaml"))) {
        if (requested.count(name)) {
            expected.push_back(name);
        }
    }
    ASSERT_EQ(expected.size(), requested.size());
    EXPECT_EQ(species_names_of(node), expected);
}

TEST(SelectSpecies, unknownSpeciesThrowsAndNamesThem) {
    // A misspelled species must not silently disappear from the phase.
    Goddard::setup_defaults();
    unordered_set<string> requested{"H2", "ZZ_MISSING", "O2", "AA_MISSING"};

    try {
        Goddard::select_species("h2o2.yaml", requested);
        FAIL() << "select_species must throw for species absent from the file";
    } catch (const std::invalid_argument& error) {
        // Missing names are listed sorted, and present ones are not listed.
        EXPECT_THAT(error.what(), ::testing::HasSubstr("AA_MISSING ZZ_MISSING"));
        EXPECT_THAT(error.what(), ::testing::Not(::testing::HasSubstr("H2")));
    }
}

TEST(SelectSpecies, ctiAndXmlFilesAreRejected) {
    // The extension check must fire before the file is looked up, so the error names the format
    // rather than reporting a missing file.
    Goddard::setup_defaults();
    for (const string& infile : {string("h2o2.cti"), string("h2o2.XML")}) {
        try {
            Goddard::select_species(infile, {"H2"});
            FAIL() << infile << " must be rejected";
        } catch (const Cantera::CanteraError& error) {
            EXPECT_THAT(error.what(), ::testing::HasSubstr("CTI and XML formats"))
                << infile;
        }
    }
}
