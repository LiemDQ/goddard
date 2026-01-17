#include "cea_loader.hpp"

#include "goddard/config.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include <string>
#include <memory>
#include <vector>
#include <array>

double FLOAT_ABSTOL = 1e-5;
double FLOAT_RELTOL = 1e-5;

using ::testing::ElementsAre;
using ::testing::ElementsAreArray;
using ::testing::UnorderedElementsAreArray;

std::string CEA_directory() {
    std::string path(DATA_DIR);
    path += "/cea_results/";

    return path;
}

std::unique_ptr<CEAResult> load_CEA(const std::string& filename) {
    CEADataLoader loader{};
    std::string path = CEA_directory() + filename;

    return loader.load_from_file(path);
}

template<typename K, typename V>
std::vector<K> get_map_keys(const std::unordered_map<K, V>& map) {
    std::vector<K> keys;
    keys.reserve(map.size());
    for (const auto& kv: map) {
        keys.push_back(kv.first);
    }
    return keys;
}


template<typename K, typename V>
std::vector<V> get_map_values(const std::unordered_map<K, V>& map) {
    std::vector<V> keys(map.size());
    for (const auto& kv: map) {
        keys.push_back(kv.second);
    }
    return keys;
}

TEST(CEALoading, JSONRead) {
    std::string file = "h2.json";

    ASSERT_NO_THROW(load_CEA(file));
    
    std::unique_ptr<CEAResult> result = load_CEA(file);

    ASSERT_TRUE(result) << "If nullptr, loading has failed.";
    CEAConditions conditions = result->conditions;
    EXPECT_EQ(conditions.fuel, "H2(L)");
    EXPECT_EQ(conditions.oxidizer, "O2(L)");
    EXPECT_DOUBLE_EQ(conditions.fuel_temp, 20.27);
    EXPECT_DOUBLE_EQ(conditions.of_ratio, 6.0);
}

TEST(CEALoading, equilibriumStates) {

    std::unique_ptr<CEAResult> result = load_CEA("h2.json");
    
    std::vector<CEAState> states = result->equilibrium_states;
    std::vector<CEAPerformance> performance = result->equilibrium_performance;

    EXPECT_EQ(states.size()-1, performance.size()); //chamber state does not have a performance value

    CEAState chamber = states[0];
    
    EXPECT_EQ(chamber.location, "CHAMBER");
    EXPECT_DOUBLE_EQ(chamber.pressure_ratio, 1.0);
    EXPECT_DOUBLE_EQ(chamber.entropy_kj_kg_k, 17.1677);

    CEAComposition chamber_comp = chamber.mass_fractions;
    std::vector<std::string> products = get_map_keys(chamber_comp);
    // std::vector<double> chamber_fracs = get_map_values(chamber_comp);
    
    const std::vector<std::string> expected_products{"H", "H2", "O", "OH", "O2"};
    
    EXPECT_EQ(chamber_comp.size(), 5);
    EXPECT_THAT(products, UnorderedElementsAreArray(expected_products));
    
    const double cstar = performance[0].cstar_m_s;

    for (std::size_t i = 0; i < performance.size(); i++){
        std::size_t eq_i = i + 1;
        CEAState state = states[eq_i];
        
        CEAComposition comp = state.mass_fractions;
        EXPECT_DOUBLE_EQ(state.entropy_kj_kg_k, chamber.entropy_kj_kg_k);
        EXPECT_THAT(get_map_keys(comp), UnorderedElementsAreArray(products));

        CEAPerformance perf = performance[i];
        EXPECT_GT(perf.ivac_m_s, perf.isp_m_s);
        EXPECT_GT(perf.isp_m_s, 0.0);
        EXPECT_DOUBLE_EQ(perf.cstar_m_s, cstar); //c star is independent of expansion
    }
}


TEST(CEALoading, frozenStates) {

    std::unique_ptr<CEAResult> result = load_CEA("h2.json");
    
    std::vector<CEAState> states = result->frozen_states;
    std::vector<CEAPerformance> performance = result->frozen_performance;

    EXPECT_EQ(states.size()-1, performance.size()); //chamber state does not have a performance value

    CEAState chamber = states[0];
    
    EXPECT_EQ(chamber.location, "CHAMBER");
    EXPECT_DOUBLE_EQ(chamber.pressure_ratio, 1.0);
    EXPECT_DOUBLE_EQ(chamber.entropy_kj_kg_k, 17.1677);

    CEAComposition chamber_comp = chamber.mass_fractions;
    std::vector<std::string> products = get_map_keys(chamber_comp);
    // std::vector<double> chamber_fracs = get_map_values(chamber_comp);
    
    const std::vector<std::string> expected_products{"H", "H2", "O", "OH", "O2"};
    
    EXPECT_EQ(chamber_comp.size(), 5);
    EXPECT_THAT(products, UnorderedElementsAreArray(expected_products));
    
    const double cstar = performance[0].cstar_m_s;

    for (std::size_t i = 0; i < performance.size(); i++){
        std::size_t eq_i = i + 1;
        CEAState state = states[eq_i];
        
        CEAComposition comp = state.mass_fractions;
        EXPECT_DOUBLE_EQ(state.entropy_kj_kg_k, chamber.entropy_kj_kg_k);
        
        //frozen composition does not change from chamber composition, by definition
        EXPECT_EQ(comp.size(), 0);
        
        //these derivatives are not calculated for frozen systems
        EXPECT_EQ(state.dlv_dlp_t, 0.0); 
        EXPECT_EQ(state.dlv_dlt_p, 0.0);

        CEAPerformance perf = performance[i];
        EXPECT_GT(perf.ivac_m_s, perf.isp_m_s);
        EXPECT_GT(perf.isp_m_s, 0.0);
        EXPECT_DOUBLE_EQ(perf.cstar_m_s, cstar); //c star is independent of expansion
    }
}