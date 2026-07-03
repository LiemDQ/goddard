cmake_minimum_required(VERSION 3.16)

project(goddard_cpp_tests CXX)

GoddardFindPackage(GTest 1.1.4 REQUIRED)
GoddardFindPackage(nlohmann_json 3.11.0 REQUIRED)
include(GoogleTest)

add_executable(
    goddardTests
    test_equilibrium.cpp
    test_thermoarray.cpp
    test_combustor.cpp
    test_cea_loader.cpp
    test_cea_integration.cpp
    cea_loader.cpp
    cea_comparer.cpp
    test_speciate.cpp
    test_nozzle.cpp
    test_prandtl_meyer.cpp
    test_prandtl_meyer_table.cpp
    test_moc_unit.cpp
    test_moc_design.cpp
    test_moc_phases.cpp
    test_moc_convergence.cpp
    test_moc_initialization.cpp
    test_kinetic_nozzle.cpp
    test_shocks.cpp
    test_gas.cpp
)

target_include_directories(goddardTests      
    PRIVATE 
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/test>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/test>
)

target_link_libraries(goddardTests PRIVATE goddard_options goddard_warnings)
target_link_libraries(goddardTests PRIVATE goddard_lib)
target_link_libraries(goddardTests PRIVATE GTest::gtest_main)
target_link_libraries(goddardTests PRIVATE GTest::gmock)
target_link_libraries(goddardTests PRIVATE nlohmann_json::nlohmann_json)
gtest_add_tests(TARGET goddardTests 
    TEST_PREFIX "goddard"
    EXTRA_ARGS ""
    TEST_LIST test_list
)
set_tests_properties(${test_list} PROPERTIES ENVIRONMENT "ASAN_OPTIONS=detect_odr_violation=0")
