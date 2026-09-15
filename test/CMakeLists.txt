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
    test_interpolation.cpp
    test_profile.cpp
    test_nozzle.cpp
    test_prandtl_meyer.cpp
    test_prandtl_meyer_table.cpp
    test_moc_unit.cpp
    test_moc_unit_processes.cpp
    test_moc_design.cpp
    test_moc_phases.cpp
    test_moc_convergence.cpp
    test_moc_initialization.cpp
    test_moc_rao.cpp
    test_moc_errors.cpp
    test_moc_invariance.cpp
    test_moc_inverse_march.cpp
    test_kinetic_nozzle.cpp
    test_shocks.cpp
    test_gas.cpp
    asan_options.cpp
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
# gtest_discover_tests lists tests by running the binary, so parameterized (TEST_P) instantiations are
# registered too; gtest_add_tests scans sources and misses them. The binary sets its own ASan options
# (asan_options.cpp) because discovery does not use the test ENVIRONMENT property.
gtest_discover_tests(goddardTests
    TEST_PREFIX "goddard"
    DISCOVERY_TIMEOUT 60
)
