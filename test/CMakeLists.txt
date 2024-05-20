cmake_minimum_required(VERSION 3.16)

project(goddard_cpp_tests CXX)

GoddardFindPackage(GTest 1.1.4 REQUIRED)

include(GoogleTest)

add_executable(
    goddardThermoTests
    test_thermo.cpp
)

target_link_libraries(goddardThermoTests PUBLIC goddard_options goddard_warnings)
target_link_libraries(goddardThermoTests PUBLIC goddard_lib)
target_link_libraries(goddardThermoTests PUBLIC GTest::gtest_main)
gtest_discover_tests(goddardThermoTests TEST_PREFIX "thermo")
