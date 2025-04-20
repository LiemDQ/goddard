cmake_minimum_required(VERSION 3.16)

project(goddard_cpp_tests CXX)

GoddardFindPackage(GTest 1.1.4 REQUIRED)

include(GoogleTest)

add_executable(
    goddardTests
    test_equilibrium.cpp
    test_thermoarray.cpp
)

target_link_libraries(goddardTests PUBLIC goddard_options goddard_warnings)
target_link_libraries(goddardTests PUBLIC goddard_lib)
target_link_libraries(goddardTests PUBLIC GTest::gtest_main)
gtest_discover_tests(goddardTests TEST_PREFIX "goddard")
