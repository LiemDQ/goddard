cmake_minimum_required(VERSION 3.16)

project(goddard_cpp_tests CXX)

GoddardFindPackage(GTest 1.1.4 REQUIRED)

include(GoogleTest)

add_executable(
    goddardEquilibriumTests
    test_equilibrium.cpp
)

target_link_libraries(goddardEquilibriumTests PUBLIC goddard_options goddard_warnings)
target_link_libraries(goddardEquilibriumTests PUBLIC goddard_lib)
target_link_libraries(goddardEquilibriumTests PUBLIC GTest::gtest_main)
gtest_discover_tests(goddardEquilibriumTests TEST_PREFIX "equilibrium")
