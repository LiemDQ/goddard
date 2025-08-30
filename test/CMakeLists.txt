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
    cea_loader.cpp
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
gtest_add_tests(TARGET goddardTests 
    TEST_PREFIX "goddard"
    EXTRA_ARGS ""
    TEST_LIST test_list
)
set_tests_properties(${test_list} PROPERTIES ENVIRONMENT "ASAN_OPTIONS=detect_odr_violation=0")
