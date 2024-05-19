cmake_minimum_required(VERSION 3.21)


if (NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
endif()

set(CMAKE_CXX_EXTENSIONS OFF)

project(
    goddard
    VERSION 0.0.1
    LANGUAGES CXX C)

# find_package(Python 3.11 Interpreter Development.Module REQUIRED)

include(ProjectOptions.cmake)

goddard_setup_options()
goddard_global_options()
include(Dependencies.cmake)
goddard_setup_dependencies()
goddard_local_options()

target_compile_features(goddard_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})

add_library(goddard::goddard_options ALIAS goddard_options)
add_library(goddard::goddard_warnings ALIAS goddard_warnings)

# add_subdirectory(extern)
add_subdirectory(src)


if(NOT PROJECT_IS_TOP_LEVEL)
  return()
endif()

include(CTest)

if(BUILD_TESTING)
  message(AUTHOR_WARNING "Building tests.")
  add_subdirectory(test)
endif()

#python bindings
include_directories("${CMAKE_SOURCE_DIR}/goddard")

# add_subdirectory(${nanobind_SOURCE_DIR} ${nanobind_BINARY_DIR})

