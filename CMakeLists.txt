cmake_minimum_required(VERSION 3.21)


if (NOT DEFINED CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 20)
endif()

set(CMAKE_CXX_EXTENSIONS OFF)

project(
    goddard
    VERSION 0.1.0
    LANGUAGES CXX C)

option(goddard_BUILD_PYTHON "Build Python bindings" ON)
option(goddard_BUILD_TOOLS "Build diagnostic tools in tools/" OFF)
option(goddard_BUILD_DOCSTRINGS
  "Extract C++ docstrings into the Python bindings via pybind11_mkdoc and \
generate nanobind stubs (requires the docs pixi environment)" OFF)
if(goddard_BUILD_PYTHON)
  find_package(Python 3.11 COMPONENTS Interpreter Development.Module REQUIRED)
endif()

include(ProjectOptions.cmake)

goddard_setup_options()
goddard_global_options()
include(Dependencies.cmake)
goddard_setup_dependencies()
goddard_local_options()

# Define where data files will be installed
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "Install path" FORCE)
endif()
# Python packages (scikit-build-core sets SKBUILD) install the data next to the extension module, at
# a path relative to the wheel root; goddard/__init__.py registers it at import. Other installs use
# <prefix>/share/goddard/data.
if(SKBUILD)
  set(DATA_INSTALL_DIR "goddard/data")
else()
  set(DATA_INSTALL_DIR "${CMAKE_INSTALL_DATADIR}/goddard/data")
endif()

# DATA_DIR is the data directory compiled into the library and registered by setup_defaults() when
# it exists. Development builds and editable installs use the source tree. A wheel must not bake in
# a path from the build machine, so it gets none and relies on goddard/__init__.py instead.
if(goddard_DEVELOPMENT_BUILD OR SKBUILD_STATE STREQUAL "editable")
  set(DATA_DIR "${CMAKE_SOURCE_DIR}/data")
elseif(SKBUILD)
  set(DATA_DIR "")
else()
  set(DATA_DIR "${CMAKE_INSTALL_PREFIX}/${DATA_INSTALL_DIR}")
endif()

configure_file(
  "${PROJECT_SOURCE_DIR}/include/goddard/config.h.in"
  "${PROJECT_BINARY_DIR}/include/goddard/config.h"
  @ONLY
)

target_compile_features(goddard_options INTERFACE cxx_std_${CMAKE_CXX_STANDARD})
add_library(goddard::goddard_options ALIAS goddard_options)
add_library(goddard::goddard_warnings ALIAS goddard_warnings)

add_subdirectory(src)


if(NOT PROJECT_IS_TOP_LEVEL)
  return()
endif()

include(CTest)

if(BUILD_TESTING)
  message(AUTHOR_WARNING "Building tests.")
  add_subdirectory(test)
  # install(DIRECTORY "${PROJECT_SOURCE_DIR}/data/cea_results"
  #         DESTINATION "${DATA_INSTALL_DIR}"
  #         FILES_MATCHING PATTERN "*.json"
  # )
endif()

if(goddard_BUILD_TOOLS)
  add_subdirectory(tools)
endif()

# Python bindings
if(goddard_BUILD_PYTHON)
  add_subdirectory(python)
endif()


# Install the thermodynamic data files. cea_results/ holds test reference data only.
install(DIRECTORY "${PROJECT_SOURCE_DIR}/data/"
        DESTINATION "${DATA_INSTALL_DIR}"
        FILES_MATCHING PATTERN "*.yaml"
        PATTERN "cea_results" EXCLUDE
)

