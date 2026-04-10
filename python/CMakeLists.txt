execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

# ---------------------------------------------------------------------------
# Generated docstrings header
#
# scripts/generate_docstrings.py either runs pybind11_mkdoc against the public
# C++ headers (when goddard_BUILD_DOCSTRINGS=ON) or emits a stub that defines
# DOC(...) as an empty string. Either way, the bindings can `#include
# "goddard_docstrings.h"` and use DOC(Goddard, Class, method) unconditionally.
# ---------------------------------------------------------------------------
set(GODDARD_DOCSTRINGS_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/src/goddard_docstrings.h")

# Public headers to scan. Keep this list explicit so we can exclude headers
# that confuse libclang (e.g. heavy template metaprogramming).
set(GODDARD_PUBLIC_HEADERS
    "${PROJECT_SOURCE_DIR}/include/goddard/gas.hpp"
    "${PROJECT_SOURCE_DIR}/include/goddard/equilibrium.hpp"
    "${PROJECT_SOURCE_DIR}/include/goddard/combustor.hpp"
    "${PROJECT_SOURCE_DIR}/include/goddard/nozzle.hpp"
    "${PROJECT_SOURCE_DIR}/include/goddard/kinetic_nozzle.hpp"
    "${PROJECT_SOURCE_DIR}/include/goddard/moc.hpp"
    "${PROJECT_SOURCE_DIR}/include/goddard/characteristics.hpp"
    "${PROJECT_SOURCE_DIR}/include/goddard/prandtlmeyer.hpp"
    "${PROJECT_SOURCE_DIR}/include/goddard/profile.hpp"
    "${PROJECT_SOURCE_DIR}/include/goddard/problem.hpp"
    "${PROJECT_SOURCE_DIR}/include/goddard/rocket_results.hpp"
    "${PROJECT_SOURCE_DIR}/include/goddard/shocks.hpp"
    "${PROJECT_SOURCE_DIR}/include/goddard/thermoarray.hpp"
    "${PROJECT_SOURCE_DIR}/include/goddard/thermo.hpp"
)

set(_mkdoc_args "--output" "${GODDARD_DOCSTRINGS_HEADER}")
if(goddard_BUILD_DOCSTRINGS)
    foreach(_h IN LISTS GODDARD_PUBLIC_HEADERS)
        list(APPEND _mkdoc_args "--header" "${_h}")
    endforeach()
    list(APPEND _mkdoc_args
        "-I" "${PROJECT_SOURCE_DIR}/include"
        "-I" "${PROJECT_BINARY_DIR}/include"
    )
    # Forward Cantera + Eigen include dirs from goddard_lib so libclang can
    # resolve transitive includes. We query the target property at configure
    # time (goddard_lib was created earlier in add_subdirectory(src)).
    get_target_property(_goddard_lib_inc goddard_lib INTERFACE_INCLUDE_DIRECTORIES)
    if(_goddard_lib_inc)
        foreach(_inc IN LISTS _goddard_lib_inc)
            # Strip generator expressions; libclang can't evaluate them and
            # the include dir resolution should work without them.
            if(NOT _inc MATCHES "^\\$<")
                list(APPEND _mkdoc_args "-I" "${_inc}")
            endif()
        endforeach()
    endif()
else()
    list(APPEND _mkdoc_args "--stub")
endif()

add_custom_command(
    OUTPUT "${GODDARD_DOCSTRINGS_HEADER}"
    COMMAND "${Python_EXECUTABLE}"
            "${PROJECT_SOURCE_DIR}/scripts/generate_docstrings.py"
            ${_mkdoc_args}
    DEPENDS
        "${PROJECT_SOURCE_DIR}/scripts/generate_docstrings.py"
        ${GODDARD_PUBLIC_HEADERS}
    COMMENT "Generating goddard_docstrings.h (BUILD_DOCSTRINGS=${goddard_BUILD_DOCSTRINGS})"
    VERBATIM
)

add_custom_target(goddard_docstrings_header
    DEPENDS "${GODDARD_DOCSTRINGS_HEADER}"
)

nanobind_add_module(
    _core
    NB_STATIC
    src/bind_main.cpp
    src/bind_enums.cpp
    src/bind_errors.cpp
    src/bind_structs.cpp
    src/bind_problem.cpp
    src/bind_thermoarray.cpp
    src/bind_combustor.cpp
    src/bind_nozzle.cpp
    src/bind_equilibrium.cpp
    src/bind_thermo.cpp
    src/bind_moc.cpp
    src/bind_kinetic_nozzle.cpp
    src/bind_gas_properties.cpp
    src/bind_shocks.cpp
    "${GODDARD_DOCSTRINGS_HEADER}"
)
add_dependencies(_core goddard_docstrings_header)

target_link_libraries(_core PRIVATE goddard_lib)
target_include_directories(_core PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# Embed RPATH so the extension can find libgoddard_lib.so at runtime
# during development (build tree) and after install
set_target_properties(_core PROPERTIES
    BUILD_RPATH "$<TARGET_FILE_DIR:goddard_lib>"
    INSTALL_RPATH "$ORIGIN"
)

# Suppress strict warning flags inherited from goddard_lib's PUBLIC linkage.
# Python C API headers and nanobind-generated code trigger many warnings
# that are not actionable.
target_compile_options(_core PRIVATE
    -Wno-error
    -Wno-old-style-cast
    -Wno-cast-align
    -Wno-conversion
    -Wno-sign-conversion
    -Wno-double-promotion
    -Wno-useless-cast
    -Wno-suggest-override
    -Wno-shadow
    -Wno-pedantic
    -Wno-non-virtual-dtor
    -Wno-missing-field-initializers
)

# Disable sanitizers for the Python extension module. ASan requires its
# runtime to be the first loaded library, which is not possible when Python
# loads extension modules dynamically.
target_compile_options(_core PRIVATE -fno-sanitize=all)
target_link_options(_core PRIVATE -fno-sanitize=all)

# Copy the built module into python/goddard/ for dev workflow (PYTHONPATH=python)
add_custom_command(TARGET _core POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:_core>
        ${PROJECT_SOURCE_DIR}/python/goddard/$<TARGET_FILE_NAME:_core>
)

# Install both the extension and the shared library into the same directory
# so that $ORIGIN RPATH resolves libgoddard_lib.so at runtime.
install(TARGETS _core goddard_lib LIBRARY DESTINATION goddard)
