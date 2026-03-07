execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")
find_package(nanobind CONFIG REQUIRED)

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
    src/bind_mixture_ratio.cpp
    src/bind_equilibrium.cpp
)

target_link_libraries(_core PRIVATE goddard_lib)
target_include_directories(_core PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/include
)

# Suppress strict warning flags inherited from goddard_lib's PUBLIC linkage
# to goddard_warnings. Python C API headers and nanobind-generated code
# trigger many warnings that are not actionable.
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

# Copy the built module into python/goddard/ for dev workflow (PYTHONPATH=python)
add_custom_command(TARGET _core POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
        $<TARGET_FILE:_core>
        ${PROJECT_SOURCE_DIR}/python/goddard/$<TARGET_FILE_NAME:_core>
)

install(TARGETS _core LIBRARY DESTINATION goddard)
