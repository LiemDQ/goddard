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
    src/bind_equilibrium.cpp
    src/bind_thermo.cpp
    src/bind_moc.cpp
    src/bind_kinetic_nozzle.cpp
    src/bind_gas_properties.cpp
    src/bind_shocks.cpp
)

target_link_libraries(_core PRIVATE goddard_lib)
target_include_directories(_core PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/include
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
