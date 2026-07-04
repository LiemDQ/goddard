cmake_minimum_required(VERSION 3.16)

#library
add_library(
  goddard_lib
  SHARED
  equilibrium.cpp
  thermoarray.cpp
  thermo.cpp
  gas.cpp
  combustor.cpp
  error.cpp
  utils.cpp
  global.cpp
  speciate.cpp
  nozzle.cpp
  gas_dynamics.cpp
  numerics.cpp
  problem.cpp
  format.cpp
  rocket_results.cpp
  prandtlmeyer.cpp
  characteristics.cpp
  characteristic_net.cpp
  moc.cpp
  moc_initialization.cpp
  moc_nozzle.cpp
  shocks.cpp
  profile.cpp
  kinetic_nozzle.cpp
  interpolation.cpp
)

target_include_directories(goddard_lib
    PUBLIC 
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>

    PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
)

target_link_libraries(goddard_lib PUBLIC cantera Eigen3::Eigen)
target_link_libraries(goddard_lib PUBLIC goddard_options)
target_link_libraries(goddard_lib PRIVATE goddard_warnings)

