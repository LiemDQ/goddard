cmake_minimum_required(VERSION 3.16)

# executable
add_executable(
    goddard_main
    main.cpp
    equilibrium.cpp
)
target_include_directories(goddard_main
    PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>

)


target_link_libraries(goddard_main PUBLIC cantera Eigen3::Eigen)
target_link_libraries(goddard_main PUBLIC goddard_options goddard_warnings)

target_include_directories(
  goddard_main 
  PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)

#library
add_library(goddard_lib SHARED equilibrium.cpp)

target_include_directories(goddard_lib
    PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
)

target_link_libraries(goddard_lib PUBLIC cantera Eigen3::Eigen)
target_link_libraries(goddard_lib PUBLIC goddard_options goddard_warnings)

target_include_directories(
  goddard_lib 
  PUBLIC
  "$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>"
)