cmake_minimum_required(VERSION 3.14)

add_executable(
    Goddard
    main.cpp
)

# target_compile_features(goddard_main PUBLIC CMAKE_CXX_STANDARD)

target_include_directories(Goddard
    PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>

)


target_link_libraries(Goddard PUBLIC cantera)