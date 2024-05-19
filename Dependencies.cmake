include(cmake/CPM.cmake)
include(GNUInstallDirs)
include(cmake/CondaAware.cmake)


# Done as a function so that updates to variables don't 
# propagate out to other targets


set(GODDARD_USE_cantera     "" CACHE PATH "Specify this option in case a specific cantera library should be used.")
set(GODDARD_USE_Threads     "" CACHE PATH "Specify this option in case a specific Threads library should be used.")

function(GoddardFindPackage name)
    if(DEFINED GODDARD_USE_${ARGV0} AND NOT GODDARD_USE_${ARGV0} STREQUAL "")
        find_package(${name} ${ARGN} QUIET PATHS ${GODDARD_USE_${name}} NO_DEFAULT_PATH)
    else()
        find_package(${name} ${ARGN} QUIET)
    endif()
    if(${name} MATCHES "Python" AND ${name}_FOUND)
        set(${name}_DIR ${Python_EXECUTABLE})  # Python_EXECUTABLE is the path to the python interpreter and not Python_DIR as usual
    endif()
    if(${name}_FOUND)
        message(STATUS "Found ${name}: ${${name}_DIR} (found version \"${${name}_VERSION}\")")
    endif()
    set(${name}_FOUND ${${name}_FOUND} PARENT_SCOPE)  # export to parent (e.g., Eigen3_FOUND, Python_FOUND)
    set(${name}_DIR ${${name}_DIR} PARENT_SCOPE)  # export to parent (e.g., Eigen3_DIR, Python_DIR)
endfunction()



function(goddard_setup_dependencies)


    GoddardFindPackage(Threads REQUIRED)
    


endfunction()