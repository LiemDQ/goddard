
# Check if a conda environment is active
if(DEFINED ENV{CONDA_PREFIX})
    # Show that conda env has been recognized
    message(STATUS "CondaAware: Conda environment detected!")
    message(STATUS "CondaAware: Found environment variable CONDA_PREFIX=$ENV{CONDA_PREFIX}")

    # Check if environment variable PYTHON is defined, and if so, set Python_EXECUTABLE to PYTHON
    if(DEFINED ENV{PYTHON})
        message(STATUS "CondaAware: Found environment variable PYTHON=$ENV{PYTHON}")
        message(STATUS "CondaAware: Setting Python_EXECUTABLE to $ENV{PYTHON}")
        set(Python_EXECUTABLE $ENV{PYTHON})
    endif()

    # Ensure the Python executable is correctly selected, either through the use
    # of CMake variables Python_EXECUTABLE or an environment variable PYTHON.
    if(NOT DEFINED Python_EXECUTABLE)
        if(UNIX)
            message(STATUS "CondaAware: Setting Python_EXECUTABLE=$ENV{CONDA_PREFIX}/bin/python")
            set(Python_EXECUTABLE "$ENV{CONDA_PREFIX}/bin/python")
        endif()

        if(WIN32)
            message(STATUS "CondaAware: Setting Python_EXECUTABLE=$ENV{CONDA_PREFIX}\\python.exe")
            set(Python_EXECUTABLE "$ENV{CONDA_PREFIX}\\python.exe")
        endif()

        if(NOT DEFINED Python_EXECUTABLE)
            message(FATAL_ERROR "CondaAware: Could not determine a value for Python_EXECUTABLE. Expecting Unix or Windows systems.")
        endif()
    endif()

    # Set auxiliary variable CONDA_AWARE_PREFIX according to the logic below.
    if(DEFINED ENV{CONDA_BUILD})
        message(STATUS "CondaAware: Detected conda build task (e.g., in a conda-forge build)!")
        if(UNIX)
            set(CONDA_AWARE_PREFIX "$ENV{PREFIX}")
        endif()
        if(WIN32)
            set(CONDA_AWARE_PREFIX "$ENV{LIBRARY_PREFIX}")
        endif()
    else()
        if(UNIX)
            set(CONDA_AWARE_PREFIX "$ENV{CONDA_PREFIX}")
        endif()

        if(WIN32)
            set(CONDA_AWARE_PREFIX "$ENV{CONDA_PREFIX}\\Library")
        endif()
    endif()

    # Check if CONDA_AWARE_PREFIX has been successfully set
    if(DEFINED CONDA_AWARE_PREFIX)
        message(STATUS "CondaAware: Setting CONDA_AWARE_PREFIX=${CONDA_AWARE_PREFIX}")
    else()
        message(FATAL_ERROR "CondaAware: Could not determine a value for CONDA_AWARE_PREFIX. Expecting Unix or Windows systems.")
    endif()

    # Set CMAKE_INSTALL_PREFIX to CONDA_AWARE_PREFIX if not specified by the user
    if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        message(STATUS "CondaAware: Setting CMAKE_INSTALL_PREFIX=CONDA_AWARE_PREFIX=${CONDA_AWARE_PREFIX}")
        set(CMAKE_INSTALL_PREFIX ${CONDA_AWARE_PREFIX})
    endif()

    # Ensure dependencies from the conda environment are used instead of those from the system.
    list(APPEND CMAKE_PREFIX_PATH ${CONDA_AWARE_PREFIX})
    message(STATUS "CondaAware: Appended ${CONDA_AWARE_PREFIX} to CMAKE_PREFIX_PATH")

    # Ensure include directory in conda environment is known to the project
    include_directories(${CONDA_AWARE_PREFIX}/include)
    message(STATUS "CondaAware: Appended ${CONDA_AWARE_PREFIX}/include to include directories")

    # Ensure library directory in conda environment is known to the project
    link_directories(${CONDA_AWARE_PREFIX}/lib)
    message(STATUS "CondaAware: Appended ${CONDA_AWARE_PREFIX}/lib to link directories")
endif()