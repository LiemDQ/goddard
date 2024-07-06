include(cmake/SystemLink.cmake)
include(cmake/LibFuzzer.cmake)
include(CMakeDependentOption)
include(CheckCXXCompilerFlag)


macro(goddard_supports_sanitizers)
  if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID MATCHES ".*GNU.*") AND NOT WIN32)
    set(SUPPORTS_UBSAN ON)
  else()
    set(SUPPORTS_UBSAN OFF)
  endif()

  if((CMAKE_CXX_COMPILER_ID MATCHES ".*Clang.*" OR CMAKE_CXX_COMPILER_ID MATCHES ".*GNU.*") AND WIN32)
    set(SUPPORTS_ASAN OFF)
  else()
    set(SUPPORTS_ASAN ON)
  endif()
endmacro()

macro(goddard_setup_options)
  option(goddard_ENABLE_HARDENING "Enable hardening" ON)
  option(goddard_ENABLE_COVERAGE "Enable coverage reporting" OFF)
  cmake_dependent_option(
    goddard_ENABLE_GLOBAL_HARDENING
    "Attempt to push hardening options to built dependencies"
    ON
    goddard_ENABLE_HARDENING
    OFF)

  goddard_supports_sanitizers()

  if(NOT PROJECT_IS_TOP_LEVEL OR goddard_PACKAGING_MAINTAINER_MODE)
    option(goddard_ENABLE_IPO "Enable IPO/LTO" OFF)
    option(goddard_WARNINGS_AS_ERRORS "Treat Warnings As Errors" OFF)
    option(goddard_ENABLE_USER_LINKER "Enable user-selected linker" OFF)
    option(goddard_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" OFF)
    option(goddard_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
    option(goddard_ENABLE_SANITIZER_UNDEFINED "Enable undefined sanitizer" OFF)
    option(goddard_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
    option(goddard_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)
    option(goddard_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
    option(goddard_ENABLE_CLANG_TIDY "Enable clang-tidy" OFF)
    option(goddard_ENABLE_CPPCHECK "Enable cpp-check analysis" OFF)
    option(goddard_ENABLE_PCH "Enable precompiled headers" OFF)
    option(goddard_ENABLE_CACHE "Enable ccache" OFF)
  else()
    option(goddard_ENABLE_IPO "Enable IPO/LTO" ON)
    option(goddard_WARNINGS_AS_ERRORS "Treat Warnings As Errors" ON)
    option(goddard_ENABLE_USER_LINKER "Enable user-selected linker" OFF)
    option(goddard_ENABLE_SANITIZER_ADDRESS "Enable address sanitizer" ${SUPPORTS_ASAN})
    option(goddard_ENABLE_SANITIZER_LEAK "Enable leak sanitizer" OFF)
    option(goddard_ENABLE_SANITIZER_UNDEFINED "Enable undefined sanitizer" ${SUPPORTS_UBSAN})
    option(goddard_ENABLE_SANITIZER_THREAD "Enable thread sanitizer" OFF)
    option(goddard_ENABLE_SANITIZER_MEMORY "Enable memory sanitizer" OFF)
    option(goddard_ENABLE_UNITY_BUILD "Enable unity builds" OFF)
    option(goddard_ENABLE_CLANG_TIDY "Enable clang-tidy" ON)
    option(goddard_ENABLE_CPPCHECK "Enable cpp-check analysis" ON)
    option(goddard_ENABLE_PCH "Enable precompiled headers" OFF)
    option(goddard_ENABLE_CACHE "Enable ccache" ON)
    option(PREFIX_SYSTEM_CANTERA "Use system-wide Cantera installation" OFF)
  endif()

  if(NOT PROJECT_IS_TOP_LEVEL)
    mark_as_advanced(
      goddard_ENABLE_IPO
      goddard_WARNINGS_AS_ERRORS
      goddard_ENABLE_USER_LINKER
      goddard_ENABLE_SANITIZER_ADDRESS
      goddard_ENABLE_SANITIZER_LEAK
      goddard_ENABLE_SANITIZER_UNDEFINED
      goddard_ENABLE_SANITIZER_THREAD
      goddard_ENABLE_SANITIZER_MEMORY
      goddard_ENABLE_UNITY_BUILD
      goddard_ENABLE_CLANG_TIDY
      goddard_ENABLE_CPPCHECK
      goddard_ENABLE_COVERAGE
      goddard_ENABLE_PCH
      goddard_ENABLE_CACHE)
  endif()

  goddard_check_libfuzzer_support(LIBFUZZER_SUPPORTED)
  if(LIBFUZZER_SUPPORTED AND (goddard_ENABLE_SANITIZER_ADDRESS OR goddard_ENABLE_SANITIZER_THREAD OR goddard_ENABLE_SANITIZER_UNDEFINED))
    set(DEFAULT_FUZZER ON)
  else()
    set(DEFAULT_FUZZER OFF)
  endif()

  option(goddard_BUILD_FUZZ_TESTS "Enable fuzz testing executable" ${DEFAULT_FUZZER})

endmacro()

macro(goddard_global_options)
  if(goddard_ENABLE_IPO)
    include(cmake/IPO.cmake)
    goddard_enable_ipo()
  endif()

  goddard_supports_sanitizers()

  if(goddard_ENABLE_HARDENING AND goddard_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN 
       OR goddard_ENABLE_SANITIZER_UNDEFINED
       OR goddard_ENABLE_SANITIZER_ADDRESS
       OR goddard_ENABLE_SANITIZER_THREAD
       OR goddard_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    message("${goddard_ENABLE_HARDENING} ${ENABLE_UBSAN_MINIMAL_RUNTIME} ${goddard_ENABLE_SANITIZER_UNDEFINED}")
    goddard_enable_hardening(goddard_options ON ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()
endmacro()

macro(goddard_local_options)
  if(PROJECT_IS_TOP_LEVEL)
    include(cmake/StandardProjectSettings.cmake)
  endif()

  add_library(goddard_warnings INTERFACE)
  add_library(goddard_options INTERFACE)

  include(cmake/CompilerWarnings.cmake)
  goddard_set_project_warnings(
    goddard_warnings
    ${goddard_WARNINGS_AS_ERRORS}
    ""
    ""
    ""
    "")

  if(goddard_ENABLE_USER_LINKER)
    include(cmake/Linker.cmake)
    configure_linker(goddard_options)
  endif()

  include(cmake/Sanitizers.cmake)
  goddard_enable_sanitizers(
    goddard_options
    ${goddard_ENABLE_SANITIZER_ADDRESS}
    ${goddard_ENABLE_SANITIZER_LEAK}
    ${goddard_ENABLE_SANITIZER_UNDEFINED}
    ${goddard_ENABLE_SANITIZER_THREAD}
    ${goddard_ENABLE_SANITIZER_MEMORY})

  set_target_properties(goddard_options PROPERTIES UNITY_BUILD ${goddard_ENABLE_UNITY_BUILD})

  if(goddard_ENABLE_PCH)
    target_precompile_headers(
      goddard_options
      INTERFACE
      <vector>
      <string>
      <utility>)
  endif()

  if(goddard_ENABLE_CACHE)
    include(cmake/Cache.cmake)
    goddard_enable_cache()
  endif()

  include(cmake/StaticAnalyzers.cmake)
  if(goddard_ENABLE_CLANG_TIDY)
    goddard_enable_clang_tidy(goddard_options ${goddard_WARNINGS_AS_ERRORS})
  endif()

  if(goddard_ENABLE_CPPCHECK)
    goddard_enable_cppcheck(${goddard_WARNINGS_AS_ERRORS} "" # override cppcheck options
    )
  endif()

  if(goddard_ENABLE_COVERAGE)
    include(cmake/Tests.cmake)
    goddard_enable_coverage(goddard_options)
  endif()

  if(goddard_WARNINGS_AS_ERRORS)
    check_cxx_compiler_flag("-Wl,--fatal-warnings" LINKER_FATAL_WARNINGS)
    if(LINKER_FATAL_WARNINGS)
      # This is not working consistently, so disabling for now
      # target_link_options(goddard_options INTERFACE -Wl,--fatal-warnings)
    endif()
  endif()

  if(goddard_ENABLE_HARDENING AND NOT goddard_ENABLE_GLOBAL_HARDENING)
    include(cmake/Hardening.cmake)
    if(NOT SUPPORTS_UBSAN 
       OR goddard_ENABLE_SANITIZER_UNDEFINED
       OR goddard_ENABLE_SANITIZER_ADDRESS
       OR goddard_ENABLE_SANITIZER_THREAD
       OR goddard_ENABLE_SANITIZER_LEAK)
      set(ENABLE_UBSAN_MINIMAL_RUNTIME FALSE)
    else()
      set(ENABLE_UBSAN_MINIMAL_RUNTIME TRUE)
    endif()
    goddard_enable_hardening(goddard_options OFF ${ENABLE_UBSAN_MINIMAL_RUNTIME})
  endif()

endmacro()