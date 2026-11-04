cmake_minimum_required(VERSION 3.25)
project(mpfix VERSION 0.1.0 LANGUAGES C CXX)

option(MPFIX_BUILD_TESTS "Build the unit tests" ${PROJECT_IS_TOP_LEVEL})
option(MPFIX_USE_LOCAL_DEPS
       "Use sibling eigen/ and, for benchmarks, highway/ checkouts when they exist"
       ON)
option(MPFIX_BUILD_BENCHMARKS "Build benchmarks; needs GMP, MPFR and Boost" OFF)
option(MPFIX_ENABLE_ORACLE
       "Build tests that need GMP and MPFR; fail if they are missing" OFF)
option(MPFIX_BUILD_FUZZERS
       "Build the libFuzzer targets in fuzz/; needs Clang and MPFIX_ENABLE_ORACLE" OFF)
set(MPFIX_SANITIZE "" CACHE STRING
    "Sanitizers for the tests and fuzzers, as -fsanitize= lists them, such as address,undefined")

include(FetchContent)

set(MPFIX_HIGHWAY_TAG 1.4.0)
set(MPFIX_EIGEN_TAG 5.0.1)
set(MPFIX_GOOGLETEST_TAG v1.17.0)

if(MPFIX_USE_LOCAL_DEPS)
  get_filename_component(mpfix_parent_dir "${PROJECT_SOURCE_DIR}" DIRECTORY)
  if(EXISTS "${mpfix_parent_dir}/highway/hwy/highway.h")
    set(FETCHCONTENT_SOURCE_DIR_HIGHWAY "${mpfix_parent_dir}/highway")
  endif()
  if(EXISTS "${mpfix_parent_dir}/eigen/Eigen/Core")
    set(FETCHCONTENT_SOURCE_DIR_EIGEN "${mpfix_parent_dir}/eigen")
  endif()
endif()

# The library is header-only and has no dependencies. Highway is needed only by
# the batch benchmark, which bench/CMakeLists.txt fetches it for.
add_library(mpfix INTERFACE)
add_library(mpfix::mpfix ALIAS mpfix)
include(GNUInstallDirs)
target_include_directories(mpfix INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
target_compile_features(mpfix INTERFACE cxx_std_20)

# cmake --install puts the headers and a package configuration in place, for
# find_package(mpfix) and the target mpfix::mpfix.
option(MPFIX_INSTALL "Generate the install rules" ${PROJECT_IS_TOP_LEVEL})
if(MPFIX_INSTALL)
  include(CMakePackageConfigHelpers)
  set(mpfix_config_dir "${CMAKE_INSTALL_DATADIR}/cmake/mpfix")
  install(DIRECTORY include/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
  install(TARGETS mpfix EXPORT mpfixTargets)
  install(EXPORT mpfixTargets NAMESPACE mpfix::
          DESTINATION "${mpfix_config_dir}")
  configure_package_config_file(cmake/mpfixConfig.cmake.in
    "${PROJECT_BINARY_DIR}/mpfixConfig.cmake"
    INSTALL_DESTINATION "${mpfix_config_dir}")
  # Header-only: any architecture, and compatible within a major version.
  write_basic_package_version_file("${PROJECT_BINARY_DIR}/mpfixConfigVersion.cmake"
    COMPATIBILITY SameMajorVersion ARCH_INDEPENDENT)
  install(FILES "${PROJECT_BINARY_DIR}/mpfixConfig.cmake"
                "${PROJECT_BINARY_DIR}/mpfixConfigVersion.cmake"
          DESTINATION "${mpfix_config_dir}")
  install(FILES LICENSE DESTINATION "${CMAKE_INSTALL_DATADIR}/doc/mpfix")
endif()

if(MPFIX_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(MPFIX_BUILD_BENCHMARKS)
  add_subdirectory(bench)
endif()
