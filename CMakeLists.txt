cmake_minimum_required(VERSION 3.25)
project(hwyMP VERSION 0.1.0 LANGUAGES C CXX)

option(HWYMP_BUILD_TESTS "Build the unit tests" ${PROJECT_IS_TOP_LEVEL})
option(HWYMP_USE_LOCAL_DEPS
       "Use sibling highway/ and eigen/ checkouts when they exist" ON)
option(HWYMP_BUILD_BENCHMARKS "Build benchmarks; needs GMP, MPFR and Boost" OFF)
option(HWYMP_ENABLE_ORACLE
       "Build tests that need GMP and MPFR; fail if they are missing" OFF)

include(FetchContent)

set(HWYMP_HIGHWAY_TAG 1.4.0)
set(HWYMP_EIGEN_TAG 5.0.1)
set(HWYMP_GOOGLETEST_TAG v1.17.0)

if(HWYMP_USE_LOCAL_DEPS)
  get_filename_component(hwymp_parent_dir "${PROJECT_SOURCE_DIR}" DIRECTORY)
  if(EXISTS "${hwymp_parent_dir}/highway/hwy/highway.h")
    set(FETCHCONTENT_SOURCE_DIR_HIGHWAY "${hwymp_parent_dir}/highway")
  endif()
  if(EXISTS "${hwymp_parent_dir}/eigen/Eigen/Core")
    set(FETCHCONTENT_SOURCE_DIR_EIGEN "${hwymp_parent_dir}/eigen")
  endif()
endif()

FetchContent_Declare(highway
  GIT_REPOSITORY https://github.com/google/highway.git
  GIT_TAG ${HWYMP_HIGHWAY_TAG}
  SYSTEM
  GIT_SHALLOW TRUE)
set(HWY_ENABLE_CONTRIB OFF CACHE BOOL "" FORCE)
set(HWY_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
set(HWY_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
set(HWY_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(highway)

add_library(hwymp INTERFACE)
add_library(hwymp::hwymp ALIAS hwymp)
target_include_directories(hwymp INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_compile_features(hwymp INTERFACE cxx_std_20)
target_link_libraries(hwymp INTERFACE hwy)

if(HWYMP_BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

if(HWYMP_BUILD_BENCHMARKS)
  add_subdirectory(bench)
endif()
