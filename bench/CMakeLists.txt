# Benchmarks compare hwyMP with MPFR and Boost.Multiprecision. Boost is
# header-only here, so only its include directory is needed.
get_filename_component(hwymp_parent_dir "${PROJECT_SOURCE_DIR}" DIRECTORY)
set(HWYMP_BOOST_DIR "${hwymp_parent_dir}/boost_1_92_0" CACHE PATH
    "Boost source tree used by the benchmarks")
if(NOT EXISTS "${HWYMP_BOOST_DIR}/boost/multiprecision/cpp_bin_float.hpp")
  message(FATAL_ERROR "Boost not found; set HWYMP_BOOST_DIR")
endif()

find_path(HWYMP_GMP_INCLUDE_DIR gmp.h REQUIRED)
find_path(HWYMP_MPFR_INCLUDE_DIR mpfr.h REQUIRED)
find_library(HWYMP_GMP_LIBRARY NAMES gmp libgmp REQUIRED)
find_library(HWYMP_MPFR_LIBRARY NAMES mpfr libmpfr REQUIRED)

add_executable(bench_scalar bench_scalar.cpp)
target_link_libraries(bench_scalar PRIVATE hwymp
  "${HWYMP_MPFR_LIBRARY}" "${HWYMP_GMP_LIBRARY}")
target_include_directories(bench_scalar SYSTEM PRIVATE
  "${HWYMP_BOOST_DIR}" "${HWYMP_GMP_INCLUDE_DIR}" "${HWYMP_MPFR_INCLUDE_DIR}")
target_compile_options(bench_scalar PRIVATE
  $<$<CXX_COMPILER_ID:Clang,GNU,AppleClang>:-march=native>)

add_executable(bench_batch bench_batch.cpp)
target_link_libraries(bench_batch PRIVATE hwymp)
target_compile_options(bench_batch PRIVATE
  $<$<CXX_COMPILER_ID:Clang,GNU,AppleClang>:-march=native>)
