# Benchmarks compare mpfix with MPFR and Boost.Multiprecision. Boost is
# header-only here, so only its include directory is needed.
get_filename_component(mpfix_parent_dir "${PROJECT_SOURCE_DIR}" DIRECTORY)
set(MPFIX_BOOST_DIR "${mpfix_parent_dir}/boost_1_92_0" CACHE PATH
    "Boost source tree used by the benchmarks")
if(NOT EXISTS "${MPFIX_BOOST_DIR}/boost/multiprecision/cpp_bin_float.hpp")
  message(FATAL_ERROR "Boost not found; set MPFIX_BOOST_DIR")
endif()

find_path(MPFIX_GMP_INCLUDE_DIR gmp.h REQUIRED)
find_path(MPFIX_MPFR_INCLUDE_DIR mpfr.h REQUIRED)
find_library(MPFIX_GMP_LIBRARY NAMES gmp libgmp REQUIRED)
find_library(MPFIX_MPFR_LIBRARY NAMES mpfr libmpfr REQUIRED)

add_executable(bench_scalar bench_scalar.cpp)
target_link_libraries(bench_scalar PRIVATE mpfix
  "${MPFIX_MPFR_LIBRARY}" "${MPFIX_GMP_LIBRARY}")
target_include_directories(bench_scalar SYSTEM PRIVATE
  "${MPFIX_BOOST_DIR}" "${MPFIX_GMP_INCLUDE_DIR}" "${MPFIX_MPFR_INCLUDE_DIR}")

# Berkeley SoftFloat 3e, optionally, for binary128: a tree with
# source/include/softfloat.h and the library built from it, as a static library
# or, on Windows, a DLL with its import library in dll/.
set(MPFIX_SOFTFLOAT_DIR "${mpfix_parent_dir}/gmp-build/softfloat" CACHE PATH
    "Berkeley SoftFloat tree for the scalar benchmark, if it exists")
find_library(MPFIX_SOFTFLOAT_LIBRARY NAMES softfloat libsoftfloat softfloat.dll
  PATHS "${MPFIX_SOFTFLOAT_DIR}/dll" "${MPFIX_SOFTFLOAT_DIR}/build"
  NO_DEFAULT_PATH)
if(MPFIX_SOFTFLOAT_LIBRARY AND
   EXISTS "${MPFIX_SOFTFLOAT_DIR}/source/include/softfloat.h")
  target_include_directories(bench_scalar SYSTEM PRIVATE
    "${MPFIX_SOFTFLOAT_DIR}/source/include")
  target_link_libraries(bench_scalar PRIVATE "${MPFIX_SOFTFLOAT_LIBRARY}")
  # SoftFloat's default build, whose binary128 functions take values.
  target_compile_definitions(bench_scalar PRIVATE MPFIX_BENCH_SOFTFLOAT
    SOFTFLOAT_FAST_INT64
    $<$<BOOL:${WIN32}>:SOFTFLOAT_DLL>)
  message(STATUS "Benchmarking SoftFloat: ${MPFIX_SOFTFLOAT_LIBRARY}")
endif()
target_compile_options(bench_scalar PRIVATE
  $<$<CXX_COMPILER_ID:Clang,GNU,AppleClang>:-march=native>)

add_executable(bench_strings bench_strings.cpp)
target_link_libraries(bench_strings PRIVATE mpfix
  "${MPFIX_MPFR_LIBRARY}" "${MPFIX_GMP_LIBRARY}")
target_include_directories(bench_strings SYSTEM PRIVATE
  "${MPFIX_BOOST_DIR}" "${MPFIX_GMP_INCLUDE_DIR}" "${MPFIX_MPFR_INCLUDE_DIR}")
target_compile_options(bench_strings PRIVATE
  $<$<CXX_COMPILER_ID:Clang,GNU,AppleClang>:-march=native>)

add_executable(bench_math bench_math.cpp)
target_link_libraries(bench_math PRIVATE mpfix
  "${MPFIX_MPFR_LIBRARY}" "${MPFIX_GMP_LIBRARY}")
target_include_directories(bench_math SYSTEM PRIVATE
  "${MPFIX_BOOST_DIR}" "${MPFIX_GMP_INCLUDE_DIR}" "${MPFIX_MPFR_INCLUDE_DIR}")
target_compile_options(bench_math PRIVATE
  $<$<CXX_COMPILER_ID:Clang,GNU,AppleClang>:-march=native>)

# The batch prototype runs on Google Highway's vectors.
FetchContent_Declare(highway
  GIT_REPOSITORY https://github.com/google/highway.git
  GIT_TAG ${MPFIX_HIGHWAY_TAG}
  SYSTEM
  GIT_SHALLOW TRUE)
set(HWY_ENABLE_CONTRIB OFF CACHE BOOL "" FORCE)
set(HWY_ENABLE_EXAMPLES OFF CACHE BOOL "" FORCE)
set(HWY_ENABLE_INSTALL OFF CACHE BOOL "" FORCE)
set(HWY_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(highway)

add_executable(bench_batch bench_batch.cpp)
target_link_libraries(bench_batch PRIVATE mpfix hwy)
target_compile_options(bench_batch PRIVATE
  $<$<CXX_COMPILER_ID:Clang,GNU,AppleClang>:-march=native>)
