set(CMAKE_CXX_EXTENSIONS OFF)

# Sanitizers for everything built here, googletest included. Undefined
# behaviour stops the test rather than printing a warning that passes.
if(MPFIX_SANITIZE)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    message(FATAL_ERROR "MPFIX_SANITIZE is for GCC and Clang")
  endif()
  add_compile_options(-fsanitize=${MPFIX_SANITIZE} -fno-sanitize-recover=all
                      -fno-omit-frame-pointer -g)
  add_link_options(-fsanitize=${MPFIX_SANITIZE})
endif()

# Eigen is header-only. SOURCE_SUBDIR points at a directory that does not exist
# so that its own CMake project, with its tests and options, is not configured.
FetchContent_Declare(eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG ${MPFIX_EIGEN_TAG}
  GIT_SHALLOW TRUE
  SOURCE_SUBDIR do-not-configure)

FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG ${MPFIX_GOOGLETEST_TAG}
  GIT_SHALLOW TRUE)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(eigen googletest)

add_library(mpfix_eigen INTERFACE)
target_include_directories(mpfix_eigen SYSTEM INTERFACE "${eigen_SOURCE_DIR}")

# Boost.Multiprecision for the test of mpfix/boost.h: a sibling Boost tree, with
# Boost.Math, when there is one, otherwise the standalone Boost.Multiprecision.
get_filename_component(mpfix_parent_dir "${PROJECT_SOURCE_DIR}" DIRECTORY)
set(MPFIX_BOOST_DIR "${mpfix_parent_dir}/boost_1_92_0" CACHE PATH
    "Boost source tree used by the benchmarks and the Boost test")
add_library(mpfix_boost INTERFACE)
if(MPFIX_USE_LOCAL_DEPS AND EXISTS "${MPFIX_BOOST_DIR}/boost/multiprecision/number.hpp")
  target_include_directories(mpfix_boost SYSTEM INTERFACE "${MPFIX_BOOST_DIR}")
else()
  FetchContent_Declare(boost_multiprecision
    GIT_REPOSITORY https://github.com/boostorg/multiprecision.git
    GIT_TAG Boost_1_92_0
    GIT_SHALLOW TRUE
    SOURCE_SUBDIR do-not-configure)
  FetchContent_MakeAvailable(boost_multiprecision)
  target_include_directories(mpfix_boost SYSTEM INTERFACE
    "${boost_multiprecision_SOURCE_DIR}/include")
endif()
include(GoogleTest)

function(mpfix_add_test name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE mpfix GTest::gtest_main ${ARGN})
  mpfix_warnings(${name})
  gtest_discover_tests(${name} DISCOVERY_MODE PRE_TEST)
endfunction()

mpfix_add_test(limb_ops_test)
mpfix_add_test(float_test)
mpfix_add_test(convert_test)
mpfix_add_test(charconv_test)
mpfix_add_test(math_test)
mpfix_add_test(complex_test)
mpfix_add_test(std_complex_test mpfix_eigen)
mpfix_add_test(windows_header_test)
mpfix_add_test(deps_smoke_test mpfix_eigen)
mpfix_add_test(eigen_test mpfix_eigen)
mpfix_add_test(boost_test mpfix_boost)

# Berkeley TestFloat: with MPFIX_TESTFLOAT_GEN set to a testfloat_gen, a test
# for each function it generates cases for, in binary16, 32, 64 and 128. The
# checker is always built, so that it compiles everywhere.
set(MPFIX_TESTFLOAT_GEN "" CACHE FILEPATH
    "testfloat_gen of Berkeley TestFloat 3e, to check mpfix against")
set(MPFIX_TESTFLOAT_LEVEL 1 CACHE STRING "testfloat_gen's -level, 1 or 2")
add_executable(testfloat_check testfloat_check.cpp)
target_link_libraries(testfloat_check PRIVATE mpfix)
mpfix_warnings(testfloat_check)
if(MPFIX_TESTFLOAT_GEN)
  set(mpfix_testfloat_formats f16 f32 f64 f128)
  set(mpfix_testfloat_functions)
  foreach(format IN LISTS mpfix_testfloat_formats)
    foreach(op add sub mul div rem sqrt mulAdd roundToInt eq le lt
               eq_signaling le_quiet lt_quiet)
      list(APPEND mpfix_testfloat_functions ${format}_${op})
    endforeach()
    foreach(integer i32 i64 ui32 ui64)
      list(APPEND mpfix_testfloat_functions
           ${format}_to_${integer} ${integer}_to_${format})
    endforeach()
    foreach(other IN LISTS mpfix_testfloat_formats)
      if(NOT other STREQUAL format)
        list(APPEND mpfix_testfloat_functions ${format}_to_${other})
      endif()
    endforeach()
  endforeach()
  foreach(function IN LISTS mpfix_testfloat_functions)
    add_test(NAME testfloat.${function}
             COMMAND testfloat_check "${MPFIX_TESTFLOAT_GEN}" ${function}
                     ${MPFIX_TESTFLOAT_LEVEL})
  endforeach()
endif()

if(MPFIX_ENABLE_ORACLE)
  find_path(MPFIX_GMP_INCLUDE_DIR gmp.h REQUIRED)
  find_path(MPFIX_MPFR_INCLUDE_DIR mpfr.h REQUIRED)
  find_library(MPFIX_GMP_LIBRARY NAMES gmp libgmp REQUIRED)
  find_library(MPFIX_MPFR_LIBRARY NAMES mpfr libmpfr REQUIRED)

  add_library(mpfix_oracle INTERFACE)
  target_include_directories(mpfix_oracle SYSTEM INTERFACE
    "${MPFIX_GMP_INCLUDE_DIR}" "${MPFIX_MPFR_INCLUDE_DIR}")
  # MPFR must come first because it depends on GMP.
  target_link_libraries(mpfix_oracle INTERFACE
    "${MPFIX_MPFR_LIBRARY}" "${MPFIX_GMP_LIBRARY}")

  mpfix_add_test(oracle_smoke_test mpfix_oracle)
  mpfix_add_test(float_oracle_test mpfix_oracle)
  mpfix_add_test(convert_oracle_test mpfix_oracle)
  mpfix_add_test(math_oracle_test mpfix_oracle)
  mpfix_add_test(pow_oracle_test mpfix_oracle)
  mpfix_add_test(trig_oracle_test mpfix_oracle)
  mpfix_add_test(inverse_trig_oracle_test mpfix_oracle)
  mpfix_add_test(special_oracle_test mpfix_oracle)

  # MPC, for the complex functions, depends on MPFR and GMP.
  find_path(MPFIX_MPC_INCLUDE_DIR mpc.h REQUIRED)
  find_library(MPFIX_MPC_LIBRARY NAMES mpc libmpc REQUIRED)
  add_library(mpfix_mpc INTERFACE)
  target_include_directories(mpfix_mpc SYSTEM INTERFACE "${MPFIX_MPC_INCLUDE_DIR}")
  target_link_libraries(mpfix_mpc INTERFACE "${MPFIX_MPC_LIBRARY}" mpfix_oracle)
  mpfix_add_test(complex_oracle_test mpfix_mpc)

  # The fuzz targets on random inputs, with any compiler, as tests; and with
  # libFuzzer, for Clang.
  foreach(name fuzz_arithmetic fuzz_strings)
    add_executable(${name}_random "${PROJECT_SOURCE_DIR}/fuzz/${name}.cpp"
                   "${PROJECT_SOURCE_DIR}/fuzz/random_main.cpp")
    target_include_directories(${name}_random PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    target_link_libraries(${name}_random PRIVATE mpfix mpfix_oracle GTest::gtest)
    mpfix_warnings(${name}_random)
    add_test(NAME ${name}_random COMMAND ${name}_random 20000)
  endforeach()
  if(MPFIX_BUILD_FUZZERS)
    add_subdirectory("${PROJECT_SOURCE_DIR}/fuzz" fuzz)
  endif()
endif()
