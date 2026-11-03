set(CMAKE_CXX_EXTENSIONS OFF)

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

include(GoogleTest)

function(mpfix_add_test name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE mpfix GTest::gtest_main ${ARGN})
  if(MSVC AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(${name} PRIVATE /W4 /permissive- /Zc:__cplusplus /bigobj)
  elseif(MSVC)
    target_compile_options(${name} PRIVATE /W4 /bigobj)
  else()
    target_compile_options(${name} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  if(WIN32)
    target_compile_definitions(${name} PRIVATE _CRT_SECURE_NO_WARNINGS)
  endif()
  gtest_discover_tests(${name} DISCOVERY_MODE PRE_TEST)
endfunction()

mpfix_add_test(limb_ops_test)
mpfix_add_test(float_test)
mpfix_add_test(convert_test)
mpfix_add_test(charconv_test)
mpfix_add_test(math_test)
mpfix_add_test(complex_test)
mpfix_add_test(windows_header_test)
mpfix_add_test(deps_smoke_test mpfix_eigen)
mpfix_add_test(eigen_test mpfix_eigen)

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
endif()
