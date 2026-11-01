set(CMAKE_CXX_EXTENSIONS OFF)

# Eigen is header-only. SOURCE_SUBDIR points at a directory that does not exist
# so that its own CMake project, with its tests and options, is not configured.
FetchContent_Declare(eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG ${HWYMP_EIGEN_TAG}
  GIT_SHALLOW TRUE
  SOURCE_SUBDIR do-not-configure)

FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG ${HWYMP_GOOGLETEST_TAG}
  GIT_SHALLOW TRUE)
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(eigen googletest)

add_library(hwymp_eigen INTERFACE)
target_include_directories(hwymp_eigen SYSTEM INTERFACE "${eigen_SOURCE_DIR}")

include(GoogleTest)

function(hwymp_add_test name)
  add_executable(${name} ${name}.cpp)
  target_link_libraries(${name} PRIVATE hwymp GTest::gtest_main ${ARGN})
  if(MSVC AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(${name} PRIVATE /W4 /permissive- /Zc:__cplusplus)
  elseif(MSVC)
    target_compile_options(${name} PRIVATE /W4)
  else()
    target_compile_options(${name} PRIVATE -Wall -Wextra -Wpedantic)
  endif()
  gtest_discover_tests(${name} DISCOVERY_MODE PRE_TEST)
endfunction()

hwymp_add_test(limb_ops_test)
hwymp_add_test(deps_smoke_test hwymp_eigen)

if(HWYMP_ENABLE_ORACLE)
  find_path(HWYMP_GMP_INCLUDE_DIR gmp.h REQUIRED)
  find_path(HWYMP_MPFR_INCLUDE_DIR mpfr.h REQUIRED)
  find_library(HWYMP_GMP_LIBRARY NAMES gmp libgmp REQUIRED)
  find_library(HWYMP_MPFR_LIBRARY NAMES mpfr libmpfr REQUIRED)

  add_library(hwymp_oracle INTERFACE)
  target_include_directories(hwymp_oracle SYSTEM INTERFACE
    "${HWYMP_GMP_INCLUDE_DIR}" "${HWYMP_MPFR_INCLUDE_DIR}")
  # MPFR must come first because it depends on GMP.
  target_link_libraries(hwymp_oracle INTERFACE
    "${HWYMP_MPFR_LIBRARY}" "${HWYMP_GMP_LIBRARY}")

  hwymp_add_test(oracle_smoke_test hwymp_oracle)
endif()
