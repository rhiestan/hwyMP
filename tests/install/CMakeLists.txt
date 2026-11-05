# Builds against mpfix as another project would: by default an installed
# mpfix, found with find_package, so that the CI checks the install rules give
# a usable package; with MPFIX_SOURCE_DIR set, the source tree through
# add_subdirectory, as FetchContent includes it.
cmake_minimum_required(VERSION 3.25)
project(mpfix_install_check LANGUAGES CXX)

set(MPFIX_SOURCE_DIR "" CACHE PATH "mpfix's source tree, to add as a subdirectory")
if(MPFIX_SOURCE_DIR)
  add_subdirectory("${MPFIX_SOURCE_DIR}" mpfix)
  # Not the top-level project, mpfix builds no tests and installs nothing.
  if(TARGET limb_ops_test OR MPFIX_INSTALL)
    message(FATAL_ERROR "mpfix as a subdirectory builds tests or installs")
  endif()
else()
  find_package(mpfix 0.1 REQUIRED)
endif()
add_executable(install_check install_check.cpp)
target_link_libraries(install_check PRIVATE mpfix::mpfix)
