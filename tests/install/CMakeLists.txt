# Builds against an installed mpfix, found with find_package: the CI checks
# that the install rules give a usable package.
cmake_minimum_required(VERSION 3.25)
project(mpfix_install_check LANGUAGES CXX)

find_package(mpfix 0.1 REQUIRED)
add_executable(install_check install_check.cpp)
target_link_libraries(install_check PRIVATE mpfix::mpfix)
