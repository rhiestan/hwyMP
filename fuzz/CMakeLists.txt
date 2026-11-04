# libFuzzer targets, for Clang. Added from tests/, whose MPFR oracle they use
# and whose MPFIX_SANITIZE flags apply here too. Run one with, for example,
#   fuzz_arithmetic -max_total_time=600 corpus/
if(NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  message(FATAL_ERROR "MPFIX_BUILD_FUZZERS needs Clang, for libFuzzer")
endif()

foreach(name fuzz_arithmetic fuzz_strings)
  add_executable(${name} ${name}.cpp)
  target_include_directories(${name} PRIVATE "${PROJECT_SOURCE_DIR}/tests")
  target_link_libraries(${name} PRIVATE mpfix mpfix_oracle GTest::gtest)
  target_compile_options(${name} PRIVATE -fsanitize=fuzzer)
  mpfix_warnings(${name})
  target_link_options(${name} PRIVATE -fsanitize=fuzzer)
endforeach()
