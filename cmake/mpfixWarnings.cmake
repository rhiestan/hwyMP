# mpfix_warnings(target): the warnings of the project's own programs (tests,
# fuzzers and benchmarks), and clang-tidy with MPFIX_CLANG_TIDY. The warnings
# that point at real bugs are errors.
function(mpfix_warnings name)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(${name} PRIVATE
      /W4 /permissive- /Zc:__cplusplus /bigobj
      # Shadowed names (4456 to 4459), uninitialized variables (4700, 4701,
      # 4703), a path without a return value (4715), the address of a local
      # returned (4172) and constant conditions (4127).
      /we4456 /we4457 /we4458 /we4459 /we4700 /we4701 /we4703 /we4715
      /we4172 /we4127)
  elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if(CMAKE_CXX_COMPILER_FRONTEND_VARIANT STREQUAL "MSVC")
      target_compile_options(${name} PRIVATE /W4 /bigobj)
    endif()
    target_compile_options(${name} PRIVATE
      -Wall -Wextra -Wshadow-all -Wpedantic
      -Werror=return-type -Werror=dangling -Werror=uninitialized
      -Werror=unreachable-code -Werror=constant-logical-operand
      -Werror=shadow)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(${name} PRIVATE
      -Wall -Wextra -Wpedantic
      -Werror=return-type -Werror=uninitialized -Werror=shadow)
  endif()
  if(MPFIX_CLANG_TIDY)
    set_target_properties(${name} PROPERTIES CXX_CLANG_TIDY "${MPFIX_CLANG_TIDY_EXE}")
  endif()
  if(WIN32)
    target_compile_definitions(${name} PRIVATE _CRT_SECURE_NO_WARNINGS)
  endif()
endfunction()
