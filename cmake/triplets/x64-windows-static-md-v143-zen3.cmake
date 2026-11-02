# For benchmarks: selects GMP's Zen 3 assembly and tuning. The default triplet
# builds for generic x86-64, which would understate GMP and MPFR. The resulting
# library only runs on CPUs that support Zen 3's instructions.
include("${CMAKE_CURRENT_LIST_DIR}/x64-windows-static-md-v143.cmake")

if(PORT STREQUAL "gmp")
  # The port's autoreconf replaces GMP's config.sub, and the stock one rejects
  # CPU names such as zen3. Seeding autoconf's cache skips that check.
  set(VCPKG_MAKE_BUILD_TRIPLET
      "--host=x86_64-pc-mingw32;--build=x86_64-pc-mingw32;ac_cv_host=zen3-pc-mingw32;ac_cv_build=zen3-pc-mingw32")
endif()
