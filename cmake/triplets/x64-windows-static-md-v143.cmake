# Static libraries against the dynamic CRT, built with Visual Studio 2022.
# VS 2026 also ships the v143 toolset, so the toolset alone is not enough to
# keep vcpkg from selecting it: the installation path is pinned as well.
set(VCPKG_TARGET_ARCHITECTURE x64)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE static)
set(VCPKG_PLATFORM_TOOLSET v143)

foreach(edition Enterprise Professional Community BuildTools)
  set(vs2022_path "C:/Program Files/Microsoft Visual Studio/2022/${edition}")
  if(EXISTS "${vs2022_path}")
    # vcpkg compares this as a string against native, backslash-separated paths.
    file(TO_NATIVE_PATH "${vs2022_path}" VCPKG_VISUAL_STUDIO_PATH)
    break()
  endif()
endforeach()

# Release only: with the debug CRT, a failing configure probe in MPFR opens an
# abort() dialog and blocks the build. Tests and benchmarks use Release anyway.
set(VCPKG_BUILD_TYPE release)
