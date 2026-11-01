# hwyMP

Fixed-size IEEE 754 floating-point types for C++, built for speed at 128 to 512
bits. Work in progress: see [docs/PLAN.md](docs/PLAN.md).

## Building

Requires CMake 3.25 or later, Ninja and a C++20 compiler. Sibling `../highway`
and `../eigen` checkouts are used when present; otherwise they are downloaded.

On Windows, run from a Visual Studio 2022 developer prompt with `VCPKG_ROOT` set.
vcpkg builds GMP and MPFR, which only the tests use:

```
cmake --preset windows-clang
cmake --build --preset windows-clang
ctest --preset windows-clang
```

On Linux and macOS, install GMP and MPFR with your package manager:

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DHWYMP_ENABLE_ORACLE=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## Licence

Apache-2.0 OR BSD-3-Clause. See [LICENSE](LICENSE). GMP and MPFR are used only by
tests and benchmarks and are not part of the library.
