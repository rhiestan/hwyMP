# mpfix

Fixed-size IEEE 754 floating-point types for C++, built for speed at 128 to 512
bits. Work in progress: see [docs/PLAN.md](docs/PLAN.md).

## Use

```cpp
#include "mpfix/charconv.h"  // Strings and streams; includes mpfix/float.h.

using mpfix::binary256;

binary256 x = 1.5;                        // Exact conversions are implicit.
binary256 third = binary256(1) / binary256(3);
double d = static_cast<double>(third);    // Rounds to nearest-even.

mpfix::Context context{mpfix::Rounding::kUp};
float f = mpfix::convert<float>(third, context);  // Any mode, with flags.
long n = mpfix::to_integer<long>(x, mpfix::Rounding::kNearestEven);

std::string text = mpfix::to_string(third);       // Shortest round trip.
binary256 y = mpfix::from_string<binary256>("0.1");
std::cout << std::setprecision(40) << y << '\n';
```

`to_chars` and `from_chars` follow `<charconv>`, with an optional `Context`.
The formats are `binary128`, `binary256` and `binary512`, or any
`mpfix::Float<P, E>`.

## Building

The library is header-only and needs nothing but a C++20 compiler: add
`include/` to the include path, or link the CMake target `mpfix::mpfix`.

The tests need CMake 3.25 or later and Ninja. A sibling `../eigen` checkout is
used when present, and `../highway` for the batch benchmark; otherwise they are
downloaded.

On Windows, run from a Visual Studio 2022 developer prompt with `VCPKG_ROOT` set.
vcpkg builds GMP and MPFR, which only the tests use:

```
cmake --preset windows-clang
cmake --build --preset windows-clang
ctest --preset windows-clang
```

On Linux and macOS, install GMP and MPFR with your package manager:

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DMPFIX_ENABLE_ORACLE=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## Licence

Apache-2.0 OR BSD-3-Clause. See [LICENSE](LICENSE). GMP and MPFR are used only by
tests and benchmarks and are not part of the library.
