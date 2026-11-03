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

`mpfix/math.h` adds `exp`, `exp2`, `exp10`, `expm1`, `exp2m1`, `exp10m1`,
`log`, `log2`, `log10`, `log1p`, `log2p1`, `log10p1`, `sinh`, `cosh`, `tanh`,
`asinh`, `acosh`, `atanh`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`,
`atan2`, `pow`, `cbrt`, `hypot`, `erf`, `erfc`, `tgamma` and `lgamma`, and the
constant `pi<F>()`, correctly rounded in every mode, with an optional
`Context` for the mode and flags. Exact results, such as `pow(9, 1.5)`, raise
no flags. The trigonometric functions reduce any argument exactly, as if with
infinitely many digits of pi:

```cpp
#include "mpfix/math.h"

binary256 e = mpfix::exp(binary256(1));
mpfix::Context down{mpfix::Rounding::kDown};
binary256 lower = mpfix::log(binary256(10), down);  // Rounded down.
binary256 r = mpfix::pow(binary256(2), binary256(0.5));
binary256 s = mpfix::sin(binary256(1e300));  // Correct to the last bit.
binary256 p = mpfix::pi<binary256>();
```

`to_chars` and `from_chars` follow `<charconv>`, with an optional `Context`.
The formats are `binary128`, `binary256` and `binary512`, or any
`mpfix::Float<P, E>`. `std::numeric_limits` is specialized for all of them.

### Complex numbers

`mpfix/complex.h` adds `mpfix::Complex<F>`, with `complex128`, `complex256` and
`complex512`, which has the interface of `std::complex` and the special values
of C99 Annex G. `std::complex<T>` is unspecified for T other than the built-in
types, and what the standard libraries make of `std::complex<binary128>`
differs; see [docs/COMPLEX.md](docs/COMPLEX.md). Addition, subtraction,
multiplication, `abs`, `arg` and `norm` are correctly rounded in each
component; division, `sqrt`, `exp`, `log`, `log10`, `pow` and the
trigonometric and hyperbolic functions and their inverses are within half an
ulp plus a tiny fraction:

```cpp
#include "mpfix/complex.h"

mpfix::complex256 z(mpfix::binary256(1), mpfix::binary256(2));
mpfix::complex256 w = mpfix::exp(z) / mpfix::sqrt(z);
```

### With Eigen

```cpp
#include <Eigen/Dense>
#include "mpfix/eigen.h"  // NumTraits for Float; includes mpfix/charconv.h.

using Matrix = Eigen::Matrix<mpfix::binary256, Eigen::Dynamic, Eigen::Dynamic>;

Matrix a = Matrix::Random(10, 10);
Matrix x = a.partialPivLu().solve(Matrix::Ones(10, 1));
std::cout << x.format(Eigen::IOFormat(Eigen::FullPrecision)) << '\n';
```

Tested with Eigen 5.0.1: the LU, QR, Cholesky, SVD and self-adjoint
eigenvalue decompositions, the real Schur form, inverses and determinants.
`mpfix/eigen.h` includes `mpfix/math.h`, so Eigen finds `exp`, `log`, `sin`,
`atan2`, `pow`, `hypot` and the rest.

`mpfix/eigen.h` also makes `mpfix::Complex<F>` a complex scalar, which mixes
with F in expressions: `Eigen::Matrix<mpfix::complex256, Eigen::Dynamic,
Eigen::Dynamic>` works with the same decompositions and with `ComplexSchur` and
`ComplexEigenSolver`. The complex eigenvalues of a real matrix, from
`EigenSolver`, are `std::complex<F>` in Eigen; see
[docs/COMPLEX.md](docs/COMPLEX.md).

### With Boost.Multiprecision

`mpfix/boost.h` makes Float a backend of `boost::multiprecision::number`, for
code written for Boost's types, and for Boost.Math:

```cpp
#include "mpfix/boost.h"  // Needs Boost 1.79 or later, or standalone Boost.Multiprecision.

mpfix::boost_binary256 x = 2;
mpfix::boost_binary256 y = sqrt(x) + boost::multiprecision::erf(x);
std::cout << std::setprecision(70) << y << '
';
```

`boost_binary128`, `boost_binary256` and `boost_binary512` give exactly Float's
results: arithmetic, `sqrt`, `fma`, `fmod`, the rounding functions and all of
the functions of `mpfix/math.h` stay correctly rounded rather than going to
Boost's generic code or Boost.Math. Expression templates are off by default;
`number<mpfix::boost_backend<P, E>, et_on>` works too. Boost is not a
dependency of the library; only this header needs it.

## Building

The library is header-only and needs nothing but a C++20 compiler: add
`include/` to the include path, or link the CMake target `mpfix::mpfix`.

The tests need CMake 3.25 or later and Ninja. A sibling `../eigen` checkout is
used when present, `../boost_1_92_0` for the Boost test and the benchmarks, and
`../highway` for the batch benchmark; otherwise Eigen, the standalone
Boost.Multiprecision and Highway are downloaded.

On Windows, run from a Visual Studio 2022 developer prompt with `VCPKG_ROOT` set.
vcpkg builds GMP, MPFR and MPC, which only the tests use:

```
cmake --preset windows-clang
cmake --build --preset windows-clang
ctest --preset windows-clang
```

On Linux and macOS, install GMP, MPFR and MPC with your package manager:

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DMPFIX_ENABLE_ORACLE=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## Licence

Apache-2.0 OR BSD-3-Clause. See [LICENSE](LICENSE). GMP, MPFR and MPC are used
only by tests and benchmarks and are not part of the library.
