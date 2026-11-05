# Changelog

## 0.1.0 (unreleased)

The first release. Until 1.0, the interface may still change between minor
versions.

### Types

- `mpfix::Float<P, E>`: IEEE 754 binary floating point with a P-bit
  significand and an E-bit exponent, header-only, trivially copyable, with the
  IEEE interchange encoding through `to_bits` and `from_bits`.
- `binary128`, `binary256` and `binary512` (P = 113, 237, 489), and any other
  format, such as `Float<53, 11>` with the values of `double`.
- `mpfix::Complex<F>`, with `complex128`, `complex256` and `complex512`.
- `std::numeric_limits` for every Float.

### Correctly rounded

In all five IEEE rounding modes, with the five exception flags through an
optional `mpfix::Context`; the operators round to nearest-even without flags.

- Addition, subtraction, multiplication, division, square root, fma,
  remainder, rounding to integers, comparisons, min and max, nextUp and
  nextDown, scalbn and logb.
- Conversions between formats, from and to `float`, `double`, `long double`
  and the integers, and from and to decimal and hexadecimal text: `to_chars`
  and `from_chars` as in `<charconv>`, `to_string` with the shortest text that
  reads back, `from_string`, and streams.
- `mpfix/math.h`: `exp`, `exp2`, `exp10`, `expm1`, `exp2m1`, `exp10m1`, `log`,
  `log2`, `log10`, `log1p`, `log2p1`, `log10p1`, `sinh`, `cosh`, `tanh`,
  `asinh`, `acosh`, `atanh`, `sin`, `cos` and `tan` with exact reduction of any
  argument, `asin`, `acos`, `atan`, `atan2`, `pow`, `cbrt`, `hypot`, `erf`,
  `erfc`, `tgamma`, `lgamma` and `pi<F>()`. Exact results raise no flags.
  These need E of at most 26.

### Complex numbers

`mpfix/complex.h` has the interface of `std::complex` and the special values of
C99 Annex G. Addition, subtraction, multiplication, `abs`, `arg` and `norm` are
correctly rounded in each component. Division is within half an ulp plus
2^-(3P) in each component; `sqrt`, `exp`, `log`, `log10`, `pow` and the
trigonometric and hyperbolic functions and their inverses are within half an
ulp plus about 2^-60.

### Integrations

- `mpfix/eigen.h`: Float and Complex as Eigen scalars, tested with Eigen 5.0.1
  in the LU, QR, Cholesky, SVD, self-adjoint and complex eigenvalue
  decompositions and the Schur forms.
- `mpfix/boost.h`: Float as a backend of `boost::multiprecision::number`, as
  `boost_binary128`, `boost_binary256` and `boost_binary512`, with Float's
  correctly rounded results also for Boost.Math's functions. Needs Boost 1.79
  or later, or the standalone Boost.Multiprecision; tested with Boost 1.92.
- CMake: the target `mpfix::mpfix`, from `add_subdirectory`, `FetchContent` or
  an installed package found with `find_package(mpfix)`.

### Speed

On an AMD Ryzen 9 5950X with Clang, against MPFR 4.3 with the steps that
IEEE-identical results need, from [docs/BENCHMARKS.md](docs/BENCHMARKS.md):

- add, mul, div and fma: 1.4 to 2.8 times faster at binary128, 1.3 to 2.6 at
  binary256, 1.0 to 1.4 at binary512;
- the elementary and special functions: from level with MPFR, for `cbrt` at
  binary512, to 25 times faster, for the C23 functions such as `exp10`;
- parsing 1.5 to 2.5 times faster; printing level with MPFR, and 0.85 times
  its speed at binary512;
- against Boost's `cpp_bin_float`, 3.5 to 46 times faster.

### Testing

Against MPFR and MPC at every precision, Berkeley TestFloat 3e for binary16,
32, 64 and 128 with all results and flags in all modes, exhaustively on small
formats, under AddressSanitizer and UndefinedBehaviorSanitizer, and by
libFuzzer against MPFR; weekly at 50 times the size. Line coverage of the
headers is 98.8%, measured with Clang.

### Supported

C++20 with MSVC and Clang on Windows x64, GCC and Clang on Linux x64, and
AppleClang on macOS ARM64, as in CI. The library needs no other library, and
uses no `__int128`.

### Known limitations

- MSVC builds are 1.2 to 3 times slower than Clang's.
- MinGW is not supported: three tests compare against its C library, whose
  `modf`, `sin` and `cos` are wrong for some arguments.
- No SIMD batch layer: AVX2 did not pay off, and AVX-512 IFMA waits for
  hardware to measure on.
