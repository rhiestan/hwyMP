# mpfix plan

Fixed-size IEEE 754 floating-point types for C++, aiming to beat Boost and match
or beat MPFR at 128 to 512 bits, with permissive licensing and Eigen support.

## Decisions

- Formats: IEEE interchange `binary128`, `binary256` and `binary512`
  (p = 113, 237, 489), through a `Float<P, E>` template.
- Speed is the main goal. A full-limb (p = 128/256/512) comparison benchmark is
  planned, using the same template where possible.
- Storage: unpacked limbs plus an 8-byte header; trivially copyable, 8-byte
  alignment. `to_bits`/`from_bits` provide the IEEE layouts.
- Operators round to nearest-even without flag tracking. Functions taking a
  context offer all 5 rounding modes and all 5 exception flags.
- No `__int128` in library code; tests may use it.
- Supported: MSVC and Clang on Windows, AppleClang on macOS ARM64, GCC and Clang
  on Linux.
- C++20. No expression templates.
- GMP and MPFR are black-box oracles, linked only into tests and benchmarks.

## Phases

0. Infrastructure: CMake, presets, vcpkg oracle, limb primitives, GitHub Actions.
1. Scalar core: add, subtract, multiply, compare, classify, round, in all modes.
   Differential tests against MPFR, and exhaustive tests on tiny formats.
2. Benchmark gate: done, see [BENCHMARKS.md](BENCHMARKS.md). Scalar beats MPFR
   everywhere with Clang. AVX2 batch failed its 1.2x target. MSVC speed is open.
3. Division, square root, FMA, remainder, roundToIntegral, min/max, nextUp,
   scalbn, logb: done. Division and square root are 1.4x to 3.5x faster than
   MPFR. FMA, tuned in phase 9, is 1.15x to 2.8x faster than MPFR IEEE.
4. Correctly rounded conversions: float, double, integers, between formats, and
   strings in both directions: done. Decimal parsing is 1.5 to 2.5 times
   faster than MPFR, printing level with it and 0.85 times its speed at
   binary512, with shortest round-trip output as well.
5. Eigen integration: NumTraits, numeric_limits, ADL math, decomposition tests:
   done, in `mpfix/eigen.h`, for Eigen 5.0.1. A Float of zero bytes is now +0,
   since Eigen clears trivially copyable scalars with memset.
6. Elementary functions, correctly rounded in all modes with flags, in
   `mpfix/math.h`. Done: exp, exp2, expm1, log, log2, log10, log1p, 1.7 to 8
   times faster than MPFR; C23's exp10, exp2m1, exp10m1, log2p1 and log10p1,
   4 to 25 times; sinh, cosh and tanh, 1.6 to 2.7 times; asinh,
   acosh and atanh, 2 to 10 times; sin, cos and tan, with exact reduction of
   any argument, 1.2 to 2.8 times; pow, 2 to 3.8 times; hypot, 1.3 to 2.9
   times; asin, acos, atan and atan2, 5 to 10 times; erf and erfc, 3 to 12
   times; tgamma and lgamma, 1.8 to 4.7 times; cbrt, level with MPFR at
   binary512 and faster below; and pi. Done.
7. Complex numbers. `mpfix/complex.h` with `mpfix::Complex<F>`,
   std::complex's interface and C99 Annex G's special values. Addition,
   subtraction, multiplication, abs, arg and norm are correctly rounded;
   division is within half an ulp plus 2^-(3P) in each component; sqrt, exp,
   log, log10, the trigonometric and hyperbolic functions and their inverses
   within half an ulp plus about 2^-60, by formulas without cancellation in
   P + 64 bits. Tested against MPC. numeric_limits is constexpr. Eigen takes
   Complex<F> as a complex scalar. What std::complex<Float> gives with
   libstdc++, libc++ and MSVC's STL is in docs/COMPLEX.md, probed by
   std_complex_test. Done.
8. Boost.Multiprecision: `mpfix/boost.h` with `mpfix::boost_backend<P, E>`, a
   backend of `boost::multiprecision::number` whose arithmetic and functions
   are Float's, correctly rounded, also those Boost would compute with
   Boost.Math, and std::numeric_limits for the numbers. Tested with Boost
   1.92, with Boost.Math, and with the standalone Boost.Multiprecision. Done.
9. Hardening and release. Done: ASan and UBSan over the whole suite in CI;
   Berkeley TestFloat 3e for binary16, 32, 64 and 128, all results and flags
   in all modes; libFuzzer targets for arithmetic and text against MPFR;
   weekly long runs; install rules with a CMake package; a benchmark summary
   and current MSVC figures. To do: tag the first release.

## Postponed

- Batch layer: AVX-512 IFMA only, and only if measured on real hardware, since
  AVX2 did not pay off. `FloatVec<P, E>`, dispatched kernels, tested through
  Intel SDE. Waiting for access to AVX-512 hardware.
