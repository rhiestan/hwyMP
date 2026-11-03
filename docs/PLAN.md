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
   MPFR. FMA is untuned and level with MPFR at binary512.
4. Correctly rounded conversions: float, double, integers, between formats, and
   strings in both directions: done. Decimal printing and parsing are faster
   than MPFR at every size, with shortest round-trip output as well.
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
7. Complex numbers. `std::complex<T>` is unspecified for T other than float,
   double and long double, and the standard libraries differ: libstdc++ (GCC
   16) works at full precision, finding mpfix's functions by
   argument-dependent lookup; MSVC's STL deprecates it (STL4037), fails to
   compile `std::log` because `numeric_limits<Float>::epsilon()` is not
   constexpr, and its generic code converts to double for exp, cos, sinh,
   atan2, isinf and signbit, silently losing precision; libc++ is untested.
   MSVC's generic code also needs constexpr arithmetic, which Float lacks.
   Done: `mpfix/complex.h` with `mpfix::Complex<F>`, std::complex's
   interface and C99 Annex G's special values. Addition, subtraction,
   multiplication, abs, arg and norm are correctly rounded; division is within
   half an ulp plus 2^-(3P) in each component; sqrt, exp, log, log10, sin,
   cos, tan, sinh, cosh and tanh within half an ulp plus about 2^-60, by
   formulas without cancellation in P + 64 bits. Tested against MPC, which
   they match in all but about 1 in 10000 components, those off by one ulp.
   numeric_limits is constexpr. To do: asin, acos, atan, asinh, acosh and
   atanh; how Eigen's complex matrices, written for std::complex, can use
   Complex; and what std::complex gives on each library, documented.
8. Hardening and release: long differential runs, TestFloat, sanitizers,
   fuzzing, and a benchmark report.

## Postponed

- Batch layer: AVX-512 IFMA only, and only if measured on real hardware, since
  AVX2 did not pay off. `FloatVec<P, E>`, dispatched kernels, tested through
  Intel SDE. Waiting for access to AVX-512 hardware.
