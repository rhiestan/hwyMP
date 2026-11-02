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
   `mpfix/math.h`. Done: exp, exp2, expm1, log, log2, log10, log1p, 1.1 to 7
   times faster than MPFR. To do: pow, the trigonometric and hyperbolic
   functions and their inverses, cbrt, hypot, and constants such as pi.
7. Batch layer: AVX-512 IFMA only, and only if measured on real hardware, since
   AVX2 did not pay off. `FloatVec<P, E>`, dispatched kernels, tested through
   Intel SDE.
8. Hardening and release: long differential runs, TestFloat, sanitizers,
   fuzzing, and a benchmark report.
