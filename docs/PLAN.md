# hwyMP plan

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
2. Benchmark gate: scalar and SoA Highway add/mul at 256 bits against MPFR and
   Boost. Scalar target: at least MPFR's speed and 5x Boost's. Batch work stays
   AVX-512-only if AVX2 throughput is below about 1.2x scalar. Full-limb
   comparison belongs here or in phase 7. Before any MPFR timing: vcpkg builds
   GMP with generic `k8` tuning and MPFR with MSVC, which understates MPFR.
   Rebuild both tuned for the CPU, for example with an overlay port.
3. Division, square root, FMA, remainder, roundToIntegral, min/max, nextUp,
   scalbn, logb.
4. Correctly rounded conversions: float, double, integers, between formats, and
   strings in both directions.
5. Eigen integration: NumTraits, numeric_limits, ADL math, decomposition tests.
6. Batch layer: `FloatVec<P, E>`, dispatched kernels, tested through Intel SDE.
7. Hardening and release: long differential runs, TestFloat, sanitizers,
   fuzzing, and a benchmark report.
