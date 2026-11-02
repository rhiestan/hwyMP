# Benchmarks

Measured on 2026-11-01 on an AMD Ryzen 9 5950X (Zen 3, AVX2, no AVX-512) under
Windows 11. One pinned thread at high priority, best of 25 runs over 1024 random
operand pairs with exponents within ±8 and random signs. Nanoseconds per
operation; lower is better. "dot step" is one operation of a multiply-add loop.

- hwyMP and Boost 1.92: Clang 22.1.2, `-O3 -march=native`.
- GMP 6.3.0: Zen 3 assembly and tuning (`cmake/triplets/...-zen3.cmake`).
- MPFR 4.2.2: its C code is compiled by MSVC, because a clang-cl build through
  vcpkg did not work. This may understate MPFR, most of all at binary128, where
  its time is spent in its own code rather than in GMP's assembly.

"MPFR" uses its native wide exponent range, which favours it. "MPFR IEEE" adds
the steps that IEEE-identical results need, matching what hwyMP delivers.

## Scalar

| Format | Library | add | mul | dot step |
|---|---|---:|---:|---:|
| binary128 | hwyMP | 8.12 | 7.55 | 10.62 |
| | MPFR | 11.83 | 15.14 | 13.39 |
| | MPFR IEEE | 16.74 | 20.13 | 18.93 |
| | Boost | 36.49 | 50.15 | 42.24 |
| binary256 | hwyMP | 12.11 | 14.20 | 15.53 |
| | MPFR | 19.56 | 36.06 | 25.69 |
| | MPFR IEEE | 24.46 | 41.27 | 30.82 |
| | Boost | 46.48 | 74.39 | 55.52 |
| binary512 | hwyMP | 16.41 | 47.92 | 32.19 |
| | MPFR | 21.92 | 57.37 | 34.20 |
| | MPFR IEEE | 27.45 | 62.62 | 39.41 |
| | Boost | 62.82 | 167.40 | 96.53 |

### Full-limb precision

Full-limb formats cost hwyMP nothing measurable: 128, 256 and 512 bits run at
the same speed as binary128, binary256 and binary512, within noise.

| Format | Library | add | mul | dot step |
|---|---|---:|---:|---:|
| full 128 | hwyMP | 8.21 | 7.56 | 10.42 |
| | MPFR | 11.12 | 27.98 | 17.30 |
| full 256 | hwyMP | 11.95 | 14.46 | 15.28 |
| | MPFR | 19.40 | 35.99 | 25.71 |
| full 512 | hwyMP | 15.64 | 47.52 | 32.02 |
| | MPFR | 21.43 | 56.43 | 34.01 |

### MSVC

MSVC 19.44 builds of hwyMP are much slower at small sizes, and at binary128 they
are slower than MPFR. This is an open issue.

| Format | add | mul | dot step |
|---|---:|---:|---:|
| binary128 | 21.41 | 24.03 | 24.13 |
| binary256 | 23.82 | 30.16 | 28.07 |
| binary512 | 28.13 | 57.80 | 40.56 |

## Batch prototype

`bench/bench_batch.cpp` runs binary256 on Highway's AVX2 target with 4 lanes, in
structure-of-arrays form with radix 2^28 limbs. It covers regular values only,
and its results match the scalar `Float` exactly.

| Operation | Scalar | Batch, per number | Speedup |
|---|---:|---:|---:|
| mul | 14.07 | 14.39 | 0.98 |
| add | 11.56 | 12.98 | 0.89 |

Converting one number from `Float` storage to this layout takes 19 ns. Eigen
stores arrays of `Float`, so it would pay that for each input and output.

## Gate decision

- **Scalar against MPFR: passed.** hwyMP is faster on every operation at every
  size, against both MPFR configurations.
- **Scalar against Boost, target 5x: not met.** Add is 3.8x to 4.5x and mul is
  3.5x to 6.6x. MPFR is the goal that matters, so this does not block anything.
- **Batch on AVX2, target 1.2x: failed.** Batch stays out of the plan for AVX2.
  AVX-512 IFMA may change the result, but it cannot be measured on this machine.
