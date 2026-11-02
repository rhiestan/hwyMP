# Benchmarks

Scalar results were measured on 2026-11-02, and the rest on 2026-11-01, on an
AMD Ryzen 9 5950X (Zen 3, AVX2, no AVX-512) under
Windows 11. One pinned thread at high priority, best of 25 runs over 1024 random
operand pairs with exponents within ±8 and random signs. Nanoseconds per
operation; lower is better. "dot step" is one operation of a multiply-add loop.
Square roots take positive operands. Boost has no correctly rounded fma for this
type, so it has no fma column.

- mpfix and Boost 1.92: Clang 22.1.2, `-O3 -march=native`.
- GMP 6.3.0: Zen 3 assembly and tuning (`cmake/triplets/...-zen3.cmake`).
- MPFR 4.2.2: its C code is compiled by MSVC, because a clang-cl build through
  vcpkg did not work. This may understate MPFR, most of all at binary128, where
  its time is spent in its own code rather than in GMP's assembly.

"MPFR" uses its native wide exponent range, which favours it. "MPFR IEEE" adds
the steps that IEEE-identical results need, matching what mpfix delivers.

## Scalar

| Format | Library | add | mul | div | sqrt | fma | dot step |
|---|---|---:|---:|---:|---:|---:|---:|
| binary128 | mpfix | 7.93 | 7.38 | 14.84 | 19.91 | 23.42 | 10.50 |
| | MPFR | 11.97 | 14.82 | 26.90 | 28.50 | 41.81 | 13.20 |
| | MPFR IEEE | 16.88 | 20.53 | 31.69 | 32.33 | 48.23 | 19.35 |
| | Boost | 35.86 | 47.99 | 178.31 | 282.44 | | 41.56 |
| binary256 | mpfix | 12.06 | 14.19 | 29.72 | 62.46 | 36.48 | 15.18 |
| | MPFR | 19.08 | 35.63 | 104.53 | 116.61 | 54.90 | 25.63 |
| | MPFR IEEE | 24.90 | 41.51 | 111.40 | 122.96 | 60.84 | 30.67 |
| | Boost | 47.11 | 72.90 | 416.67 | 536.17 | | 54.30 |
| binary512 | mpfix | 15.74 | 51.80 | 81.12 | 140.59 | 78.96 | 32.02 |
| | MPFR | 21.41 | 56.05 | 157.43 | 251.44 | 77.15 | 33.71 |
| | MPFR IEEE | 27.01 | 60.80 | 166.07 | 258.67 | 83.38 | 38.58 |
| | Boost | 61.77 | 166.68 | 1069.02 | 1086.34 | | 94.27 |

mpfix's fma is a straightforward double-width addition and has not been tuned.
At binary512 it is level with MPFR rather than ahead.

### Full-limb precision

Full-limb formats cost mpfix little: 128, 256 and 512 bits run at the same speed
as binary128, binary256 and binary512, except for division at 128 and 256 bits,
which needs one more quotient limb.

| Format | Library | add | mul | div | sqrt | fma | dot step |
|---|---|---:|---:|---:|---:|---:|---:|
| full 128 | mpfix | 8.02 | 7.42 | 17.55 | 19.71 | 24.03 | 10.25 |
| | MPFR | 11.06 | 27.09 | 42.47 | 113.70 | 42.79 | 17.41 |
| full 256 | mpfix | 11.87 | 14.22 | 34.97 | 63.46 | 36.14 | 14.96 |
| | MPFR | 18.88 | 35.68 | 103.95 | 169.48 | 54.97 | 25.48 |
| full 512 | mpfix | 16.34 | 48.62 | 84.12 | 138.23 | 79.02 | 32.27 |
| | MPFR | 21.58 | 56.52 | 157.96 | 268.47 | 79.67 | 33.75 |

### MSVC

MSVC 19.44 builds of mpfix are much slower at small sizes, and at binary128 they
are slower than MPFR. This is an open issue. These figures predate division and
square root.

| Format | add | mul | dot step |
|---|---:|---:|---:|
| binary128 | 21.41 | 24.03 | 24.13 |
| binary256 | 23.82 | 30.16 | 28.07 |
| binary512 | 28.13 | 57.80 | 40.56 |

## Strings

`bench/bench_strings.cpp`, measured on 2026-11-02 with the same setup: best of
15 runs over 512 random values, either near 1 (binary exponents within ±8) or
"wide" (binary exponents within ±1000). Nanoseconds per conversion.

- **print:** enough significant digits to read back (36, 73 and 149), in
  scientific notation: `to_chars` with a precision, `mpfr_get_str`, and Boost's
  `str`.
- **shortest:** the shortest string that reads back. Only mpfix offers this.
- **parse:** reading the shortest strings back: `from_chars`, `mpfr_strtofr`,
  and Boost's constructor from a string.

All mpfix conversions are correctly rounded and raise the IEEE flags.

| Format | Values | Library | print | shortest | parse |
|---|---|---|---:|---:|---:|
| binary128 | near 1 | mpfix | 341 | 966 | 144 |
| | | MPFR | 480 | | 347 |
| | | Boost | 784 | | 824 |
| binary256 | near 1 | mpfix | 452 | 1387 | 212 |
| | | MPFR | 598 | | 472 |
| | | Boost | 1501 | | 1554 |
| binary512 | near 1 | mpfix | 779 | 2110 | 479 |
| | | MPFR | 857 | | 667 |
| | | Boost | 4089 | | 3644 |
| binary128 | wide | mpfix | 350 | 929 | 150 |
| | | MPFR | 577 | | 397 |
| | | Boost | 1316 | | 976 |
| binary256 | wide | mpfix | 477 | 1251 | 216 |
| | | MPFR | 688 | | 508 |
| | | Boost | 2137 | | 1645 |
| binary512 | wide | mpfix | 805 | 2279 | 474 |
| | | MPFR | 938 | | 684 |
| | | Boost | 4627 | | 3768 |

Conversions are computed with interval arithmetic in a slightly wider `Float`,
which decides nearly every case; the rest fall back to exact big-integer
arithmetic. That fallback is what costs time near the far ends of binary512's
exponent range, where its integers have millions of bits. There, a decimal
string within a relative 10^-170 of a rounding boundary takes tens of
milliseconds, and so does output with more than about 780 significant digits.
Fixed notation of such values has up to 1.3 million digits and is slow in
proportion.

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

- **Scalar against MPFR: passed.** mpfix is faster on every operation at every
  size, against both MPFR configurations.
- **Scalar against Boost, target 5x: not met.** Add is 3.8x to 4.5x and mul is
  3.5x to 6.6x. MPFR is the goal that matters, so this does not block anything.
- **Batch on AVX2, target 1.2x: failed.** Batch stays out of the plan for AVX2.
  AVX-512 IFMA may change the result, but it cannot be measured on this machine.
