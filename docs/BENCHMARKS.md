# Benchmarks

Measured on 2026-11-02 on an AMD Ryzen 9 5950X (Zen 3, AVX2, no AVX-512) under
Windows 11. One pinned thread at high priority, best of 25 runs over 1024 random
operand pairs with exponents within ±8 and random signs. Nanoseconds per
operation; lower is better. "dot step" is one operation of a multiply-add loop.
Square roots take positive operands. Boost has no correctly rounded fma for this
type, so it has no fma column.

- mpfix and Boost 1.92: Clang 22.1.2, `-O3 -march=native`.
- GMP 6.3.0 and MPFR 4.3.0-dev: built as DLLs with MSYS2's GCC, GMP as a fat
  binary that picks its assembly for the running CPU.

"MPFR" uses its native wide exponent range, which favours it. "MPFR IEEE" adds
the steps that IEEE-identical results need, matching what mpfix delivers.

An earlier MPFR 4.2.2 build through vcpkg had its C code compiled by MSVC, and
was much slower: at binary128 its add, mul and sqrt took 12.0, 14.8 and 28.5 ns,
against 8.4, 7.3 and 14.0 here. Benchmarks against MPFR on Windows should use a
GCC build.

## Scalar

| Format | Library | add | mul | div | sqrt | fma | dot step |
|---|---|---:|---:|---:|---:|---:|---:|
| binary128 | mpfix | 8.20 | 7.53 | 14.93 | 20.49 | 24.23 | 10.92 |
| | MPFR | 8.37 | 7.31 | 16.17 | 13.95 | 30.27 | 7.78 |
| | MPFR IEEE | 13.36 | 13.20 | 20.32 | 19.47 | 36.45 | 13.10 |
| | Boost | 36.66 | 48.75 | 182.81 | 282.39 | | 41.92 |
| binary256 | mpfix | 12.30 | 14.76 | 30.39 | 63.76 | 37.06 | 15.53 |
| | MPFR | 15.60 | 31.92 | 73.56 | 104.82 | 42.13 | 20.00 |
| | MPFR IEEE | 20.45 | 33.82 | 79.30 | 111.16 | 47.46 | 26.06 |
| | Boost | 46.59 | 75.00 | 424.64 | 547.42 | | 55.02 |
| binary512 | mpfix | 16.73 | 48.93 | 84.06 | 143.84 | 80.27 | 32.56 |
| | MPFR | 17.91 | 57.29 | 112.74 | 216.89 | 75.54 | 31.31 |
| | MPFR IEEE | 22.70 | 62.49 | 117.64 | 220.71 | 81.25 | 36.39 |
| | Boost | 63.36 | 169.23 | 1087.66 | 1121.13 | | 96.05 |

Against MPFR IEEE, which gives the same results, mpfix is faster everywhere
except binary128 square root, and level on binary512 fma. Native MPFR is ahead
on binary128 square root (14.0 against 20.5 ns), on the binary128 dot step, and
on binary512 fma and dot step. These are the targets for tuning: mpfix's fma is
a straightforward double-width addition, and its square root has no special
case for two limbs.

### Full-limb precision

Full-limb formats cost mpfix little: 128, 256 and 512 bits run at the same speed
as binary128, binary256 and binary512, except for division at 128 and 256 bits,
which needs one more quotient limb. MPFR is slower at full limbs, most of all
for square root.

| Format | Library | add | mul | div | sqrt | fma | dot step |
|---|---|---:|---:|---:|---:|---:|---:|
| full 128 | mpfix | 8.32 | 7.58 | 17.57 | 20.22 | 23.91 | 10.66 |
| | MPFR | 7.88 | 13.95 | 36.02 | 101.39 | 32.60 | 9.61 |
| full 256 | mpfix | 12.11 | 14.65 | 35.82 | 64.25 | 37.08 | 15.45 |
| | MPFR | 15.38 | 28.40 | 74.20 | 143.93 | 44.45 | 20.21 |
| full 512 | mpfix | 16.12 | 48.94 | 85.88 | 143.78 | 80.42 | 32.43 |
| | MPFR | 17.65 | 57.43 | 111.38 | 231.01 | 77.16 | 31.23 |

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

`bench/bench_strings.cpp`, with the same setup: best of 15 runs over 512 random
values, either near 1 (binary exponents within ±8) or "wide" (binary exponents
within ±1000). Nanoseconds per conversion.

- **print:** enough significant digits to read back (36, 73 and 149), in
  scientific notation: `to_chars` with a precision, `mpfr_get_str`, and Boost's
  `str`.
- **shortest:** the shortest string that reads back. Only mpfix offers this.
- **parse:** reading the shortest strings back: `from_chars`, `mpfr_strtofr`,
  and Boost's constructor from a string.

All mpfix conversions are correctly rounded and raise the IEEE flags.

| Format | Values | Library | print | shortest | parse |
|---|---|---|---:|---:|---:|
| binary128 | near 1 | mpfix | 342 | 965 | 151 |
| | | MPFR | 356 | | 371 |
| | | Boost | 794 | | 825 |
| binary256 | near 1 | mpfix | 446 | 1261 | 229 |
| | | MPFR | 455 | | 494 |
| | | Boost | 1538 | | 1559 |
| binary512 | near 1 | mpfix | 780 | 2089 | 441 |
| | | MPFR | 665 | | 668 |
| | | Boost | 3271 | | 3398 |
| binary128 | wide | mpfix | 353 | 931 | 156 |
| | | MPFR | 431 | | 429 |
| | | Boost | 1319 | | 974 |
| binary256 | wide | mpfix | 474 | 1252 | 231 |
| | | MPFR | 508 | | 525 |
| | | Boost | 2120 | | 1604 |
| binary512 | wide | mpfix | 804 | 2199 | 438 |
| | | MPFR | 704 | | 713 |
| | | Boost | 3855 | | 3335 |

Parsing is 1.5 to 2.5 times faster than MPFR. Printing is level with MPFR at
binary128 and binary256, and about 15% slower at binary512.

Conversions are computed with interval arithmetic in a slightly wider `Float`,
which decides nearly every case; the rest fall back to exact big-integer
arithmetic. That fallback is what costs time near the far ends of binary512's
exponent range, where its integers have millions of bits. There, a decimal
string within a relative 10^-170 of a rounding boundary takes tens of
milliseconds, and so does output with more than about 780 significant digits.
Fixed notation of such values has up to 1.3 million digits and is slow in
proportion.

## Elementary functions

### Exponentials and logarithms

`bench/bench_math.cpp`, same machine and libraries, not pinned. Best of 9 runs
over 256 arguments: exponentials in (-20, 20), logarithms of positive values
with binary exponents in [-20, 20), expm1 and log1p below 1/4 in magnitude.
Nanoseconds per call, all rounding to nearest. Both libraries round correctly.
Boost's `cpp_bin_float` is shown for exp and log only.

| Format | Library | exp | exp2 | expm1 | log | log2 | log10 | log1p |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| binary128 | mpfix | 600 | 533 | 422 | 392 | 414 | 406 | 281 |
| | MPFR | 1032 | 1277 | 1007 | 1741 | 1929 | 3433 | 1863 |
| | Boost | 5034 | | | 6533 | | | |
| binary256 | mpfix | 1047 | 1008 | 846 | 777 | 819 | 793 | 605 |
| | MPFR | 1624 | 1901 | 1550 | 2777 | 3254 | 5793 | 4697 |
| | Boost | 17296 | | | 18331 | | | |
| binary512 | mpfix | 2493 | 2368 | 2055 | 2152 | 2273 | 2200 | 1803 |
| | MPFR | 2892 | 3237 | 2799 | 4687 | 4962 | 8978 | 7386 |
| | Boost | 75066 | | | 61405 | | | |

mpfix is 1.2 to 2.4 times faster than MPFR for the exponentials and 2.2 to 8
times faster for the logarithms. The margin for exp narrows with size: at
binary512 its series still costs about 20 full multiplications and as many
additions, where MPFR works on fixed-point limbs.

Most results are decided by the bits of the working value alone: rounding
boundaries all have at most P + 1 bits, so when the bits below those are
neither nearly all zeros nor nearly all ones, the value rounds directly.

### Powers and roots

The same benchmark. pow takes x as for log and y as for exp; cbrt takes values
as for log, and hypot pairs of them.

| Format | Library | pow | cbrt | hypot |
|---|---|---:|---:|---:|
| binary128 | mpfix | 978 | 273 | 94 |
| | MPFR | 3069 | 425 | 214 |
| | Boost | 12715 | | |
| binary256 | mpfix | 1818 | 461 | 161 |
| | MPFR | 4964 | 513 | 318 |
| | Boost | 37547 | | |
| binary512 | mpfix | 4798 | 668 | 345 |
| | MPFR | 8092 | 662 | 462 |
| | Boost | 143011 | | |

pow is e^(y log x) with exact and tied results detected beforehand. cbrt runs
Newton's iteration on the inverse cube root, which needs no division, from a
third-order step at 128 bits; at binary512 it is level with MPFR's integer
root. hypot works in P + 64 bits and only falls back to an exact square root
of 4P bits for results that may be exact or ties.

The first call in each format builds constants and reduction tables, once:
0.4 ms for exp and 0.5 ms for log at binary128, 2.2 ms and 3.7 ms at binary512.

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

- **Scalar against MPFR: passed when measured** against an MSVC build of MPFR,
  with mpfix faster on every operation. Against the faster GCC build above it
  still beats MPFR IEEE nearly everywhere, but native MPFR wins some binary128
  operations and binary512 fma.
- **Scalar against Boost, target 5x: not met.** Add is 3.8x to 4.5x and mul is
  3.5x to 6.6x. MPFR is the goal that matters, so this does not block anything.
- **Batch on AVX2, target 1.2x: failed.** Batch stays out of the plan for AVX2.
  AVX-512 IFMA may change the result, but it cannot be measured on this machine.
