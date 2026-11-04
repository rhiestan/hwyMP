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
- Berkeley SoftFloat 3e, for binary128 only: a DLL built with MSYS2's GCC,
  `-O2`, with its default 8086-SSE configuration and
  `SOFTFLOAT_FAST_INT64`. Measured on 2026-11-04; the benchmark first checks
  that its results equal mpfix's on all inputs. Like mpfix, it rounds in five
  modes and raises the IEEE flags.

"MPFR" uses its native wide exponent range, which favours it. "MPFR IEEE" adds
the steps that IEEE-identical results need, matching what mpfix delivers.

An earlier MPFR 4.2.2 build through vcpkg had its C code compiled by MSVC, and
was much slower: at binary128 its add, mul and sqrt took 12.0, 14.8 and 28.5 ns,
against 8.4, 7.3 and 14.0 here. Benchmarks against MPFR on Windows should use a
GCC build.

## Summary

Speedups of mpfix over MPFR and Boost's `cpp_bin_float`, from the tables
below, with Clang. All three give correctly rounded results, mpfix and "MPFR
IEEE" the same ones; Boost has no subnormals or flags.

| | binary128 | binary256 | binary512 |
|---|---|---|---|
| add, mul, div, fma over MPFR IEEE | 1.4x to 2.8x | 1.3x to 2.6x | 1.0x to 1.4x |
| add, mul, div over SoftFloat | 1.6x, 1.7x, 1.2x | | |
| sqrt, fma over SoftFloat | 1.1x, 1.6x | | |
| sqrt over MPFR IEEE | 0.95x | 1.7x | 1.5x |
| add and mul over Boost | 4.5x, 6.5x | 3.8x, 5.1x | 3.8x, 3.5x |
| div and sqrt over Boost | 12x, 14x | 14x, 8.6x | 13x, 7.8x |
| parsing over MPFR | 2.5x | 2.2x | 1.5x |
| printing over MPFR | 1.0x | 1.0x | 0.85x |
| exp and log over MPFR | 2.5x, 4.5x | 2.5x, 3.6x | 1.7x, 2.1x |
| sin, cos, tan over MPFR | 2.5x, 1.8x, 2.8x | 2.3x, 1.7x, 2.5x | 1.6x, 1.2x, 1.8x |
| asin, acos, atan over MPFR | 6.9x, 6.5x, 5.8x | 8.1x, 7.5x, 10x | 5.7x, 5.3x, 7.8x |
| erf, erfc, tgamma over MPFR | 4.5x, 12x, 4.7x | 4.1x, 9.0x, 3.7x | 3.2x, 5.1x, 2.2x |
| exp and log over Boost | 12x, 16x | 26x, 24x | 46x, 28x |

The exceptions are binary128 square root, where MPFR IEEE is 5% faster and
native MPFR 1.5 times; native MPFR's binary128 and binary512 dot step; and
printing at binary512. MSVC builds of mpfix are 1.2 to 3 times slower than
Clang's.

## Scalar

| Format | Library | add | mul | div | sqrt | fma | dot step |
|---|---|---:|---:|---:|---:|---:|---:|
| binary128 | mpfix | 8.20 | 7.53 | 14.93 | 20.49 | 13.07 | 10.92 |
| | MPFR | 8.37 | 7.31 | 16.17 | 13.95 | 30.27 | 7.78 |
| | MPFR IEEE | 13.36 | 13.20 | 20.32 | 19.47 | 36.45 | 13.10 |
| | Boost | 36.66 | 48.75 | 182.81 | 282.39 | | 41.92 |
| | SoftFloat | 12.69 | 13.09 | 17.24 | 22.67 | 20.94 | 12.54 |
| binary256 | mpfix | 12.30 | 14.76 | 30.39 | 63.76 | 33.64 | 15.53 |
| | MPFR | 15.60 | 31.92 | 73.56 | 104.82 | 42.13 | 20.00 |
| | MPFR IEEE | 20.45 | 33.82 | 79.30 | 111.16 | 47.46 | 26.06 |
| | Boost | 46.59 | 75.00 | 424.64 | 547.42 | | 55.02 |
| binary512 | mpfix | 16.73 | 48.93 | 84.06 | 143.84 | 70.45 | 32.56 |
| | MPFR | 17.91 | 57.29 | 112.74 | 216.89 | 75.54 | 31.31 |
| | MPFR IEEE | 22.70 | 62.49 | 117.64 | 220.71 | 81.25 | 36.39 |
| | Boost | 63.36 | 169.23 | 1087.66 | 1121.13 | | 96.05 |

Against MPFR IEEE, which gives the same results, mpfix is faster everywhere
except binary128 square root. Native MPFR is ahead on binary128 square root
(14.0 against 20.5 ns) and on the binary128 and binary512 dot step. The square
root has no special case for two limbs yet.

The fma figures are from 2026-11-04, after tuning. It adds c to the exact
product in a frame of 2L + 1 limbs, set by whichever of the two is larger, so
that it needs no comparison or swap first, takes the sticky bit from the limbs
shifted out, and shifts within a limb, the usual case, without the general
shifts. That took binary128 from 24.2 to 13.1 ns, binary256 from 37.1 to 33.6
and binary512 from 80.3 to 70.5, ahead of native MPFR at every size.

SoftFloat, the reference software implementation of IEEE binary128, is 1.6
and 1.7 times slower than mpfix for add and mul, and 1.1 to 1.2 times for
division, square root and the dot step, and 1.6 times for fma.

The compiler matters: built with the Clang 19.1.5 that Visual Studio 2022
ships, mpfix's binary128 square root took 30 ns instead of 20, the other
operations about the same.

### Full-limb precision

Full-limb formats cost mpfix little: 128, 256 and 512 bits run at the same speed
as binary128, binary256 and binary512, except for division at 128 and 256 bits,
which needs one more quotient limb. MPFR is slower at full limbs, most of all
for square root.

| Format | Library | add | mul | div | sqrt | fma | dot step |
|---|---|---:|---:|---:|---:|---:|---:|
| full 128 | mpfix | 8.32 | 7.58 | 17.57 | 20.22 | 13.29 | 10.66 |
| | MPFR | 7.88 | 13.95 | 36.02 | 101.39 | 32.60 | 9.61 |
| full 256 | mpfix | 12.11 | 14.65 | 35.82 | 64.25 | 33.50 | 15.45 |
| | MPFR | 15.38 | 28.40 | 74.20 | 143.93 | 44.45 | 20.21 |
| full 512 | mpfix | 16.12 | 48.94 | 85.88 | 143.78 | 69.10 | 32.43 |
| | MPFR | 17.65 | 57.43 | 111.38 | 231.01 | 77.16 | 31.23 |

### MSVC

MSVC 19.44 builds of mpfix, `/O2`, remeasured on 2026-11-04, against the Clang
figures above. They are 1.2 to 3 times slower than Clang's, most at binary128,
where they are also slower than MPFR (here MPFR 4.2.2 from vcpkg, built with
MSVC) except in division. Why is still open; Clang is the compiler to use for
speed on Windows.

| Format | Compiler | add | mul | div | sqrt | fma | dot step |
|---|---|---:|---:|---:|---:|---:|---:|
| binary128 | Clang | 8.20 | 7.53 | 14.93 | 20.49 | 24.23 | 10.92 |
| | MSVC | 21.59 | 23.91 | 22.53 | 34.43 | 49.04 | 24.06 |
| binary256 | Clang | 12.30 | 14.76 | 30.39 | 63.76 | 37.06 | 15.53 |
| | MSVC | 23.71 | 29.31 | 45.04 | 77.64 | 58.72 | 27.40 |
| binary512 | Clang | 16.73 | 48.93 | 84.06 | 143.84 | 80.27 | 32.56 |
| | MSVC | 27.58 | 56.58 | 120.51 | 176.40 | 104.75 | 40.29 |

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
| binary128 | mpfix | 416 | 383 | 273 | 393 | 440 | 431 | 300 |
| | MPFR | 1033 | 1305 | 1019 | 1749 | 2087 | 3628 | 2064 |
| | Boost | 5078 | | | 6464 | | | |
| binary256 | mpfix | 698 | 661 | 490 | 826 | 865 | 846 | 652 |
| | MPFR | 1720 | 2045 | 1665 | 2973 | 3498 | 6230 | 5117 |
| | Boost | 18316 | | | 19515 | | | |
| binary512 | mpfix | 1856 | 1663 | 1337 | 2363 | 2475 | 2397 | 1963 |
| | MPFR | 3129 | 3498 | 3057 | 5048 | 5328 | 9730 | 8059 |
| | Boost | 84623 | | | 66530 | | | |

mpfix is 1.7 to 3.7 times faster than MPFR for the exponentials and 2.1 to 8
times faster for the logarithms. The exponential's series runs in fixed point
on raw limbs, (e^r - 1)/r being near 1, which skips the normalization and
rounding of every Float operation and made exp 1.4 to 1.7 times faster.

Most results are decided by the bits of the working value alone: rounding
boundaries all have at most P + 1 bits, so when the bits below those are
neither nearly all zeros nor nearly all ones, the value rounds directly.

C23's exp10, exp2m1, exp10m1, log2p1 and log10p1, with arguments as for exp
and, for the last four, below 1/4 in magnitude:

| Format | Library | exp10 | exp2m1 | exp10m1 | log2p1 | log10p1 |
|---|---|---:|---:|---:|---:|---:|
| binary128 | mpfix | 427 | 258 | 339 | 319 | 320 |
| | MPFR | 2681 | 2648 | 8411 | 2488 | 6813 |
| binary256 | mpfix | 747 | 469 | 586 | 676 | 680 |
| | MPFR | 4741 | 3578 | 12615 | 5469 | 13465 |
| binary512 | mpfix | 1901 | 1566 | 1578 | 2022 | 2029 |
| | MPFR | 9362 | 6741 | 21076 | 8211 | 24434 |

They share the kernels of exp and log1p; MPFR's versions are newer and 4 to
25 times slower. exp10 reduces x log2(10) with log2(10) split in a 62-bit
part, whose product with x is exact, and the rest. The functions are exact,
tied or just beside a rounding boundary at integers and powers of the base,
which they detect.

### Hyperbolic functions

The same benchmark. sinh, cosh, tanh and asinh take arguments as for exp,
acosh one plus values as for log, and atanh arguments below 1/2 in magnitude.

| Format | Library | sinh | cosh | tanh | asinh | acosh | atanh |
|---|---|---:|---:|---:|---:|---:|---:|
| binary128 | mpfix | 509 | 492 | 520 | 516 | 518 | 436 |
| | MPFR | 1307 | 1311 | 1389 | 2195 | 3390 | 4005 |
| | Boost | 5838 | 5910 | 6098 | 7357 | 7297 | 7958 |
| binary256 | mpfix | 812 | 783 | 814 | 998 | 1015 | 853 |
| | MPFR | 1983 | 1911 | 1999 | 3618 | 5236 | 9088 |
| | Boost | 18932 | 18882 | 19395 | 19098 | 18162 | 20202 |
| binary512 | mpfix | 2039 | 2014 | 2119 | 2719 | 2749 | 2531 |
| | MPFR | 3446 | 3438 | 3469 | 5649 | 8161 | 12631 |
| | Boost | 84599 | 83566 | 84817 | 66423 | 70306 | 78130 |

sinh, cosh and tanh take e^|x| from the exponential's kernel and cost about as
much as exp, with one division more; they are 1.6 to 2.7 times faster than
MPFR.
asinh, acosh and atanh are logarithms of arguments formed without
cancellation, through the same tables as log, and 2 to 10 times faster.

### Trigonometric functions

The same benchmark, with arguments as for exp.

| Format | Library | sin | cos | tan |
|---|---|---:|---:|---:|
| binary128 | mpfix | 551 | 551 | 613 |
| | MPFR | 1393 | 1005 | 1711 |
| | Boost | 4556 | 4630 | 9412 |
| binary256 | mpfix | 895 | 897 | 964 |
| | MPFR | 2016 | 1552 | 2371 |
| | Boost | 9200 | 9408 | 19068 |
| binary512 | mpfix | 2168 | 2209 | 2325 |
| | MPFR | 3385 | 2661 | 4150 |
| | Boost | 32546 | 32529 | 66790 |

mpfix is 1.2 to 2.8 times faster than MPFR. Arguments below 2^27 are reduced
by Cody and Waite's method, and larger ones, or those extremely close to
a multiple of pi/2, by Payne and Hanek's, exactly. Two tables of sin and
cos at multiples of 2^-7 and 2^-14 leave an argument below 2^-15, whose sine
series runs in fixed point on raw limbs, without the rounding steps of Float
arithmetic.

Arguments with binary exponents above about 32000 need more of 2/pi than the
first 32768 bits computed on first use. They are computed once, when first
needed, by the Chudnovsky series: 0.3 s for 2^19 bits and 8 s for 2^22,
which binary512's largest arguments need.

### Inverse trigonometric functions

The same benchmark. asin and acos take arguments below 1 in magnitude, atan
arguments as for exp, and atan2(y, x) y as for log and x as for exp.

| Format | Library | asin | acos | atan | atan2 |
|---|---|---:|---:|---:|---:|
| binary128 | mpfix | 411 | 471 | 429 | 365 |
| | MPFR | 2844 | 3052 | 2503 | 2859 |
| | Boost | 10402 | 10517 | 23210 | 7868 |
| binary256 | mpfix | 679 | 767 | 703 | 585 |
| | MPFR | 5514 | 5746 | 7062 | 5226 |
| | Boost | 32886 | 33354 | 61194 | 21616 |
| binary512 | mpfix | 1663 | 1817 | 1603 | 1403 |
| | MPFR | 9402 | 9568 | 12580 | 8754 |
| | Boost | 177122 | 177012 | 263992 | 89099 |

mpfix is 5 to 10 times faster than MPFR. All four are atan of an argument
formed without cancellation: two steps of t = (a - c) / (1 + a c) with tables
of atan at multiples of 2^-7 and 2^-14 leave t below 2^-15, whose series in
t^2 runs in fixed point. Their tables take 0.25 ms at binary128 and 2 ms at
binary512 on the first call.

### Error and gamma functions

The same benchmark, with arguments in (-20, 20).

| Format | Library | erf | erfc | tgamma | lgamma |
|---|---|---:|---:|---:|---:|
| binary128 | mpfix | 2390 | 3737 | 3379 | 4197 |
| | MPFR | 10827 | 44325 | 15791 | 15476 |
| binary256 | mpfix | 6341 | 9472 | 7427 | 9172 |
| | MPFR | 26257 | 85608 | 27752 | 27273 |
| binary512 | mpfix | 25085 | 35817 | 28973 | 35168 |
| | MPFR | 81382 | 182547 | 63952 | 62975 |

erf sums e^-(x^2) times a series of positive terms below 1, and is 1 - erfc
above. erfc takes e^(x^2) erfc x from a table at multiples of 1/16, built
from the asymptotic series downwards, and a Taylor series whose coefficients
follow from its differential equation; above, the asymptotic series alone.
The first erfc builds the table: 0.9 ms at binary128, 10 ms at binary512.
lgamma and tgamma use Stirling's series from x = P/4, with coefficients from
the tangent numbers, shifting smaller arguments up and reflecting negative
ones. Their errors are tracked as they go, so that the cancellation near the
zeros of lgamma passes the result to a wider working format.

### Powers and roots

The same benchmark. pow takes x as for log and y as for exp; cbrt takes values
as for log, and hypot pairs of them.

| Format | Library | pow | cbrt | hypot |
|---|---|---:|---:|---:|
| binary128 | mpfix | 887 | 300 | 101 |
| | MPFR | 3377 | 419 | 229 |
| | Boost | 13649 | | |
| binary256 | mpfix | 1550 | 510 | 174 |
| | MPFR | 5606 | 564 | 501 |
| | Boost | 40898 | | |
| binary512 | mpfix | 4382 | 724 | 372 |
| | MPFR | 8755 | 752 | 503 |
| | Boost | 155526 | | |

pow is e^(y log x) with exact and tied results detected beforehand. cbrt runs
Newton's iteration on the inverse cube root, which needs no division, from a
third-order step at 128 bits; at binary512 it is level with MPFR's integer
root. hypot works in P + 64 bits and only falls back to an exact square root
of 4P bits for results that may be exact or ties.

The first call in each format builds constants and reduction tables, once:
0.15 ms for exp and 0.5 ms for log at binary128, 0.8 ms and 3.7 ms at
binary512.
The first sin also computes pi and 2/pi, for all formats, in 3 ms; its own
tables take 0.4 ms at binary256 and 1 ms at binary512.

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
  still beats MPFR IEEE nearly everywhere, but native MPFR wins binary128
  square root and the dot step.
- **Scalar against Boost, target 5x: not met.** Add is 3.8x to 4.5x and mul is
  3.5x to 6.6x. MPFR is the goal that matters, so this does not block anything.
- **Batch on AVX2, target 1.2x: failed.** Batch stays out of the plan for AVX2.
  AVX-512 IFMA may change the result, but it cannot be measured on this machine.
