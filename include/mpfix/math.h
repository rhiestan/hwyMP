// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Exponentials, logarithms, hyperbolic and trigonometric functions, powers and
// roots: exp, exp2, exp10, expm1, exp2m1, exp10m1, log, log2, log10, log1p,
// log2p1, log10p1, sinh, cosh, tanh, asinh, acosh, atanh, sin, cos, tan, asin,
// acos, atan, atan2, pow, cbrt and hypot, and the constant pi, correctly
// rounded in all five modes, with the IEEE 754 flags when a Context is given.
//
// Each function evaluates in a wider working format with a proven bound on
// its error, then rounds both ends of the interval that bound gives. If they
// round alike, that is the correctly rounded result; if not, the true value
// lies very close to a rounding boundary, and the evaluation is repeated with
// about twice and then four times the precision (Ziv's strategy). Results that
// are exact, such as log2 of a power of two, are detected separately, since
// no interval can decide them. Correct rounding at the last step would only
// fail for an input whose result lies within about 2^-(4P) of a rounding
// boundary, relative to its size; no such input is known, and for random
// inputs the first step almost always decides.
//
// The first working format reduces arguments with tables of about 100
// entries, built on the first call in each format along with the constants:
// about 0.5 ms for binary128 and 4 ms for binary512. Formats with more than 26
// exponent bits are not supported.

#ifndef MPFIX_MATH_H_
#define MPFIX_MATH_H_

#include <algorithm>
#include <array>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <utility>

#include "mpfix/bigint.h"
#include "mpfix/charconv.h"
#include "mpfix/float.h"
#include "mpfix/pi.h"

namespace mpfix {
namespace detail {

// The working formats, with `scale` times one limb more than F.
template <class F, int scale>
using MathWorking =
    Float<static_cast<int>(kLimbBits * scale * (F::kLimbs + 1)), 29>;

constexpr int IntegerSqrt(int n) {
  int root = 0;
  while ((root + 1) * (root + 1) <= n) ++root;
  return root;
}

// 2 atanh(1 / q) * 2^bits, rounded down, with an error below 2^-bits times
// the number of terms.
inline BigUint TwiceAtanhOfInverse(limb_t q, std::size_t bits) {
  BigUint power(1);
  power.ShiftLeft(bits + 1);
  power.DivSmall(q);  // 2 / q^(2k + 1), scaled.
  const limb_t q2 = q * q;
  BigUint sum;
  for (limb_t k = 1; !power.is_zero(); k += 2) {
    BigUint term = power;
    term.DivSmall(k);
    sum.Add(term);
    power.DivSmall(q2);
  }
  return sum;
}

template <class W>
W FromFixedPoint(const BigUint& value, std::size_t bits) {
  if (value.is_zero()) return W::Zero();
  return RoundInteger<W>(false, value, -static_cast<std::int64_t>(bits), false,
                         SilentEnv{Rounding::kNearestEven});
}

// The constants and series coefficients for a working format W. Every value
// is within one rounding, 2^-W::kPrecision relative, of the exact one, except
// as noted.
template <class W>
struct MathConstants {
  static constexpr int kPrecision = W::kPrecision;
  // The exponential series runs on arguments below 2^-kExpScale, reached by
  // halving; each halving costs one squaring afterwards. The series costs
  // about 2 sqrt(n) multiplications for n terms, so few halvings pay.
  static constexpr int kExpScale = (std::max)(2, IntegerSqrt(kPrecision / 8));
  // The logarithm series runs on arguments below 2^-kLogScale, reached by
  // square roots, which cost more than series terms.
  static constexpr int kLogScale = (std::max)(2, IntegerSqrt(kPrecision / 7));
  static constexpr int kLogTerms =
      (kPrecision + 2 + 2 * kLogScale - 1) / (2 * kLogScale);

  W ln2;
  // ln2_hi has kPrecision - 32 bits, so k * ln2_hi is exact for |k| < 2^32,
  // and ln2_hi + ln2_lo approximates ln 2 to about 2^-(2 kPrecision - 32).
  W ln2_hi;
  W ln2_lo;
  W inv_ln2;   // Within two roundings.
  W ln10;
  W inv_ln10;  // Within two roundings.
  // log2(10) = log2_10_hi + log2_10_lo, with a 62-bit high part, so that
  // x log2_10_hi is exact for x of up to kPrecision - 62 bits; the sum is
  // within 2^-(2 kPrecision) or so.
  W log2_10_hi;
  W log2_10_lo;
  std::array<W, kLogTerms> inv_odd;  // 1 / (2i + 1).

  MathConstants() {
    // ln 2 = 2 atanh(1/3), to twice the precision for the split.
    const std::size_t bits = 2 * static_cast<std::size_t>(kPrecision) + 64;
    const BigUint ln2_fixed = TwiceAtanhOfInverse(3, bits);
    ln2 = FromFixedPoint<W>(ln2_fixed, bits);
    BigUint hi = ln2_fixed;
    const std::size_t drop = ln2_fixed.bit_length() - (kPrecision - 32);
    hi.ShiftRight(drop);
    hi.ShiftLeft(drop);
    BigUint lo = ln2_fixed;
    lo.Sub(hi);
    ln2_hi = FromFixedPoint<W>(hi, bits);
    ln2_lo = FromFixedPoint<W>(lo, bits);
    const W one(1);
    inv_ln2 = one / ln2;

    // ln 10 = 3 ln 2 + ln(5/4), and ln(5/4) = 2 atanh(1/9).
    BigUint ln10_fixed = ln2_fixed;
    ln10_fixed.MulAdd(3, 0);
    ln10_fixed.Add(TwiceAtanhOfInverse(9, bits));
    ln10 = FromFixedPoint<W>(ln10_fixed, bits);
    inv_ln10 = one / ln10;

    BigUint scaled = ln10_fixed;
    scaled.ShiftLeft(bits);
    BigUint log2_10_fixed;
    BigUint remainder;
    BigUint::DivRem(scaled, ln2_fixed, log2_10_fixed, remainder);
    BigUint high = log2_10_fixed;
    const std::size_t high_drop = high.bit_length() - 62;
    high.ShiftRight(high_drop);
    high.ShiftLeft(high_drop);
    BigUint low = log2_10_fixed;
    low.Sub(high);
    log2_10_hi = FromFixedPoint<W>(high, bits);
    log2_10_lo = FromFixedPoint<W>(low, bits);

    for (int i = 0; i < kLogTerms; ++i) {
      inv_odd[static_cast<std::size_t>(i)] = one / W(2 * i + 1);
    }
  }
};

template <class W>
const MathConstants<W>& GetMathConstants() {
  static const MathConstants<W> constants;
  return constants;
}

// --- Kernels, in a working format W with unit roundoff u = 2^-kPrecision ---

// x * m, rounded to nearest, for 0 < m < 2^63. Takes linear time.
template <class W>
W MulSmall(const W& x, limb_t m) {
  if (x.is_zero() || m == 1) return x;
  constexpr std::size_t kCount = W::kLimbs + 1;
  std::array<limb_t, kCount> limbs;
  limb_t carry = 0;
  for (std::size_t i = 0; i < W::kLimbs; ++i) {
    limb_t high;
    limb_t low = mul_wide(x.raw_limbs()[i], m, high);
    carry_t c = 0;
    limbs[i] = adc(low, carry, c);
    carry = high + c;
  }
  limbs[kCount - 1] = carry;
  // x * m < 2^(64 kCount), and its top limb is nonzero for m > 1.
  const int shift = countl_zero(carry);
  if (shift != 0) {
    for (std::size_t i = kCount; i-- > 1;) {
      limbs[i] = (limbs[i] << shift) | (limbs[i - 1] >> (kLimbBits - shift));
    }
    limbs[0] <<= shift;
  }
  return FloatAccess::Round<W>(x.signbit(),
                               std::int64_t{x.raw_exponent()} + kLimbBits - shift,
                               limbs.data(), kCount, false, NearestEnv{});
}

// x / m, rounded to nearest, for 0 < m < 2^63. Takes linear time.
template <class W>
W DivSmall(const W& x, limb_t m) {
  if (x.is_zero() || m == 1) return x;
  constexpr std::size_t kCount = W::kLimbs + 1;
  std::array<limb_t, kCount> limbs;
  limbs[0] = 0;
  for (std::size_t i = 0; i < W::kLimbs; ++i) limbs[i + 1] = x.raw_limbs()[i];
  const bool sticky = DivideLimbs(limbs.data(), kCount, m) != 0;
  // The quotient is at least 2^(64 kCount - 1) / m >= 2^(64 kCount - 64).
  const int shift = countl_zero(limbs[kCount - 1]);
  if (shift != 0) {
    for (std::size_t i = kCount; i-- > 1;) {
      limbs[i] = (limbs[i] << shift) | (limbs[i - 1] >> (kLimbBits - shift));
    }
    limbs[0] <<= shift;
  }
  return FloatAccess::Round<W>(x.signbit(),
                               std::int64_t{x.raw_exponent()} - shift,
                               limbs.data(), kCount, sticky, NearestEnv{});
}

// Fixed-point numbers v 2^(64 kN) in kN + 1 limbs, least significant first,
// in two's complement, for the series of the exponential and the sine: their
// sums are near 1, so that an absolute error is a relative one, and raw limbs
// skip the normalization and rounding of every step.
template <std::size_t kN>
using FixedLimbs = std::array<limb_t, kN + 1>;

// a b for nonnegative a and b, truncated: within 1.01 2^-64kN below. The
// partial products below limb kN - 2 of the full product are left out: each
// is below 2^(64 (kN - 1)), and there are fewer than 2^64 of them.
template <std::size_t kN>
FixedLimbs<kN> FixedMul(const FixedLimbs<kN>& a, const FixedLimbs<kN>& b) {
  constexpr std::size_t kCount = kN + 1;
  constexpr std::size_t kFirstColumn = kN >= 2 ? kN - 2 : 0;
  std::array<limb_t, 2 * kCount> full{};
  for (std::size_t i = 0; i < kCount; ++i) {
    if (a[i] == 0) continue;
    limb_t carry = 0;
    for (std::size_t j = kFirstColumn > i ? kFirstColumn - i : 0; j < kCount;
         ++j) {
      limb_t hi;
      limb_t lo = mul_wide(a[i], b[j], hi);
      carry_t c = 0;
      lo = adc(lo, full[i + j], c);
      hi += c;
      c = 0;
      lo = adc(lo, carry, c);
      full[i + j] = lo;
      carry = hi + c;
    }
    full[i + kCount] = carry;
  }
  FixedLimbs<kN> result;
  for (std::size_t i = 0; i < kCount; ++i) result[i] = full[i + kN];
  return result;
}

// a m for a nonnegative a and a product that fits.
template <std::size_t kN>
FixedLimbs<kN> FixedMulSmall(const FixedLimbs<kN>& a, limb_t m) {
  FixedLimbs<kN> result;
  limb_t carry = 0;
  for (std::size_t i = 0; i <= kN; ++i) {
    limb_t hi;
    limb_t lo = mul_wide(a[i], m, hi);
    carry_t c = 0;
    result[i] = adc(lo, carry, c);
    carry = hi + c;
  }
  return result;
}

// acc += x, or acc -= x.
template <std::size_t kN>
void FixedAccumulate(FixedLimbs<kN>& acc, const FixedLimbs<kN>& x,
                     bool subtract) {
  carry_t carry = 0;
  for (std::size_t i = 0; i <= kN; ++i) {
    acc[i] = subtract ? sbb(acc[i], x[i], carry) : adc(acc[i], x[i], carry);
  }
}

template <std::size_t kN>
void FixedNegate(FixedLimbs<kN>& a) {
  carry_t carry = 1;
  for (std::size_t i = 0; i <= kN; ++i) a[i] = adc(~a[i], 0, carry);
}

// A nonnegative v < 2^63 of the working format, truncated: within 2^-64kN.
template <std::size_t kN, class W>
FixedLimbs<kN> ToFixed(const W& v) {
  static_assert(W::kLimbs == kN);
  FixedLimbs<kN> result{};
  if (v.is_zero()) return result;
  // v = limbs 2^(e + 1 - 64 kN), so v 2^(64 kN) = limbs 2^(e + 1).
  for (std::size_t i = 0; i < kN; ++i) result[i] = v.raw_limbs()[i];
  const int shift = v.raw_exponent() + 1;
  if (shift > 0) {
    for (std::size_t i = kN; i > 0; --i) {
      result[i] = (result[i] << shift) | (result[i - 1] >> (kLimbBits - shift));
    }
    result[0] <<= shift;
  } else if (shift < 0) {
    const int right = -shift;
    const std::size_t words = static_cast<std::size_t>(right / kLimbBits);
    const int bits = right % kLimbBits;
    for (std::size_t i = 0; i <= kN; ++i) {
      const std::size_t from = i + words;
      const limb_t low = from <= kN ? result[from] : 0;
      const limb_t high = from + 1 <= kN ? result[from + 1] : 0;
      result[i] = bits == 0 ? low : (low >> bits) | (high << (kLimbBits - bits));
    }
  }
  return result;
}

// A positive fixed-point value, rounded to the working format.
template <class W, std::size_t kN>
W FromFixed(FixedLimbs<kN> a) {
  std::size_t top = kN;
  while (a[top] == 0) --top;
  const int shift = countl_zero(a[top]);
  // The value's top bit has weight 2^(64 (top - kN) + 63 - shift).
  const std::int64_t exponent =
      static_cast<std::int64_t>(kLimbBits) *
          (static_cast<std::int64_t>(top) - static_cast<std::int64_t>(kN)) +
      (kLimbBits - 1 - shift);
  if (shift != 0) {
    for (std::size_t i = top; i > 0; --i) {
      a[i] = (a[i] << shift) | (a[i - 1] >> (kLimbBits - shift));
    }
    a[0] <<= shift;
  }
  return FloatAccess::Round<W>(false, exponent, a.data(), top + 1, false,
                               NearestEnv{});
}

// Paterson and Stockmeyer's evaluation of the sum over j < terms of c_j y^j,
// for c_0 = 1 and c_j = c_(j-1) / d(j), with (-1)^j more if `alternate`, in
// fixed point, from powers[i] = y^i for i <= k and 0 <= y < 1/4. In blocks of
// k terms, the coefficients are small integers over one common denominator,
// so that only the powers and one product per block are full
// multiplications. With acc_b the sum from term bk on over its first
// coefficient, acc_b = 1 + T / D for T = sum over 0 < i < k of +-M_i y^i
// +- y^k acc_(b+1), D = d(bk + 1) ... d(bk + k) and
// M_i = d(bk + i + 1) ... d(bk + k); k factors d(j) must fit in 63 bits.
//
// Each power, product and quotient truncates by below 2^-64kN, and the
// growth of the multipliers is undone by the division, so that the sum is
// within 4 2^-64kN as long as it stays above 1/2.
template <std::size_t kN, class D>
FixedLimbs<kN> FixedSeries(const FixedLimbs<kN>* powers, int terms, int k,
                           const D& d, bool alternate) {
  FixedLimbs<kN> one{};
  one[kN] = 1;
  const int blocks = (terms + k - 1) / k;
  FixedLimbs<kN> acc = one;
  for (int block = blocks - 1; block >= 0; --block) {
    const int base = block * k;
    const int size = (std::min)(k, terms - base);
    const bool has_next = block < blocks - 1;
    if (size < 2 && !has_next) {
      acc = one;
      continue;
    }
    FixedLimbs<kN> t{};
    if (has_next) {
      FixedAccumulate<kN>(t, FixedMul<kN>(powers[size], acc),
                          alternate && (size & 1) != 0);
    }
    limb_t multiplier = 1;
    for (int i = size - 1; i >= 1; --i) {
      if (has_next || i < size - 1) multiplier *= d(base + i + 1);
      FixedAccumulate<kN>(t, FixedMulSmall<kN>(powers[i], multiplier),
                          alternate && (i & 1) != 0);
    }
    // T may be negative; divide its magnitude.
    const bool negative = (t[kN] >> (kLimbBits - 1)) != 0;
    if (negative) FixedNegate<kN>(t);
    DivideLimbs(t.data(), kN + 1, multiplier * d(base + 1));
    acc = one;
    FixedAccumulate<kN>(acc, t, negative);
  }
  return acc;
}

// The block size for `terms` terms whose largest factor d has d_bits bits.
inline int SeriesBlock(int terms, int d_bits) {
  constexpr int kMaxBlock = 16;
  return (std::max)(
      1, (std::min)({IntegerSqrt(terms), 63 / d_bits, kMaxBlock}));
}

// e^r - 1 for |r| < 1/2, with a relative error below `units` * u. The error
// of r itself is not included.
template <class W>
W ExpM1Kernel(const W& r, int& units) {
  using Constants = MathConstants<W>;
  if (r.is_zero()) {
    units = 0;
    return r;
  }
  // Halve r below 2^-kExpScale; each halving is undone by one squaring.
  const int halvings =
      (std::max)(0, r.raw_exponent() + 1 + Constants::kExpScale);
  const W x = scalbn(r, -halvings);
  // The first omitted term, relative to x, is below |x|^n / (n + 1)!.
  const int scale = -(x.raw_exponent() + 1);
  int n = 1;
  int log2_factorial = 1;  // A lower bound on log2((n + 1)!).
  while (n * scale + log2_factorial < Constants::kPrecision + 3) {
    ++n;
    log2_factorial += std::bit_width(static_cast<unsigned>(n + 1)) - 1;
  }

  // e^x - 1 = x H(x) for H(x) = sum over i < n of x^i / (i + 1)!, which is
  // within 2^-(kPrecision + 3) of the true value, and between 0.88 and 1.14.
  // In fixed point, it is within 4 u absolute, and so 5 u relative; rounding
  // it and the product add 2 u.
  constexpr std::size_t kN = W::kLimbs;
  const auto d = [](int j) { return static_cast<limb_t>(j + 1); };
  const int k = SeriesBlock(n, std::bit_width(d(n)));
  std::array<FixedLimbs<kN>, 17> powers;
  powers[1] = ToFixed<kN>(abs(x));
  for (int i = 2; i <= k; ++i) {
    powers[static_cast<std::size_t>(i)] =
        FixedMul<kN>(powers[static_cast<std::size_t>(i - 1)], powers[1]);
  }
  const W acc = FromFixed<W, kN>(
      FixedSeries<kN>(powers.data(), n, k, d, x.signbit()));
  W q = x * acc;
  // e^(2a) - 1 = q (q + 2) for q = e^a - 1, which keeps the relative error
  // of a small result.
  const W two(2);
  for (int i = 0; i < halvings; ++i) q = q * (q + two);
  // The series stays within 8 u, and each squaring adds 2 u; for |r| <= 1/2
  // the squarings at most double the earlier error.
  units = 2 * (10 + 2 * halvings);
  return q;
}

// log(1 + d) for -0.3 < d < 0.42, with a relative error below `units` * u.
// d is exact.
template <class W>
W Log1pKernel(W d, int& units) {
  using Constants = MathConstants<W>;
  if (d.is_zero()) {
    units = 0;
    return d;
  }
  const Constants& c = GetMathConstants<W>();
  const W one(1);
  // log(1 + d) = 2 log(1 + d') with d' = d / (1 + sqrt(1 + d)), which keeps
  // the relative error of d' within 3 u more than that of d.
  int roots = 0;
  while (d.raw_exponent() >= -Constants::kLogScale) {
    d = d / (one + sqrt(one + d));
    ++roots;
  }
  // log(1 + d) = 2 atanh(t) for t = d / (2 + d).
  const W t = d / (W(2) + d);
  const W t2 = t * t;
  const int scale = -(t.raw_exponent() + 1);
  const int n = (Constants::kPrecision + 2 + 2 * scale - 1) / (2 * scale);
  W p = c.inv_odd[static_cast<std::size_t>(n - 1)];
  for (int i = n - 2; i >= 0; --i) {
    p = p * t2 + c.inv_odd[static_cast<std::size_t>(i)];
  }
  units = 4 * roots + 2 * n + 16;
  return scalbn(t * p, roots + 1);
}

// Tables for reducing r in [-0.36, 0.36] to below 2^-15 in two steps, as
// r = i/2^7 + j/2^14 + r2. They hold e^(i/2^7) - 1 and e^(j/2^14) - 1, which
// keeps the relative error of small results, computed in twice the precision
// and rounded once.
template <class W>
struct ExpTables {
  static constexpr int kOffset1 = 47;  // i in [-46, 46].
  static constexpr int kOffset2 = 65;  // j in [-64, 64].
  std::array<W, 2 * kOffset1 + 1> expm1_1;
  std::array<W, 2 * kOffset2 + 1> expm1_2;

  ExpTables() {
    using Wide = Float<2 * W::kPrecision, 29>;
    for (int i = -kOffset1; i <= kOffset1; ++i) {
      int units;
      expm1_1[static_cast<std::size_t>(i + kOffset1)] =
          FloatAccess::Convert<W>(ExpM1Kernel(scalbn(Wide(i), -7), units),
                                  NearestEnv{});
    }
    for (int j = -kOffset2; j <= kOffset2; ++j) {
      int units;
      expm1_2[static_cast<std::size_t>(j + kOffset2)] =
          FloatAccess::Convert<W>(ExpM1Kernel(scalbn(Wide(j), -14), units),
                                  NearestEnv{});
    }
  }
};

template <class W>
const ExpTables<W>& GetExpTables() {
  static const ExpTables<W> tables;
  return tables;
}

// e^r - 1 for |r| <= 0.36, with a relative error below `units` * u. The
// error of r itself is not included.
template <class W, bool kUseTables>
W ExpM1Reduced(const W& r, int& units) {
  if constexpr (kUseTables) {
    const ExpTables<W>& t = GetExpTables<W>();
    // Both subtractions are exact.
    const int i = to_integer<int>(scalbn(r, 7), Rounding::kNearestEven);
    const W r1 = i != 0 ? r - scalbn(W(i), -7) : r;
    const int j = to_integer<int>(scalbn(r1, 14), Rounding::kNearestEven);
    const W r2 = j != 0 ? r1 - scalbn(W(j), -14) : r1;
    W q = ExpM1Kernel(r2, units);
    // e^(a + b) - 1 = A + B + A B for A = e^a - 1 and B = e^b - 1. Each
    // step cancels at most half of its larger term, since the reduced
    // arguments are each at most half the size of the table step.
    const auto combine = [&q](const W& a) { q = (a + q) + a * q; };
    if (j != 0) combine(t.expm1_2[static_cast<std::size_t>(j + t.kOffset2)]);
    if (i != 0) combine(t.expm1_1[static_cast<std::size_t>(i + t.kOffset1)]);
    if (i != 0 || j != 0) units = 3 * units + 16;
    return q;
  } else {
    return ExpM1Kernel(r, units);
  }
}

// Tables for reducing m in [0.705, 1.41) to near 1 in two steps, as
// m c1 c2 = 1 + d2 with |d2| < 2^-13.9: c1 = 1 / (1 + i/2^7) to 12 bits and
// c2 = 1 / (1 + j/2^13) to 20 bits, so that for an m of at most
// W::kPrecision - 34 bits both products are exact. The logarithms are
// computed in twice the precision and rounded once.
template <class W>
struct LogTables {
  static constexpr int kOffset1 = 53;  // i in [-38, 52].
  static constexpr int kOffset2 = 64;  // j in [-47, 47].
  std::array<W, 2 * kOffset1 + 1> c1;
  std::array<W, 2 * kOffset1 + 1> minus_log_c1;
  std::array<W, 2 * kOffset2 + 1> c2;
  std::array<W, 2 * kOffset2 + 1> minus_log_c2;

  LogTables();
};

template <class W, bool kUseTables>
W LogOfReduced(const W& m, int& units);

template <class W>
LogTables<W>::LogTables() {
  using Wide = Float<2 * W::kPrecision, 29>;
  const auto fill = [](int offset, int index_bits, int value_bits, W* c,
                       W* minus_log_c) {
    for (int i = -offset; i <= offset; ++i) {
      const double scale = std::ldexp(1.0, value_bits);
      const double value =
          std::ldexp(std::nearbyint(scale / (1.0 + std::ldexp(i, -index_bits))),
                     -value_bits);
      c[i + offset] = W(value);
      int units;
      minus_log_c[i + offset] = FloatAccess::Convert<W>(
          -LogOfReduced<Wide, false>(Wide(value), units), NearestEnv{});
    }
  };
  fill(kOffset1, 7, 12, c1.data(), minus_log_c1.data());
  fill(kOffset2, 13, 20, c2.data(), minus_log_c2.data());
}

template <class W>
const LogTables<W>& GetLogTables() {
  static const LogTables<W> tables;
  return tables;
}

// log(m) for m in [0.705, 1.41), with a relative error below `units` * u.
// With tables, m must have at most W::kPrecision - 34 bits.
template <class W, bool kUseTables>
W LogOfReduced(const W& m, int& units) {
  const W one(1);
  W d = m - one;  // Exact.
  if constexpr (kUseTables) {
    const LogTables<W>& t = GetLogTables<W>();
    W product = m;
    W tail = W::Zero();
    const int i = to_integer<int>(scalbn(d, 7), Rounding::kNearestEven);
    if (i != 0) {
      const auto index = static_cast<std::size_t>(i + t.kOffset1);
      product = product * t.c1[index];
      tail = t.minus_log_c1[index];
      d = product - one;
    }
    const int j = to_integer<int>(scalbn(d, 13), Rounding::kNearestEven);
    if (j != 0) {
      const auto index = static_cast<std::size_t>(j + t.kOffset2);
      product = product * t.c2[index];
      tail = tail + t.minus_log_c2[index];
      d = product - one;
    }
    W result = Log1pKernel(d, units);
    if (i != 0 || j != 0) {
      // Each table entry is at most 2.4 times |log m|, and log(1 + d) at
      // most 1.3 times.
      result = result + tail;
      units = 2 * units + 12;
    }
    return result;
  } else {
    return Log1pKernel(d, units);
  }
}

// log(z) for a finite positive z, with a relative error below `units` * u
// plus that of z. Sets `binary_exponent` and `log_m` to e and log(m) with
// z = 2^e m, m in [0.705, 1.41). With tables, z must have at most
// W::kPrecision - 34 bits.
template <class W, bool kUseTables = false>
W LogKernel(const W& z, int& units, std::int64_t& binary_exponent, W& log_m) {
  const MathConstants<W>& c = GetMathConstants<W>();
  std::int32_t e = z.raw_exponent();
  W m = scalbn(z, -e);
  if (m >= W(1.41)) {
    m = scalbn(m, -1);
    ++e;
  }
  log_m = LogOfReduced<W, kUseTables>(m, units);
  binary_exponent = e;
  if (e == 0) return log_m;
  // |e ln 2 + log(m)| >= |e ln 2| / 2, so the sum adds at most 6 u.
  units += 6;
  return W(e) * c.ln2 + log_m;
}

// x = k ln 2 + r with |r| < 0.35; the error of r is below 0.36 u absolute.
template <class W, class F>
W ReduceByLn2(const F& x, std::int64_t& k) {
  const MathConstants<W>& c = GetMathConstants<W>();
  const W wide = FloatAccess::Convert<W>(x, NearestEnv{});
  const W multiple = round(wide * c.inv_ln2);
  k = to_integer<std::int64_t>(multiple, Rounding::kNearestEven);
  // x - k ln2_hi is exact: both are multiples of 2^(32 - kPrecision) and
  // their difference is below 1.
  return (wide - multiple * c.ln2_hi) - multiple * c.ln2_lo;
}

// --- Rounding ---

// Whether bits [low, high) of the limbs, counted from bit 0 of limbs[0], are
// all equal to `ones`.
template <std::size_t kCount>
bool BitsAllEqual(const std::array<limb_t, kCount>& limbs, int low, int high,
                  bool ones) {
  const limb_t fill = ones ? ~limb_t{0} : 0;
  for (int bit = low; bit < high;) {
    const int word = bit / kLimbBits;
    const int offset = bit % kLimbBits;
    const int take = (std::min)(kLimbBits - offset, high - bit);
    const limb_t mask = (take == kLimbBits ? ~limb_t{0}
                                           : ((limb_t{1} << take) - 1))
                        << offset;
    if ((limbs[static_cast<std::size_t>(word)] & mask) != (fill & mask)) {
      return false;
    }
    bit += take;
  }
  return true;
}

// Rounds a value known to lie within 2^-bits relative of y and to be
// inexact, if no number of at most P + 1 bits lies that close to y. Every
// rounding boundary has at most P + 1 bits, in every mode and on the
// subnormal grid too, so y then rounds as the value does. Returns false if
// it cannot tell.
template <class F, class W, class Env>
bool RoundClearOfBoundaries(const W& y, int bits, const Env& env, F& result) {
  constexpr int kTotal = static_cast<int>(W::kLimbs) * kLimbBits;
  // The bits below the first P + 1, and the error in units of the lowest
  // limb bit, below 2^slack.
  constexpr int kBelow = kTotal - (F::kPrecision + 1);
  const int slack = kTotal - bits;
  if (slack + 1 >= kBelow) return false;
  const auto& limbs = y.raw_limbs();
  if (BitsAllEqual(limbs, slack, kBelow, false) ||
      BitsAllEqual(limbs, slack, kBelow, true)) {
    return false;
  }
  result = FloatAccess::Round<F>(y.signbit(), y.raw_exponent(), limbs.data(),
                                 W::kLimbs, true, env);
  return true;
}

// Rounds a value known to lie within 2^-bits relative of y, and known to be
// inexact. Returns false if the interval does not decide the result, unless
// `last` is set, in which case y is rounded.
template <class F, class W, class Env>
bool RoundTranscendental(const W& y, int bits, bool last, const Env& env,
                         F& result) {
  if (RoundClearOfBoundaries(y, bits, env, result)) return true;
  Context down{Rounding::kDown, 0};
  Context up{Rounding::kUp, 0};
  const W delta = scalbn(abs(y), -bits);
  const W low = last ? y : sub(y, delta, down);
  const W high = last ? y : add(y, delta, up);
  const SilentEnv rounding{env.rounding()};

  // Rounded to P bits with a wider exponent range, which decides overflow
  // and tininess after rounding, as IEEE 754 defines them.
  using Unbounded = Float<F::kPrecision, F::kExponentBits + 1>;
  const Unbounded wide = FloatAccess::Convert<Unbounded>(low, rounding);
  if (wide.to_bits() !=
      FloatAccess::Convert<Unbounded>(high, rounding).to_bits()) {
    return false;
  }
  const Unbounded magnitude = abs(wide);
  if (magnitude > Unbounded(F::MaxFinite())) {
    result = FloatAccess::Convert<F>(wide, env);  // Overflow and inexact.
    return true;
  }
  if (magnitude >= Unbounded(F::MinNormal())) {
    result = FloatAccess::Convert<F>(wide, rounding);
    env.raise(kFlagInexact);
    return true;
  }
  // Tiny: rounding to the subnormal grid must be decided on its own.
  const F a = FloatAccess::Convert<F>(low, rounding);
  if (a.to_bits() != FloatAccess::Convert<F>(high, rounding).to_bits()) {
    return false;
  }
  env.raise(kFlagInexact | kFlagUnderflow);
  result = a;
  return true;
}

// Evaluates in up to three working formats. `level(tag, last, result)` works
// in the format of the std::type_identity tag.
// Tests set kFirstLevel to start at the second or third format.
template <class F, int kFirstLevel, class Level>
F ZivEvaluate(const Level& level) {
  F result;
  if constexpr (kFirstLevel <= 1) {
    if (level(std::type_identity<MathWorking<F, 1>>{}, false, result)) {
      return result;
    }
  }
  if constexpr (kFirstLevel <= 2) {
    if (level(std::type_identity<MathWorking<F, 2>>{}, false, result)) {
      return result;
    }
  }
  level(std::type_identity<MathWorking<F, 4>>{}, true, result);
  return result;
}

// The lost bits for an error of `units` u, with a margin.
constexpr int LostBits(int units) {
  return std::bit_width(static_cast<unsigned>(units)) + 2;
}

// Rounds v + t for a nonzero finite v and a t with sign `t_positive` and
// |t| < |v| 2^-(P + 2): t only matters in directed modes, but it always makes
// the result inexact.
template <class F, class Env>
F RoundBeside(const F& v, bool t_positive, const Env& env) {
  using W = MathWorking<F, 1>;
  const W wide = abs(FloatAccess::Convert<W>(v, NearestEnv{}));
  const bool negative = v.signbit();
  // Rounding with a sticky bit adds an infinitesimal to the magnitude.
  W magnitude = wide;
  if (t_positive == negative) {
    // The magnitude shrinks: start below it instead. This is exact.
    magnitude = wide - scalbn(wide, -(F::kPrecision + 2));
  }
  return FloatAccess::Round<F>(negative, magnitude.raw_exponent(),
                               magnitude.raw_limbs().data(), W::kLimbs, true,
                               env);
}

// A value far below the smallest subnormal.
template <class F, class Env>
F RoundTiny(bool negative, const Env& env) {
  const limb_t top_bit = limb_t{1} << (kLimbBits - 1);
  return FloatAccess::Round<F>(negative,
                               std::int64_t{F::kEmin} - F::kPrecision - 4,
                               &top_bit, 1, true, env);
}

// The lost bits for an error of `units` u, with a margin.
inline int LostBitsWide(std::int64_t units) {
  return std::bit_width(static_cast<std::uint64_t>(units)) + 2;
}

// Rounds w, a value truncated towards zero, to F with the given sign. `sticky`
// tells whether the truncation lost anything. W must have at least P + 2 bits.
template <class F, class W, class Env>
F RoundTruncated(bool negative, const W& w, bool sticky, const Env& env) {
  return FloatAccess::Round<F>(negative, w.raw_exponent(),
                               w.raw_limbs().data(), W::kLimbs, sticky, env);
}

// Whether a finite x is an integer.
template <class F>
bool IsInteger(const F& x) {
  return x.is_zero() || LowestBitExponent(x) >= 0;
}

// The number of significant bits of a finite nonzero x.
template <class F>
std::int64_t SignificantBits(const F& x) {
  return std::int64_t{x.raw_exponent()} - LowestBitExponent(x) + 1;
}

template <class F, class Env>
F RoundOverflow(bool negative, const Env& env) {
  const limb_t top_bit = limb_t{1} << (kLimbBits - 1);
  return FloatAccess::Round<F>(negative, std::int64_t{F::kEmax} + 1, &top_bit, 1,
                               false, env);
}

template <class F, class Env>
F InvalidResult(const Env& env) {
  env.raise(kFlagInvalid);
  return F::QuietNaN();
}

template <class F>
using Wide1 = MathWorking<F, 1>;

// --- The functions ---

template <class F, class Env, int kFirstLevel = 1>
F ExpImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_inf()) return x.signbit() ? F::Zero() : x;
  if (x.is_zero()) return F(1);
  if (x.raw_exponent() <= -(F::kPrecision + 4)) {
    return RoundBeside(F(1), !x.signbit(), env);
  }
  const Wide1<F> wide = FloatAccess::Convert<Wide1<F>>(x, NearestEnv{});
  // e^(kEmax + 2) overflows, and e^-(kEmax + P + 4) is below a quarter of the
  // smallest subnormal.
  if (wide >= Wide1<F>(F::kEmax + 2)) return RoundOverflow<F>(false, env);
  if (wide <= Wide1<F>(-(F::kEmax + F::kPrecision + 4))) {
    return RoundTiny<F>(false, env);
  }
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    std::int64_t k;
    const W r = ReduceByLn2<W>(x, k);
    int units;
    const W q =
        ExpM1Reduced<W, std::is_same_v<W, MathWorking<F, 1>>>(r, units);
    const W y = scalbn(q + W(1), static_cast<int>(k));
    return RoundTranscendental(y, W::kPrecision - LostBits(units + 4), last,
                               env, result);
  });
}

template <class F, class Env, int kFirstLevel = 1>
F Exp2Impl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_inf()) return x.signbit() ? F::Zero() : x;
  if (x.raw_exponent() <= -(F::kPrecision + 4) && !x.is_zero()) {
    return RoundBeside(F(1), !x.signbit(), env);
  }
  const Wide1<F> wide = FloatAccess::Convert<Wide1<F>>(x, NearestEnv{});
  if (wide >= Wide1<F>(F::kEmax + 1)) return RoundOverflow<F>(false, env);
  // Below 2^(kEmin - P - 1), a quarter of the smallest subnormal.
  if (wide < Wide1<F>(F::kEmin - F::kPrecision - 1)) {
    return RoundTiny<F>(false, env);
  }
  if (trunc(x) == x) {
    // Exact, or a subnormal that rounds away.
    const limb_t top_bit = limb_t{1} << (kLimbBits - 1);
    return FloatAccess::Round<F>(
        false, to_integer<std::int64_t>(x, Rounding::kTowardZero), &top_bit, 1,
        false, env);
  }
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    const MathConstants<W>& c = GetMathConstants<W>();
    const W wx = FloatAccess::Convert<W>(x, NearestEnv{});
    const W nearest = round(wx);
    const std::int64_t k =
        to_integer<std::int64_t>(nearest, Rounding::kNearestEven);
    const W r = (wx - nearest) * c.ln2;  // Exact difference, |r| <= 0.35.
    int units;
    const W q =
        ExpM1Reduced<W, std::is_same_v<W, MathWorking<F, 1>>>(r, units);
    const W y = scalbn(q + W(1), static_cast<int>(k));
    return RoundTranscendental(y, W::kPrecision - LostBits(units + 6), last,
                               env, result);
  });
}

template <class F, class Env, int kFirstLevel = 1>
F Expm1Impl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_inf()) return x.signbit() ? F(-1) : x;
  if (x.is_zero()) return x;
  // e^x - 1 = x + x^2/2 + ..., always above x.
  if (x.raw_exponent() <= -(F::kPrecision + 3)) {
    return RoundBeside(x, true, env);
  }
  const Wide1<F> wide = FloatAccess::Convert<Wide1<F>>(x, NearestEnv{});
  if (wide >= Wide1<F>(F::kEmax + 2)) return RoundOverflow<F>(false, env);
  // Then e^x < 2^-(P + 3), and e^x - 1 is just above -1.
  if (wide <= Wide1<F>(-(F::kPrecision + 3))) {
    return RoundBeside(F(-1), true, env);
  }
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    std::int64_t k;
    const W r = ReduceByLn2<W>(x, k);
    int units;
    W y = ExpM1Reduced<W, std::is_same_v<W, MathWorking<F, 1>>>(r, units);
    if (k != 0) {
      // |2^k e^r| <= 3.4 |2^k e^r - 1| for |x| >= ln 2 / 2.
      y = scalbn(y + W(1), static_cast<int>(k)) - W(1);
      units = 4 * units + 8;
    }
    return RoundTranscendental(y, W::kPrecision - LostBits(units + 4), last,
                               env, result);
  });
}

// The special cases shared by the logarithms of x. Returns true if they
// decide the result.
template <class F, class Env>
bool LogSpecialCase(const F& x, const Env& env, F& result) {
  if (x.is_nan()) {
    result = FloatAccess::PropagateNaN(x, env);
  } else if (x.is_zero()) {
    env.raise(kFlagDivByZero);
    result = F::Infinity(true);
  } else if (x.signbit()) {
    result = InvalidResult<F>(env);
  } else if (x.is_inf()) {
    result = x;
  } else if (x == F(1)) {
    result = F::Zero();
  } else {
    return false;
  }
  return true;
}

template <class F, class Env, int kFirstLevel = 1>
F LogImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  F special;
  if (LogSpecialCase(x, env, special)) return special;
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    int units;
    std::int64_t e;
    W log_m;
    const W y = LogKernel<W, std::is_same_v<W, MathWorking<F, 1>>>(
        FloatAccess::Convert<W>(x, NearestEnv{}), units, e, log_m);
    return RoundTranscendental(y, W::kPrecision - LostBits(units), last, env,
                               result);
  });
}

template <class F, class Env, int kFirstLevel = 1>
F Log2Impl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  F special;
  if (LogSpecialCase(x, env, special)) return special;
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    const MathConstants<W>& c = GetMathConstants<W>();
    int units;
    std::int64_t e;
    W log_m;
    LogKernel<W, std::is_same_v<W, MathWorking<F, 1>>>(
        FloatAccess::Convert<W>(x, NearestEnv{}), units, e, log_m);
    if (log_m.is_zero()) {
      // A power of two: exact, unless the format cannot hold e.
      result = FloatAccess::Convert<F>(e, env);
      return true;
    }
    // |log2(m)| <= 1/2, so e + log2(m) adds at most 8 u.
    W y = log_m * c.inv_ln2;
    if (e != 0) y = W(e) + y;
    return RoundTranscendental(y, W::kPrecision - LostBits(units + 8), last,
                               env, result);
  });
}

// Whether x is 10^t for some t >= 1, which it returns in t.
template <class F>
bool IsPowerOfTen(const F& x, std::int64_t& t) {
  if (x.raw_exponent() < 3) return false;
  t = LowestBitExponent(x);  // 10^t = 5^t 2^t.
  if (t < 1) return false;
  // 10^t lies in [2^floor(t log2 10), 2^(floor(t log2 10) + 1)).
  const double estimate = static_cast<double>(t) * 3.32192809488736234787;
  if (static_cast<double>(x.raw_exponent()) < estimate - 2 ||
      static_cast<double>(x.raw_exponent()) > estimate + 2) {
    return false;
  }
  return EqualsDecimal(x, BigUint(1), t);
}

template <class F, class Env, int kFirstLevel = 1>
F Log10Impl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  F special;
  if (LogSpecialCase(x, env, special)) return special;
  std::int64_t t;
  if (IsPowerOfTen(x, t)) return FloatAccess::Convert<F>(t, env);
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    const MathConstants<W>& c = GetMathConstants<W>();
    int units;
    std::int64_t e;
    W log_m;
    const W y = LogKernel<W, std::is_same_v<W, MathWorking<F, 1>>>(
                    FloatAccess::Convert<W>(x, NearestEnv{}), units, e,
                    log_m) *
                c.inv_ln10;
    return RoundTranscendental(y, W::kPrecision - LostBits(units + 4), last,
                               env, result);
  });
}

// The special cases shared by log1p, log2p1 and log10p1. Returns true if they
// decide the result.
template <class F, class Env>
bool Log1pSpecialCase(const F& x, const Env& env, F& result) {
  if (x.is_nan()) {
    result = FloatAccess::PropagateNaN(x, env);
  } else if (x.is_zero() || (x.is_inf() && !x.signbit())) {
    result = x;
  } else if (x == F(-1)) {
    env.raise(kFlagDivByZero);
    result = F::Infinity(true);
  } else if (x < F(-1)) {
    result = InvalidResult<F>(env);
  } else {
    return false;
  }
  return true;
}

// log(1 + x) in base e, 2 or 10, for x outside the special cases.
template <int kBase, class F, int kFirstLevel, class Env>
F Log1pEvaluate(const F& x, const Env& env) {
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    const W wx = FloatAccess::Convert<W>(x, NearestEnv{});
    int units;
    W y;
    // The first format reduces 1 + x with tables, which needs 1 + x exact
    // and short; the others only take that path from 1/4 up.
    constexpr bool kFirst = std::is_same_v<W, MathWorking<F, 1>>;
    constexpr int kKernelBelow = kFirst ? -13 : -2;
    const std::int32_t ex = x.raw_exponent();
    std::int64_t e;
    W log_m;
    if (ex < kKernelBelow) {
      y = Log1pKernel(wx, units);
    } else if (ex < F::kPrecision) {
      // 1 + x is exact, with at most P + 14 bits.
      y = LogKernel<W, kFirst>(W(1) + wx, units, e, log_m);
    } else {
      // 1 + x is exact unless x is large, where log(1 + x) >= 0.22 and the
      // rounding of 1 + x adds below 5 u.
      y = LogKernel(W(1) + wx, units, e, log_m);
      units += 5;
    }
    if constexpr (kBase != 0) {
      // The inverse logarithm is within 2 u, the product 1 u more.
      const MathConstants<W>& c = GetMathConstants<W>();
      y = y * (kBase == 2 ? c.inv_ln2 : c.inv_ln10);
      units += 3;
    }
    return RoundTranscendental(y, W::kPrecision - LostBits(units), last, env,
                               result);
  });
}

template <class F, class Env, int kFirstLevel = 1>
F Log1pImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  F special;
  if (Log1pSpecialCase(x, env, special)) return special;
  // log(1 + x) = x - x^2/2 + ..., always below x.
  if (x.raw_exponent() <= -(F::kPrecision + 3)) {
    return RoundBeside(x, false, env);
  }
  return Log1pEvaluate<0, F, kFirstLevel>(x, env);
}

// --- Hyperbolic functions ---
//
// sinh, cosh and tanh take e^a for a = |x| from the exponential's kernel, and
// none has an exact result besides sinh(0), cosh(0) and tanh(0), since e^a is
// transcendental for every rational a != 0. Their inverses are logarithms.

// Whether x^2 < 2^-(P + 2), so that sinh, cosh and tanh differ from x, 1 and x
// by less than their value times 2^-(P + 2).
template <class F>
bool HyperbolicIsLinear(const F& x) {
  return 2 * std::int64_t{x.raw_exponent()} <= -(F::kPrecision + 4);
}

// Whether a = |x| makes sinh and cosh overflow: e^a / 2 is far above 2^(kEmax
// + 1) there. Below, they are evaluated, and RoundTranscendental decides.
template <class F>
bool HyperbolicOverflows(const F& a) {
  return FloatAccess::Convert<Wide1<F>>(a, NearestEnv{}) >=
         Wide1<F>(F::kEmax + 2);
}

// sinh a, cosh a and tanh a for a finite a > 0, below the overflow threshold
// and, for tanh, below P/2 + 2, with a relative error below `units` * u. The
// first working format reduces through tables.
template <class W, bool kUseTables>
W SinhKernel(const W& a, int& units) {
  std::int64_t k;
  const W r = ReduceByLn2<W>(a, k);  // Exact if k = 0.
  const W q = ExpM1Reduced<W, kUseTables>(r, units);
  const W one(1);
  if (k == 0) {
    // sinh a = q (q + 2) / (2 (q + 1)) for q = e^a - 1 in [0, 0.42], which
    // cancels nothing.
    units = 2 * units + 8;
    return scalbn(q * (q + W(2)) / (q + one), -1);
  }
  // a >= 0.34, so e^(2a) >= 2, and e^a - e^-a loses at most a factor of
  // (e^(2a) + 1) / (e^(2a) - 1) <= 3.
  const W e = scalbn(q + one, static_cast<int>(k));
  units = 4 * units + 16;
  return scalbn(e - one / e, -1);
}

template <class W, bool kUseTables>
W CoshKernel(const W& a, int& units) {
  std::int64_t k;
  const W r = ReduceByLn2<W>(a, k);
  const W q = ExpM1Reduced<W, kUseTables>(r, units);
  // (e^a + e^-a) / 2, a sum of positive terms.
  const W one(1);
  const W e = scalbn(q + one, static_cast<int>(k));
  units += 8;
  return scalbn(e + one / e, -1);
}

template <class W, bool kUseTables>
W TanhKernel(const W& a, int& units) {
  // tanh a = q / (q + 2) for q = e^(2a) - 1, computed as for expm1.
  std::int64_t k;
  const W r = ReduceByLn2<W>(scalbn(a, 1), k);
  W q = ExpM1Reduced<W, kUseTables>(r, units);
  if (k != 0) {
    // Within 4 u more for the reduction, as for expm1.
    q = scalbn(q + W(1), static_cast<int>(k)) - W(1);
    units = 4 * units + 12;
  }
  // Dividing by q + 2 at most doubles the error, and adds 1 u.
  units = 2 * units + 8;
  return q / (q + W(2));
}

// log z for z = 1 + d > 1 known to within 4 u relative, with a relative
// error below `units` * u. `small_d()` gives d to within 4 u relative; it is
// only called for d below about 2^-8, where log z is evaluated as log1p(d).
template <class W, bool kUseTables, class SmallD>
W LogOfInexact(const W& z, const SmallD& small_d, std::int64_t& units) {
  constexpr std::int64_t kInputUnits = 4;
  const W one(1);
  int kernel_units;
  if ((z - one).raw_exponent() < -8) {
    // A relative error of d carries over to log1p(d) at most as large.
    const W y = Log1pKernel(small_d(), kernel_units);
    units = kernel_units + kInputUnits + 1;
    return y;
  }
  // z has more bits than the tables' reduction keeps exact, so its two
  // products each round, by below 1 u relative: an absolute error below 2 u
  // in log z, besides the kernel's.
  std::int64_t e;
  W log_m;
  const W y = LogKernel<W, kUseTables>(z, kernel_units, e, log_m);
  constexpr std::int64_t kReductionUnits = kUseTables ? 2 : 0;
  // A relative error below delta in z is an absolute one below delta in
  // log z, and log z >= 2^raw_exponent, above 2^-10 here.
  const int shift = (std::max)(0, -y.raw_exponent());
  units = kernel_units + ((kInputUnits + kReductionUnits + 1) << shift);
  return y;
}

// asinh a, acosh a and atanh a, for a finite a > 0, a > 1 for acosh and
// a < 1 for atanh, with a relative error below `units` * u. The first working
// format reduces through tables.
template <class W, bool kUseTables>
W AsinhKernel(const W& a, std::int64_t& units) {
  const W one(1);
  const W a2 = a * a;
  if (a < one) {
    // asinh a = log1p(d) for d = a + a^2 / (1 + sqrt(1 + a^2)), within 2 u,
    // and 1 + d within 3 u.
    const W d = a + a2 / (one + sqrt(one + a2));
    return LogOfInexact<W, kUseTables>(one + d, [&] { return d; }, units);
  }
  // log(a + sqrt(a^2 + 1)), a sum of positive terms within 3 u.
  const W z = a + sqrt(a2 + one);
  return LogOfInexact<W, kUseTables>(z, [&] { return z - one; }, units);
}

template <class W, bool kUseTables>
W AcoshKernel(const W& a, std::int64_t& units) {
  const W one(1);
  if (a < W(2)) {
    // acosh a = log1p(d) for d = t + sqrt(t (t + 2)), t = a - 1 exact, which
    // cancels nothing: d is within 2 u, and 1 + d within 3 u.
    const W t = a - one;
    const W d = t + sqrt(t * (t + W(2)));
    return LogOfInexact<W, kUseTables>(one + d, [&] { return d; }, units);
  }
  // log(a + sqrt(a^2 - 1)): the subtraction loses at most a factor of 4/3,
  // and the sum is within 3 u.
  const W z = a + sqrt(a * a - one);
  return LogOfInexact<W, kUseTables>(z, [&] { return z - one; }, units);
}

template <class W, bool kUseTables>
W AtanhKernel(const W& a, std::int64_t& units) {
  // atanh a = log((1 + a) / (1 - a)) / 2 = log1p(2a / (1 - a)) / 2, with the
  // quotients within 2 u.
  const W one(1);
  const W below = one - a;
  const W z = (one + a) / below;
  return scalbn(LogOfInexact<W, kUseTables>(
                    z, [&] { return scalbn(a, 1) / below; }, units),
                -1);
}

enum class Hyperbolic { kSinh, kCosh, kTanh, kAsinh, kAcosh, kAtanh };

// Evaluates a hyperbolic function of |x| and rounds it, with the sign of x
// for the odd functions.
template <Hyperbolic kKind, class F, int kFirstLevel, class Env>
F HyperbolicEvaluate(const F& x, const Env& env) {
  const bool negative = kKind != Hyperbolic::kCosh &&
                        kKind != Hyperbolic::kAcosh && x.signbit();
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    constexpr bool kFirst = std::is_same_v<W, Wide1<F>>;
    const W a = abs(FloatAccess::Convert<W>(x, NearestEnv{}));
    std::int64_t units;
    W y;
    if constexpr (kKind == Hyperbolic::kAsinh) {
      y = AsinhKernel<W, kFirst>(a, units);
    } else if constexpr (kKind == Hyperbolic::kAcosh) {
      y = AcoshKernel<W, kFirst>(a, units);
    } else if constexpr (kKind == Hyperbolic::kAtanh) {
      y = AtanhKernel<W, kFirst>(a, units);
    } else {
      int small_units;
      if constexpr (kKind == Hyperbolic::kSinh) {
        y = SinhKernel<W, kFirst>(a, small_units);
      } else if constexpr (kKind == Hyperbolic::kCosh) {
        y = CoshKernel<W, kFirst>(a, small_units);
      } else {
        y = TanhKernel<W, kFirst>(a, small_units);
      }
      units = small_units;
    }
    return RoundTranscendental(negative ? -y : y,
                               W::kPrecision - LostBitsWide(units), last, env,
                               result);
  });
}

template <class F, class Env, int kFirstLevel = 1>
F SinhImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_inf() || x.is_zero()) return x;
  const bool negative = x.signbit();
  // sinh x = x + x^3/6 + ..., always beyond x.
  if (HyperbolicIsLinear(x)) return RoundBeside(x, !negative, env);
  if (HyperbolicOverflows(abs(x))) return RoundOverflow<F>(negative, env);
  return HyperbolicEvaluate<Hyperbolic::kSinh, F, kFirstLevel>(x, env);
}

template <class F, class Env, int kFirstLevel = 1>
F CoshImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_inf()) return F::Infinity();
  if (x.is_zero()) return F(1);
  // cosh x = 1 + x^2/2 + ..., always above 1.
  if (HyperbolicIsLinear(x)) return RoundBeside(F(1), true, env);
  if (HyperbolicOverflows(abs(x))) return RoundOverflow<F>(false, env);
  return HyperbolicEvaluate<Hyperbolic::kCosh, F, kFirstLevel>(x, env);
}

template <class F, class Env, int kFirstLevel = 1>
F TanhImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_zero()) return x;
  const bool negative = x.signbit();
  if (x.is_inf()) return negative ? F(-1) : F(1);
  // tanh x = x - x^3/3 + ..., always closer to zero than x.
  if (HyperbolicIsLinear(x)) return RoundBeside(x, negative, env);
  // 1 - tanh a = 2 / (e^(2a) + 1) < 2 e^-(P + 4) < 2^-(P + 2) from a >= P/2 + 2.
  if (FloatAccess::Convert<Wide1<F>>(abs(x), NearestEnv{}) >=
      Wide1<F>(F::kPrecision / 2 + 2)) {
    return RoundBeside(negative ? F(-1) : F(1), negative, env);
  }
  return HyperbolicEvaluate<Hyperbolic::kTanh, F, kFirstLevel>(x, env);
}

// asinh, acosh and atanh have no exact results besides asinh(0), acosh(1)
// and atanh(0): their values at algebraic arguments are transcendental.

template <class F, class Env, int kFirstLevel = 1>
F AsinhImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_inf() || x.is_zero()) return x;
  // asinh x = x - x^3/6 + ..., always closer to zero than x.
  if (HyperbolicIsLinear(x)) return RoundBeside(x, x.signbit(), env);
  return HyperbolicEvaluate<Hyperbolic::kAsinh, F, kFirstLevel>(x, env);
}

template <class F, class Env, int kFirstLevel = 1>
F AcoshImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  const F one(1);
  if (x < one) return InvalidResult<F>(env);
  if (x == one) return F::Zero();
  if (x.is_inf()) return x;
  return HyperbolicEvaluate<Hyperbolic::kAcosh, F, kFirstLevel>(x, env);
}

template <class F, class Env, int kFirstLevel = 1>
F AtanhImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_zero()) return x;
  const F a = abs(x);
  const F one(1);
  if (a == one) {
    env.raise(kFlagDivByZero);
    return F::Infinity(x.signbit());
  }
  if (a > one) return InvalidResult<F>(env);
  // atanh x = x + x^3/3 + ..., always beyond x.
  if (HyperbolicIsLinear(x)) return RoundBeside(x, !x.signbit(), env);
  return HyperbolicEvaluate<Hyperbolic::kAtanh, F, kFirstLevel>(x, env);
}

// --- Trigonometric functions ---

// pi/2 and 2/pi in a working format, rounded to nearest, and pi/2 split as
// half_pi_hi + half_pi_lo: half_pi_hi has kPrecision - 32 bits, so that
// k half_pi_hi is exact for |k| < 2^32, and the sum is within
// 2^-(2 kPrecision - 34) of pi/2.
template <class W>
struct TrigConstants {
  W half_pi;
  W half_pi_hi;
  W half_pi_lo;
  W two_over_pi;

  TrigConstants() {
    constexpr std::size_t kBits = 2 * W::kPrecision + 64;
    WithPiBits(kBits, [this](const PiFixedPoint& p) {
      // pi/2 2^kBits, and 2/pi 2^kBits, within 1.
      BigUint pi_fixed = p.pi;
      pi_fixed.ShiftRight(p.bits - kBits + 1);
      BigUint inverse = p.two_over_pi;
      inverse.ShiftRight(p.bits - kBits);
      half_pi = FromFixedPoint<W>(pi_fixed, kBits);
      two_over_pi = FromFixedPoint<W>(inverse, kBits);
      BigUint hi = pi_fixed;
      const std::size_t drop = hi.bit_length() - (W::kPrecision - 32);
      hi.ShiftRight(drop);
      hi.ShiftLeft(drop);
      BigUint lo = pi_fixed;
      lo.Sub(hi);
      half_pi_hi = FromFixedPoint<W>(hi, kBits);
      half_pi_lo = FromFixedPoint<W>(lo, kBits);
      return 0;
    });
  }
};

template <class W>
const TrigConstants<W>& GetTrigConstants() {
  static const TrigConstants<W> constants;
  return constants;
}

// |x| = (k + f) pi/2 for an integer k and |f| <= 1/2, for a finite nonzero x:
// returns r = f pi/2, with a relative error below 3.1 u, and k modulo 4 in
// `quadrant`. This is Payne and Hanek's reduction: with |x| = M 2^q for an
// integer M, only the bits of 2/pi from 2^-(q - 1) on matter modulo 4, and
// the product of M with a window of them is exact.
// The first attempt keeps kFirstGuard bits beyond W's precision, and each
// further one 64 more; tests start lower to exercise the later attempts.
template <class W, std::int64_t kFirstGuard = 64, class F>
W ReduceByHalfPi(const F& x, int& quadrant) {
  constexpr std::int64_t kSignificandBits =
      static_cast<std::int64_t>(F::kLimbs) * kLimbBits;
  const BigUint m = BigUint::FromLimbs(x.raw_limbs().data(), F::kLimbs);
  const std::int64_t q = std::int64_t{x.raw_exponent()} - kSignificandBits + 1;
  // 2/pi = sum of t_j 2^-j over j >= 1; M 2^q t_j 2^-j is a multiple of 4
  // for j <= q - 2.
  const std::int64_t first = (std::max)(std::int64_t{1}, q - 1);
  for (std::int64_t guard = kFirstGuard;; guard += 64) {
    // The bits of 2/pi after 2^-last, and the table's error, below
    // 2^-(last + 8), move x 2/pi by less than 1.01 M 2^(q - last), which
    // is 1.01 2^-(W::kPrecision + guard).
    const std::int64_t last = q + kSignificandBits + W::kPrecision + guard;
    const auto table_bits = static_cast<std::size_t>(last + 8);
    const BigUint product =
        WithPiBits(table_bits, [&](const PiFixedPoint& p) {
          return m * BitRange(p.two_over_pi,
                                 p.bits - static_cast<std::size_t>(last),
                                 static_cast<std::size_t>(last - first + 1));
        });
    // product = x 2/pi 2^fraction_bits, modulo 2^(fraction_bits + 2).
    const auto fraction_bits = static_cast<std::size_t>(last - q);
    const bool above_half = product.bit(fraction_bits - 1);
    BigUint fraction = BitRange(product, 0, fraction_bits);
    if (above_half) {
      BigUint whole(1);
      whole.ShiftLeft(fraction_bits);
      whole.Sub(fraction);
      fraction = std::move(whole);
    }
    // With |f| >= 2^-(guard - 8), the error is below 2^-(W::kPrecision + 7)
    // relative. Otherwise x lies very close to a multiple of pi/2, and more
    // bits are needed.
    if (static_cast<std::int64_t>(fraction.bit_length()) + guard <
        static_cast<std::int64_t>(fraction_bits) + 9) {
      continue;
    }
    quadrant = ((product.bit(fraction_bits) ? 1 : 0) +
                (product.bit(fraction_bits + 1) ? 2 : 0) +
                (above_half ? 1 : 0)) &
               3;
    // Rounding f, pi/2 and their product add 3 u.
    const W f = RoundInteger<W>(above_half, fraction,
                                -static_cast<std::int64_t>(fraction_bits),
                                false, SilentEnv{Rounding::kNearestEven});
    return f * GetTrigConstants<W>().half_pi;
  }
}

// The reduction of ReduceByHalfPi for |x| < 2^27, by Cody and Waite's method:
// r = (|x| - k half_pi_hi) - k half_pi_lo, where the first difference is
// exact. Then r is within 2^-(2 kPrecision - 61) of the true value, from
// k half_pi_lo and the error of the split, plus one rounding: within 1.01 u
// for |r| >= 2^-(kPrecision - 68). Returns false for a smaller r, which needs
// the exact reduction.
template <class W, class F>
bool ReduceSmallByHalfPi(const F& x, int& quadrant, W& r) {
  const TrigConstants<W>& c = GetTrigConstants<W>();
  const W wide = FloatAccess::Convert<W>(x, NearestEnv{});
  const W multiple = round(wide * c.two_over_pi);
  const auto k = to_integer<std::int64_t>(multiple, Rounding::kNearestEven);
  quadrant = static_cast<int>(k & 3);
  if (k == 0) {
    r = wide;
    return true;
  }
  // x and k half_pi_hi are multiples of 2^-(kPrecision - 34) below 2^29,
  // and their difference is below 1.
  r = (wide - multiple * c.half_pi_hi) - multiple * c.half_pi_lo;
  return r.raw_exponent() >= 68 - W::kPrecision;
}

// sin b for an exact b with 0 < |b| < 1/2, with a relative error below
// 6 u: b G(b^2) for G(y) = sum over j of (-y)^j / (2j + 1)!, in fixed point.
template <class W>
W SinSeries(const W& b) {
  constexpr std::size_t kN = W::kLimbs;
  // The first omitted term, j = n, relative to b, is below
  // |b|^2n / (2n + 1)! < 2^-(kPrecision + 4).
  const int scale = -(b.raw_exponent() + 1);
  int n = 1;
  int log2_factorial = 2;  // A lower bound on log2((2n + 1)!).
  while (2 * n * scale + log2_factorial < W::kPrecision + 4) {
    ++n;
    log2_factorial += std::bit_width(static_cast<unsigned>(2 * n)) - 1 +
                      std::bit_width(static_cast<unsigned>(2 * n + 1)) - 1;
  }
  const auto d = [](int j) {
    return static_cast<limb_t>(2 * j) * static_cast<limb_t>(2 * j + 1);
  };
  const int k = SeriesBlock(n, std::bit_width(d(n)));
  std::array<FixedLimbs<kN>, 17> powers;
  powers[1] = ToFixed<kN>(b * b);
  for (int i = 2; i <= k; ++i) {
    powers[static_cast<std::size_t>(i)] =
        FixedMul<kN>(powers[static_cast<std::size_t>(i - 1)], powers[1]);
  }
  // G is above 0.98 and within 4 u absolute; rounding it and b G add 2 u.
  return b * FromFixed<W, kN>(FixedSeries<kN>(powers.data(), n, k, d, true));
}

// sin r and cos r - 1 for an exact r with |r| < 0.8, with relative errors
// below `sin_units` and `cos_units` u, by halving r and doubling back. For
// working formats without tables, and for building them.
template <class W>
void SinCosByHalving(const W& r, W& sin_r, W& cos_r_minus_1, double& sin_units,
                     double& cos_units) {
  const W one(1);
  if (r.is_zero()) {
    sin_r = r;
    cos_r_minus_1 = W::Zero();
    sin_units = 0;
    cos_units = 0;
    return;
  }
  // Halve r below 2^-kScale; each halving is undone by two multiplications
  // afterwards.
  constexpr int kScale = (std::max)(2, IntegerSqrt(W::kPrecision) / 3);
  const int halvings = (std::max)(0, r.raw_exponent() + 1 + kScale);
  W s = SinSeries(scalbn(r, -halvings));
  // Each rounding adds up to 1 u.
  double ds = 6;
  // cos b - 1 = -s^2 / (1 + sqrt(1 - s^2)), with no cancellation: s^2 within
  // 2 ds + 1, the denominator within 1.8 u and the quotient 1 u more.
  W c = -(s * s) / (one + sqrt(one - s * s));
  double dc = 2 * ds + 4;
  // sin 2a = 2 sin a (1 + (cos a - 1)) and cos 2a - 1 = -2 sin^2 a. A relative
  // error e of cos a - 1 makes one of |c| e / (1 + c) <= 2 |c| e in 1 + c,
  // for c >= -1/2.
  for (int i = 0; i < halvings; ++i) {
    const double ratio = std::ldexp(1.0, c.raw_exponent() + 2);
    const W next_c = -scalbn(s * s, 1);
    s = scalbn(s * (one + c), 1);
    const double next_ds = ds + ratio * dc + 2;
    dc = 2 * ds + 1;
    ds = next_ds;
    c = next_c;
  }
  sin_r = s;
  cos_r_minus_1 = c;
  sin_units = ds;
  cos_units = dc;
}

// A value and a bound on its relative error, in units of u.
template <class W>
struct Tracked {
  W value;
  double units;
};

template <class W>
Tracked<W> TrackedMul(const Tracked<W>& a, const Tracked<W>& b) {
  return {a.value * b.value, a.units + b.units + 1};
}

// a + b, whose error is that of each term relative to the sum: with
// |a| < 2^(ea + 1) and |a + b| >= 2^es, 2^(ea + 1 - es) times a's.
template <class W>
Tracked<W> TrackedAdd(const Tracked<W>& a, const Tracked<W>& b) {
  if (a.value.is_zero()) return b;
  if (b.value.is_zero()) return a;
  const W sum = a.value + b.value;
  if (sum.is_zero()) return {sum, 0};  // Only for exactly opposite terms.
  const int es = sum.raw_exponent();
  // The computed sum is within 1 u of the exact one, hence the factor.
  return {sum, (std::ldexp(a.units, a.value.raw_exponent() + 1 - es) +
                std::ldexp(b.units, b.value.raw_exponent() + 1 - es)) *
                       1.001 +
                   1};
}

// sin and cos of i/2^7 for i in [0, 102] and of j/2^14 for j in [0, 64],
// computed in twice the precision and rounded once: within 1.1 u.
template <class W>
struct TrigTables {
  static constexpr int kCount1 = 103;
  static constexpr int kCount2 = 65;
  std::array<W, kCount1> sin1;
  std::array<W, kCount1> cos1;
  std::array<W, kCount2> sin2;
  std::array<W, kCount2> cos2;

  TrigTables() {
    using Wide = Float<2 * W::kPrecision, 29>;
    const auto fill = [](int count, int shift, W* sin_out, W* cos_out) {
      for (int i = 0; i < count; ++i) {
        Wide s;
        Wide c;
        double s_units;
        double c_units;
        SinCosByHalving(scalbn(Wide(i), -shift), s, c, s_units, c_units);
        sin_out[i] = FloatAccess::Convert<W>(s, NearestEnv{});
        cos_out[i] = FloatAccess::Convert<W>(Wide(1) + c, NearestEnv{});
      }
    };
    fill(kCount1, 7, sin1.data(), cos1.data());
    fill(kCount2, 14, sin2.data(), cos2.data());
  }
};

template <class W>
const TrigTables<W>& GetTrigTables() {
  static const TrigTables<W> tables;
  return tables;
}

// sin r and cos r for an exact r with |r| < 0.8, with relative errors below
// `sin_units` and `cos_units` u. With tables, r = i/2^7 + j/2^14 + r2 with
// |r2| <= 2^-15, and sin(a + b) = sin a cos b + cos a sin b and
// cos(a + b) = cos a cos b - sin a sin b combine the parts; the errors are
// tracked as they go, since the parts may have opposite signs.
template <class W, bool kUseTables>
void SinCosReduced(const W& r, W& sin_r, W& cos_r, double& sin_units,
                   double& cos_units) {
  const W one(1);
  if constexpr (kUseTables) {
    const TrigTables<W>& t = GetTrigTables<W>();
    const W a = abs(r);
    // Both subtractions are exact.
    const int i = to_integer<int>(scalbn(a, 7), Rounding::kNearestEven);
    const W r1 = i != 0 ? a - scalbn(W(i), -7) : a;
    const int j = to_integer<int>(scalbn(r1, 14), Rounding::kNearestEven);
    const W r2 = j != 0 ? r1 - scalbn(W(j), -14) : r1;
    Tracked<W> s{W::Zero(), 0};
    Tracked<W> c{one, 0};
    if (!r2.is_zero()) {
      s = {SinSeries(r2), 6};
      // cos r2 = sqrt(1 - s^2), where s^2 < 2^-30 adds nothing to speak of.
      c = {sqrt(one - s.value * s.value), 1.6};
    }
    const auto combine = [&s, &c](const W& sin_part, const W& cos_part,
                                  bool negative) {
      const Tracked<W> sp{negative ? -sin_part : sin_part, 1.1};
      const Tracked<W> cp{cos_part, 1.1};
      const Tracked<W> next_s =
          TrackedAdd(TrackedMul(sp, c), TrackedMul(cp, s));
      const Tracked<W> sin_sin = TrackedMul(sp, s);
      c = TrackedAdd(TrackedMul(cp, c),
                     Tracked<W>{-sin_sin.value, sin_sin.units});
      s = next_s;
    };
    if (j != 0) {
      const auto index = static_cast<std::size_t>(j < 0 ? -j : j);
      combine(t.sin2[index], t.cos2[index], j < 0);
    }
    if (i != 0) {
      const auto index = static_cast<std::size_t>(i);
      combine(t.sin1[index], t.cos1[index], false);
    }
    sin_r = r.signbit() ? -s.value : s.value;
    cos_r = c.value;
    sin_units = s.units;
    cos_units = c.units;
  } else {
    W c;
    double c_units;
    SinCosByHalving(r, sin_r, c, sin_units, c_units);
    // 1 + c, with |c| / (1 + c) <= 2 |c| for c >= -1/2.
    cos_r = one + c;
    cos_units =
        std::ldexp(c_units, c.is_zero() ? -64 : c.raw_exponent() + 2) + 1;
  }
}

enum class Trig { kSin, kCos, kTan };

// The units of an error bound, rounded up.
inline std::int64_t CeilUnits(double units) {
  return static_cast<std::int64_t>(units) + 1;
}

template <Trig kKind, class F, int kFirstLevel, class Env>
F TrigEvaluate(const F& x, const Env& env) {
  const bool negative = x.signbit();
  const F a = abs(x);
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    // Below 1/2, r = |x| is exact.
    int quadrant = 0;
    W r;
    double reduction_units = 0;
    if (a.raw_exponent() < -1) {
      r = FloatAccess::Convert<W>(a, NearestEnv{});
    } else if (a.raw_exponent() < 27 && ReduceSmallByHalfPi<W>(a, quadrant, r)) {
      reduction_units = 1.5;
    } else {
      r = ReduceByHalfPi<W>(a, quadrant);
      reduction_units = 3.1;
    }
    W s;
    W cos_r;
    double s_units;
    double c_units;
    SinCosReduced<W, std::is_same_v<W, Wide1<F>>>(r, s, cos_r, s_units,
                                                  c_units);
    // A relative error e of r makes at most e in sin r and 0.8 e in cos r
    // for |r| <= pi/4.
    s_units += reduction_units;
    c_units += reduction_units;
    W y;
    double units;
    if constexpr (kKind == Trig::kTan) {
      // tan x = tan r, or -1 / tan r in odd quadrants.
      units = s_units + c_units + 1;
      y = (quadrant & 1) == 0 ? s / cos_r : -(cos_r / s);
      if (negative) y = -y;
    } else {
      // sin(r + k pi/2) is sin r, cos r, -sin r and -cos r for k = 0..3, and
      // cos(r + k pi/2) is sin(r + (k + 1) pi/2).
      constexpr bool kOdd = kKind == Trig::kSin;
      const int k = kOdd ? quadrant : quadrant + 1;
      const bool use_cos = (k & 1) != 0;
      y = use_cos ? cos_r : s;
      units = use_cos ? c_units : s_units;
      if ((k & 2) != 0) y = -y;
      if constexpr (kOdd) {
        if (negative) y = -y;
      }
    }
    return RoundTranscendental(y,
                               W::kPrecision - LostBitsWide(CeilUnits(units)),
                               last, env, result);
  });
}

// sin, cos and tan have no exact results besides sin(0), cos(0) and tan(0):
// their values at algebraic arguments are transcendental.

template <class F, class Env, int kFirstLevel = 1>
F SinImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_inf()) return InvalidResult<F>(env);
  if (x.is_zero()) return x;
  // sin x = x - x^3/6 + ..., always closer to zero than x.
  if (HyperbolicIsLinear(x)) return RoundBeside(x, x.signbit(), env);
  return TrigEvaluate<Trig::kSin, F, kFirstLevel>(x, env);
}

template <class F, class Env, int kFirstLevel = 1>
F CosImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_inf()) return InvalidResult<F>(env);
  if (x.is_zero()) return F(1);
  // cos x = 1 - x^2/2 + ..., always below 1.
  if (HyperbolicIsLinear(x)) return RoundBeside(F(1), false, env);
  return TrigEvaluate<Trig::kCos, F, kFirstLevel>(x, env);
}

template <class F, class Env, int kFirstLevel = 1>
F TanImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_inf()) return InvalidResult<F>(env);
  if (x.is_zero()) return x;
  // tan x = x + x^3/3 + ..., always beyond x.
  if (HyperbolicIsLinear(x)) return RoundBeside(x, !x.signbit(), env);
  return TrigEvaluate<Trig::kTan, F, kFirstLevel>(x, env);
}

// +-pi m 2^shift rounded to F, for m = 1 or 3, from P + 64 bits whose last 62
// are far from a rounding boundary, which every value of pi's bits so far
// has been.
template <class F, class Env>
F PiMultiple(limb_t m, int shift, bool negative, const Env& env) {
  for (std::size_t bits = F::kPrecision + 64;; bits *= 2) {
    bool decided = false;
    const F result = WithPiBits(bits + 2, [&](const PiFixedPoint& p) {
      BigUint top = p.pi;
      top.ShiftRight(p.bits - bits);  // Within 2 of pi 2^bits.
      top.MulAdd(m, 0);               // Within 6 of pi m 2^bits.
      const limb_t low = top.limbs()[0] & ((limb_t{1} << 62) - 1);
      decided = low >= 8 && low <= (limb_t{1} << 62) - 8;
      return RoundInteger<F>(negative, top,
                             shift - static_cast<std::int64_t>(bits), true,
                             env);
    });
    if (decided) return result;
  }
}

template <class F, class Env>
F PiImpl(const Env& env) {
  return PiMultiple<F>(1, 0, false, env);
}

// --- Inverse trigonometric functions ---
//
// All four reduce to atan of a nonnegative argument: asin x =
// atan(x / sqrt((1 - x)(1 + x))), acos x = 2 atan(sqrt((1 - x) / (1 + x)))
// and atan2(y, x) = atan(|y / x|) in the right quadrant. None has an exact
// result besides at 0 and, for acos, at 1: their values at other algebraic
// arguments are transcendental.

// atan t for an exact t with |t| < 1/2, with a relative error below 10 u:
// t G(t^2) for G(y) = sum over j of (-y)^j / (2j + 1), in fixed point by
// Paterson and Stockmeyer. In blocks of k terms, acc_b = S / D + (-y)^k
// acc_(b+1) with D = (2bk + 1) ... (2bk + 2k - 1) and S the sum over i < k
// of (-1)^i (D / (2bk + 2i + 1)) y^i, whose terms decrease.
template <class W>
W AtanSeries(const W& t) {
  constexpr std::size_t kN = W::kLimbs;
  // The first omitted term, j = n, is below t^2n < 2^-(kPrecision + 4)
  // relative to t.
  const int scale = -(t.raw_exponent() + 1);
  const int n = (W::kPrecision + 4 + 2 * scale - 1) / (2 * scale);
  const int k = SeriesBlock(n, std::bit_width(static_cast<unsigned>(2 * n + 1)));
  std::array<FixedLimbs<kN>, 17> powers{};
  powers[0][kN] = 1;
  powers[1] = ToFixed<kN>(t * t);
  for (int i = 2; i <= k; ++i) {
    powers[static_cast<std::size_t>(i)] =
        FixedMul<kN>(powers[static_cast<std::size_t>(i - 1)], powers[1]);
  }
  const int blocks = (n + k - 1) / k;
  FixedLimbs<kN> acc{};
  for (int block = blocks - 1; block >= 0; --block) {
    const int base = block * k;
    const int size = (std::min)(k, n - base);
    limb_t denominator = 1;
    for (int i = 0; i < size; ++i) {
      denominator *= static_cast<limb_t>(2 * (base + i) + 1);
    }
    FixedLimbs<kN> sum{};
    for (int i = 0; i < size; ++i) {
      FixedAccumulate<kN>(
          sum,
          FixedMulSmall<kN>(powers[static_cast<std::size_t>(i)],
                            denominator / static_cast<limb_t>(2 * (base + i) + 1)),
          (i & 1) != 0);
    }
    DivideLimbs(sum.data(), kN + 1, denominator);
    if (block < blocks - 1) {
      FixedAccumulate<kN>(
          sum, FixedMul<kN>(powers[static_cast<std::size_t>(size)], acc),
          (size & 1) != 0);
    }
    acc = sum;
  }
  // G is above 0.9 and within 8 u absolute; rounding it and t G add 2 u.
  return t * FromFixed<W, kN>(acc);
}

// atan a for a >= 0, by halving: atan a = 2 atan(a / (1 + sqrt(1 + a^2))),
// each step within 3.5 u more. For working formats without tables, and for
// building them.
template <class W>
Tracked<W> AtanByHalving(W a) {
  if (a.is_zero()) return {a, 0};
  const W one(1);
  constexpr int kScale = (std::max)(2, IntegerSqrt(W::kPrecision) / 4);
  int halvings = 0;
  double units = 0;
  while (a.raw_exponent() >= -kScale) {
    a = a / (one + sqrt(one + a * a));
    ++halvings;
    units += 3.5;
  }
  return {scalbn(AtanSeries(a), halvings), units + 10};
}

// atan(i/2^7) for i in [0, 128] and atan(j/2^14) for j in [0, 64], computed in
// twice the precision and rounded once: within 1.1 u.
template <class W>
struct AtanTables {
  static constexpr int kCount1 = 129;
  static constexpr int kCount2 = 65;
  std::array<W, kCount1> atan1;
  std::array<W, kCount2> atan2;

  AtanTables() {
    using Wide = Float<2 * W::kPrecision, 29>;
    for (int i = 0; i < kCount1; ++i) {
      atan1[static_cast<std::size_t>(i)] = FloatAccess::Convert<W>(
          AtanByHalving(scalbn(Wide(i), -7)).value, NearestEnv{});
    }
    for (int j = 0; j < kCount2; ++j) {
      atan2[static_cast<std::size_t>(j)] = FloatAccess::Convert<W>(
          AtanByHalving(scalbn(Wide(j), -14)).value, NearestEnv{});
    }
  }
};

template <class W>
const AtanTables<W>& GetAtanTables() {
  static const AtanTables<W> tables;
  return tables;
}

// (t - c) / (1 + t c) for a table point c, where t - c is exact, and t has a
// relative error of `t_units`: the subtraction turns it into an absolute one.
template <class W>
Tracked<W> AtanStep(const Tracked<W>& t, const W& c) {
  const W one(1);
  const W numerator = t.value - c;
  if (numerator.is_zero()) return {numerator, 0};
  const double numerator_units =
      t.value.is_zero()
          ? 0
          : std::ldexp(t.units,
                       t.value.raw_exponent() + 1 - numerator.raw_exponent());
  // The denominator is within 2 u, the quotient 1 u more.
  return {numerator / (one + t.value * c), numerator_units + 3};
}

// atan a for an a >= 0 known within `a_units` u relative. With tables,
// a <= 1 is reduced as a = c1 + ... with atan a = atan c1 + atan c2 + atan t,
// for c1 = i/2^7, c2 = j/2^14 and t = (t1 - c2) / (1 + t1 c2), where
// t1 = (a - c1) / (1 + a c1) and |t| <= 2^-15 or so; a > 1 goes through
// pi/2 - atan(1/a).
template <class W, bool kUseTables>
Tracked<W> AtanNonnegative(const W& a, double a_units) {
  const W one(1);
  if (a.is_zero()) return {a, 0};
  if (a > one) {
    const Tracked<W> inner = AtanNonnegative<W, kUseTables>(one / a, a_units + 1);
    return TrackedAdd(Tracked<W>{GetTrigConstants<W>().half_pi, 1},
                      Tracked<W>{-inner.value, inner.units});
  }
  if constexpr (kUseTables) {
    const AtanTables<W>& tables = GetAtanTables<W>();
    // a - c1 is exact, and so is t1 - c2, both being at most half a step.
    const int i = to_integer<int>(scalbn(a, 7), Rounding::kNearestEven);
    Tracked<W> t{a, 0};
    if (i != 0) t = AtanStep(t, scalbn(W(i), -7));
    const int j = to_integer<int>(scalbn(t.value, 14), Rounding::kNearestEven);
    if (j != 0) t = AtanStep(t, scalbn(W(j), -14));
    Tracked<W> result{W::Zero(), 0};
    if (!t.value.is_zero()) result = {AtanSeries(t.value), t.units + 10};
    if (j != 0) {
      const W part = tables.atan2[static_cast<std::size_t>(j < 0 ? -j : j)];
      result = TrackedAdd(Tracked<W>{j < 0 ? -part : part, 1.1}, result);
    }
    if (i != 0) {
      result = TrackedAdd(
          Tracked<W>{tables.atan1[static_cast<std::size_t>(i)], 1.1}, result);
    }
    // A relative error e of a makes at most e in atan a.
    result.units += a_units;
    return result;
  } else {
    Tracked<W> result = AtanByHalving(a);
    result.units += a_units;
    return result;
  }
}

// Rounds a value just below the exact positive v, of at most P + 1 bits, as
// RoundBeside does for P bits.
template <class F, class W, class Env>
F RoundJustBelow(bool negative, const W& v, const Env& env) {
  const W magnitude = v - scalbn(v, -(F::kPrecision + 4));  // Exact.
  return FloatAccess::Round<F>(negative, magnitude.raw_exponent(),
                               magnitude.raw_limbs().data(), W::kLimbs, true,
                               env);
}

enum class InverseTrig { kAtan, kAsin, kAcos, kAtan2 };

// Evaluates atan, asin, acos or atan2(y, x) for finite nonzero arguments
// outside their special cases.
template <InverseTrig kKind, class F, int kFirstLevel, class Env>
F InverseTrigEvaluate(const F& x, const F& y, const Env& env) {
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    constexpr bool kFirst = std::is_same_v<W, Wide1<F>>;
    const W one(1);
    const W wx = FloatAccess::Convert<W>(x, NearestEnv{});
    Tracked<W> value{W::Zero(), 0};
    bool negative = false;
    if constexpr (kKind == InverseTrig::kAtan) {
      negative = x.signbit();
      value = AtanNonnegative<W, kFirst>(abs(wx), 0);
    } else if constexpr (kKind == InverseTrig::kAsin) {
      // x / sqrt((1 - |x|)(1 + |x|)): 1 - |x| is exact from 1/2 on, so the
      // quotient is within 4 u.
      negative = x.signbit();
      const W a = abs(wx);
      value = AtanNonnegative<W, kFirst>(a / sqrt((one - a) * (one + a)), 4);
    } else if constexpr (kKind == InverseTrig::kAcos) {
      // 2 atan(sqrt((1 - x) / (1 + x))), within 3 u before the atan.
      value = AtanNonnegative<W, kFirst>(sqrt((one - wx) / (one + wx)), 3);
      value.value = scalbn(value.value, 1);
    } else {
      // atan2(y, x) = +-atan(|y / x|), or +-(pi - atan(|y / x|)) for x < 0.
      negative = y.signbit();
      const W wy = FloatAccess::Convert<W>(y, NearestEnv{});
      value = AtanNonnegative<W, kFirst>(abs(wy) / abs(wx), 1);
      if (x.signbit()) {
        value = TrackedAdd(
            Tracked<W>{scalbn(GetTrigConstants<W>().half_pi, 1), 1},
            Tracked<W>{-value.value, value.units});
      }
    }
    const W signed_value = negative ? -value.value : value.value;
    return RoundTranscendental(
        signed_value, W::kPrecision - LostBitsWide(CeilUnits(value.units)),
        last, env, result);
  });
}

template <class F, class Env, int kFirstLevel = 1>
F AtanImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_zero()) return x;
  if (x.is_inf()) return PiMultiple<F>(1, -1, x.signbit(), env);
  // atan x = x - x^3/3 + ..., always closer to zero than x.
  if (HyperbolicIsLinear(x)) return RoundBeside(x, x.signbit(), env);
  return InverseTrigEvaluate<InverseTrig::kAtan, F, kFirstLevel>(x, x, env);
}

template <class F, class Env, int kFirstLevel = 1>
F AsinImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_zero()) return x;
  const F a = abs(x);
  if (a > F(1)) return InvalidResult<F>(env);
  if (a == F(1)) return PiMultiple<F>(1, -1, x.signbit(), env);
  // asin x = x + x^3/6 + ..., always beyond x.
  if (HyperbolicIsLinear(x)) return RoundBeside(x, !x.signbit(), env);
  return InverseTrigEvaluate<InverseTrig::kAsin, F, kFirstLevel>(x, x, env);
}

template <class F, class Env, int kFirstLevel = 1>
F AcosImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (abs(x) > F(1)) return InvalidResult<F>(env);
  if (x == F(1)) return F::Zero();
  if (x == F(-1)) return PiMultiple<F>(1, 0, false, env);
  return InverseTrigEvaluate<InverseTrig::kAcos, F, kFirstLevel>(x, x, env);
}

template <class F, class Env, int kFirstLevel = 1>
F Atan2Impl(const F& y, const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  // The special cases of IEEE 754 and C99 Annex F.
  if (y.is_nan() || x.is_nan()) return FloatAccess::PropagateNaN(y, x, env);
  const bool negative = y.signbit();
  if (y.is_zero()) {
    return x.signbit() ? PiMultiple<F>(1, 0, negative, env) : y;
  }
  if (x.is_zero()) return PiMultiple<F>(1, -1, negative, env);
  if (y.is_inf()) {
    if (!x.is_inf()) return PiMultiple<F>(1, -1, negative, env);
    return x.signbit() ? PiMultiple<F>(3, -2, negative, env)
                       : PiMultiple<F>(1, -2, negative, env);
  }
  if (x.is_inf()) {
    return x.signbit() ? PiMultiple<F>(1, 0, negative, env) : F::Zero(negative);
  }
  if (!x.signbit()) {
    // For |y / x| < 2^-(1.5 P + 8), atan2 is just below the exact quotient q,
    // by less than q 2^-(3P + 14). Truncated to P + 2 bits, q either is exact,
    // or lies more than that above its truncation: the quotient of two P-bit
    // numbers is farther than q 2^-(2P + 4) from any other P + 2-bit number.
    using Q = Float<F::kPrecision + 2, 29>;
    const std::int64_t estimate =
        std::int64_t{y.raw_exponent()} - x.raw_exponent();
    if (estimate < -(3 * F::kPrecision / 2 + 10)) {
      Context context{Rounding::kTowardZero};
      const Q q = div(abs(FloatAccess::Convert<Q>(y, NearestEnv{})),
                      abs(FloatAccess::Convert<Q>(x, NearestEnv{})), context);
      if ((context.flags & kFlagInexact) != 0) {
        return FloatAccess::Round<F>(negative, q.raw_exponent(),
                                     q.raw_limbs().data(), Q::kLimbs, true,
                                     env);
      }
      return RoundJustBelow<F>(
          negative, FloatAccess::Convert<Wide1<F>>(q, NearestEnv{}), env);
    }
  }
  return InverseTrigEvaluate<InverseTrig::kAtan2, F, kFirstLevel>(x, y, env);
}



// --- exp10, exp2m1, exp10m1, log2p1 and log10p1 ---

// x log2(10) = k + f with |f| <= 1/2 or so: returns r = f ln 2, within 1.1 u
// absolute. x log2_10_hi is exact, and so is its difference with k, both
// being multiples of 2^-(kPrecision - 2) or more below 2^28.
template <class W, class F>
W ReduceByLog2Of10(const F& x, std::int64_t& k) {
  const MathConstants<W>& c = GetMathConstants<W>();
  const W wx = FloatAccess::Convert<W>(x, NearestEnv{});
  const W high = wx * c.log2_10_hi;
  const W low = wx * c.log2_10_lo;
  const W multiple = round(high + low);
  k = to_integer<std::int64_t>(multiple, Rounding::kNearestEven);
  const W f = multiple.is_zero() ? high + low : (high - multiple) + low;
  return f * c.ln2;
}

// x log2(10), to check the thresholds of overflow and underflow.
template <class F>
Wide1<F> TimesLog2Of10(const F& x) {
  const MathConstants<Wide1<F>>& c = GetMathConstants<Wide1<F>>();
  return FloatAccess::Convert<Wide1<F>>(x, NearestEnv{}) *
         (c.log2_10_hi + c.log2_10_lo);
}

// Rounds a value just above an exact v of the first working format.
template <class F, class W, class Env>
F RoundJustAbove(const W& v, const Env& env) {
  return FloatAccess::Round<F>(v.signbit(), v.raw_exponent(),
                               v.raw_limbs().data(), W::kLimbs, true, env);
}

template <class F, class Env, int kFirstLevel = 1>
F Exp10Impl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_inf()) return x.signbit() ? F::Zero() : x;
  if (x.is_zero()) return F(1);
  // 10^x = 1 + x ln 10 + ..., with |x ln 10| < 2^-(P + 2).
  if (x.raw_exponent() <= -(F::kPrecision + 5)) {
    return RoundBeside(F(1), !x.signbit(), env);
  }
  const Wide1<F> t = TimesLog2Of10(x);
  if (t >= Wide1<F>(F::kEmax + 2)) return RoundOverflow<F>(false, env);
  if (t <= Wide1<F>(-(F::kEmax + F::kPrecision + 4))) {
    return RoundTiny<F>(false, env);
  }
  // 10^n = 5^n 2^n for an integer n >= 0 is exact or a tie if 5^n has at
  // most P + 1 bits, and otherwise as far from both as 10^x at other x.
  if (!x.signbit() && IsInteger(x) && x <= F(F::kPrecision / 2 + 1)) {
    const auto n = to_integer<std::int64_t>(x, Rounding::kTowardZero);
    const BigUint power = BigUint::Pow5(static_cast<std::uint64_t>(n));
    if (power.bit_length() <= static_cast<std::size_t>(F::kPrecision) + 1) {
      return RoundInteger<F>(false, power, n, false, env);
    }
  }
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    std::int64_t k;
    const W r = ReduceByLog2Of10<W>(x, k);
    int units;
    const W q = ExpM1Reduced<W, std::is_same_v<W, Wide1<F>>>(r, units);
    const W y = scalbn(q + W(1), static_cast<int>(k));
    return RoundTranscendental(y, W::kPrecision - LostBits(units + 4), last,
                               env, result);
  });
}

template <class F, class Env, int kFirstLevel = 1>
F Exp2m1Impl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_inf()) return x.signbit() ? F(-1) : x;
  if (x.is_zero()) return x;
  const Wide1<F> wide = FloatAccess::Convert<Wide1<F>>(x, NearestEnv{});
  // 2^x - 1 >= 2^(kEmax + 1) for x > kEmax + 1, while at kEmax + 1 itself it
  // is just below, which rounding decides.
  if (wide > Wide1<F>(F::kEmax + 1)) return RoundOverflow<F>(false, env);
  // Then 2^x < 2^-(P + 3), and 2^x - 1 is just above -1.
  if (wide <= Wide1<F>(-(F::kPrecision + 3))) {
    return RoundBeside(F(-1), true, env);
  }
  if (IsInteger(x)) {
    // 2^n - 1 is exact in the working format up to n = P + 3, and rounds as
    // just below 2^n beyond, with no rounding boundary between them.
    const auto n = to_integer<std::int64_t>(x, Rounding::kTowardZero);
    const Wide1<F> one(1);
    const Wide1<F> power = scalbn(one, static_cast<int>(n));
    if (n <= F::kPrecision + 3) {
      return FloatAccess::Convert<F>(power - one, env);
    }
    return RoundJustBelow<F>(false, power, env);
  }
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    const MathConstants<W>& c = GetMathConstants<W>();
    const W wx = FloatAccess::Convert<W>(x, NearestEnv{});
    const W nearest = round(wx);
    const auto k = to_integer<std::int64_t>(nearest, Rounding::kNearestEven);
    // The difference is exact, and r is within 2 u.
    const W r = (wx - nearest) * c.ln2;
    int units;
    W y = ExpM1Reduced<W, std::is_same_v<W, Wide1<F>>>(r, units);
    if (k != 0) {
      // |2^x| <= 3.5 |2^x - 1| for |x| >= 1/2.
      y = scalbn(y + W(1), static_cast<int>(k)) - W(1);
      units = 4 * units + 8;
    }
    return RoundTranscendental(y, W::kPrecision - LostBits(units + 4), last,
                               env, result);
  });
}

template <class F, class Env, int kFirstLevel = 1>
F Exp10m1Impl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_inf()) return x.signbit() ? F(-1) : x;
  if (x.is_zero()) return x;
  const Wide1<F> t = TimesLog2Of10(x);
  if (t >= Wide1<F>(F::kEmax + 2)) return RoundOverflow<F>(false, env);
  // Then 10^x < 2^-(P + 3), and 10^x - 1 is just above -1.
  if (t <= Wide1<F>(-(F::kPrecision + 3))) {
    return RoundBeside(F(-1), true, env);
  }
  if (!x.signbit() && IsInteger(x) && x <= F(F::kPrecision / 2 + 1)) {
    // 10^n - 1 is exact or a tie if it has at most P + 1 bits. Beyond, if
    // 10^n = 5^n 2^n has at most P + 1 bits itself, 10^n - 1 lies just below
    // it, with no rounding boundary between them, which are 2 or more apart
    // there. Otherwise 10^n - 1 rounds as 10^n, which is as far from a
    // boundary as 10^x at other x.
    const auto n = to_integer<std::int64_t>(x, Rounding::kTowardZero);
    const BigUint power = BigUint::Pow5(static_cast<std::uint64_t>(n));
    BigUint value = power;
    value.ShiftLeft(static_cast<std::size_t>(n));
    value.Sub(BigUint(1));
    if (value.bit_length() <= static_cast<std::size_t>(F::kPrecision) + 1) {
      return RoundInteger<F>(false, value, 0, false, env);
    }
    if (power.bit_length() <= static_cast<std::size_t>(F::kPrecision) + 1) {
      return RoundJustBelow<F>(
          false,
          RoundInteger<Wide1<F>>(false, power, n, false,
                                 SilentEnv{Rounding::kNearestEven}),
          env);
    }
  }
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    std::int64_t k;
    const W r = ReduceByLog2Of10<W>(x, k);
    int units;
    W y = ExpM1Reduced<W, std::is_same_v<W, Wide1<F>>>(r, units);
    if (k != 0) {
      // |10^x| <= 3.5 |10^x - 1| for |x log2(10)| >= 1/2.
      y = scalbn(y + W(1), static_cast<int>(k)) - W(1);
      units = 4 * units + 8;
    }
    return RoundTranscendental(y, W::kPrecision - LostBits(units + 4), last,
                               env, result);
  });
}

// log2(1 + x) and log10(1 + x). They are exact where 1 + x is a power of the
// base, and for x a power of the base above 2^(P + 4), just above its
// logarithm, an integer n, by less than any working format resolves:
// log(1 + x) - log x < 2^-(P + 4) / ln 2, while rounding boundaries lie
// n 2^-(P + 1) or more apart there.
template <int kBase, class F, class Env, int kFirstLevel = 1>
F LogPlus1Impl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  F special;
  if (Log1pSpecialCase(x, env, special)) return special;
  using W1 = Wide1<F>;
  std::int64_t n;
  // 1 + x is a power of the base only for x >= -1/2, from where it is exact
  // below 2^P.
  if (x.raw_exponent() >= -1 && x.raw_exponent() < F::kPrecision) {
    const W1 z = W1(1) + FloatAccess::Convert<W1>(x, NearestEnv{});
    if constexpr (kBase == 2) {
      if (SignificantBits(z) == 1) {
        return FloatAccess::Convert<F>(std::int64_t{z.raw_exponent()}, env);
      }
    } else {
      if (IsPowerOfTen(z, n)) return FloatAccess::Convert<F>(n, env);
    }
  } else if (x.raw_exponent() >= F::kPrecision + 4) {
    bool power = false;
    if constexpr (kBase == 2) {
      power = SignificantBits(x) == 1;
      n = x.raw_exponent();
    } else {
      power = IsPowerOfTen(x, n);
    }
    if (power) return RoundJustAbove<F>(W1(static_cast<double>(n)), env);
  }
  return Log1pEvaluate<kBase, F, kFirstLevel>(x, env);
}

// --- Powers and roots ---

// The exact x^2 + y^2 always fits: the exponents differ by at most P + 2.
template <class F>
using HypotWorking = Float<4 * F::kPrecision + 8, 29>;

template <class F, class Env>
F HypotImpl(const F& x, const F& y, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  // An infinity wins over a quiet NaN, but not over a signaling one.
  if (x.is_signaling() || y.is_signaling()) {
    return FloatAccess::PropagateNaN(x, y, env);
  }
  if (x.is_inf() || y.is_inf()) return F::Infinity();
  if (x.is_nan() || y.is_nan()) return FloatAccess::PropagateNaN(x, y, env);
  F a = abs(x);
  F b = abs(y);
  if (a < b) std::swap(a, b);
  if (b.is_zero()) return a;
  // Then b / a < 2^-(P + 2), and hypot(a, b) - a < a 2^-(2P + 4).
  if (std::int64_t{a.raw_exponent()} - b.raw_exponent() > F::kPrecision + 2) {
    return RoundBeside(a, true, env);
  }
  // Within 2 u in the first working format: a^2 + b^2 within 1.5 u, its root
  // within 0.75 u, and 0.5 u more for rounding the root.
  using W = Wide1<F>;
  const W wa = FloatAccess::Convert<W>(a, NearestEnv{});
  const W wb = FloatAccess::Convert<W>(b, NearestEnv{});
  const W root = sqrt(wa * wa + wb * wb);
  // Exact results have at most P + 1 bits too, so this never decides one.
  F result;
  if (RoundClearOfBoundaries(root, W::kPrecision - 4, env, result)) {
    return result;
  }
  // Exactly then, with x^2 + y^2 exact.
  using Exact = HypotWorking<F>;
  const Exact ea = FloatAccess::Convert<Exact>(a, NearestEnv{});
  const Exact eb = FloatAccess::Convert<Exact>(b, NearestEnv{});
  Context context{Rounding::kTowardZero};
  const Exact exact_root = sqrt(ea * ea + eb * eb, context);
  return RoundTruncated<F>(false, exact_root,
                           (context.flags & kFlagInexact) != 0, env);
}

// The precision in which a W-precision Newton step can take its start and
// compute its correction.
template <class W>
using HalfWorking = Float<64 * ((W::kPrecision / 2 + 16 + 63) / 64), 29>;

// m^(-1/3) for m in [1, 8), with a relative error below 2 u. Newton's
// iteration z <- z + z (1 - m z^3) / 3 takes the relative error from e to
// about 2 e^2, without a division, and adds below 1.4 u of rounding, which the
// next step removes. Wide formats start from a root in about half their
// precision, and compute the correction in that precision.
template <class W>
W InverseCbrtOfReduced(const W& m) {
  const W one(1);
  if constexpr (W::kPrecision <= 128) {
    // The double estimate is within e = 2^-50. With r = 1 - m z^3, the root is
    // z (1 - r)^(-1/3) = z (1 + r/3 + 2r^2/9 + 14r^3/81 + ...), so one step
    // of z <- z (1 + r (3 + 2r) / 9) leaves about 14 (3e)^3 / 81 < 2^-147.
    W z(1.0 / std::cbrt(static_cast<double>(m)));
    const W residual = one - m * ((z * z) * z);
    if (residual.is_zero()) return z;
    return z + DivSmall(z * (residual * (W(3) + scalbn(residual, 1))), 9);
  } else if constexpr (W::kPrecision <= 360) {
    // The same third-order step from a 128-bit root, within 2^-126.6, leaves
    // below 2^-377.
    using S = Float<128, 29>;
    const W z = FloatAccess::Convert<W>(
        InverseCbrtOfReduced(FloatAccess::Convert<S>(m, NearestEnv{})),
        NearestEnv{});
    const W residual = one - m * ((z * z) * z);
    if (residual.is_zero()) return z;
    return z + DivSmall(z * (residual * (W(3) + scalbn(residual, 1))), 9);
  } else {
    // The half-precision root is within 2^-(H - 2), so one step leaves
    // 2^-(2H - 3) < 2^-(kPrecision + 28) besides the rounding.
    using H = HalfWorking<W>;
    const H half =
        InverseCbrtOfReduced(FloatAccess::Convert<H>(m, NearestEnv{}));
    const W z = FloatAccess::Convert<W>(half, NearestEnv{});
    const W residual = one - m * ((z * z) * z);
    if (residual.is_zero()) return z;
    const H correction =
        DivSmall(half * FloatAccess::Convert<H>(residual, NearestEnv{}), 3);
    return z + FloatAccess::Convert<W>(correction, NearestEnv{});
  }
}

// cbrt(m) for m in [1, 8), with a relative error below 5 u.
template <class W>
W CbrtOfReduced(const W& m) {
  if constexpr (W::kPrecision <= 128) {
    // From z within 2 u, m z^2 is within 5 u.
    const W z = InverseCbrtOfReduced(m);
    return m * (z * z);
  } else {
    // One Newton step y <- y - (y^3 - m) z^2 / 3 on y = m z^2 from half the
    // precision, where 1 / y^2 = z^2 to that precision. That leaves the
    // rounding of y^3, at most 0.7 u, and of the step, 0.5 u.
    using H = HalfWorking<W>;
    const H mh = FloatAccess::Convert<H>(m, NearestEnv{});
    const H z = InverseCbrtOfReduced(mh);
    const H z2 = z * z;
    const W y = FloatAccess::Convert<W>(mh * z2, NearestEnv{});
    const W residual = (y * y) * y - m;
    if (residual.is_zero()) return y;
    const H correction =
        DivSmall(FloatAccess::Convert<H>(residual, NearestEnv{}) * z2, 3);
    return y - FloatAccess::Convert<W>(correction, NearestEnv{});
  }
}

template <class F, class Env, int kFirstLevel = 1>
F CbrtImpl(const F& x, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  if (x.is_nan()) return FloatAccess::PropagateNaN(x, env);
  if (x.is_inf() || x.is_zero()) return x;
  const bool negative = x.signbit();
  const F magnitude = abs(x);
  // |x| = 2^(3k) m with m in [1, 8).
  const std::int64_t e = x.raw_exponent();
  const std::int64_t k = (e >= 0 ? e : e - 2) / 3;
  // The root is exact only if x = a 2^b, a odd, has 3 | b and a a cube.
  const std::int64_t low = LowestBitExponent(x);
  const std::int64_t bits = SignificantBits(x);
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    const W m = scalbn(FloatAccess::Convert<W>(magnitude, NearestEnv{}),
                       static_cast<int>(-3 * k));
    const W y = scalbn(CbrtOfReduced(m), static_cast<int>(k));
    const W signed_y = negative ? -y : y;
    constexpr int kBits = W::kPrecision - LostBits(5);
    // An exact root has at most P bits, so this never decides one.
    if (RoundClearOfBoundaries(signed_y, kBits, env, result)) return true;
    if (low % 3 == 0) {
      const F root = FloatAccess::Convert<F>(y, NearestEnv{});
      const std::int64_t root_bits = SignificantBits(root);
      if (3 * root_bits - 2 <= bits && bits <= 3 * root_bits) {
        using Cube = Float<3 * F::kPrecision, 29>;
        const Cube c = FloatAccess::Convert<Cube>(root, NearestEnv{});
        if (c * c * c == FloatAccess::Convert<Cube>(magnitude, NearestEnv{})) {
          result = negative ? -root : root;
          return true;
        }
      }
    }
    return RoundTranscendental(signed_y, kBits, last, env, result);
  });
}

// Rounds x^y for a finite positive x != 1 and a finite nonzero y if it has at
// most P + 1 bits, which covers every exact result and every tie. Returns
// false if it has more. |log2(x^y)| must be below 2^27.
template <class F, class Env>
bool ExactPow(const F& x, const F& y, bool negative, const Env& env,
              F& result) {
  const std::int64_t low = LowestBitExponent(x);
  if (SignificantBits(x) == 1) {
    // x = 2^low, and x^y is exact if low y is an integer.
    using W = Wide1<F>;  // low * y is exact.
    const W product = W(low) * FloatAccess::Convert<W>(y, NearestEnv{});
    if (!IsInteger(product)) return false;
    const limb_t top_bit = limb_t{1} << (kLimbBits - 1);
    result = FloatAccess::Round<F>(
        negative, to_integer<std::int64_t>(product, Rounding::kTowardZero),
        &top_bit, 1, false, env);
    return true;
  }
  // x = a 2^low with a odd and above 1, so x^y needs y = n / 2^s with n > 0
  // and x a perfect 2^s-th power, z^(2^s) = x; then x^y = z^n.
  const std::int64_t s = (std::max)(std::int64_t{0}, -LowestBitExponent(y));
  // a >= 3 has at least 1.5 * 2^s bits.
  if (y.signbit() || (s > 0 && (s > 20 || (std::int64_t{1} << s) >=
                                              SignificantBits(x)))) {
    return false;
  }
  F z = x;
  for (std::int64_t i = 0; i < s; ++i) {
    Context context;
    z = sqrt(z, context);
    if (context.flags != 0) return false;
  }
  const F n_float = scalbn(y, static_cast<int>(s));
  if (n_float >= F(static_cast<double>(F::kPrecision) + 1)) return false;
  const auto n = to_integer<std::int64_t>(n_float, Rounding::kTowardZero);
  // a^n, for z = a' 2^z_low, has at least n (bits - 1) + 1 bits.
  const std::int64_t z_low = LowestBitExponent(z);
  if (n * (SignificantBits(z) - 1) > F::kPrecision) return false;
  // So z^n has at most 2P bits, and every partial product is exact.
  using W = MathWorking<F, 2>;
  W base = scalbn(FloatAccess::Convert<W>(z, NearestEnv{}),
                  static_cast<int>(-z_low));
  W power(1);
  for (std::int64_t rest = n;;) {
    if ((rest & 1) != 0) power = power * base;
    rest >>= 1;
    if (rest == 0) break;
    base = base * base;
  }
  if (SignificantBits(power) > F::kPrecision + 1) return false;
  result = FloatAccess::Round<F>(negative, power.raw_exponent() + z_low * n,
                                 power.raw_limbs().data(), W::kLimbs, false,
                                 env);
  return true;
}

template <class F, class Env, int kFirstLevel = 1>
F PowImpl(const F& x, const F& y, const Env& env) {
  static_assert(F::kExponentBits <= 26, "exponent range too wide");
  // The special cases of IEEE 754 and C99 Annex F.
  if (x.is_signaling() || y.is_signaling()) {
    return FloatAccess::PropagateNaN(x, y, env);
  }
  const F one(1);
  if (y.is_zero() || x == one) return one;
  if (x.is_nan() || y.is_nan()) return FloatAccess::PropagateNaN(x, y, env);
  const F magnitude = abs(x);
  if (y.is_inf()) {
    if (magnitude == one) return one;
    return (magnitude > one) == y.signbit() ? F::Zero() : F::Infinity();
  }
  const bool odd = LowestBitExponent(y) == 0;
  const bool negative = x.signbit() && odd;
  if (x.is_zero()) {
    if (!y.signbit()) return F::Zero(negative);
    env.raise(kFlagDivByZero);
    return F::Infinity(negative);
  }
  if (x.is_inf()) return y.signbit() ? F::Zero(negative) : F::Infinity(negative);
  if (x.signbit() && !IsInteger(y)) return InvalidResult<F>(env);
  if (magnitude == one) return negative ? -one : one;

  // t = y log|x|, with |x^y| = e^t.
  using W1 = Wide1<F>;
  int log_units;
  std::int64_t log_exponent;
  W1 log_m;
  const W1 y1 = FloatAccess::Convert<W1>(y, NearestEnv{});
  const W1 t1 =
      y1 * LogKernel<W1, true>(FloatAccess::Convert<W1>(magnitude, NearestEnv{}),
                               log_units, log_exponent, log_m);
  // As for exp, with a relative error of t far below these margins.
  if (t1 >= W1(F::kEmax + 2)) return RoundOverflow<F>(negative, env);
  if (t1 <= W1(-(F::kEmax + F::kPrecision + 4))) {
    return RoundTiny<F>(negative, env);
  }
  if (t1.raw_exponent() <= -(F::kPrecision + 4)) {
    return RoundBeside(negative ? -one : one, (t1 > W1::Zero()) != negative,
                       env);
  }
  {
    F exact;
    if (ExactPow(magnitude, y, negative, env, exact)) return exact;
  }
  return ZivEvaluate<F, kFirstLevel>([&](auto tag, bool last, F& result) {
    using W = typename decltype(tag)::type;
    constexpr bool kFirst = std::is_same_v<W, W1>;
    W t;
    int units_log;
    if constexpr (kFirst) {
      t = t1;
      units_log = log_units;
    } else {
      std::int64_t e;
      W m;
      t = FloatAccess::Convert<W>(y, NearestEnv{}) *
          LogKernel<W>(FloatAccess::Convert<W>(magnitude, NearestEnv{}),
                       units_log, e, m);
    }
    // t has a relative error below (units_log + 1) u, so an absolute one
    // below |t| (units_log + 1) u, which carries over to e^t relatively.
    // The reduction is exact up to 0.36 u since |t| < 2^26.
    std::int64_t k;
    const W r = ReduceByLn2<W>(t, k);
    int units_exp;
    const W q = ExpM1Reduced<W, kFirst>(r, units_exp);
    const W value = scalbn(q + W(1), static_cast<int>(k));
    const std::int64_t size =
        t.raw_exponent() >= 0 ? std::int64_t{2} << t.raw_exponent() : 1;
    const std::int64_t units = units_exp + 4 + size * (units_log + 2);
    return RoundTranscendental(negative ? -value : value,
                               W::kPrecision - LostBitsWide(units), last, env,
                               result);
  });
}

}  // namespace detail

// e^x.
template <int P, int E>
Float<P, E> exp(const Float<P, E>& x) {
  return detail::ExpImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> exp(const Float<P, E>& x, Context& context) {
  return detail::ExpImpl(x, detail::ContextEnv{context});
}

// 2^x.
template <int P, int E>
Float<P, E> exp2(const Float<P, E>& x) {
  return detail::Exp2Impl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> exp2(const Float<P, E>& x, Context& context) {
  return detail::Exp2Impl(x, detail::ContextEnv{context});
}

// e^x - 1, accurate for small x.
template <int P, int E>
Float<P, E> expm1(const Float<P, E>& x) {
  return detail::Expm1Impl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> expm1(const Float<P, E>& x, Context& context) {
  return detail::Expm1Impl(x, detail::ContextEnv{context});
}

// The natural logarithm. log(+-0) is -infinity and raises divide by zero, and
// a negative x gives NaN and raises invalid.
template <int P, int E>
Float<P, E> log(const Float<P, E>& x) {
  return detail::LogImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> log(const Float<P, E>& x, Context& context) {
  return detail::LogImpl(x, detail::ContextEnv{context});
}

// The base-2 logarithm, exact for powers of two.
template <int P, int E>
Float<P, E> log2(const Float<P, E>& x) {
  return detail::Log2Impl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> log2(const Float<P, E>& x, Context& context) {
  return detail::Log2Impl(x, detail::ContextEnv{context});
}

// The base-10 logarithm, exact for powers of ten.
template <int P, int E>
Float<P, E> log10(const Float<P, E>& x) {
  return detail::Log10Impl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> log10(const Float<P, E>& x, Context& context) {
  return detail::Log10Impl(x, detail::ContextEnv{context});
}

// log(1 + x), accurate for small x. log1p(-1) is -infinity and raises divide
// by zero, and x < -1 gives NaN and raises invalid.
template <int P, int E>
Float<P, E> log1p(const Float<P, E>& x) {
  return detail::Log1pImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> log1p(const Float<P, E>& x, Context& context) {
  return detail::Log1pImpl(x, detail::ContextEnv{context});
}

// The hyperbolic sine.
template <int P, int E>
Float<P, E> sinh(const Float<P, E>& x) {
  return detail::SinhImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> sinh(const Float<P, E>& x, Context& context) {
  return detail::SinhImpl(x, detail::ContextEnv{context});
}

// The hyperbolic cosine.
template <int P, int E>
Float<P, E> cosh(const Float<P, E>& x) {
  return detail::CoshImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> cosh(const Float<P, E>& x, Context& context) {
  return detail::CoshImpl(x, detail::ContextEnv{context});
}

// The hyperbolic tangent.
template <int P, int E>
Float<P, E> tanh(const Float<P, E>& x) {
  return detail::TanhImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> tanh(const Float<P, E>& x, Context& context) {
  return detail::TanhImpl(x, detail::ContextEnv{context});
}

// The inverse hyperbolic sine.
template <int P, int E>
Float<P, E> asinh(const Float<P, E>& x) {
  return detail::AsinhImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> asinh(const Float<P, E>& x, Context& context) {
  return detail::AsinhImpl(x, detail::ContextEnv{context});
}

// The inverse hyperbolic cosine. x < 1 gives NaN and raises invalid.
template <int P, int E>
Float<P, E> acosh(const Float<P, E>& x) {
  return detail::AcoshImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> acosh(const Float<P, E>& x, Context& context) {
  return detail::AcoshImpl(x, detail::ContextEnv{context});
}

// The inverse hyperbolic tangent. atanh(+-1) is +-infinity and raises divide
// by zero, and |x| > 1 gives NaN and raises invalid.
template <int P, int E>
Float<P, E> atanh(const Float<P, E>& x) {
  return detail::AtanhImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> atanh(const Float<P, E>& x, Context& context) {
  return detail::AtanhImpl(x, detail::ContextEnv{context});
}

// sin x, for x in radians. sin(+-infinity) is NaN and raises invalid. Large
// arguments are reduced exactly, as if with infinitely many digits of pi.
template <int P, int E>
Float<P, E> sin(const Float<P, E>& x) {
  return detail::SinImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> sin(const Float<P, E>& x, Context& context) {
  return detail::SinImpl(x, detail::ContextEnv{context});
}

// cos x, for x in radians. cos(+-infinity) is NaN and raises invalid.
template <int P, int E>
Float<P, E> cos(const Float<P, E>& x) {
  return detail::CosImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> cos(const Float<P, E>& x, Context& context) {
  return detail::CosImpl(x, detail::ContextEnv{context});
}

// tan x, for x in radians. tan(+-infinity) is NaN and raises invalid. No
// finite x is a pole, so the result is always finite.
template <int P, int E>
Float<P, E> tan(const Float<P, E>& x) {
  return detail::TanImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> tan(const Float<P, E>& x, Context& context) {
  return detail::TanImpl(x, detail::ContextEnv{context});
}

// pi, correctly rounded; mpfix::pi<binary256>(). It raises inexact.
template <class F>
F pi() {
  return detail::PiImpl<F>(detail::NearestEnv{});
}
template <class F>
F pi(Context& context) {
  return detail::PiImpl<F>(detail::ContextEnv{context});
}

// 10^x, exact for 10^n that F holds.
template <int P, int E>
Float<P, E> exp10(const Float<P, E>& x) {
  return detail::Exp10Impl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> exp10(const Float<P, E>& x, Context& context) {
  return detail::Exp10Impl(x, detail::ContextEnv{context});
}

// 2^x - 1, accurate for small x.
template <int P, int E>
Float<P, E> exp2m1(const Float<P, E>& x) {
  return detail::Exp2m1Impl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> exp2m1(const Float<P, E>& x, Context& context) {
  return detail::Exp2m1Impl(x, detail::ContextEnv{context});
}

// 10^x - 1, accurate for small x.
template <int P, int E>
Float<P, E> exp10m1(const Float<P, E>& x) {
  return detail::Exp10m1Impl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> exp10m1(const Float<P, E>& x, Context& context) {
  return detail::Exp10m1Impl(x, detail::ContextEnv{context});
}

// log2(1 + x), accurate for small x. log2p1(-1) is -infinity and raises divide
// by zero, and x < -1 gives NaN and raises invalid.
template <int P, int E>
Float<P, E> log2p1(const Float<P, E>& x) {
  return detail::LogPlus1Impl<2>(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> log2p1(const Float<P, E>& x, Context& context) {
  return detail::LogPlus1Impl<2>(x, detail::ContextEnv{context});
}

// log10(1 + x), accurate for small x, with the special cases of log2p1.
template <int P, int E>
Float<P, E> log10p1(const Float<P, E>& x) {
  return detail::LogPlus1Impl<10>(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> log10p1(const Float<P, E>& x, Context& context) {
  return detail::LogPlus1Impl<10>(x, detail::ContextEnv{context});
}

// The inverse tangent, in (-pi/2, pi/2); atan(+-infinity) is +-pi/2 rounded.
template <int P, int E>
Float<P, E> atan(const Float<P, E>& x) {
  return detail::AtanImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> atan(const Float<P, E>& x, Context& context) {
  return detail::AtanImpl(x, detail::ContextEnv{context});
}

// The inverse sine, in [-pi/2, pi/2]. |x| > 1 gives NaN and raises invalid.
template <int P, int E>
Float<P, E> asin(const Float<P, E>& x) {
  return detail::AsinImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> asin(const Float<P, E>& x, Context& context) {
  return detail::AsinImpl(x, detail::ContextEnv{context});
}

// The inverse cosine, in [0, pi]. |x| > 1 gives NaN and raises invalid.
template <int P, int E>
Float<P, E> acos(const Float<P, E>& x) {
  return detail::AcosImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> acos(const Float<P, E>& x, Context& context) {
  return detail::AcosImpl(x, detail::ContextEnv{context});
}

// The angle of the point (x, y), in [-pi, pi], with the special cases of
// IEEE 754 for zeros and infinities: atan2(+-0, -0) is +-pi and
// atan2(+-infinity, -infinity) is +-3pi/4, for example. It raises no flag for
// zeros, and underflow for results below the normal range, as usual.
template <int P, int E>
Float<P, E> atan2(const Float<P, E>& y, const Float<P, E>& x) {
  return detail::Atan2Impl(y, x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> atan2(const Float<P, E>& y, const Float<P, E>& x,
                  Context& context) {
  return detail::Atan2Impl(y, x, detail::ContextEnv{context});
}

// x^y, with the special cases of IEEE 754: pow(x, +-0) and pow(1, y) are 1
// even for a quiet NaN, a negative x needs an integer y, and pow(+-0, y) for
// y < 0 is an infinity and raises divide by zero. Exact results, such as
// pow(9, 1.5), raise no flags.
template <int P, int E>
Float<P, E> pow(const Float<P, E>& x, const Float<P, E>& y) {
  return detail::PowImpl(x, y, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> pow(const Float<P, E>& x, const Float<P, E>& y, Context& context) {
  return detail::PowImpl(x, y, detail::ContextEnv{context});
}

// The real cube root, exact for perfect cubes.
template <int P, int E>
Float<P, E> cbrt(const Float<P, E>& x) {
  return detail::CbrtImpl(x, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> cbrt(const Float<P, E>& x, Context& context) {
  return detail::CbrtImpl(x, detail::ContextEnv{context});
}

// sqrt(x^2 + y^2) without intermediate overflow or underflow. An infinity
// gives +infinity even with a quiet NaN.
template <int P, int E>
Float<P, E> hypot(const Float<P, E>& x, const Float<P, E>& y) {
  return detail::HypotImpl(x, y, detail::NearestEnv{});
}
template <int P, int E>
Float<P, E> hypot(const Float<P, E>& x, const Float<P, E>& y,
                  Context& context) {
  return detail::HypotImpl(x, y, detail::ContextEnv{context});
}

}  // namespace mpfix

#endif  // MPFIX_MATH_H_
