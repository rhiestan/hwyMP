// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// pi and 2/pi in fixed point to any number of bits, for the trigonometric
// functions and their argument reduction. They come from the Chudnovsky series
// by binary splitting, with Newton's method for the large quotients and the
// square root, so that the cost grows like that of Karatsuba multiplication.
// The first 32768 bits are computed once, on first use, in about 3 ms; more,
// which only arguments with binary exponents above about 32000 need, are
// computed and kept when first asked for: about 0.3 s for 2^19 bits and 8 s
// for 2^22, which binary512's largest arguments need.

#ifndef MPFIX_PI_H_
#define MPFIX_PI_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>

#include "mpfix/bigint.h"

namespace mpfix::detail {

// Bits [low, low + count) of v, as an integer.
inline BigUint BitRange(const BigUint& v, std::size_t low,
                           std::size_t count) {
  const auto& limbs = v.limbs();
  const std::size_t first = low / kLimbBits;
  const int shift = static_cast<int>(low % kLimbBits);
  const std::size_t words = (count + kLimbBits - 1) / kLimbBits;
  const auto limb = [&](std::size_t i) -> limb_t {
    return i < limbs.size() ? limbs[i] : 0;
  };
  const auto fill = [&](limb_t* out) {
    for (std::size_t i = 0; i < words; ++i) {
      const limb_t low_part = limb(first + i);
      out[i] = shift == 0 ? low_part
                          : (low_part >> shift) |
                                (limb(first + i + 1) << (kLimbBits - shift));
    }
    const int top = static_cast<int>(count % kLimbBits);
    if (top != 0) out[words - 1] &= (limb_t{1} << top) - 1;
  };
  constexpr std::size_t kInline = 32;
  if (words <= kInline) {
    std::array<limb_t, kInline> out;
    fill(out.data());
    return BigUint::FromLimbs(out.data(), words);
  }
  std::vector<limb_t> out(words);
  fill(out.data());
  return BigUint::FromLimbs(out.data(), words);
}

// Divisors of at least this many limbs divide by Newton's method.
inline constexpr std::size_t kNewtonLimbs = 64;

// About 2^(n + m) / d for n = bit_length(d) and m >= 1, within a relative
// 2^-(m - 2). Only the top m + 64 bits of d matter.
inline BigUint ApproximateReciprocal(const BigUint& d, std::size_t m) {
  const std::size_t n = d.bit_length();
  if (m <= kLimbBits * kNewtonLimbs) {
    // Directly, from the top m + 64 bits of d: within 2^-m from the floor and
    // 2^-(m + 63) from the truncation.
    const std::size_t drop = n > m + 64 ? n - (m + 64) : 0;
    BigUint top = d;
    top.ShiftRight(drop);
    BigUint numerator(1);
    numerator.ShiftLeft(n - drop + m);
    BigUint quotient;
    BigUint remainder;
    BigUint::DivRem(numerator, top, quotient, remainder);
    return quotient;
  }
  // From half the precision, within e = 2^-(m/2 + 30), one Newton step
  // x += x (1 - d x) leaves e^2 < 2^-(m + 60), plus truncations below
  // 2^-(m - 1).
  const std::size_t h = m / 2 + 32;
  BigUint x = ApproximateReciprocal(d, h);
  x.ShiftLeft(m - h);
  const std::size_t drop = n > m + 32 ? n - (m + 32) : 0;
  BigUint top = d;
  top.ShiftRight(drop);
  const std::size_t scale = n - drop + m;  // x ~ 2^scale / top.
  BigUint residual = top * x;
  BigUint one(1);
  one.ShiftLeft(scale);
  const bool above = Compare(residual, one) > 0;
  if (above) {
    residual.Sub(one);
  } else {
    one.Sub(residual);
    residual = std::move(one);
  }
  BigUint delta = x * residual;
  delta.ShiftRight(scale);
  if (above) {
    x.Sub(delta);
  } else {
    x.Add(delta);
  }
  return x;
}

// floor(a / d), for d != 0. Divisors of at least `newton_limbs` limbs divide
// through a reciprocal.
inline BigUint Quotient(const BigUint& a, const BigUint& d,
                        std::size_t newton_limbs = kNewtonLimbs) {
  if (Compare(a, d) < 0) return BigUint();
  if (d.limbs().size() < newton_limbs) {
    BigUint quotient;
    BigUint remainder;
    BigUint::DivRem(a, d, quotient, remainder);
    return quotient;
  }
  const std::size_t n = d.bit_length();
  const std::size_t quotient_bits = a.bit_length() - n + 1;
  const std::size_t m = quotient_bits + 16;
  const BigUint reciprocal = ApproximateReciprocal(d, m);  // ~ 2^(n + m) / d.
  // a / d ~ a r / 2^(n + m), from the top bits of a; the estimate is within a
  // few units, and the remainder corrects it.
  const std::size_t keep = quotient_bits + 32;
  const std::size_t drop = a.bit_length() > keep ? a.bit_length() - keep : 0;
  BigUint top = a;
  top.ShiftRight(drop);
  BigUint q = top * reciprocal;
  q.ShiftRight(n + m - drop);
  BigUint product = q * d;
  const BigUint one(1);
  while (Compare(product, a) > 0) {
    q.Sub(one);
    product.Sub(d);
  }
  BigUint remainder = a;
  remainder.Sub(product);
  while (Compare(remainder, d) >= 0) {
    q.Add(one);
    remainder.Sub(d);
  }
  return q;
}

// floor(sqrt(a)).
inline BigUint SqrtFloor(const BigUint& a) {
  const std::size_t bits = a.bit_length();
  if (bits == 0) return BigUint();
  BigUint x;
  if (bits <= 128) {
    x = BigUint(1);
    x.ShiftLeft((bits + 1) / 2);  // Above sqrt(a).
  } else {
    // From the root of the top half of the bits, which is within a relative
    // 2^-(bits/4), so that two or three Newton steps follow.
    const std::size_t shift = bits / 4;
    BigUint high = a;
    high.ShiftRight(2 * shift);
    x = SqrtFloor(high);
    x.Add(BigUint(1));
    x.ShiftLeft(shift);  // Above sqrt(a).
  }
  // Newton's iteration from above decreases until it reaches the floor.
  for (;;) {
    BigUint y = Quotient(a, x);
    y.Add(x);
    y.ShiftRight(1);
    if (Compare(y, x) >= 0) return x;
    x = std::move(y);
  }
}

// A signed big integer, for the binary splitting.
struct SignedBig {
  BigUint magnitude;
  bool negative = false;
};

inline SignedBig AddSigned(SignedBig a, const SignedBig& b) {
  if (a.negative == b.negative) {
    a.magnitude.Add(b.magnitude);
    return a;
  }
  if (Compare(a.magnitude, b.magnitude) >= 0) {
    a.magnitude.Sub(b.magnitude);
    if (a.magnitude.is_zero()) a.negative = false;
    return a;
  }
  SignedBig result{b.magnitude, b.negative};
  result.magnitude.Sub(a.magnitude);
  return result;
}

// P, Q and R of the Chudnovsky series over the terms [a, b), by binary
// splitting: P(a, a + 1) = -(6a - 5)(2a - 1)(6a - 1),
// Q(a, a + 1) = a^3 640320^3 / 24 and R(a, a + 1) = P (13591409 + 545140134 a),
// with P(a, c) = P(a, b) P(b, c), Q likewise and
// R(a, c) = Q(b, c) R(a, b) + P(a, b) R(b, c). P is negative for an odd
// number of terms; `p` holds its magnitude, and only if `need_p`.
struct ChudnovskyTerms {
  BigUint p;
  BigUint q;
  SignedBig r;
};

inline void ChudnovskySplit(limb_t a, limb_t b, bool need_p,
                            ChudnovskyTerms& out) {
  if (b - a == 1) {
    out.p = BigUint(6 * a - 5);
    out.p.MulAdd(2 * a - 1, 0);
    out.p.MulAdd(6 * a - 1, 0);
    out.q = BigUint(a);
    out.q.MulAdd(a, 0);
    out.q.MulAdd(a, 0);
    out.q.MulAdd(10939058860032000u, 0);  // 640320^3 / 24.
    out.r.magnitude = out.p;
    out.r.magnitude.MulAdd(13591409 + 545140134 * a, 0);
    out.r.negative = true;
    return;
  }
  const limb_t middle = a + (b - a) / 2;
  ChudnovskyTerms left;
  ChudnovskyTerms right;
  ChudnovskySplit(a, middle, true, left);
  ChudnovskySplit(middle, b, need_p, right);
  const bool left_p_negative = ((middle - a) & 1) != 0;
  out.r = AddSigned(SignedBig{right.q * left.r.magnitude, left.r.negative},
                    SignedBig{left.p * right.r.magnitude,
                              right.r.negative != left_p_negative});
  out.q = left.q * right.q;
  if (need_p) out.p = left.p * right.p;
}

// pi and 2/pi as integers within 1 of pi 2^bits and 2/pi 2^bits.
struct PiFixedPoint {
  std::size_t bits = 0;
  BigUint pi;
  BigUint two_over_pi;
};

inline PiFixedPoint ComputePi(std::size_t bits) {
  // With f fraction bits and terms below 2^-(f + 64): each term gains 47.1.
  const std::size_t f = bits + 64;
  const auto terms = static_cast<limb_t>((f + 64) / 47 + 2);
  ChudnovskyTerms s;
  ChudnovskySplit(1, terms, false, s);
  // pi = 426880 sqrt(10005) Q / D for D = 13591409 Q + R, and D / Q is
  // about 2^23.7.
  SignedBig d{s.q, false};
  d.magnitude.MulAdd(13591409, 0);
  d = AddSigned(std::move(d), s.r);
  // ratio = D / Q 2^f, from D and Q cut to their top f + 96 bits, which
  // moves it by a relative 2^-(f + 94): within 2 units.
  const std::size_t keep = f + 96;
  BigUint dt = d.magnitude;
  BigUint qt = s.q;
  const std::size_t d_drop = dt.bit_length() > keep ? dt.bit_length() - keep : 0;
  const std::size_t q_drop = qt.bit_length() > keep ? qt.bit_length() - keep : 0;
  dt.ShiftRight(d_drop);
  qt.ShiftRight(q_drop);
  dt.ShiftLeft(f + d_drop - q_drop);  // D's top bits exceed Q's by about 24.
  const BigUint ratio = Quotient(dt, qt);
  // sqrt(10005) 2^f, within 1.
  BigUint radicand(10005);
  radicand.ShiftLeft(2 * f);
  const BigUint root = SqrtFloor(radicand);

  PiFixedPoint result;
  result.bits = bits;
  // 2/pi = (D / Q) sqrt(10005) / 2135467200, within 3 units of 2^-f.
  result.two_over_pi = ratio * root;
  result.two_over_pi.ShiftRight(f);
  result.two_over_pi.DivSmall(2135467200);
  result.two_over_pi.ShiftRight(f - bits);
  // pi = 426880 sqrt(10005) / (D / Q), within 2 units of 2^-f.
  BigUint numerator = root;
  numerator.MulAdd(426880, 0);
  numerator.ShiftLeft(f);
  result.pi = Quotient(numerator, ratio);
  result.pi.ShiftRight(f - bits);
  return result;
}

inline constexpr std::size_t kPiBaseBits = std::size_t{1} << 15;

inline const PiFixedPoint& PiBase() {
  static const PiFixedPoint base = ComputePi(kPiBaseBits);
  return base;
}

// pi and 2/pi to at least `bits` bits. The largest size asked for so far is
// kept, rounded up to a multiple of 1/8 of its power of two, so that a growing
// series of requests recomputes only a few times.
inline std::shared_ptr<const PiFixedPoint> PiAtLeast(std::size_t bits) {
  static std::mutex mutex;
  static std::shared_ptr<const PiFixedPoint> large;
  const std::lock_guard<std::mutex> lock(mutex);
  if (!large || large->bits < bits) {
    std::size_t step = 1;
    while (step * 16 <= bits) step *= 2;
    std::size_t target = (bits + step - 1) / step * step;
    if (large) target = (std::max)(target, large->bits + large->bits / 8);
    large = std::make_shared<const PiFixedPoint>(
        ComputePi((std::max)(target, 2 * kPiBaseBits)));
  }
  return large;
}

// Calls fn(const PiFixedPoint&) with at least `bits` bits.
template <class Fn>
decltype(auto) WithPiBits(std::size_t bits, const Fn& fn) {
  if (bits <= kPiBaseBits) return fn(PiBase());
  const std::shared_ptr<const PiFixedPoint> large = PiAtLeast(bits);
  return fn(*large);
}

}  // namespace mpfix::detail

#endif  // MPFIX_PI_H_
