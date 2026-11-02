// Copyright 2026 The hwyMP Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Fixed-size IEEE 754 binary floating point.
//
// Float<P, E> has P bits of precision, including the hidden bit, and E exponent
// bits. Values are stored unpacked: a sign, an unbiased exponent and a
// significand of L = ceil(P / 64) limbs. The significand is left-aligned, so
// its top bit is always set, and the 64L - P spare bits below it are zero.
//
// Subnormals are also kept normalized, with an exponent below kEmin and extra
// zero bits at the bottom. Operands therefore need no subnormal handling;
// only the final rounding step and the conversion to the interchange format
// do. to_bits() and from_bits() give the exact IEEE interchange encoding.
//
// Operators and functions without a Context round to nearest-even and do not
// track flags, which is the fast path. Overloads taking a Context offer all five
// rounding modes and all five exception flags. Tininess is detected after
// rounding.
//
// NaN rules: an operation with NaN operands returns the first NaN operand,
// quieted. An invalid operation returns the positive default quiet NaN.
// fma(0, inf, qNaN) returns the NaN without raising invalid.

#ifndef HWYMP_FLOAT_H_
#define HWYMP_FLOAT_H_

#include <array>
#include <cmath>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <limits>

#include "hwymp/limb_ops.h"

namespace hwymp {

enum class Rounding : std::uint8_t {
  kNearestEven,
  kNearestAway,
  kTowardZero,
  kUp,    // Toward positive infinity.
  kDown,  // Toward negative infinity.
};

// IEEE 754 exception flags. They are sticky: operations only ever set them.
enum Flag : unsigned {
  kFlagInvalid = 1u << 0,
  kFlagDivByZero = 1u << 1,
  kFlagOverflow = 1u << 2,
  kFlagUnderflow = 1u << 3,
  kFlagInexact = 1u << 4,
};

struct Context {
  Rounding rounding = Rounding::kNearestEven;
  unsigned flags = 0;
};

enum class FloatClass : std::uint8_t {
  kNaN,
  kInfinite,
  kZero,
  kSubnormal,
  kNormal,
};

namespace detail {

// Rounds to nearest-even and discards flags, for the operators.
struct NearestEnv {
  static constexpr Rounding rounding() noexcept {
    return Rounding::kNearestEven;
  }
  static constexpr void raise(unsigned) noexcept {}
};

// Rounds in a given mode and discards flags, for tentative rounding steps.
struct SilentEnv {
  Rounding mode;
  constexpr Rounding rounding() const noexcept { return mode; }
  static constexpr void raise(unsigned) noexcept {}
};

struct ContextEnv {
  Context& context;
  constexpr Rounding rounding() const noexcept { return context.rounding; }
  constexpr void raise(unsigned flags) const noexcept {
    context.flags |= flags;
  }
};

// Like ContextEnv, but never raises inexact, for roundToIntegral.
struct NoInexactEnv {
  Context& context;
  constexpr Rounding rounding() const noexcept { return context.rounding; }
  constexpr void raise(unsigned flags) const noexcept {
    context.flags |= flags & ~unsigned{kFlagInexact};
  }
};

constexpr int LimbsFor(int bits) { return (bits + kLimbBits - 1) / kLimbBits; }

// Decides whether a discarded fraction rounds the magnitude up. `round_bit` is
// the most significant discarded bit and `sticky` is the OR of the rest.
HWYMP_ALWAYS_INLINE constexpr bool ShouldRoundUp(Rounding mode, bool negative,
                                                 bool is_odd, bool round_bit,
                                                 bool sticky) noexcept {
  switch (mode) {
    case Rounding::kNearestEven:
      return round_bit & (sticky | is_odd);
    case Rounding::kNearestAway:
      return round_bit;
    case Rounding::kTowardZero:
      return false;
    case Rounding::kUp:
      return !negative & (round_bit | sticky);
    case Rounding::kDown:
      return negative & (round_bit | sticky);
  }
  return false;
}

// Shifts the little-endian limb array left by 0 <= count < 64 * N bits.
template <std::size_t N>
constexpr void ShiftLeft(std::array<limb_t, N>& x, int count) noexcept {
  const int words = count / kLimbBits;
  const int bits = count % kLimbBits;
  for (int i = static_cast<int>(N) - 1; i >= 0; --i) {
    const int source = i - words;
    limb_t value = 0;
    if (source >= 0) {
      value = x[source];
      if (bits != 0) {
        value <<= bits;
        if (source >= 1) value |= x[source - 1] >> (kLimbBits - bits);
      }
    }
    x[i] = value;
  }
}

// Shifts the little-endian limb array right by 0 <= count < 64 * N bits.
template <std::size_t N>
constexpr void ShiftRight(std::array<limb_t, N>& x, int count) noexcept {
  const int words = count / kLimbBits;
  const int bits = count % kLimbBits;
  const int size = static_cast<int>(N);
  for (int i = 0; i < size; ++i) {
    const int source = i + words;
    limb_t value = 0;
    if (source < size) {
      value = x[source];
      if (bits != 0) {
        value >>= bits;
        if (source + 1 < size) value |= x[source + 1] << (kLimbBits - bits);
      }
    }
    x[i] = value;
  }
}

template <std::size_t N>
constexpr bool IsZero(const std::array<limb_t, N>& x) noexcept {
  limb_t combined = 0;
  for (const limb_t limb : x) combined |= limb;
  return combined == 0;
}

// Returns <0, 0 or >0 for x < y, x == y or x > y.
template <std::size_t N>
constexpr int CompareLimbs(const std::array<limb_t, N>& x,
                           const std::array<limb_t, N>& y) noexcept {
  for (std::size_t i = N; i-- > 0;) {
    if (x[i] != y[i]) return x[i] < y[i] ? -1 : 1;
  }
  return 0;
}

// Returns floor(sqrt(hi:lo)) and stores hi:lo minus its square in
// remainder_hi:remainder_lo. Requires hi >= 2^62, so that the root fills 64
// bits. The hardware estimate is only a starting point, so the host's rounding
// mode does not affect the result.
HWYMP_ALWAYS_INLINE limb_t IntegerSqrt(limb_t hi, limb_t lo,
                                       limb_t& remainder_hi,
                                       limb_t& remainder_lo) noexcept {
  // root^2 is about hi / 2^10, so root * 2^37 estimates sqrt(hi:lo).
  const double root =
      std::sqrt(static_cast<double>(static_cast<std::int64_t>(hi >> 10)));
  limb_t estimate = static_cast<limb_t>(
      static_cast<std::int64_t>(root * 34359738368.0));  // 2^35
  constexpr limb_t kLimit = limb_t{1} << 62;
  if (estimate >= kLimit) estimate = kLimit - 1;
  estimate <<= 2;

  // One Newton step takes the estimate from about 50 bits to nearly all 64.
  limb_t square_hi;
  limb_t square_lo = mul_wide(estimate, estimate, square_hi);
  carry_t borrow = 0;
  limb_t delta_lo = sbb(lo, square_lo, borrow);
  limb_t delta_hi = sbb(hi, square_hi, borrow);
  const bool below = borrow != 0;
  if (below) {
    carry_t negate = 0;
    delta_lo = sbb(0, delta_lo, negate);
    delta_hi = sbb(0, delta_hi, negate);
  }
  limb_t unused;
  const limb_t step = div_2by1(
      delta_hi >> 1, (delta_hi << 63) | (delta_lo >> 1), estimate, unused);
  limb_t result = below ? estimate - step : estimate + step;
  if (!below && result < estimate) result = ~limb_t{0};

  // A Newton step never lands below the true root, and truncating the step
  // loses less than one, so the result is at or above the floor.
  square_lo = mul_wide(result, result, square_hi);
  borrow = 0;
  limb_t rest_lo = sbb(lo, square_lo, borrow);
  limb_t rest_hi = sbb(hi, square_hi, borrow);
  bool negative = borrow != 0;
  while (negative) {
    // (s - 1)^2 = s^2 - s - (s - 1).
    carry_t carry = 0;
    rest_lo = adc(rest_lo, result, carry);
    rest_hi = adc(rest_hi, 0, carry);
    bool crossed = carry != 0;
    --result;
    carry = 0;
    rest_lo = adc(rest_lo, result, carry);
    rest_hi = adc(rest_hi, 0, carry);
    crossed = crossed || carry != 0;
    negative = !crossed;
  }
  remainder_hi = rest_hi;
  remainder_lo = rest_lo;
  return result;
}

}  // namespace detail

template <int P, int E>
class Float {
  static_assert(P >= 3, "NaN encodings need at least two trailing bits");
  static_assert(E >= 2 && E <= 29,
                "Exponent width must keep sums of exponents within int32_t");

 public:
  static constexpr int kPrecision = P;
  static constexpr int kExponentBits = E;
  static constexpr int kEmax = (1 << (E - 1)) - 1;
  static constexpr int kEmin = 1 - kEmax;
  static constexpr std::size_t kLimbs = detail::LimbsFor(P);
  static constexpr std::size_t kBitsLimbs = detail::LimbsFor(P + E);

  using Limbs = std::array<limb_t, kLimbs>;
  // The interchange encoding, least significant limb first.
  using Bits = std::array<limb_t, kBitsLimbs>;

  // Leaves the value uninitialized, like the built-in floating-point types.
  Float() noexcept = default;

  // --- Interchange format ---

  static Float from_bits(const Bits& bits) noexcept {
    const bool negative = ExtractBits(bits, P + E - 1, 1) != 0;
    const limb_t biased = ExtractBits(bits, P - 1, E);

    Limbs trailing;
    for (std::size_t i = 0; i < kLimbs; ++i) trailing[i] = bits[i];
    constexpr int top_bits = (P - 1) % kLimbBits;
    if constexpr (top_bits != 0) {
      trailing[(P - 1) / kLimbBits] &= (limb_t{1} << top_bits) - 1;
    }
    for (std::size_t i = (P - 1) / kLimbBits + (top_bits != 0 ? 1 : 0);
         i < kLimbs; ++i) {
      trailing[i] = 0;
    }
    detail::ShiftLeft(trailing, kSpareBits);

    const bool all_trailing_zero = detail::IsZero(trailing);
    if (biased == (limb_t{1} << E) - 1) {
      if (all_trailing_zero) return Infinity(negative);
      trailing[kLimbs - 1] |= kTopBit;
      return FromParts(negative, kExpNaN, trailing);
    }
    if (biased == 0) {
      if (all_trailing_zero) return Zero(negative);
      const int leading = CountLeadingZeros(trailing);
      detail::ShiftLeft(trailing, leading);
      return FromParts(negative, kEmin - leading, trailing);
    }
    trailing[kLimbs - 1] |= kTopBit;
    return FromParts(negative, static_cast<std::int32_t>(biased) - kEmax,
                     trailing);
  }

  Bits to_bits() const noexcept {
    Limbs trailing = limbs_;
    limb_t biased;
    if (exponent_ == kExpNaN || exponent_ == kExpInfinity) {
      biased = (limb_t{1} << E) - 1;
      trailing[kLimbs - 1] &= ~kTopBit;
    } else if (exponent_ == kExpZero) {
      biased = 0;
    } else if (exponent_ < kEmin) {
      biased = 0;
      detail::ShiftRight(trailing, kEmin - exponent_);
    } else {
      biased = static_cast<limb_t>(exponent_ + kEmax);
      trailing[kLimbs - 1] &= ~kTopBit;
    }
    detail::ShiftRight(trailing, kSpareBits);

    Bits bits{};
    for (std::size_t i = 0; i < kLimbs; ++i) bits[i] = trailing[i];
    InsertBits(bits, P - 1, E, biased);
    InsertBits(bits, P + E - 1, 1, negative_ ? 1 : 0);
    return bits;
  }

  // --- Special values ---

  static Float Zero(bool negative = false) noexcept {
    return FromParts(negative, kExpZero, Limbs{});
  }

  static Float Infinity(bool negative = false) noexcept {
    return FromParts(negative, kExpInfinity, TopBitOnly());
  }

  static Float QuietNaN(bool negative = false) noexcept {
    Limbs limbs = TopBitOnly();
    SetQuietBit(limbs);
    return FromParts(negative, kExpNaN, limbs);
  }

  // The signaling NaN whose payload has only its lowest bit set.
  static Float SignalingNaN(bool negative = false) noexcept {
    Limbs limbs = TopBitOnly();
    limbs[0] |= limb_t{1} << kSpareBits;
    return FromParts(negative, kExpNaN, limbs);
  }

  static Float MaxFinite(bool negative = false) noexcept {
    return FromParts(negative, kEmax, AllOnes());
  }

  static Float MinNormal(bool negative = false) noexcept {
    return FromParts(negative, kEmin, TopBitOnly());
  }

  static Float MinSubnormal(bool negative = false) noexcept {
    return FromParts(negative, kEmin - (P - 1), TopBitOnly());
  }

  // --- Classification ---

  bool signbit() const noexcept { return negative_ != 0; }
  bool is_nan() const noexcept { return exponent_ == kExpNaN; }
  bool is_inf() const noexcept { return exponent_ == kExpInfinity; }
  bool is_zero() const noexcept { return exponent_ == kExpZero; }
  bool is_finite() const noexcept { return exponent_ < kExpInfinity; }
  bool is_subnormal() const noexcept {
    return !is_zero() && exponent_ < kEmin;
  }
  bool is_normal() const noexcept {
    return exponent_ >= kEmin && exponent_ < kExpInfinity;
  }
  bool is_signaling() const noexcept {
    return is_nan() && (limbs_[kLimbs - 1] & kQuietBit) == 0;
  }

  FloatClass classify() const noexcept {
    if (is_nan()) return FloatClass::kNaN;
    if (is_inf()) return FloatClass::kInfinite;
    if (is_zero()) return FloatClass::kZero;
    return exponent_ < kEmin ? FloatClass::kSubnormal : FloatClass::kNormal;
  }

  // --- Sign operations: these never raise flags, even for signaling NaNs ---

  Float operator+() const noexcept { return *this; }

  Float operator-() const noexcept {
    Float result = *this;
    result.negative_ ^= 1;
    return result;
  }

  friend Float abs(const Float& x) noexcept {
    Float result = x;
    result.negative_ = 0;
    return result;
  }

  friend Float copysign(const Float& magnitude, const Float& sign) noexcept {
    Float result = magnitude;
    result.negative_ = sign.negative_;
    return result;
  }

  // --- Arithmetic ---

  friend Float operator+(const Float& a, const Float& b) noexcept {
    return AddImpl(a, b, b.negative_, detail::NearestEnv{});
  }
  friend Float operator-(const Float& a, const Float& b) noexcept {
    return AddImpl(a, b, b.negative_ ^ 1, detail::NearestEnv{});
  }
  friend Float operator*(const Float& a, const Float& b) noexcept {
    return MulImpl(a, b, detail::NearestEnv{});
  }

  Float& operator+=(const Float& other) noexcept {
    return *this = *this + other;
  }
  Float& operator-=(const Float& other) noexcept {
    return *this = *this - other;
  }
  friend Float operator/(const Float& a, const Float& b) noexcept {
    return DivImpl(a, b, detail::NearestEnv{});
  }

  Float& operator*=(const Float& other) noexcept {
    return *this = *this * other;
  }
  Float& operator/=(const Float& other) noexcept {
    return *this = *this / other;
  }

  friend Float add(const Float& a, const Float& b, Context& context) noexcept {
    return AddImpl(a, b, b.negative_, detail::ContextEnv{context});
  }
  friend Float sub(const Float& a, const Float& b, Context& context) noexcept {
    return AddImpl(a, b, b.negative_ ^ 1, detail::ContextEnv{context});
  }
  friend Float mul(const Float& a, const Float& b, Context& context) noexcept {
    return MulImpl(a, b, detail::ContextEnv{context});
  }
  friend Float div(const Float& a, const Float& b, Context& context) noexcept {
    return DivImpl(a, b, detail::ContextEnv{context});
  }

  friend Float sqrt(const Float& x) noexcept {
    return SqrtImpl(x, detail::NearestEnv{});
  }
  friend Float sqrt(const Float& x, Context& context) noexcept {
    return SqrtImpl(x, detail::ContextEnv{context});
  }

  // a * b + c with a single rounding.
  friend Float fma(const Float& a, const Float& b, const Float& c) noexcept {
    return FmaImpl(a, b, c, detail::NearestEnv{});
  }
  friend Float fma(const Float& a, const Float& b, const Float& c,
                   Context& context) noexcept {
    return FmaImpl(a, b, c, detail::ContextEnv{context});
  }

  // The IEEE remainder x - n * y, where n is x / y rounded to the nearest
  // integer, ties to even. It is always exact. The time grows linearly with the
  // difference between the operands' exponents.
  friend Float remainder(const Float& x, const Float& y) noexcept {
    return RemainderImpl(x, y, detail::NearestEnv{});
  }
  friend Float remainder(const Float& x, const Float& y,
                         Context& context) noexcept {
    return RemainderImpl(x, y, detail::ContextEnv{context});
  }

  // --- Rounding to integers ---

  // Rounds in context.rounding. Inexact results raise no flag.
  friend Float round_to_integral(const Float& x, Context& context) noexcept {
    return RoundToIntegralImpl(x, detail::NoInexactEnv{context});
  }
  // Rounds in context.rounding and raises inexact for inexact results.
  friend Float round_to_integral_exact(const Float& x,
                                       Context& context) noexcept {
    return RoundToIntegralImpl(x, detail::ContextEnv{context});
  }

  friend Float floor(const Float& x) noexcept {
    return RoundToIntegralImpl(x, detail::SilentEnv{Rounding::kDown});
  }
  friend Float ceil(const Float& x) noexcept {
    return RoundToIntegralImpl(x, detail::SilentEnv{Rounding::kUp});
  }
  friend Float trunc(const Float& x) noexcept {
    return RoundToIntegralImpl(x, detail::SilentEnv{Rounding::kTowardZero});
  }
  // Ties round away from zero, as in C.
  friend Float round(const Float& x) noexcept {
    return RoundToIntegralImpl(x, detail::SilentEnv{Rounding::kNearestAway});
  }

  // --- Minimum and maximum, as in IEEE 754-2019 ---
  //
  // -0 is below +0. minimum and maximum return NaN if an operand is NaN, while
  // the _number forms prefer a number over a NaN. A signaling NaN raises invalid
  // in both forms.

  friend Float minimum(const Float& a, const Float& b) noexcept {
    return MinMaxImpl<false, false>(a, b, detail::NearestEnv{});
  }
  friend Float minimum(const Float& a, const Float& b,
                       Context& context) noexcept {
    return MinMaxImpl<false, false>(a, b, detail::ContextEnv{context});
  }
  friend Float maximum(const Float& a, const Float& b) noexcept {
    return MinMaxImpl<true, false>(a, b, detail::NearestEnv{});
  }
  friend Float maximum(const Float& a, const Float& b,
                       Context& context) noexcept {
    return MinMaxImpl<true, false>(a, b, detail::ContextEnv{context});
  }
  friend Float minimum_number(const Float& a, const Float& b) noexcept {
    return MinMaxImpl<false, true>(a, b, detail::NearestEnv{});
  }
  friend Float minimum_number(const Float& a, const Float& b,
                              Context& context) noexcept {
    return MinMaxImpl<false, true>(a, b, detail::ContextEnv{context});
  }
  friend Float maximum_number(const Float& a, const Float& b) noexcept {
    return MinMaxImpl<true, true>(a, b, detail::NearestEnv{});
  }
  friend Float maximum_number(const Float& a, const Float& b,
                              Context& context) noexcept {
    return MinMaxImpl<true, true>(a, b, detail::ContextEnv{context});
  }

  // --- Neighbours and exponents ---

  // The next representable value above or below x. These never overflow or
  // underflow: next_up(MaxFinite) is infinity, and zeros of either sign step to
  // the smallest subnormal.
  friend Float next_up(const Float& x) noexcept {
    return NextImpl(x, true, detail::NearestEnv{});
  }
  friend Float next_up(const Float& x, Context& context) noexcept {
    return NextImpl(x, true, detail::ContextEnv{context});
  }
  friend Float next_down(const Float& x) noexcept {
    return NextImpl(x, false, detail::NearestEnv{});
  }
  friend Float next_down(const Float& x, Context& context) noexcept {
    return NextImpl(x, false, detail::ContextEnv{context});
  }

  // x * 2^n.
  friend Float scalbn(const Float& x, int n) noexcept {
    return ScalbnImpl(x, n, detail::NearestEnv{});
  }
  friend Float scalbn(const Float& x, int n, Context& context) noexcept {
    return ScalbnImpl(x, n, detail::ContextEnv{context});
  }

  static constexpr int kIlogbZero = std::numeric_limits<int>::min();
  static constexpr int kIlogbInfinity = std::numeric_limits<int>::max();
  static constexpr int kIlogbNaN = std::numeric_limits<int>::max();

  // The exponent of x as an integer, with subnormals treated as normalized.
  // Zero, infinity and NaN return the constants above and raise invalid.
  friend int ilogb(const Float& x) noexcept {
    return IlogbImpl(x, detail::NearestEnv{});
  }
  friend int ilogb(const Float& x, Context& context) noexcept {
    return IlogbImpl(x, detail::ContextEnv{context});
  }

  // The exponent of x as a Float. logb(0) is -infinity and raises divide by
  // zero, and logb(+-infinity) is +infinity.
  friend Float logb(const Float& x) noexcept {
    return LogbImpl(x, detail::NearestEnv{});
  }
  friend Float logb(const Float& x, Context& context) noexcept {
    return LogbImpl(x, detail::ContextEnv{context});
  }

  // --- Comparison ---
  //
  // The operators are quiet: they raise nothing. NaN compares unordered and
  // -0 equals +0.

  friend bool operator==(const Float& a, const Float& b) noexcept {
    return Order(a, b) == std::partial_ordering::equivalent;
  }

  friend std::partial_ordering operator<=>(const Float& a,
                                           const Float& b) noexcept {
    return Order(a, b);
  }

  // Raises invalid for signaling NaN operands only.
  friend std::partial_ordering compare_quiet(const Float& a, const Float& b,
                                             Context& context) noexcept {
    if (a.is_signaling() || b.is_signaling()) context.flags |= kFlagInvalid;
    return Order(a, b);
  }

  // Raises invalid for any NaN operand, as IEEE 754 requires for <, <=, > and
  // >= on quiet NaNs.
  friend std::partial_ordering compare_signaling(const Float& a, const Float& b,
                                                 Context& context) noexcept {
    if (a.is_nan() || b.is_nan()) context.flags |= kFlagInvalid;
    return Order(a, b);
  }

  // --- Raw access, for tests and for interoperating with other libraries ---

  std::int32_t raw_exponent() const noexcept { return exponent_; }
  const Limbs& raw_limbs() const noexcept { return limbs_; }

  // `limbs` must be normalized and `exponent` within the format's range.
  static Float FromRaw(bool negative, std::int32_t exponent,
                       const Limbs& limbs) noexcept {
    return FromParts(negative, exponent, limbs);
  }

 private:
  static constexpr int kSpareBits = static_cast<int>(kLimbs) * kLimbBits - P;
  static constexpr limb_t kTopBit = limb_t{1} << (kLimbBits - 1);
  static constexpr limb_t kQuietBit = kTopBit >> 1;

  // Exponent values that mark non-numbers. They lie far outside the range of
  // any regular value and keep the ordering zero < regular < infinity.
  static constexpr std::int32_t kExpZero =
      std::numeric_limits<std::int32_t>::min();
  static constexpr std::int32_t kExpInfinity =
      std::numeric_limits<std::int32_t>::max() - 1;
  static constexpr std::int32_t kExpNaN =
      std::numeric_limits<std::int32_t>::max();

  // A significand with the extra limb that holds the bits below the last limb:
  // `rest` has the highest of those bits at its top, and its lowest bit is
  // sticky, meaning that it also stands for everything shifted out below it.
  struct Extended {
    Limbs limbs;
    limb_t rest;
  };

  static Float FromParts(bool negative, std::int32_t exponent,
                         const Limbs& limbs) noexcept {
    Float result;
    result.limbs_ = limbs;
    result.exponent_ = exponent;
    result.negative_ = negative ? 1 : 0;
    return result;
  }

  static Limbs TopBitOnly() noexcept {
    Limbs limbs{};
    limbs[kLimbs - 1] = kTopBit;
    return limbs;
  }

  static Limbs AllOnes() noexcept {
    Limbs limbs;
    limbs.fill(~limb_t{0});
    if constexpr (kSpareBits != 0) limbs[0] &= ~LowMask(kSpareBits);
    return limbs;
  }

  static void SetQuietBit(Limbs& limbs) noexcept {
    limbs[kLimbs - 1] |= kQuietBit;
  }

  // Requires 0 < count < 64.
  static constexpr limb_t LowMask(int count) noexcept {
    return (limb_t{1} << count) - 1;
  }

  // Bit fields are at most E + 1 bits wide and span at most two limbs.
  static limb_t ExtractBits(const Bits& bits, int position,
                            int width) noexcept {
    const std::size_t word = static_cast<std::size_t>(position) / kLimbBits;
    const int offset = position % kLimbBits;
    limb_t value = bits[word] >> offset;
    if (offset + width > kLimbBits) {
      value |= bits[word + 1] << (kLimbBits - offset);
    }
    return value & LowMask(width);
  }

  static void InsertBits(Bits& bits, int position, int width,
                         limb_t value) noexcept {
    const std::size_t word = static_cast<std::size_t>(position) / kLimbBits;
    const int offset = position % kLimbBits;
    bits[word] |= value << offset;
    if (offset + width > kLimbBits) {
      bits[word + 1] |= value >> (kLimbBits - offset);
    }
  }

  static int CountLeadingZeros(const Limbs& limbs) noexcept {
    for (std::size_t i = kLimbs; i-- > 0;) {
      if (limbs[i] != 0) {
        return static_cast<int>(kLimbs - 1 - i) * kLimbBits +
               countl_zero(limbs[i]);
      }
    }
    return static_cast<int>(kLimbs) * kLimbBits;
  }

  static bool IsRegular(const Float& x) noexcept {
    return x.exponent_ != kExpZero && x.exponent_ < kExpInfinity;
  }

  // --- Rounding ---

  // Rounds to P bits. Returns true when the significand overflowed, in which
  // case it is left as 2^(64L - 1) and the exponent must grow by one.
  template <class Env>
  HWYMP_ALWAYS_INLINE static bool RoundToPrecision(bool negative,
                                                   Extended& value,
                                                   const Env& env) noexcept {
    bool is_odd;
    bool round_bit;
    bool sticky;
    if constexpr (kSpareBits > 0) {
      const limb_t half = limb_t{1} << (kSpareBits - 1);
      const limb_t low = value.limbs[0] & LowMask(kSpareBits);
      is_odd = ((value.limbs[0] >> kSpareBits) & 1) != 0;
      round_bit = (low & half) != 0;
      sticky = ((low & (half - 1)) | value.rest) != 0;
      value.limbs[0] &= ~LowMask(kSpareBits);
    } else {
      is_odd = (value.limbs[0] & 1) != 0;
      round_bit = (value.rest >> (kLimbBits - 1)) != 0;
      sticky = (value.rest << 1) != 0;
    }
    if (round_bit || sticky) env.raise(kFlagInexact);
    // Always adding the increment avoids a branch on what is a coin flip.
    const limb_t increment = detail::ShouldRoundUp(env.rounding(), negative,
                                                   is_odd, round_bit, sticky)
                                 ? limb_t{1} << kSpareBits
                                 : 0;
    carry_t carry = 0;
    value.limbs[0] = adc(value.limbs[0], increment, carry);
    for (std::size_t i = 1; i < kLimbs; ++i) {
      value.limbs[i] = adc(value.limbs[i], 0, carry);
    }
    if (carry != 0) value.limbs[kLimbs - 1] = kTopBit;
    return carry != 0;
  }

  // Produces the final result from a normalized value of any exponent.
  template <class Env>
  HWYMP_ALWAYS_INLINE static Float Finish(bool negative, std::int32_t exponent,
                                          Extended& value,
                                          const Env& env) noexcept {
    if (exponent < kEmin || exponent > kEmax) [[unlikely]] {
      if (exponent > kEmax) return Overflow(negative, env);
      return FinishTiny(negative, exponent, value, env);
    }
    if (RoundToPrecision(negative, value, env)) {
      ++exponent;
      if (exponent > kEmax) [[unlikely]] return Overflow(negative, env);
    }
    return FromParts(negative, exponent, value.limbs);
  }

  template <class Env>
  static Float Overflow(bool negative, const Env& env) noexcept {
    env.raise(kFlagOverflow | kFlagInexact);
    const Rounding mode = env.rounding();
    const bool to_infinity = mode == Rounding::kNearestEven ||
                             mode == Rounding::kNearestAway ||
                             (mode == Rounding::kUp && !negative) ||
                             (mode == Rounding::kDown && negative);
    return to_infinity ? Infinity(negative) : MaxFinite(negative);
  }

  // Rounds a value below 2^kEmin onto the subnormal grid.
  template <class Env>
  static Float FinishTiny(bool negative, std::int32_t exponent, Extended& value,
                          const Env& env) noexcept {
    constexpr int kTotalBits = static_cast<int>(kLimbs) * kLimbBits;

    // A value just below 2^kEmin is not tiny if rounding it with an unbounded
    // exponent reaches 2^kEmin.
    bool tiny = true;
    if (exponent == kEmin - 1) {
      Extended copy = value;
      tiny = !RoundToPrecision(negative, copy,
                               detail::SilentEnv{env.rounding()});
    }

    const std::int64_t discard =
        std::int64_t{kSpareBits} + (std::int64_t{kEmin} - exponent);
    Limbs& limbs = value.limbs;

    bool is_odd = false;
    bool round_bit = false;
    bool sticky = value.rest != 0;
    if (discard > kTotalBits) {
      sticky = true;  // The value is nonzero and entirely below the round bit.
    } else {
      InspectLowBits(limbs, static_cast<int>(discard), is_odd, round_bit,
                     sticky);
    }

    if (round_bit || sticky) {
      env.raise(kFlagInexact);
      if (tiny) env.raise(kFlagUnderflow);
    }
    const bool round_up = detail::ShouldRoundUp(env.rounding(), negative,
                                                is_odd, round_bit, sticky);

    if (discard >= kTotalBits) {
      return round_up ? MinSubnormal(negative) : Zero(negative);
    }

    const int count = static_cast<int>(discard);
    detail::ShiftRight(limbs, count);
    detail::ShiftLeft(limbs, count);
    if (round_up && IncrementAtBit(limbs, count)) ++exponent;
    return FromParts(negative, exponent, limbs);
  }

  // Finds the rounding inputs for discarding the low 0 < count <= 64L bits.
  // Bits already known to be lost are passed in through `sticky`.
  static void InspectLowBits(const Limbs& limbs, int count, bool& is_odd,
                             bool& round_bit, bool& sticky) noexcept {
    constexpr int kTotalBits = static_cast<int>(kLimbs) * kLimbBits;
    const auto bit_at = [&limbs](int position) {
      return ((limbs[static_cast<std::size_t>(position) / kLimbBits] >>
               (position % kLimbBits)) &
              1) != 0;
    };
    round_bit = bit_at(count - 1);
    if (count < kTotalBits) is_odd = bit_at(count);
    Limbs below = limbs;
    detail::ShiftLeft(below, kTotalBits - count);
    // `below` now starts with the round bit; drop it to keep the rest.
    below[kLimbs - 1] &= ~kTopBit;
    sticky = sticky || !detail::IsZero(below);
  }

  // Adds 2^bit. On a carry out of the top, leaves the significand as 2^(64L-1)
  // and returns true: the exponent must then grow by one.
  static bool IncrementAtBit(Limbs& limbs, int bit) noexcept {
    const std::size_t word = static_cast<std::size_t>(bit) / kLimbBits;
    carry_t carry = 0;
    limbs[word] = adc(limbs[word], limb_t{1} << (bit % kLimbBits), carry);
    for (std::size_t i = word + 1; i < kLimbs; ++i) {
      limbs[i] = adc(limbs[i], 0, carry);
    }
    if (carry != 0) limbs[kLimbs - 1] = kTopBit;
    return carry != 0;
  }

  // --- NaN handling ---

  template <class Env>
  static Float PropagateNaN(const Float& a, const Float& b,
                            const Env& env) noexcept {
    if (a.is_signaling() || b.is_signaling()) env.raise(kFlagInvalid);
    Float result = a.is_nan() ? a : b;
    SetQuietBit(result.limbs_);
    return result;
  }

  template <class Env>
  static Float InvalidOperation(const Env& env) noexcept {
    env.raise(kFlagInvalid);
    return QuietNaN();
  }

  // --- Addition and subtraction ---

  // Shifts `small` right by distance >= 0, jamming lost bits into `rest`.
  HWYMP_ALWAYS_INLINE static Extended AlignRight(const Limbs& small,
                                                 std::int32_t distance) noexcept {
    constexpr int kWords = static_cast<int>(kLimbs);
    if (distance >= (kWords + 1) * kLimbBits) [[unlikely]] {
      return Extended{Limbs{}, 1};
    }
    const int words = distance / kLimbBits;
    const int bits = distance % kLimbBits;

    if (words == 0) {
      // By far the most common case, with no limbs to move.
      const int complement = kLimbBits - 1 - bits;
      Extended result;
      result.rest = (small[0] << 1) << complement;
      for (std::size_t i = 0; i + 1 < kLimbs; ++i) {
        result.limbs[i] =
            (small[i] >> bits) | ((small[i + 1] << 1) << complement);
      }
      result.limbs[kLimbs - 1] = small[kLimbs - 1] >> bits;
      return result;
    }

    // Element 0 is `rest`; elements 1 to L are the shifted limbs.
    std::array<limb_t, kLimbs + 1> wide{};
    for (int i = 0; i + words <= kWords; ++i) {
      const int source = i + words - 1;
      wide[static_cast<std::size_t>(i)] =
          source >= 0 ? small[static_cast<std::size_t>(source)] : 0;
    }
    limb_t lost = 0;
    for (int i = 0; i + 1 < words; ++i) {
      lost |= small[static_cast<std::size_t>(i)];
    }
    // Shifting in two steps keeps a bit count of zero well defined, so no
    // branch is needed for it.
    const int complement = kLimbBits - 1 - bits;
    lost |= (wide[0] << 1) << complement;
    Extended result;
    result.rest = ((wide[0] >> bits) | ((wide[1] << 1) << complement)) |
                  (lost != 0 ? 1 : 0);
    for (std::size_t i = 0; i < kLimbs; ++i) {
      const limb_t above = i + 1 < kLimbs ? wide[i + 2] : 0;
      result.limbs[i] = (wide[i + 1] >> bits) | ((above << 1) << complement);
    }
    return result;
  }

  template <class Env>
  static Float AddSpecial(const Float& a, const Float& b, std::uint32_t b_sign,
                          const Env& env) noexcept {
    if (a.is_nan() || b.is_nan()) return PropagateNaN(a, b, env);
    if (a.is_inf()) {
      if (b.is_inf() && a.negative_ != b_sign) return InvalidOperation(env);
      return a;
    }
    if (b.is_inf()) return Infinity(b_sign != 0);
    if (a.is_zero() && b.is_zero()) {
      if (a.negative_ == b_sign) return a;
      return Zero(env.rounding() == Rounding::kDown);
    }
    if (a.is_zero()) {
      Float result = b;
      result.negative_ = b_sign;
      return result;
    }
    return a;
  }

  // Computes a + b, using b_sign in place of b's sign so that subtraction
  // shares this code.
  template <class Env>
  HWYMP_ALWAYS_INLINE static Float AddImpl(const Float& a, const Float& b,
                                           std::uint32_t b_sign,
                                           const Env& env) noexcept {
    if (!IsRegular(a) || !IsRegular(b)) [[unlikely]] {
      return AddSpecial(a, b, b_sign, env);
    }

    const bool swap = a.exponent_ < b.exponent_;
    const Limbs* big = swap ? &b.limbs_ : &a.limbs_;
    const Limbs* small = swap ? &a.limbs_ : &b.limbs_;
    std::int32_t exponent = swap ? b.exponent_ : a.exponent_;
    std::int32_t distance = swap ? b.exponent_ - a.exponent_
                                 : a.exponent_ - b.exponent_;
    bool negative = (swap ? b_sign : a.negative_) != 0;

    if (a.negative_ == b_sign) {
      Extended sum = AlignRight(*small, distance);
      carry_t carry = 0;
      for (std::size_t i = 0; i < kLimbs; ++i) {
        sum.limbs[i] = adc((*big)[i], sum.limbs[i], carry);
      }
      // A carry happens about half the time, so shift by it without a branch.
      const int shift = carry;
      const int complement = kLimbBits - 1 - shift;
      sum.rest = ((sum.limbs[0] << 1) << complement) | (sum.rest >> shift) |
                 (sum.rest & 1);
      for (std::size_t i = 0; i + 1 < kLimbs; ++i) {
        sum.limbs[i] = (sum.limbs[i] >> shift) |
                       ((sum.limbs[i + 1] << 1) << complement);
      }
      sum.limbs[kLimbs - 1] = (sum.limbs[kLimbs - 1] >> shift) |
                              ((limb_t{carry} << 1) << complement);
      exponent += shift;
      return Finish(negative, exponent, sum, env);
    }

    if (distance == 0) {
      const int order = detail::CompareLimbs(*big, *small);
      if (order == 0) return Zero(env.rounding() == Rounding::kDown);
      if (order < 0) {
        const Limbs* temporary = big;
        big = small;
        small = temporary;
        negative = !negative;
      }
    }

    Extended difference = AlignRight(*small, distance);
    carry_t borrow = 0;
    difference.rest = sbb(0, difference.rest, borrow);
    for (std::size_t i = 0; i < kLimbs; ++i) {
      difference.limbs[i] = sbb((*big)[i], difference.limbs[i], borrow);
    }
    if (difference.limbs[kLimbs - 1] >> (kLimbBits - 2) == 0) [[unlikely]] {
      exponent -= NormalizeDifference(difference);
    } else {
      // At most one bit of cancellation, which is a coin flip.
      const int shift = static_cast<int>(difference.limbs[kLimbs - 1] >>
                                         (kLimbBits - 1)) ^
                        1;
      ShiftLeftByBit(difference, shift);
      exponent -= shift;
    }
    return Finish(negative, exponent, difference, env);
  }

  // Shifts left by `shift`, which is 0 or 1, without branching on it.
  HWYMP_ALWAYS_INLINE static void ShiftLeftByBit(Extended& value,
                                                 int shift) noexcept {
    const int complement = kLimbBits - 1 - shift;
    for (std::size_t i = kLimbs - 1; i > 0; --i) {
      value.limbs[i] = (value.limbs[i] << shift) |
                       ((value.limbs[i - 1] >> 1) >> complement);
    }
    value.limbs[0] =
        (value.limbs[0] << shift) | ((value.rest >> 1) >> complement);
    value.rest <<= shift;
  }

  // Shifts a nonzero difference left until its top bit is set, returning the
  // distance. Only the exact case with distance <= 1 shifts by more than one,
  // and there `rest` holds no jammed bits.
  static int NormalizeDifference(Extended& value) noexcept {
    std::array<limb_t, kLimbs + 1> wide;
    wide[0] = value.rest;
    for (std::size_t i = 0; i < kLimbs; ++i) wide[i + 1] = value.limbs[i];

    int leading = 0;
    for (std::size_t i = kLimbs + 1; i-- > 0;) {
      if (wide[i] != 0) {
        leading += countl_zero(wide[i]);
        break;
      }
      leading += kLimbBits;
    }
    detail::ShiftLeft(wide, leading);
    value.rest = wide[0];
    for (std::size_t i = 0; i < kLimbs; ++i) value.limbs[i] = wide[i + 1];
    return leading;
  }

  // --- Multiplication ---

  template <class Env>
  static Float MulSpecial(const Float& a, const Float& b, bool negative,
                          const Env& env) noexcept {
    if (a.is_nan() || b.is_nan()) return PropagateNaN(a, b, env);
    if (a.is_inf() || b.is_inf()) {
      if (a.is_zero() || b.is_zero()) return InvalidOperation(env);
      return Infinity(negative);
    }
    return Zero(negative);
  }

  template <class Env>
  HWYMP_ALWAYS_INLINE static Float MulImpl(const Float& a, const Float& b,
                                           const Env& env) noexcept {
    const bool negative = (a.negative_ ^ b.negative_) != 0;
    if (!IsRegular(a) || !IsRegular(b)) [[unlikely]] {
      return MulSpecial(a, b, negative, env);
    }

    const std::array<limb_t, 2 * kLimbs> product =
        MultiplySignificands(a.limbs_, b.limbs_);

    Extended value;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      value.limbs[i] = product[i + kLimbs];
    }
    value.rest = product[kLimbs - 1];
    limb_t lost = 0;
    for (std::size_t i = 0; i + 1 < kLimbs; ++i) lost |= product[i];

    // Both significands are in [1, 2), so the product is in [1, 4). Either
    // case is equally likely, so this avoids a branch.
    const int shift =
        static_cast<int>(value.limbs[kLimbs - 1] >> (kLimbBits - 1)) ^ 1;
    ShiftLeftByBit(value, shift);
    value.rest |= lost != 0 ? 1 : 0;
    const std::int32_t exponent = a.exponent_ + b.exponent_ + 1 - shift;
    return Finish(negative, exponent, value, env);
  }

  HWYMP_ALWAYS_INLINE static std::array<limb_t, 2 * kLimbs>
  MultiplySignificands(const Limbs& a, const Limbs& b) noexcept {
    std::array<limb_t, 2 * kLimbs> product;
    if constexpr (kLimbs > 4) {
      // Row by row, with one unbroken carry chain for each of the low and high
      // halves. Compilers unroll these rectangular loops at this size.
      product.fill(0);
      for (std::size_t i = 0; i < kLimbs; ++i) {
        Limbs lows;
        Limbs highs;
        for (std::size_t j = 0; j < kLimbs; ++j) {
          lows[j] = mul_wide(a[j], b[i], highs[j]);
        }
        carry_t carry = 0;
        for (std::size_t j = 0; j < kLimbs; ++j) {
          product[i + j] = adc(product[i + j], lows[j], carry);
        }
        product[i + kLimbs] = carry;
        carry = 0;
        for (std::size_t j = 0; j < kLimbs; ++j) {
          product[i + j + 1] = adc(product[i + j + 1], highs[j], carry);
        }
      }
      return product;
    }
    // Product scanning: each column is summed into a three-limb accumulator,
    // which needs fewer carry operations than adding row by row.
    limb_t column_low = 0;
    limb_t column_middle = 0;
    limb_t column_high = 0;
    for (std::size_t k = 0; k + 1 < 2 * kLimbs; ++k) {
      const std::size_t first = k < kLimbs ? 0 : k - kLimbs + 1;
      const std::size_t last = k < kLimbs ? k : kLimbs - 1;
      for (std::size_t i = first; i <= last; ++i) {
        limb_t high;
        const limb_t low = mul_wide(a[i], b[k - i], high);
        carry_t carry = 0;
        column_low = adc(column_low, low, carry);
        column_middle = adc(column_middle, high, carry);
        column_high += carry;
      }
      product[k] = column_low;
      column_low = column_middle;
      column_middle = column_high;
      column_high = 0;
    }
    product[2 * kLimbs - 1] = column_low;
    return product;
  }

  // --- Division ---

  // Sets r to r * 2^64 mod divisor and returns the quotient limb. Requires
  // r < divisor, a normalized divisor, and `inverse` from reciprocal_3by2 on the
  // divisor's top two limbs.
  HWYMP_ALWAYS_INLINE static limb_t DivideStep(Limbs& r, const Limbs& divisor,
                                               limb_t inverse) noexcept {
    const limb_t d1 = divisor[kLimbs - 1];
    limb_t d0 = 0;
    limb_t u1 = 0;
    limb_t u0 = 0;
    if constexpr (kLimbs >= 2) {
      d0 = divisor[kLimbs - 2];
      u1 = r[kLimbs - 2];
    }
    if constexpr (kLimbs >= 3) u0 = r[kLimbs - 3];
    const limb_t u2 = r[kLimbs - 1];

    limb_t q;
    limb_t r1;
    limb_t r0;
    int over = 0;  // The remainder's bits above the L limbs of `r`.
    if (u2 == d1 && u1 == d0) [[unlikely]] {
      // The top limbs alone would give a quotient limb of 2^64.
      q = ~limb_t{0};
      carry_t carry = 0;
      r0 = adc(d0, u0, carry);
      r1 = adc(d1, 0, carry);
      over = carry;
    } else {
      q = div_3by2(u2, u1, u0, d1, d0, inverse, r1, r0);
    }

    if constexpr (kLimbs >= 3) {
      for (std::size_t i = kLimbs - 3; i > 0; --i) r[i] = r[i - 1];
      r[0] = 0;
    }
    r[kLimbs - 1] = r1;
    if constexpr (kLimbs >= 2) r[kLimbs - 2] = r0;

    if constexpr (kLimbs >= 3) {
      // Subtract q times the divisor limbs that the estimate did not see.
      limb_t incoming = 0;
      carry_t borrow = 0;
      for (std::size_t i = 0; i + 2 < kLimbs; ++i) {
        limb_t high;
        limb_t low = mul_wide(q, divisor[i], high);
        carry_t carry = 0;
        low = adc(low, incoming, carry);
        incoming = high + carry;
        r[i] = sbb(r[i], low, borrow);
      }
      r[kLimbs - 2] = sbb(r[kLimbs - 2], incoming, borrow);
      r[kLimbs - 1] = sbb(r[kLimbs - 1], 0, borrow);
      over -= borrow;
      while (over < 0) [[unlikely]] {
        --q;
        carry_t carry = 0;
        for (std::size_t i = 0; i < kLimbs; ++i) {
          r[i] = adc(r[i], divisor[i], carry);
        }
        over += carry;
      }
    }
    return q;
  }

  template <class Env>
  static Float DivSpecial(const Float& a, const Float& b, bool negative,
                          const Env& env) noexcept {
    if (a.is_nan() || b.is_nan()) return PropagateNaN(a, b, env);
    if (a.is_inf()) {
      if (b.is_inf()) return InvalidOperation(env);
      return Infinity(negative);
    }
    if (b.is_inf()) return Zero(negative);
    if (b.is_zero()) {
      if (a.is_zero()) return InvalidOperation(env);
      env.raise(kFlagDivByZero);
      return Infinity(negative);
    }
    return Zero(negative);
  }

  template <class Env>
  HWYMP_ALWAYS_INLINE static Float DivImpl(const Float& a, const Float& b,
                                           const Env& env) noexcept {
    const bool negative = (a.negative_ ^ b.negative_) != 0;
    if (!IsRegular(a) || !IsRegular(b)) [[unlikely]] {
      return DivSpecial(a, b, negative, env);
    }

    // The integer part of a's significand over b's is 0 or 1. Which one is a
    // coin flip, so both of the steps that depend on it avoid a branch.
    Limbs r;
    carry_t borrow = 0;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      r[i] = sbb(a.limbs_[i], b.limbs_[i], borrow);
    }
    const limb_t mask = limb_t{0} - static_cast<limb_t>(borrow);
    carry_t carry = 0;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      r[i] = adc(r[i], b.limbs_[i] & mask, carry);
    }
    const int integer_part = borrow ^ 1;

    // With a spare bit, L fraction limbs reach the round bit even when the
    // integer part is 0, so the last and slowest step can be skipped.
    constexpr std::size_t kFractionLimbs = kSpareBits >= 1 ? kLimbs : kLimbs + 1;
    std::array<limb_t, kFractionLimbs + 1> quotient;
    quotient[kFractionLimbs] = static_cast<limb_t>(integer_part);
    const limb_t d1 = b.limbs_[kLimbs - 1];
    limb_t d0 = 0;
    if constexpr (kLimbs >= 2) d0 = b.limbs_[kLimbs - 2];
    const limb_t inverse = reciprocal_3by2(d1, d0);
    for (std::size_t i = kFractionLimbs; i-- > 0;) {
      quotient[i] = DivideStep(r, b.limbs_, inverse);
    }

    // Normalize: the integer bit, if set, moves into the top of the limbs.
    const int complement = kLimbBits - 1 - integer_part;
    constexpr std::size_t kOffset = kFractionLimbs - kLimbs;
    Extended value;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      value.limbs[i] = (quotient[i + kOffset] >> integer_part) |
                       ((quotient[i + kOffset + 1] << 1) << complement);
    }
    if constexpr (kOffset == 0) {
      value.rest = (quotient[0] << 1) << complement;
    } else {
      // A bit shifted out here needs no sticky: an exact quotient has at most
      // P significant bits, and any other sets the sticky bit below.
      value.rest = (quotient[0] >> integer_part) |
                   ((quotient[1] << 1) << complement);
    }
    value.rest |= detail::IsZero(r) ? 0 : 1;
    const std::int32_t exponent = a.exponent_ - b.exponent_ - 1 + integer_part;
    return Finish(negative, exponent, value, env);
  }

  // --- Square root ---

  template <class Env>
  static Float SqrtImpl(const Float& x, const Env& env) noexcept {
    if (!IsRegular(x) || x.negative_ != 0) [[unlikely]] {
      if (x.is_nan()) return PropagateNaN(x, x, env);
      if (x.is_zero()) return x;
      if (x.negative_ != 0) return InvalidOperation(env);
      return x;  // Positive infinity.
    }

    // Limb by limb, the way long division works. With S the root so far, N the
    // radicand and R = N - S^2, the next limb q of the root is the largest one
    // that keeps R - 2 S q - q^2 non-negative. R / 2S bounds it from above, and
    // the remainder stays exact, so rounding needs no final squaring.
    //
    // The radicand has 2L limbs and its top two bits are not both zero, so the
    // root fills L limbs. An odd exponent becomes even by doubling the
    // significand.
    constexpr std::size_t kWide = 2 * kLimbs;
    const int shift = (x.exponent_ & 1) ^ 1;
    const int complement = kLimbBits - 1 - shift;
    std::array<limb_t, kWide> rem{};
    rem[kLimbs - 1] = (x.limbs_[0] << 1) << complement;
    for (std::size_t i = 0; i + 1 < kLimbs; ++i) {
      rem[kLimbs + i] = (x.limbs_[i] >> shift) |
                        ((x.limbs_[i + 1] << 1) << complement);
    }
    rem[kWide - 1] = x.limbs_[kLimbs - 1] >> shift;

    Limbs root{};
    root[kLimbs - 1] = detail::IntegerSqrt(rem[kWide - 1], rem[kWide - 2],
                                           rem[kWide - 1], rem[kWide - 2]);
    if constexpr (kLimbs >= 2) {
      // The second limb of the root is not known yet, so it counts as zero.
      limb_t inverse = reciprocal_3by2(root[kLimbs - 1], 0);
      for (std::size_t j = kLimbs - 1; j-- > 0;) {
        if (j + 3 == kLimbs) {
          inverse = reciprocal_3by2(root[kLimbs - 1], root[kLimbs - 2]);
        }
        NextRootLimb(root, rem, j, inverse);
      }
    }

    // There is never a tie: R = S would need a root ending in exactly one half,
    // and its square is not an integer at this scale.
    bool above_half = false;
    for (std::size_t i = kLimbs; i < kWide; ++i) {
      above_half = above_half || rem[i] != 0;
    }
    Limbs low;
    for (std::size_t i = 0; i < kLimbs; ++i) low[i] = rem[i];
    if (!above_half) above_half = detail::CompareLimbs(low, root) > 0;
    Extended value;
    value.limbs = root;
    value.rest = (above_half ? kTopBit : 0) | (detail::IsZero(rem) ? 0 : 1);
    return Finish(false, x.exponent_ >> 1, value, env);
  }

  // Finds limb j of the root and updates the remainder. Limbs above j are final,
  // and `inverse` matches the root's top two limbs.
  HWYMP_ALWAYS_INLINE static void NextRootLimb(
      Limbs& root, std::array<limb_t, 2 * kLimbs>& rem, std::size_t j,
      limb_t inverse) noexcept {
    constexpr std::size_t kWide = 2 * kLimbs;
    const auto halved = [&rem](std::size_t i) {
      const limb_t above = i + 1 < kWide ? rem[i + 1] : 0;
      return (rem[i] >> 1) | (above << (kLimbBits - 1));
    };
    const limb_t d1 = root[kLimbs - 1];
    const limb_t d0 = root[kLimbs - 2];
    const limb_t u2 = halved(kLimbs + j);
    const limb_t u1 = halved(kLimbs + j - 1);
    const limb_t u0 = halved(kLimbs + j - 2);

    // The estimate can only be too high.
    limb_t q = ~limb_t{0};
    if (u2 < d1 || (u2 == d1 && u1 < d0)) {
      limb_t unused_high;
      limb_t unused_low;
      q = div_3by2(u2, u1, u0, d1, d0, inverse, unused_high, unused_low);
    }

    // Subtract (2 root + q 2^(64j)) q 2^(64j), which starts at limb 2j.
    std::array<limb_t, kWide> subtrahend{};
    subtrahend[2 * j] = mul_wide(q, q, subtrahend[2 * j + 1]);
    limb_t incoming = 0;
    limb_t top_of_previous = 0;
    carry_t add_carry = 0;
    for (std::size_t i = j + 1; i < kLimbs; ++i) {
      limb_t high;
      limb_t low = mul_wide(root[i], q, high);
      carry_t carry = 0;
      low = adc(low, incoming, carry);
      incoming = high + carry;
      const limb_t doubled = (low << 1) | top_of_previous;
      top_of_previous = low >> (kLimbBits - 1);
      subtrahend[i + j] = adc(subtrahend[i + j], doubled, add_carry);
    }
    const limb_t doubled = (incoming << 1) | top_of_previous;
    subtrahend[kLimbs + j] = adc(subtrahend[kLimbs + j], doubled, add_carry);
    subtrahend[kLimbs + j + 1] =
        adc(subtrahend[kLimbs + j + 1], incoming >> (kLimbBits - 1), add_carry);
    for (std::size_t i = kLimbs + j + 2; i < kWide; ++i) {
      subtrahend[i] = adc(subtrahend[i], 0, add_carry);
    }

    carry_t borrow = 0;
    for (std::size_t i = 2 * j; i < kWide; ++i) {
      rem[i] = sbb(rem[i], subtrahend[i], borrow);
    }
    root[j] = q;

    // While negative, step q down: (s - e)^2 = s^2 - (2s - e)e for e = 2^(64j).
    int negative = borrow;
    while (negative > 0) [[unlikely]] {
      limb_t previous_top = 0;
      carry_t carry = 0;
      for (std::size_t i = j; i < kLimbs; ++i) {
        const limb_t doubled_limb = (root[i] << 1) | previous_top;
        previous_top = root[i] >> (kLimbBits - 1);
        rem[i + j] = adc(rem[i + j], doubled_limb, carry);
      }
      rem[kLimbs + j] = adc(rem[kLimbs + j], previous_top, carry);
      for (std::size_t i = kLimbs + j + 1; i < kWide; ++i) {
        rem[i] = adc(rem[i], 0, carry);
      }
      carry_t step_borrow = 0;
      rem[2 * j] = sbb(rem[2 * j], 1, step_borrow);
      for (std::size_t i = 2 * j + 1; i < kWide; ++i) {
        rem[i] = sbb(rem[i], 0, step_borrow);
      }
      --root[j];
      negative += step_borrow;
      negative -= carry;
    }
  }

  // --- Fused multiply-add ---

  template <class Env>
  static Float FmaSpecial(const Float& a, const Float& b, const Float& c,
                          bool product_negative, const Env& env) noexcept {
    if (a.is_nan() || b.is_nan() || c.is_nan()) {
      if (a.is_signaling() || b.is_signaling() || c.is_signaling()) {
        env.raise(kFlagInvalid);
      }
      Float result = a.is_nan() ? a : (b.is_nan() ? b : c);
      SetQuietBit(result.limbs_);
      return result;
    }
    if (a.is_inf() || b.is_inf()) {
      if (a.is_zero() || b.is_zero()) return InvalidOperation(env);
      if (c.is_inf() && (c.negative_ != 0) != product_negative) {
        return InvalidOperation(env);
      }
      return Infinity(product_negative);
    }
    if (c.is_inf()) return c;
    if (a.is_zero() || b.is_zero()) {
      if (!c.is_zero()) return c;
      if ((c.negative_ != 0) == product_negative) return c;
      return Zero(env.rounding() == Rounding::kDown);
    }
    return MulImpl(a, b, env);  // c is zero, and the product is not.
  }

  template <class Env>
  static Float FmaImpl(const Float& a, const Float& b, const Float& c,
                       const Env& env) noexcept {
    const bool product_negative = (a.negative_ ^ b.negative_) != 0;
    if (!IsRegular(a) || !IsRegular(b) || !IsRegular(c)) [[unlikely]] {
      return FmaSpecial(a, b, c, product_negative, env);
    }

    // Addition at double width, on the exact product. Limb 0 is a guard limb
    // whose lowest bit is sticky, as in `Extended`.
    constexpr std::size_t kWide = 2 * kLimbs + 1;
    using Wide = std::array<limb_t, kWide>;

    const std::array<limb_t, 2 * kLimbs> exact =
        MultiplySignificands(a.limbs_, b.limbs_);
    Wide big{};
    Wide small{};
    for (std::size_t i = 0; i < 2 * kLimbs; ++i) big[i + 1] = exact[i];
    std::int32_t big_exponent = a.exponent_ + b.exponent_ + 1;
    if ((big[kWide - 1] & kTopBit) == 0) {
      detail::ShiftLeft(big, 1);
      --big_exponent;
    }
    for (std::size_t i = 0; i < kLimbs; ++i) {
      small[i + kLimbs + 1] = c.limbs_[i];
    }
    std::int32_t small_exponent = c.exponent_;
    bool big_negative = product_negative;
    bool small_negative = c.negative_ != 0;

    if (big_exponent < small_exponent ||
        (big_exponent == small_exponent &&
         detail::CompareLimbs(big, small) < 0)) {
      const Wide temporary = big;
      big = small;
      small = temporary;
      const std::int32_t exponent = big_exponent;
      big_exponent = small_exponent;
      small_exponent = exponent;
      const bool negative = big_negative;
      big_negative = small_negative;
      small_negative = negative;
    }

    const std::int32_t distance = big_exponent - small_exponent;
    constexpr int kWideBits = static_cast<int>(kWide) * kLimbBits;
    if (distance >= kWideBits) {
      small.fill(0);
      small[0] = 1;
    } else if (distance > 0) {
      const int count = distance;
      Wide lost = small;
      detail::ShiftLeft(lost, kWideBits - count);
      detail::ShiftRight(small, count);
      small[0] |= detail::IsZero(lost) ? 0 : 1;
    }

    std::int32_t exponent = big_exponent;
    Wide total;
    if (big_negative == small_negative) {
      carry_t carry = 0;
      for (std::size_t i = 0; i < kWide; ++i) {
        total[i] = adc(big[i], small[i], carry);
      }
      if (carry != 0) {
        const limb_t sticky = total[0] & 1;
        detail::ShiftRight(total, 1);
        total[0] |= sticky;
        total[kWide - 1] |= kTopBit;
        ++exponent;
      }
    } else {
      carry_t borrow = 0;
      for (std::size_t i = 0; i < kWide; ++i) {
        total[i] = sbb(big[i], small[i], borrow);
      }
      if (detail::IsZero(total)) return Zero(env.rounding() == Rounding::kDown);
      int leading = 0;
      for (std::size_t i = kWide; i-- > 0;) {
        if (total[i] != 0) {
          leading += countl_zero(total[i]);
          break;
        }
        leading += kLimbBits;
      }
      if (leading != 0) detail::ShiftLeft(total, leading);
      exponent -= leading;
    }

    Extended value;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      value.limbs[i] = total[i + kLimbs + 1];
    }
    limb_t lost = 0;
    for (std::size_t i = 0; i < kLimbs; ++i) lost |= total[i];
    value.rest = total[kLimbs] | (lost != 0 ? 1 : 0);
    return Finish(big_negative, exponent, value, env);
  }

  // --- Remainder ---

  // Builds an exact nonzero result whose significand may have leading zeros.
  static Float NormalizedExact(bool negative, std::int32_t exponent,
                               Limbs limbs) noexcept {
    const int leading = CountLeadingZeros(limbs);
    detail::ShiftLeft(limbs, leading);
    return FromParts(negative, exponent - leading, limbs);
  }

  template <class Env>
  static Float RemainderImpl(const Float& x, const Float& y,
                             const Env& env) noexcept {
    if (x.is_nan() || y.is_nan()) return PropagateNaN(x, y, env);
    if (x.is_inf() || y.is_zero()) return InvalidOperation(env);
    if (y.is_inf() || x.is_zero()) return x;

    const bool negative = x.negative_ != 0;
    const Limbs& divisor = y.limbs_;
    const std::int32_t distance = x.exponent_ - y.exponent_;
    if (distance < -1) return x;  // |x| < |y| / 2, so the quotient rounds to 0.
    if (distance == -1) {
      // |x| < |y|, and the quotient rounds to 1 exactly when 2|x| > |y|.
      if (detail::CompareLimbs(x.limbs_, divisor) <= 0) return x;
      // 2|y| - |x| fits in L limbs, so wrapping past the top bit is harmless.
      Limbs doubled = divisor;
      detail::ShiftLeft(doubled, 1);
      carry_t borrow = 0;
      Limbs r;
      for (std::size_t i = 0; i < kLimbs; ++i) {
        r[i] = sbb(doubled[i], x.limbs_[i], borrow);
      }
      return NormalizedExact(!negative, y.exponent_ - 1, r);
    }

    // Long division, keeping the remainder r and the quotient's lowest bit.
    Limbs r = x.limbs_;
    bool odd = false;
    const auto subtract_divisor = [&r, &divisor]() {
      carry_t borrow = 0;
      for (std::size_t i = 0; i < kLimbs; ++i) {
        r[i] = sbb(r[i], divisor[i], borrow);
      }
    };
    if (detail::CompareLimbs(r, divisor) >= 0) {
      subtract_divisor();
      odd = true;
    }
    for (int i = distance % kLimbBits; i > 0; --i) {
      const bool overflow = (r[kLimbs - 1] & kTopBit) != 0;
      detail::ShiftLeft(r, 1);
      odd = overflow || detail::CompareLimbs(r, divisor) >= 0;
      if (odd) subtract_divisor();
    }
    const std::int32_t words = distance / kLimbBits;
    if (words != 0) {
      const limb_t d1 = divisor[kLimbs - 1];
      limb_t d0 = 0;
      if constexpr (kLimbs >= 2) d0 = divisor[kLimbs - 2];
      const limb_t inverse = reciprocal_3by2(d1, d0);
      limb_t q = 0;
      for (std::int32_t i = 0; i < words; ++i) {
        q = DivideStep(r, divisor, inverse);
      }
      odd = (q & 1) != 0;
    }
    if (detail::IsZero(r)) return Zero(negative);

    // Round the quotient to nearest by comparing 2r with the divisor.
    int order = 1;
    if ((r[kLimbs - 1] & kTopBit) == 0) {
      Limbs doubled = r;
      detail::ShiftLeft(doubled, 1);
      order = detail::CompareLimbs(doubled, divisor);
    }
    bool result_negative = negative;
    if (order > 0 || (order == 0 && odd)) {
      carry_t borrow = 0;
      for (std::size_t i = 0; i < kLimbs; ++i) {
        r[i] = sbb(divisor[i], r[i], borrow);
      }
      result_negative = !negative;
    }
    return NormalizedExact(result_negative, y.exponent_, r);
  }

  // --- Rounding to integers ---

  template <class Env>
  static Float RoundToIntegralImpl(const Float& x, const Env& env) noexcept {
    if (x.is_nan()) return PropagateNaN(x, x, env);
    if (!IsRegular(x) || x.exponent_ >= P - 1) return x;

    const bool negative = x.negative_ != 0;
    if (x.exponent_ < 0) {
      const bool round_bit = x.exponent_ == -1;
      const bool sticky = x.exponent_ < -1 || x.limbs_ != TopBitOnly();
      env.raise(kFlagInexact);
      if (detail::ShouldRoundUp(env.rounding(), negative, false, round_bit,
                                sticky)) {
        return FromParts(negative, 0, TopBitOnly());
      }
      return Zero(negative);
    }

    constexpr int kTotalBits = static_cast<int>(kLimbs) * kLimbBits;
    const int count = kTotalBits - 1 - x.exponent_;
    Limbs limbs = x.limbs_;
    bool is_odd = false;
    bool round_bit = false;
    bool sticky = false;
    InspectLowBits(limbs, count, is_odd, round_bit, sticky);
    if (!round_bit && !sticky) return x;
    env.raise(kFlagInexact);
    detail::ShiftRight(limbs, count);
    detail::ShiftLeft(limbs, count);
    std::int32_t exponent = x.exponent_;
    if (detail::ShouldRoundUp(env.rounding(), negative, is_odd, round_bit,
                              sticky) &&
        IncrementAtBit(limbs, count)) {
      ++exponent;
      // Only formats with kEmax < P - 1 can get here.
      if (exponent > kEmax) return Overflow(negative, env);
    }
    return FromParts(negative, exponent, limbs);
  }

  // --- Minimum and maximum ---

  // Orders -0 below +0. Neither operand may be NaN.
  static bool IsLess(const Float& a, const Float& b) noexcept {
    if (a.is_zero() && b.is_zero()) return a.negative_ > b.negative_;
    return Order(a, b) == std::partial_ordering::less;
  }

  template <bool kMaximum, bool kPreferNumber, class Env>
  static Float MinMaxImpl(const Float& a, const Float& b,
                          const Env& env) noexcept {
    if (a.is_nan() || b.is_nan()) {
      if constexpr (kPreferNumber) {
        if (!a.is_nan() || !b.is_nan()) {
          if (a.is_signaling() || b.is_signaling()) env.raise(kFlagInvalid);
          return a.is_nan() ? b : a;
        }
      }
      return PropagateNaN(a, b, env);
    }
    if constexpr (kMaximum) {
      return IsLess(a, b) ? b : a;
    } else {
      return IsLess(b, a) ? b : a;
    }
  }

  // --- Neighbours and exponents ---

  template <class Env>
  static Float NextImpl(const Float& x, bool up, const Env& env) noexcept {
    if (x.is_nan()) return PropagateNaN(x, x, env);
    if (x.is_zero()) return MinSubnormal(!up);
    const bool negative = x.negative_ != 0;
    if (x.is_inf()) return negative == up ? MaxFinite(negative) : x;

    // In the interchange encoding, the magnitude's neighbours are the adjacent
    // integers, across subnormals and up to infinity alike.
    Bits bits = x.to_bits();
    if (up != negative) {
      carry_t carry = 1;
      for (limb_t& limb : bits) limb = adc(limb, 0, carry);
    } else {
      carry_t borrow = 1;
      for (limb_t& limb : bits) limb = sbb(limb, 0, borrow);
    }
    return from_bits(bits);
  }

  template <class Env>
  static Float ScalbnImpl(const Float& x, int n, const Env& env) noexcept {
    if (x.is_nan()) return PropagateNaN(x, x, env);
    if (!IsRegular(x)) return x;
    // Clamping keeps the sum in range without changing any result, since a
    // value this far out overflows or vanishes whatever its significand.
    const std::int64_t limit = std::int64_t{1} << (E + 1);
    std::int64_t exponent = std::int64_t{x.exponent_} + n;
    if (exponent > limit + P) exponent = limit + P;
    if (exponent < -limit - P) exponent = -limit - P;
    Extended value;
    value.limbs = x.limbs_;
    value.rest = 0;
    return Finish(x.negative_ != 0, static_cast<std::int32_t>(exponent), value,
                  env);
  }

  template <class Env>
  static int IlogbImpl(const Float& x, const Env& env) noexcept {
    if (IsRegular(x)) return x.exponent_;
    env.raise(kFlagInvalid);
    if (x.is_nan()) return kIlogbNaN;
    return x.is_zero() ? kIlogbZero : kIlogbInfinity;
  }

  template <class Env>
  static Float LogbImpl(const Float& x, const Env& env) noexcept {
    if (x.is_nan()) return PropagateNaN(x, x, env);
    if (x.is_inf()) return Infinity(false);
    if (x.is_zero()) {
      env.raise(kFlagDivByZero);
      return Infinity(true);
    }
    if (x.exponent_ == 0) return Zero(false);
    const bool negative = x.exponent_ < 0;
    const limb_t magnitude = static_cast<limb_t>(
        negative ? -std::int64_t{x.exponent_} : std::int64_t{x.exponent_});
    const int leading = countl_zero(magnitude);
    Extended value;
    value.limbs = Limbs{};
    value.limbs[kLimbs - 1] = magnitude << leading;
    value.rest = 0;
    // An exponent always fits the significand in practice. Tiny test formats
    // may need rounding here, which must not disturb the caller's mode or flags.
    return Finish(negative, kLimbBits - 1 - leading, value,
                  detail::SilentEnv{Rounding::kNearestEven});
  }

  // --- Comparison ---

  static std::partial_ordering Order(const Float& a, const Float& b) noexcept {
    if (a.is_nan() || b.is_nan()) return std::partial_ordering::unordered;
    if (a.is_zero() && b.is_zero()) return std::partial_ordering::equivalent;
    if (a.negative_ != b.negative_) {
      return a.negative_ != 0 ? std::partial_ordering::less
                              : std::partial_ordering::greater;
    }
    int order = 0;
    if (a.exponent_ != b.exponent_) {
      order = a.exponent_ < b.exponent_ ? -1 : 1;
    } else {
      order = detail::CompareLimbs(a.limbs_, b.limbs_);
    }
    if (a.negative_ != 0) order = -order;
    if (order == 0) return std::partial_ordering::equivalent;
    return order < 0 ? std::partial_ordering::less
                     : std::partial_ordering::greater;
  }

  Limbs limbs_;
  std::int32_t exponent_;
  std::uint32_t negative_;
};

using binary128 = Float<113, 15>;
using binary256 = Float<237, 19>;
using binary512 = Float<489, 23>;

}  // namespace hwymp

#endif  // HWYMP_FLOAT_H_
