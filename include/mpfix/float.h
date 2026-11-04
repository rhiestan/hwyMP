// Copyright 2026 The mpfix Authors
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

#ifndef MPFIX_FLOAT_H_
#define MPFIX_FLOAT_H_

#include <array>
#include <bit>
#include <cmath>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "mpfix/limb_ops.h"

namespace mpfix {

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
  static constexpr void raise(unsigned /*flags*/) noexcept {}
};

// Rounds in a given mode and discards flags, for tentative rounding steps.
struct SilentEnv {
  Rounding mode;
  constexpr Rounding rounding() const noexcept { return mode; }
  static constexpr void raise(unsigned /*flags*/) noexcept {}
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
MPFIX_ALWAYS_INLINE constexpr bool ShouldRoundUp(Rounding mode, bool negative,
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
      return (!negative) & (round_bit | sticky);
    case Rounding::kDown:
      return negative & (round_bit | sticky);
  }
  return false;
}

// Shifts the little-endian limb array left by 0 <= count < 64 * N bits.
template <std::size_t N>
constexpr void ShiftLeft(std::array<limb_t, N>& x, int count) noexcept {
  // Unsigned, so that compilers see 0 <= words and 0 <= bits < 64.
  const auto shift = static_cast<unsigned>(count);
  const int words = static_cast<int>(shift / kLimbBits);
  const int bits = static_cast<int>(shift % kLimbBits);
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
  // Unsigned, so that compilers see 0 <= words and 0 <= bits < 64.
  const auto shift = static_cast<unsigned>(count);
  const int words = static_cast<int>(shift / kLimbBits);
  const int bits = static_cast<int>(shift % kLimbBits);
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
MPFIX_ALWAYS_INLINE limb_t IntegerSqrt(limb_t hi, limb_t lo,
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
  // A division here is faster than computing the step in floating point on
  // Zen 3: measured 17.5 against 19.7 ns for binary128's square root.
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
class Float;

namespace detail {

template <class T>
inline constexpr bool kIsFloat = false;
template <int P, int E>
inline constexpr bool kIsFloat<Float<P, E>> = true;

// Built-in types that convert to and from Float. bool is left out, so that a
// Float never silently becomes 0 or 1.
template <class T>
concept Arithmetic = std::is_arithmetic_v<T> && !std::is_same_v<T, bool>;

template <class T>
concept Integer = std::is_integral_v<T> && !std::is_same_v<T, bool>;

template <class T>
concept BinaryFloatingPoint = std::is_floating_point_v<T> &&
                              std::numeric_limits<T>::radix == 2 &&
                              std::numeric_limits<T>::has_infinity;

// The Float with the same values as a built-in floating-point type: float is
// Float<24, 8>, double is Float<53, 11>, and the x87 long double is
// Float<64, 15>.
template <BinaryFloatingPoint T>
struct NativeFormatOf {
  static constexpr int kDigits = std::numeric_limits<T>::digits;
  static constexpr int kExponentBits =
      std::bit_width(
          static_cast<unsigned>(std::numeric_limits<T>::max_exponent));
  using type = Float<kDigits, kExponentBits>;
  // Whether T is laid out as the interchange encoding of `type`, which lets
  // its bits be copied directly.
  static constexpr bool kInterchange =
      (sizeof(T) == 4 || sizeof(T) == 8) &&
      sizeof(T) * 8 == static_cast<std::size_t>(kDigits + kExponentBits);
};

template <class T>
using NativeFormat = typename NativeFormatOf<T>::type;

template <class T>
using SameSizeUnsigned =
    std::conditional_t<sizeof(T) == 4, std::uint32_t, std::uint64_t>;

// Whether every value of T converts exactly to Float<P, E>.
template <int P, int E, class T>
constexpr bool ExactlyConvertible() {
  constexpr int kEmax = (1 << (E - 1)) - 1;
  if constexpr (kIsFloat<T>) {
    return T::kPrecision <= P && T::kExponentBits <= E;
  } else if constexpr (std::is_integral_v<T>) {
    return P >= std::numeric_limits<T>::digits &&
           kEmax >= std::numeric_limits<T>::digits;
  } else {
    return NativeFormatOf<T>::kDigits <= P &&
           NativeFormatOf<T>::kExponentBits <= E;
  }
}

// Gives the free conversion functions and the string conversions access to the
// internals of every Float.
struct FloatAccess;

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
    if (Exponent() == kExpNaN || Exponent() == kExpInfinity) {
      biased = (limb_t{1} << E) - 1;
      trailing[kLimbs - 1] &= ~kTopBit;
    } else if (Exponent() == kExpZero) {
      biased = 0;
    } else if (Exponent() < kEmin) {
      biased = 0;
      detail::ShiftRight(trailing, kEmin - Exponent());
    } else {
      biased = static_cast<limb_t>(Exponent() + kEmax);
      trailing[kLimbs - 1] &= ~kTopBit;
    }
    detail::ShiftRight(trailing, kSpareBits);

    Bits bits{};
    for (std::size_t i = 0; i < kLimbs; ++i) bits[i] = trailing[i];
    InsertBits(bits, P - 1, E, biased);
    InsertBits(bits, P + E - 1, 1, negative_ ? 1 : 0);
    return bits;
  }

  // --- Conversions ---
  //
  // Constructing a Float from a built-in arithmetic type or from another format
  // rounds to nearest-even. The constructors are implicit when they are always
  // exact, as from double or int64_t to binary128. The free functions convert()
  // and to_integer() below take a Context for the other modes and for flags.
  //
  // NaN payloads keep their most significant bits, so a NaN that is widened and
  // narrowed again comes back unchanged. Signaling NaNs are quieted and raise
  // invalid.

  // Whether every value of T converts exactly.
  template <class T>
  static constexpr bool kExactFrom = detail::ExactlyConvertible<P, E, T>();

  template <detail::Arithmetic T>
  explicit(!kExactFrom<T>) Float(T value) noexcept
      : Float(ConvertFrom(value, detail::NearestEnv{})) {}

  template <int P2, int E2>
    requires(P2 != P || E2 != E)
  explicit(!kExactFrom<Float<P2, E2>>)
      Float(const Float<P2, E2>& value) noexcept
      : Float(ConvertFormat(value, detail::NearestEnv{})) {}

  // Rounds to nearest-even.
  template <detail::BinaryFloatingPoint T>
  explicit operator T() const noexcept {
    return ConvertTo<T>(detail::NearestEnv{});
  }

  // Truncates toward zero and saturates: NaN converts to 0, and values beyond
  // the range of I to its minimum or maximum.
  template <detail::Integer I>
  explicit operator I() const noexcept {
    return ToInteger<I>(Rounding::kTowardZero, false,
                        detail::NearestEnv{});
  }

  // True for anything but zero, NaN included, as for the built-in types.
  explicit operator bool() const noexcept { return !is_zero(); }

  // --- Special values ---

  static constexpr Float Zero(bool negative = false) noexcept {
    return FromParts(negative, kExpZero, Limbs{});
  }

  static constexpr Float Infinity(bool negative = false) noexcept {
    return FromParts(negative, kExpInfinity, TopBitOnly());
  }

  static constexpr Float QuietNaN(bool negative = false) noexcept {
    Limbs limbs = TopBitOnly();
    SetQuietBit(limbs);
    return FromParts(negative, kExpNaN, limbs);
  }

  // The signaling NaN whose payload has only its lowest bit set.
  static constexpr Float SignalingNaN(bool negative = false) noexcept {
    Limbs limbs = TopBitOnly();
    limbs[0] |= limb_t{1} << kSpareBits;
    return FromParts(negative, kExpNaN, limbs);
  }

  static constexpr Float MaxFinite(bool negative = false) noexcept {
    return FromParts(negative, kEmax, AllOnes());
  }

  static constexpr Float MinNormal(bool negative = false) noexcept {
    return FromParts(negative, kEmin, TopBitOnly());
  }

  static constexpr Float MinSubnormal(bool negative = false) noexcept {
    return FromParts(negative, kEmin - (P - 1), TopBitOnly());
  }

  // 2^e, for e from kEmin - (P - 1), the smallest subnormal, to kEmax.
  static constexpr Float PowerOfTwo(int e) noexcept {
    return FromParts(false, e, TopBitOnly());
  }

  // --- Classification ---

  constexpr bool signbit() const noexcept { return negative_ != 0; }
  constexpr bool is_nan() const noexcept { return Exponent() == kExpNaN; }
  constexpr bool is_inf() const noexcept {
    return Exponent() == kExpInfinity;
  }
  constexpr bool is_zero() const noexcept { return Exponent() == kExpZero; }
  constexpr bool is_finite() const noexcept {
    return Exponent() < kExpInfinity;
  }
  constexpr bool is_subnormal() const noexcept {
    return !is_zero() && Exponent() < kEmin;
  }
  constexpr bool is_normal() const noexcept {
    return Exponent() >= kEmin && Exponent() < kExpInfinity;
  }
  constexpr bool is_signaling() const noexcept {
    return is_nan() && (limbs_[kLimbs - 1] & kQuietBit) == 0;
  }

  // Free forms, as <cmath> has them for the built-in types, so that generic
  // code calling them unqualified finds these.
  friend constexpr bool isnan(const Float& x) noexcept { return x.is_nan(); }
  friend constexpr bool isinf(const Float& x) noexcept { return x.is_inf(); }
  friend constexpr bool isfinite(const Float& x) noexcept {
    return x.is_finite();
  }
  friend constexpr bool isnormal(const Float& x) noexcept {
    return x.is_normal();
  }
  friend constexpr bool signbit(const Float& x) noexcept {
    return x.signbit();
  }

  FloatClass classify() const noexcept {
    if (is_nan()) return FloatClass::kNaN;
    if (is_inf()) return FloatClass::kInfinite;
    if (is_zero()) return FloatClass::kZero;
    return Exponent() < kEmin ? FloatClass::kSubnormal : FloatClass::kNormal;
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

  // The neighbour of x in the direction of y, as in C: y itself if x == y,
  // which gives the sign of y to a zero. Unlike next_up, a finite x stepping
  // to infinity raises overflow and inexact, and a result that is subnormal or
  // zero raises underflow and inexact.
  friend Float nextafter(const Float& x, const Float& y) noexcept {
    return NextAfterImpl(x, y, detail::NearestEnv{});
  }
  friend Float nextafter(const Float& x, const Float& y,
                         Context& context) noexcept {
    return NextAfterImpl(x, y, detail::ContextEnv{context});
  }
  // As nextafter, toward a long double, compared exactly.
  friend Float nexttoward(const Float& x, long double y) noexcept {
    return NextTowardImpl(x, y, detail::NearestEnv{});
  }
  friend Float nexttoward(const Float& x, long double y,
                          Context& context) noexcept {
    return NextTowardImpl(x, y, detail::ContextEnv{context});
  }

  // x * 2^n.
  friend Float scalbn(const Float& x, int n) noexcept {
    return ScalbnImpl(x, n, detail::NearestEnv{});
  }
  friend Float scalbn(const Float& x, int n, Context& context) noexcept {
    return ScalbnImpl(x, n, detail::ContextEnv{context});
  }
  friend Float scalbln(const Float& x, long n) noexcept {
    return ScalbnImpl(x, n, detail::NearestEnv{});
  }
  friend Float scalbln(const Float& x, long n, Context& context) noexcept {
    return ScalbnImpl(x, n, detail::ContextEnv{context});
  }
  // The same as scalbn, since the radix is 2.
  friend Float ldexp(const Float& x, int n) noexcept {
    return ScalbnImpl(x, n, detail::NearestEnv{});
  }
  friend Float ldexp(const Float& x, int n, Context& context) noexcept {
    return ScalbnImpl(x, n, detail::ContextEnv{context});
  }

  // Splits x into a significand in [0.5, 1) with the sign of x, which it
  // returns, and a power of two, stored in *exponent. Subnormals return a
  // normalized significand. Zeros, infinities and NaNs return x, quieted, and
  // store 0. Always exact; only a signaling NaN raises a flag, invalid.
  friend Float frexp(const Float& x, int* exponent) noexcept {
    return FrexpImpl(x, exponent, detail::NearestEnv{});
  }
  friend Float frexp(const Float& x, int* exponent, Context& context) noexcept {
    return FrexpImpl(x, exponent, detail::ContextEnv{context});
  }

  // Splits x into its integer part, stored in *integral, and its fractional
  // part, which it returns. Both have the sign of x, so modf(-2) gives -0 and
  // -2. Infinities give a zero fraction, and NaNs give NaN for both. Always
  // exact; only a signaling NaN raises a flag, invalid.
  friend Float modf(const Float& x, Float* integral) noexcept {
    return ModfImpl(x, integral, detail::NearestEnv{});
  }
  friend Float modf(const Float& x, Float* integral,
                    Context& context) noexcept {
    return ModfImpl(x, integral, detail::ContextEnv{context});
  }

  static constexpr int kIlogbZero = (std::numeric_limits<int>::min)();
  static constexpr int kIlogbInfinity = (std::numeric_limits<int>::max)();
  static constexpr int kIlogbNaN = (std::numeric_limits<int>::max)();

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

  std::int32_t raw_exponent() const noexcept { return Exponent(); }
  const Limbs& raw_limbs() const noexcept { return limbs_; }

  // `limbs` must be normalized and `exponent` within the format's range.
  static Float FromRaw(bool negative, std::int32_t exponent,
                       const Limbs& limbs) noexcept {
    return FromParts(negative, exponent, limbs);
  }

 private:
  template <int, int>
  friend class Float;
  friend struct detail::FloatAccess;

  static constexpr int kSpareBits = static_cast<int>(kLimbs) * kLimbBits - P;
  static constexpr limb_t kTopBit = limb_t{1} << (kLimbBits - 1);
  static constexpr limb_t kQuietBit = kTopBit >> 1;

  // Exponent values that mark non-numbers. They lie far outside the range of
  // any regular value and keep the ordering zero < regular < infinity.
  static constexpr std::int32_t kExpZero =
      (std::numeric_limits<std::int32_t>::min)();
  static constexpr std::int32_t kExpInfinity =
      (std::numeric_limits<std::int32_t>::max)() - 1;
  static constexpr std::int32_t kExpNaN =
      (std::numeric_limits<std::int32_t>::max)();

  // A significand with the extra limb that holds the bits below the last limb:
  // `rest` has the highest of those bits at its top, and its lowest bit is
  // sticky, meaning that it also stands for everything shifted out below it.
  struct Extended {
    Limbs limbs;
    limb_t rest;
  };

  static constexpr std::uint32_t Bias(std::int32_t exponent) noexcept {
    return static_cast<std::uint32_t>(exponent) ^ 0x80000000u;
  }
  constexpr std::int32_t Exponent() const noexcept {
    return static_cast<std::int32_t>(biased_exponent_ ^ 0x80000000u);
  }

  static constexpr Float FromParts(bool negative, std::int32_t exponent,
                                   const Limbs& limbs) noexcept {
    Float result;
    result.limbs_ = limbs;
    result.biased_exponent_ = Bias(exponent);
    result.negative_ = negative ? 1 : 0;
    return result;
  }

  static constexpr Limbs TopBitOnly() noexcept {
    Limbs limbs{};
    limbs[kLimbs - 1] = kTopBit;
    return limbs;
  }

  static constexpr Limbs AllOnes() noexcept {
    Limbs limbs{};
    limbs.fill(~limb_t{0});
    if constexpr (kSpareBits != 0) limbs[0] &= ~LowMask(kSpareBits);
    return limbs;
  }

  static constexpr void SetQuietBit(Limbs& limbs) noexcept {
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
    return x.Exponent() != kExpZero && x.Exponent() < kExpInfinity;
  }

  // --- Rounding ---

  // Rounds to P bits. Returns true when the significand overflowed, in which
  // case it is left as 2^(64L - 1) and the exponent must grow by one.
  template <class Env>
  MPFIX_ALWAYS_INLINE static bool RoundToPrecision(bool negative,
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
  MPFIX_ALWAYS_INLINE static Float Finish(bool negative, std::int32_t exponent,
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
    constexpr int kStorageBits = static_cast<int>(kLimbs) * kLimbBits;

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
    if (discard > kStorageBits) {
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

    if (discard >= kStorageBits) {
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
    constexpr int kStorageBits = static_cast<int>(kLimbs) * kLimbBits;
    const auto bit_at = [&limbs](int position) {
      return ((limbs[static_cast<std::size_t>(position) / kLimbBits] >>
               (position % kLimbBits)) &
              1) != 0;
    };
    round_bit = bit_at(count - 1);
    if (count < kStorageBits) is_odd = bit_at(count);
    Limbs below = limbs;
    detail::ShiftLeft(below, kStorageBits - count);
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

  // Shifts `lesser` right by distance >= 0, jamming lost bits into `rest`.
  MPFIX_ALWAYS_INLINE static Extended AlignRight(const Limbs& lesser,
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
      result.rest = (lesser[0] << 1) << complement;
      for (std::size_t i = 0; i + 1 < kLimbs; ++i) {
        result.limbs[i] =
            (lesser[i] >> bits) | ((lesser[i + 1] << 1) << complement);
      }
      result.limbs[kLimbs - 1] = lesser[kLimbs - 1] >> bits;
      return result;
    }

    // Element 0 is `rest`; elements 1 to L are the shifted limbs.
    std::array<limb_t, kLimbs + 1> wide{};
    for (int i = 0; i + words <= kWords; ++i) {
      const int source = i + words - 1;
      wide[static_cast<std::size_t>(i)] =
          source >= 0 ? lesser[static_cast<std::size_t>(source)] : 0;
    }
    limb_t lost = 0;
    for (int i = 0; i + 1 < words; ++i) {
      lost |= lesser[static_cast<std::size_t>(i)];
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
  MPFIX_ALWAYS_INLINE static Float AddImpl(const Float& a, const Float& b,
                                           std::uint32_t b_sign,
                                           const Env& env) noexcept {
    if (!IsRegular(a) || !IsRegular(b)) [[unlikely]] {
      return AddSpecial(a, b, b_sign, env);
    }

    const bool swap = a.Exponent() < b.Exponent();
    const Limbs* big = swap ? &b.limbs_ : &a.limbs_;
    const Limbs* lesser = swap ? &a.limbs_ : &b.limbs_;
    std::int32_t exponent = swap ? b.Exponent() : a.Exponent();
    std::int32_t distance = swap ? b.Exponent() - a.Exponent()
                                 : a.Exponent() - b.Exponent();
    bool negative = (swap ? b_sign : a.negative_) != 0;

    if (a.negative_ == b_sign) {
      Extended sum = AlignRight(*lesser, distance);
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
      const int order = detail::CompareLimbs(*big, *lesser);
      if (order == 0) return Zero(env.rounding() == Rounding::kDown);
      if (order < 0) {
        const Limbs* temporary = big;
        big = lesser;
        lesser = temporary;
        negative = !negative;
      }
    }

    Extended difference = AlignRight(*lesser, distance);
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
  MPFIX_ALWAYS_INLINE static void ShiftLeftByBit(Extended& value,
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
  MPFIX_ALWAYS_INLINE static Float MulImpl(const Float& a, const Float& b,
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
    const std::int32_t exponent = a.Exponent() + b.Exponent() + 1 - shift;
    return Finish(negative, exponent, value, env);
  }

  MPFIX_ALWAYS_INLINE static std::array<limb_t, 2 * kLimbs>
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
  MPFIX_ALWAYS_INLINE static limb_t DivideStep(Limbs& r, const Limbs& divisor,
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
  MPFIX_ALWAYS_INLINE static Float DivImpl(const Float& a, const Float& b,
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
    const std::int32_t exponent = a.Exponent() - b.Exponent() - 1 + integer_part;
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
    const int shift = (x.Exponent() & 1) ^ 1;
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
    if constexpr (kLimbs == 2) {
      // One more limb, whose estimate is (R / 2) / S by a single division:
      // S has one limb so far. Too high by at most 2, which the update
      // corrects.
      NextRootLimbOfTwo(root, rem);
    } else if constexpr (kLimbs > 2) {
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
    return Finish(false, x.Exponent() >> 1, value, env);
  }

  // Finds limb j of the root and updates the remainder. Limbs above j are final,
  // and `inverse` matches the root's top two limbs.
  MPFIX_ALWAYS_INLINE static void NextRootLimb(
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
    ApplyRootLimb(root, rem, j, q);
  }

  // The second and last limb of a two-limb root, from the first: (R / 2) / S
  // by one division, as S has a single limb. Like the general estimate, it
  // can only be too high, and ApplyRootLimb corrects it.
  MPFIX_ALWAYS_INLINE static void NextRootLimbOfTwo(
      Limbs& root, std::array<limb_t, 2 * kLimbs>& rem) noexcept {
    const limb_t u2 = (rem[2] >> 1) | (rem[3] << (kLimbBits - 1));
    const limb_t u1 = (rem[1] >> 1) | (rem[2] << (kLimbBits - 1));
    limb_t q = ~limb_t{0};
    if (u2 < root[1]) {
      limb_t unused;
      q = div_2by1(u2, u1, root[1], unused);
    }
    ApplyRootLimb(root, rem, 0, q);
  }

  // Sets limb j of the root to the estimate q, which is at or above the true
  // limb, and updates the remainder, stepping q down while it is negative.
  MPFIX_ALWAYS_INLINE static void ApplyRootLimb(
      Limbs& root, std::array<limb_t, 2 * kLimbs>& rem, std::size_t j,
      limb_t q) noexcept {
    constexpr std::size_t kWide = 2 * kLimbs;
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

  // Shifts right by count >= 0 and folds the bits shifted out into the lowest
  // bit, which is sticky.
  template <std::size_t N>
  static void ShiftRightSticky(std::array<limb_t, N>& x,
                               std::int32_t count) noexcept {
    if (count == 0) return;
    if (count < kLimbBits) {
      // The common case, within a limb.
      const int bits = static_cast<int>(count);
      const limb_t lost = x[0] << (kLimbBits - bits);
      for (std::size_t i = 0; i + 1 < N; ++i) {
        x[i] = (x[i] >> bits) | (x[i + 1] << (kLimbBits - bits));
      }
      x[N - 1] >>= bits;
      x[0] |= lost != 0 ? 1 : 0;
      return;
    }
    constexpr std::int32_t kBits = static_cast<std::int32_t>(N) * kLimbBits;
    if (count >= kBits) {
      const bool any = !detail::IsZero(x);
      x.fill(0);
      x[0] = any ? 1 : 0;
      return;
    }
    const int words = static_cast<int>(count / kLimbBits);
    const int bits = static_cast<int>(count % kLimbBits);
    limb_t lost = 0;
    for (int i = 0; i < words; ++i) lost |= x[static_cast<std::size_t>(i)];
    if (bits != 0) lost |= x[static_cast<std::size_t>(words)] << (kLimbBits - bits);
    detail::ShiftRight(x, static_cast<int>(count));
    x[0] |= lost != 0 ? 1 : 0;
  }

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

    // Addition at double width, on the exact product. Both operands sit in a
    // frame of 2L + 1 limbs whose top bit has weight 2^exponent; limb 0 is a
    // guard limb whose lowest bit is sticky, as in `Extended`. The product,
    // in [1, 4), fills limbs 1 to 2L with its top bit possibly clear; c fills
    // the top L limbs. Whichever has the higher frame keeps it, and the other
    // shifts right into it.
    constexpr std::size_t kWide = 2 * kLimbs + 1;
    using Wide = std::array<limb_t, kWide>;

    const std::array<limb_t, 2 * kLimbs> exact =
        MultiplySignificands(a.limbs_, b.limbs_);
    Wide product;
    product[0] = 0;
    for (std::size_t i = 0; i < 2 * kLimbs; ++i) product[i + 1] = exact[i];
    Wide addend;
    for (std::size_t i = 0; i <= kLimbs; ++i) addend[i] = 0;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      addend[i + kLimbs + 1] = c.limbs_[i];
    }
    const std::int32_t product_exponent = a.Exponent() + b.Exponent() + 1;
    const std::int32_t distance = product_exponent - c.Exponent();
    const bool addend_negative = c.negative_ != 0;

    std::int32_t exponent;
    if (distance >= 0) {
      exponent = product_exponent;
      ShiftRightSticky(addend, distance);
    } else {
      exponent = c.Exponent();
      ShiftRightSticky(product, -distance);
    }

    Wide total;
    bool negative = product_negative;
    if (product_negative == addend_negative) {
      carry_t carry = 0;
      for (std::size_t i = 0; i < kWide; ++i) {
        total[i] = adc(product[i], addend[i], carry);
      }
      if (carry != 0) {
        const limb_t sticky = total[0] & 1;
        detail::ShiftRight(total, 1);
        total[0] |= sticky;
        total[kWide - 1] |= kTopBit;
        ++exponent;
      }
    } else if (distance < 0) {
      // c's frame: the shifted product lies below c.
      carry_t borrow = 0;
      for (std::size_t i = 0; i < kWide; ++i) {
        total[i] = sbb(addend[i], product[i], borrow);
      }
      negative = addend_negative;
    } else {
      carry_t borrow = 0;
      for (std::size_t i = 0; i < kWide; ++i) {
        total[i] = sbb(product[i], addend[i], borrow);
      }
      if (borrow != 0) {
        // |c| exceeded the product. In the product's frame that needs c
        // shifted by at most one bit, which loses nothing, so the difference
        // is exact: negate it.
        carry_t carry = 1;
        for (std::size_t i = 0; i < kWide; ++i) {
          total[i] = adc(~total[i], 0, carry);
        }
        negative = addend_negative;
      }
    }
    // Leading zeros from cancellation, or the one of a product below 2.
    int leading = 0;
    for (std::size_t i = kWide; i-- > 0;) {
      if (total[i] != 0) {
        leading += countl_zero(total[i]);
        break;
      }
      leading += kLimbBits;
    }
    if (leading == static_cast<int>(kWide) * kLimbBits) {
      return Zero(env.rounding() == Rounding::kDown);
    }
    if (leading >= kLimbBits) {
      detail::ShiftLeft(total, leading);
    } else if (leading != 0) {
      // Within a limb: a product below 2, or cancellation of a few bits.
      for (std::size_t i = kWide; i-- > 1;) {
        total[i] = (total[i] << leading) | (total[i - 1] >> (kLimbBits - leading));
      }
      total[0] <<= leading;
    }
    exponent -= leading;

    Extended value;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      value.limbs[i] = total[i + kLimbs + 1];
    }
    limb_t lost = 0;
    for (std::size_t i = 0; i < kLimbs; ++i) lost |= total[i];
    value.rest = total[kLimbs] | (lost != 0 ? 1 : 0);
    return Finish(negative, exponent, value, env);
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
    const std::int32_t distance = x.Exponent() - y.Exponent();
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
      return NormalizedExact(!negative, y.Exponent() - 1, r);
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
    return NormalizedExact(result_negative, y.Exponent(), r);
  }

  // --- Rounding to integers ---

  template <class Env>
  static Float RoundToIntegralImpl(const Float& x, const Env& env) noexcept {
    if (x.is_nan()) return PropagateNaN(x, x, env);
    if (!IsRegular(x) || x.Exponent() >= P - 1) return x;

    const bool negative = x.negative_ != 0;
    if (x.Exponent() < 0) {
      const bool round_bit = x.Exponent() == -1;
      const bool sticky = x.Exponent() < -1 || x.limbs_ != TopBitOnly();
      env.raise(kFlagInexact);
      if (detail::ShouldRoundUp(env.rounding(), negative, false, round_bit,
                                sticky)) {
        return FromParts(negative, 0, TopBitOnly());
      }
      return Zero(negative);
    }

    constexpr int kStorageBits = static_cast<int>(kLimbs) * kLimbBits;
    const int count = kStorageBits - 1 - x.Exponent();
    Limbs limbs = x.limbs_;
    bool is_odd = false;
    bool round_bit = false;
    bool sticky = false;
    InspectLowBits(limbs, count, is_odd, round_bit, sticky);
    if (!round_bit && !sticky) return x;
    env.raise(kFlagInexact);
    detail::ShiftRight(limbs, count);
    detail::ShiftLeft(limbs, count);
    std::int32_t exponent = x.Exponent();
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

  // Steps from x toward y, given how they compare. Neither is NaN.
  template <class Env>
  static Float StepToward(const Float& x, std::partial_ordering order,
                          bool y_negative, const Env& env) noexcept {
    if (order == std::partial_ordering::equivalent) {
      Float result = x;
      result.negative_ = y_negative;
      return result;
    }
    const Float result =
        NextImpl(x, order == std::partial_ordering::less, env);
    if (result.is_inf()) {
      env.raise(kFlagOverflow | kFlagInexact);
    } else if (result.is_zero() || result.is_subnormal()) {
      env.raise(kFlagUnderflow | kFlagInexact);
    }
    return result;
  }

  template <class Env>
  static Float NextAfterImpl(const Float& x, const Float& y,
                             const Env& env) noexcept {
    if (x.is_nan() || y.is_nan()) return PropagateNaN(x, y, env);
    return StepToward(x, Order(x, y), y.negative_ != 0, env);
  }

  template <class Env>
  static Float NextTowardImpl(const Float& x, long double y,
                              const Env& env) noexcept {
    if (x.is_nan()) return PropagateNaN(x, x, env);
    if (std::isnan(y)) return QuietNaN();
    // A format holding both exactly, in which they compare.
    using Native = detail::NativeFormatOf<long double>;
    using Wide = Float<(std::max)(P, Native::kDigits),
                       (std::max)(E, Native::kExponentBits)>;
    const Wide wide_x = Wide::ConvertFormat(x, detail::NearestEnv{});
    const Wide wide_y(y);
    return StepToward(x, wide_x <=> wide_y, std::signbit(y), env);
  }

  template <class Env>
  static Float FrexpImpl(const Float& x, int* exponent,
                         const Env& env) noexcept {
    *exponent = 0;
    if (x.is_nan()) return PropagateNaN(x, x, env);
    if (!IsRegular(x)) return x;
    *exponent = x.Exponent() + 1;
    Float result = x;
    result.biased_exponent_ = Bias(-1);
    return result;
  }

  template <class Env>
  static Float ModfImpl(const Float& x, Float* integral,
                        const Env& env) noexcept {
    if (x.is_nan()) return *integral = PropagateNaN(x, x, env);
    *integral = RoundToIntegralImpl(x, detail::SilentEnv{Rounding::kTowardZero});
    if (x.is_inf()) return Zero(x.negative_ != 0);
    // Exact: below 1 the integer part is zero, and above it both operands
    // lie on the grid of x's last place.
    Float fraction = AddImpl(x, *integral, integral->negative_ ^ 1,
                             detail::NearestEnv{});
    fraction.negative_ = x.negative_;
    return fraction;
  }

  template <class Env>
  static Float ScalbnImpl(const Float& x, std::int64_t n,
                          const Env& env) noexcept {
    if (x.is_nan()) return PropagateNaN(x, x, env);
    if (!IsRegular(x)) return x;
    // Clamping keeps the sum in range without changing any result, since a
    // value this far out overflows or vanishes whatever its significand.
    const std::int64_t limit = std::int64_t{1} << (E + 1);
    // Clamping n first keeps the sum from overflowing for any long.
    if (n > 2 * (limit + P)) n = 2 * (limit + P);
    if (n < -2 * (limit + P)) n = -2 * (limit + P);
    std::int64_t exponent = std::int64_t{x.Exponent()} + n;
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
    if (IsRegular(x)) return x.Exponent();
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
    if (x.Exponent() == 0) return Zero(false);
    const bool negative = x.Exponent() < 0;
    const limb_t magnitude = static_cast<limb_t>(
        negative ? -std::int64_t{x.Exponent()} : std::int64_t{x.Exponent()});
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

  // --- Conversions ---

  // Rounds the value limbs / 2^(64 count - 1) * 2^exponent, where `limbs` has
  // `count` limbs, least significant first, and its top bit set. `sticky` adds
  // an infinitesimal below the last limb: it marks a value known to be inexact.
  template <class Env>
  static Float RoundSignificand(bool negative, std::int64_t exponent,
                                const limb_t* limbs, std::size_t count,
                                bool sticky, const Env& env) noexcept {
    // Initialized, although the loop sets every limb, for GCC's
    // -Wmaybe-uninitialized.
    Extended value{};
    for (std::size_t i = 0; i < kLimbs; ++i) {
      const std::size_t from_top = kLimbs - i;
      value.limbs[i] = from_top <= count ? limbs[count - from_top] : 0;
    }
    value.rest = 0;
    if (count > kLimbs) {
      value.rest = limbs[count - kLimbs - 1];
      for (std::size_t i = 0; i + kLimbs + 1 < count; ++i) {
        sticky = sticky || limbs[i] != 0;
      }
    }
    if (sticky) value.rest |= 1;
    // At or above 2^(kEmax + 1) every mode overflows. Far below the subnormals
    // the exponent is clamped, which keeps it within int32_t and changes
    // nothing, as with scalbn.
    if (exponent > kEmax) return Overflow(negative, env);
    const std::int64_t limit = -(std::int64_t{1} << (E + 1)) - P;
    if (exponent < limit) exponent = limit;
    return Finish(negative, static_cast<std::int32_t>(exponent), value, env);
  }

  template <int P2, int E2, class Env>
  static Float ConvertFormat(const Float<P2, E2>& x, const Env& env) noexcept {
    using Source = Float<P2, E2>;
    const bool negative = x.negative_ != 0;
    if (x.is_nan()) {
      if (x.is_signaling()) env.raise(kFlagInvalid);
      // Both significands are left-aligned, so the payload keeps its top bits.
      Limbs limbs;
      for (std::size_t i = 0; i < kLimbs; ++i) {
        const std::size_t from_top = kLimbs - i;
        limbs[i] = from_top <= Source::kLimbs
                       ? x.limbs_[Source::kLimbs - from_top]
                       : 0;
      }
      if constexpr (kSpareBits != 0) limbs[0] &= ~LowMask(kSpareBits);
      SetQuietBit(limbs);
      return FromParts(negative, kExpNaN, limbs);
    }
    if (x.is_inf()) return Infinity(negative);
    if (x.is_zero()) return Zero(negative);
    return RoundSignificand(negative, x.Exponent(), x.limbs_.data(),
                            Source::kLimbs, false, env);
  }

  // From a built-in arithmetic type.
  template <detail::Arithmetic T, class Env>
  static Float ConvertFrom(T value, const Env& env) noexcept {
    if constexpr (std::is_integral_v<T>) {
      static_assert(sizeof(T) <= sizeof(limb_t), "Integers up to 64 bits");
      bool negative = false;
      limb_t magnitude;
      if constexpr (std::is_signed_v<T>) {
        negative = value < 0;
        // Negating modulo 2^64 handles the most negative value too.
        magnitude = static_cast<limb_t>(static_cast<std::int64_t>(value));
        if (negative) magnitude = limb_t{0} - magnitude;
      } else {
        magnitude = static_cast<limb_t>(value);
      }
      if (magnitude == 0) return Zero(false);
      const int leading = countl_zero(magnitude);
      const limb_t normalized = magnitude << leading;
      return RoundSignificand(negative, kLimbBits - 1 - leading, &normalized,
                              1, false, env);
    } else {
      return ConvertFormat(DecodeNative(value), env);
    }
  }

  // The exact Float equivalent of a built-in floating-point value.
  template <detail::BinaryFloatingPoint T>
  static detail::NativeFormat<T> DecodeNative(T value) noexcept {
    using Native = detail::NativeFormat<T>;
    if constexpr (detail::NativeFormatOf<T>::kInterchange) {
      using Unsigned = detail::SameSizeUnsigned<T>;
      return Native::from_bits(
          {static_cast<limb_t>(std::bit_cast<Unsigned>(value))});
    } else {
      // Such as the x87 extended type, with its explicit leading bit, or a
      // long double that is binary128. Payloads are not portable here, so a
      // NaN becomes the default NaN.
      const bool negative = std::signbit(value);
      if (std::isnan(value)) return Native::QuietNaN(negative);
      if (std::isinf(value)) return Native::Infinity(negative);
      if (value == 0) return Native::Zero(negative);
      int exponent = 0;
      // Exact: m is in [1/2, 1), and splitting off whole limbs loses nothing.
      T m = std::frexp(std::fabs(value), &exponent);
      typename Native::Limbs limbs;
      for (std::size_t i = Native::kLimbs; i-- > 0;) {
        m = std::ldexp(m, kLimbBits);
        const limb_t chunk = static_cast<limb_t>(m);
        m -= static_cast<T>(chunk);
        limbs[i] = chunk;
      }
      return Native::FromParts(negative, exponent - 1, limbs);
    }
  }

  // To a built-in floating-point type.
  template <detail::BinaryFloatingPoint T, class Env>
  T ConvertTo(const Env& env) const noexcept {
    using Native = detail::NativeFormat<T>;
    const Native x = Native::ConvertFormat(*this, env);
    if constexpr (detail::NativeFormatOf<T>::kInterchange) {
      using Unsigned = detail::SameSizeUnsigned<T>;
      return std::bit_cast<T>(static_cast<Unsigned>(x.to_bits()[0]));
    } else {
      const bool negative = x.signbit();
      T result;
      if (x.is_nan()) {
        result = std::numeric_limits<T>::quiet_NaN();
      } else if (x.is_inf()) {
        result = std::numeric_limits<T>::infinity();
      } else if (x.is_zero()) {
        result = 0;
      } else {
        // Every partial sum lies on T's grid, so the additions are exact.
        result = 0;
        for (std::size_t i = Native::kLimbs; i-- > 0;) {
          const int scale =
              x.Exponent() - (kLimbBits - 1) -
              static_cast<int>(Native::kLimbs - 1 - i) * kLimbBits;
          result += std::ldexp(static_cast<T>(x.limbs_[i]), scale);
        }
      }
      return negative ? -result : result;
    }
  }

  // Rounds in `mode`, then converts to I. NaN and values outside I's range
  // raise invalid and saturate, with NaN giving 0. Inexact is raised only when
  // `exact` is set.
  template <detail::Integer I, class Env>
  I ToInteger(Rounding mode, bool exact, const Env& env) const noexcept {
    static_assert(sizeof(I) <= sizeof(limb_t), "Integers up to 64 bits");
    using Limits = std::numeric_limits<I>;
    const bool negative = negative_ != 0;
    const I saturated = negative ? (Limits::min)() : (Limits::max)();
    if (is_nan()) {
      env.raise(kFlagInvalid);
      return 0;
    }
    Context local{mode, 0};
    const Float rounded = RoundToIntegralImpl(*this, detail::ContextEnv{local});
    if (rounded.is_inf()) {
      env.raise(kFlagInvalid);
      return saturated;
    }
    limb_t magnitude = 0;
    if (!rounded.is_zero()) {
      if (rounded.Exponent() >= kLimbBits) {
        env.raise(kFlagInvalid);
        return saturated;
      }
      // An integer below 2^64 lies entirely in the top limb. A nonzero
      // integer has an exponent of at least 0, so the shift is below 64.
      // NOLINTBEGIN(clang-analyzer-core.BitwiseShift)
      magnitude =
          rounded.limbs_[kLimbs - 1] >> (kLimbBits - 1 - rounded.Exponent());
      // NOLINTEND(clang-analyzer-core.BitwiseShift)
      const limb_t limit =
          static_cast<limb_t>((Limits::max)()) + (negative ? 1 : 0);
      if ((negative && !Limits::is_signed) || magnitude > limit) {
        env.raise(kFlagInvalid);
        return saturated;
      }
    }
    if (exact) env.raise(local.flags & kFlagInexact);
    return static_cast<I>(negative ? limb_t{0} - magnitude : magnitude);
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
    if (a.Exponent() != b.Exponent()) {
      order = a.Exponent() < b.Exponent() ? -1 : 1;
    } else {
      order = detail::CompareLimbs(a.limbs_, b.limbs_);
    }
    if (a.negative_ != 0) order = -order;
    if (order == 0) return std::partial_ordering::equivalent;
    return order < 0 ? std::partial_ordering::less
                     : std::partial_ordering::greater;
  }

  Limbs limbs_;
  // The exponent is stored offset by 2^31, so that kExpZero is stored as 0
  // and a Float whose bytes are all zero is +0: memset and calloc give zeros,
  // as libraries such as Eigen assume of trivially copyable scalars.
  std::uint32_t biased_exponent_;
  std::uint32_t negative_;
};

using binary128 = Float<113, 15>;
using binary256 = Float<237, 19>;
using binary512 = Float<489, 23>;

namespace detail {

struct FloatAccess {
  template <class To, class From, class Env>
  static To Convert(const From& x, const Env& env) noexcept {
    if constexpr (kIsFloat<To>) {
      if constexpr (kIsFloat<From>) {
        return To::ConvertFormat(x, env);
      } else {
        return To::ConvertFrom(x, env);
      }
    } else {
      return x.template ConvertTo<To>(env);
    }
  }

  template <class I, int P, int E, class Env>
  static I ToInteger(const Float<P, E>& x, Rounding mode, bool exact,
                     const Env& env) noexcept {
    return x.template ToInteger<I>(mode, exact, env);
  }

  // See Float::RoundSignificand.
  template <class F, class Env>
  static F Round(bool negative, std::int64_t exponent, const limb_t* limbs,
                 std::size_t count, bool sticky, const Env& env) noexcept {
    return F::RoundSignificand(negative, exponent, limbs, count, sticky, env);
  }

  // A NaN operand, quieted, raising invalid if it was signaling.
  template <class F, class Env>
  static F PropagateNaN(const F& x, const Env& env) noexcept {
    return F::PropagateNaN(x, x, env);
  }

  // The first NaN of two operands, at least one a NaN, quieted, raising
  // invalid if either was signaling.
  template <class F, class Env>
  static F PropagateNaN(const F& x, const F& y, const Env& env) noexcept {
    return F::PropagateNaN(x, y, env);
  }
};

template <class To, class From>
concept Convertible =
    (kIsFloat<To> && (kIsFloat<From> || Arithmetic<From>)) ||
    (BinaryFloatingPoint<To> && kIsFloat<From>);

}  // namespace detail

// Converts between formats, and between a Float and a built-in arithmetic
// type, rounding in context.rounding and raising flags. Built-in targets must
// be floating-point; integers go through to_integer().
template <class To, class From>
  requires detail::Convertible<To, From>
To convert(const From& x, Context& context) noexcept {
  return detail::FloatAccess::Convert<To>(x, detail::ContextEnv{context});
}

// As above, rounding to nearest-even without flags.
template <class To, class From>
  requires detail::Convertible<To, From>
To convert(const From& x) noexcept {
  return detail::FloatAccess::Convert<To>(x, detail::NearestEnv{});
}

// IEEE 754 convertToInteger: rounds x in `mode`, then converts to I. NaN and
// values outside I's range raise invalid and saturate, with NaN giving 0. A
// rounded result raises no inexact flag.
template <detail::Integer I, int P, int E>
I to_integer(const Float<P, E>& x, Rounding mode, Context& context) noexcept {
  return detail::FloatAccess::ToInteger<I>(x, mode, false,
                                           detail::ContextEnv{context});
}

// As above, without flags.
template <detail::Integer I, int P, int E>
I to_integer(const Float<P, E>& x, Rounding mode) noexcept {
  return detail::FloatAccess::ToInteger<I>(x, mode, false,
                                           detail::NearestEnv{});
}

// IEEE 754 convertToIntegerExact: like to_integer, but an inexact result
// raises inexact.
template <detail::Integer I, int P, int E>
I to_integer_exact(const Float<P, E>& x, Rounding mode,
                   Context& context) noexcept {
  return detail::FloatAccess::ToInteger<I>(x, mode, true,
                                           detail::ContextEnv{context});
}

}  // namespace mpfix

// The limits follow the conventions of the built-in types: min() is the
// smallest normal value, and min_exponent and max_exponent are one more than
// the IEEE exponents, as for a significand in [0.5, 1). The values are not
// constexpr, since Float's operations are not.
namespace std {

// Specializing numeric_limits for a program-defined type is allowed.
template <int P, int E>
// NOLINTNEXTLINE(bugprone-std-namespace-modification)
class numeric_limits<mpfix::Float<P, E>> {
  using F = mpfix::Float<P, E>;

 public:
  static constexpr bool is_specialized = true;
  static constexpr bool is_signed = true;
  static constexpr bool is_integer = false;
  static constexpr bool is_exact = false;
  static constexpr bool has_infinity = true;
  static constexpr bool has_quiet_NaN = true;
  static constexpr bool has_signaling_NaN = true;
  static constexpr bool has_denorm_loss = false;
  static constexpr std::float_round_style round_style = std::round_to_nearest;
  static constexpr bool is_iec559 = true;
  static constexpr bool is_bounded = true;
  static constexpr bool is_modulo = false;
  static constexpr int radix = 2;
  static constexpr int digits = P;
  // floor((P - 1) log10 2) and ceil(P log10 2) + 1, from log10 2 < 0.30103.
  static constexpr int digits10 = static_cast<int>((P - 1) * 30103LL / 100000);
  static constexpr int max_digits10 =
      static_cast<int>((P * 30103LL + 99999) / 100000) + 1;
  static constexpr int min_exponent = F::kEmin + 1;
  static constexpr int max_exponent = F::kEmax + 1;
  // The smallest n with 10^n normal, and the largest with 10^n finite: for
  // double, -307 and 308.
  static constexpr int min_exponent10 =
      -static_cast<int>((-F::kEmin) * 30103LL / 100000);
  static constexpr int max_exponent10 =
      static_cast<int>((F::kEmax + 1) * 30103LL / 100000);
  static constexpr bool traps = false;
  // Tininess is detected after rounding.
  static constexpr bool tinyness_before = false;

  // Parenthesized against the min and max macros of <windows.h>. All are
  // constexpr, as for the built-in types, which some standard library code
  // relies on.
  static constexpr F(min)() noexcept { return F::MinNormal(); }
  static constexpr F(max)() noexcept { return F::MaxFinite(); }
  static constexpr F lowest() noexcept { return F::MaxFinite(true); }
  // 2^(1 - P) and 1/2.
  static constexpr F epsilon() noexcept { return F::PowerOfTwo(1 - P); }
  static constexpr F round_error() noexcept { return F::PowerOfTwo(-1); }
  static constexpr F infinity() noexcept { return F::Infinity(); }
  static constexpr F quiet_NaN() noexcept { return F::QuietNaN(); }
  static constexpr F signaling_NaN() noexcept { return F::SignalingNaN(); }
  static constexpr F denorm_min() noexcept { return F::MinSubnormal(); }
};

}  // namespace std

#endif  // MPFIX_FLOAT_H_
