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
// Operators round to nearest-even and do not track flags, which is the fast
// path. add(), sub() and mul() take a Context for all five rounding modes and
// all five exception flags. Tininess is detected after rounding.
//
// NaN rules: an operation with NaN operands returns the first NaN operand,
// quieted. An invalid operation returns the positive default quiet NaN.

#ifndef HWYMP_FLOAT_H_
#define HWYMP_FLOAT_H_

#include <array>
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

constexpr int LimbsFor(int bits) { return (bits + kLimbBits - 1) / kLimbBits; }

// Decides whether a discarded fraction rounds the magnitude up. `round_bit` is
// the most significant discarded bit and `sticky` is the OR of the rest.
HWYMP_ALWAYS_INLINE constexpr bool ShouldRoundUp(Rounding mode, bool negative,
                                                 bool is_odd, bool round_bit,
                                                 bool sticky) noexcept {
  switch (mode) {
    case Rounding::kNearestEven:
      return round_bit && (sticky || is_odd);
    case Rounding::kNearestAway:
      return round_bit;
    case Rounding::kTowardZero:
      return false;
    case Rounding::kUp:
      return !negative && (round_bit || sticky);
    case Rounding::kDown:
      return negative && (round_bit || sticky);
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
  Float& operator*=(const Float& other) noexcept {
    return *this = *this * other;
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
    if (!detail::ShouldRoundUp(env.rounding(), negative, is_odd, round_bit,
                               sticky)) {
      return false;
    }
    carry_t carry = 0;
    value.limbs[0] = adc(value.limbs[0], limb_t{1} << kSpareBits, carry);
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
      const int count = static_cast<int>(discard);
      const auto bit_at = [&limbs](int position) {
        return ((limbs[static_cast<std::size_t>(position) / kLimbBits] >>
                 (position % kLimbBits)) &
                1) != 0;
      };
      round_bit = bit_at(count - 1);
      if (count < kTotalBits) is_odd = bit_at(count);
      Limbs below = limbs;
      detail::ShiftLeft(below, kTotalBits - (count - 1) - 1);
      // `below` now starts with the round bit; drop it to keep the rest.
      below[kLimbs - 1] &= ~kTopBit;
      sticky = sticky || !detail::IsZero(below);
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
    if (round_up) {
      const std::size_t word = static_cast<std::size_t>(count) / kLimbBits;
      const limb_t increment = limb_t{1} << (count % kLimbBits);
      carry_t carry = 0;
      limbs[word] = adc(limbs[word], increment, carry);
      for (std::size_t i = word + 1; i < kLimbs; ++i) {
        limbs[i] = adc(limbs[i], 0, carry);
      }
      if (carry != 0) {
        limbs[kLimbs - 1] = kTopBit;
        ++exponent;
      }
    }
    return FromParts(negative, exponent, limbs);
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
    if (bits != 0) {
      lost |= wide[0] << (kLimbBits - bits);
      detail::ShiftRight(wide, bits);
    }
    Extended result;
    for (std::size_t i = 0; i < kLimbs; ++i) result.limbs[i] = wide[i + 1];
    result.rest = wide[0] | (lost != 0 ? 1 : 0);
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
      if (carry != 0) {
        sum.rest = (sum.limbs[0] << (kLimbBits - 1)) | (sum.rest >> 1) |
                   (sum.rest & 1);
        detail::ShiftRight(sum.limbs, 1);
        sum.limbs[kLimbs - 1] |= kTopBit;
        ++exponent;
      }
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
    if ((difference.limbs[kLimbs - 1] & kTopBit) == 0) {
      exponent -= NormalizeDifference(difference);
    }
    return Finish(negative, exponent, difference, env);
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

    std::array<limb_t, 2 * kLimbs> product{};
    for (std::size_t i = 0; i < kLimbs; ++i) {
      limb_t previous_high = 0;
      for (std::size_t j = 0; j < kLimbs; ++j) {
        limb_t high;
        limb_t low = mul_wide(a.limbs_[j], b.limbs_[i], high);
        carry_t carry = 0;
        low = adc(low, product[i + j], carry);
        high = adc(high, 0, carry);
        carry = 0;
        low = adc(low, previous_high, carry);
        high = adc(high, 0, carry);
        product[i + j] = low;
        previous_high = high;
      }
      product[i + kLimbs] = previous_high;
    }

    Extended value;
    for (std::size_t i = 0; i < kLimbs; ++i) {
      value.limbs[i] = product[i + kLimbs];
    }
    value.rest = product[kLimbs - 1];
    limb_t lost = 0;
    for (std::size_t i = 0; i + 1 < kLimbs; ++i) lost |= product[i];

    // Both significands are in [1, 2), so the product is in [1, 4).
    std::int32_t exponent = a.exponent_ + b.exponent_ + 1;
    if ((value.limbs[kLimbs - 1] & kTopBit) == 0) {
      --exponent;
      detail::ShiftLeft(value.limbs, 1);
      value.limbs[0] |= value.rest >> (kLimbBits - 1);
      value.rest <<= 1;
    }
    value.rest |= lost != 0 ? 1 : 0;
    return Finish(negative, exponent, value, env);
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
