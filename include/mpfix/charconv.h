// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Correctly rounded conversions between Float and character strings, in
// decimal and hexadecimal, in the manner of <charconv>, plus stream operators.
//
// All conversions are exact or correctly rounded, for any number of digits:
// they work with exact big integers. A conversion near the ends of a wide
// exponent range handles numbers of millions of bits; binary512's largest
// value takes a few milliseconds.
//
// Shortest output is the shortest decimal string that reads back to the same
// value under round-to-nearest-even, and the nearest such string if there are
// several. Output with a precision rounds in the Context's mode, if one is
// given, and to nearest-even otherwise.
//
// Differences from <charconv>: from_chars stores infinity, zero or a subnormal
// for values out of range, as IEEE 754 requires, and reports success; the
// Context flags tell overflow and underflow apart. Hexadecimal output always
// has a leading 1, subnormals included.

#ifndef MPFIX_CHARCONV_H_
#define MPFIX_CHARCONV_H_

#include <algorithm>
#include <bit>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <istream>
#include <ostream>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>
#include <vector>

#include "mpfix/bigint.h"
#include "mpfix/float.h"

namespace mpfix {
namespace detail {

inline constexpr double kLog10Of2 = 0.301029995663981195;

inline bool HasFormat(std::chars_format fmt, std::chars_format part) {
  return (fmt & part) == part;
}

inline char ToLower(char c) {
  return c >= 'A' && c <= 'Z' ? static_cast<char>(c - 'A' + 'a') : c;
}

// Whether [p, last) starts with `word`, ignoring case.
inline bool StartsWithIgnoringCase(const char* p, const char* last,
                                   std::string_view word) {
  if (static_cast<std::size_t>(last - p) < word.size()) return false;
  for (std::size_t i = 0; i < word.size(); ++i) {
    if (ToLower(p[i]) != word[i]) return false;
  }
  return true;
}

inline int DigitValue(char c, bool hex) {
  if (c >= '0' && c <= '9') return c - '0';
  if (hex) {
    const char lower = ToLower(c);
    if (lower >= 'a' && lower <= 'f') return lower - 'a' + 10;
  }
  return -1;
}

// Saturates far beyond any exponent that matters, so that sums of a few such
// values cannot overflow.
inline constexpr std::int64_t kExponentLimit = std::int64_t{1} << 60;

inline std::int64_t SaturatingAdd(std::int64_t a, std::int64_t b) {
  const std::int64_t sum = a + b;
  return std::clamp(sum, -kExponentLimit, kExponentLimit);
}

// Collects the significant digits of a number into a big integer, dropping
// leading and trailing zeros. The value is digits * base^adjust.
class DigitAccumulator {
 public:
  explicit DigitAccumulator(bool hex) : hex_(hex) {}

  void Push(int digit, bool fractional) {
    if (fractional) --adjust_;
    if (digit == 0) {
      if (count_ > 0) ++pending_zeros_;
      return;
    }
    for (; pending_zeros_ > 0; --pending_zeros_) Append(0);
    Append(digit);
  }

  // Finishes and returns the digits as an integer.
  BigUint& Finish() {
    if (chunk_digits_ > 0) digits_.MulAdd(Power(chunk_digits_), chunk_);
    chunk_digits_ = 0;
    chunk_ = 0;
    adjust_ += pending_zeros_;
    pending_zeros_ = 0;
    return digits_;
  }

  std::int64_t adjust() const { return adjust_; }
  // The number of significant digits, internal zeros included.
  std::int64_t count() const { return count_; }

 private:
  // 10^19 and 16^15 are the largest powers that fit a limb.
  int ChunkSize() const { return hex_ ? 15 : 19; }

  limb_t Power(int digits) const {
    limb_t power = 1;
    for (int i = 0; i < digits; ++i) power *= hex_ ? 16 : 10;
    return power;
  }

  void Append(int digit) {
    chunk_ = chunk_ * (hex_ ? 16 : 10) + static_cast<limb_t>(digit);
    ++count_;
    if (++chunk_digits_ == ChunkSize()) {
      digits_.MulAdd(Power(chunk_digits_), chunk_);
      chunk_digits_ = 0;
      chunk_ = 0;
    }
  }

  bool hex_;
  BigUint digits_;
  limb_t chunk_ = 0;
  int chunk_digits_ = 0;
  std::int64_t count_ = 0;
  std::int64_t pending_zeros_ = 0;
  std::int64_t adjust_ = 0;
};

// The value of eight decimal digits, the first in the lowest byte of `chars`.
// Pairs, then quadruples, of digits combine in parallel within the word.
inline limb_t EightDigits(limb_t chars) {
  chars -= 0x3030303030303030u;
  // Byte 2i becomes the two-digit number of digits 2i and 2i + 1.
  chars = chars * 10 + (chars >> 8);
  constexpr limb_t kPairs = 0x000000FF000000FFu;
  // The high halves of these products sum to the eight-digit number.
  const limb_t outer = (chars & kPairs) * (100 + (limb_t{1000000} << 32));
  const limb_t inner = ((chars >> 16) & kPairs) * (1 + (limb_t{10000} << 32));
  return (outer + inner) >> 32;
}

// A decimal number with leading and trailing zeros dropped: value *
// 10^exponent, where value has `count` digits.
struct DecimalDigits {
  BigUint value;
  std::int64_t count;
  std::int64_t exponent;
};

// Reads the integer digits [integer_begin, integer_end) and the fraction digits
// [fraction_begin, fraction_end), which must all be decimal digits.
inline DecimalDigits ReadDecimalDigits(const char* integer_begin,
                                       const char* integer_end,
                                       const char* fraction_begin,
                                       const char* fraction_end) {
  // The digits form one sequence across the point, worth that integer times
  // 10^-(fraction digits). Only its part from the first to the last nonzero
  // digit is kept.
  const char* begin[2] = {integer_begin, fraction_begin};
  const char* end[2] = {integer_end, fraction_end};
  DecimalDigits result{BigUint(), 0, -(fraction_end - fraction_begin)};
  while (begin[0] != end[0] && *begin[0] == '0') ++begin[0];
  if (begin[0] == end[0]) {
    while (begin[1] != end[1] && *begin[1] == '0') ++begin[1];
  }
  while (end[1] != begin[1] && end[1][-1] == '0') {
    --end[1];
    ++result.exponent;
  }
  if (begin[1] == end[1]) {
    while (end[0] != begin[0] && end[0][-1] == '0') {
      --end[0];
      ++result.exponent;
    }
  }

  constexpr limb_t kPowers[20] = {1u,
                                  10u,
                                  100u,
                                  1000u,
                                  10000u,
                                  100000u,
                                  1000000u,
                                  10000000u,
                                  100000000u,
                                  1000000000u,
                                  10000000000u,
                                  100000000000u,
                                  1000000000000u,
                                  10000000000000u,
                                  100000000000000u,
                                  1000000000000000u,
                                  10000000000000000u,
                                  100000000000000000u,
                                  1000000000000000000u,
                                  10000000000000000000u};
  limb_t chunk = 0;
  int chunk_digits = 0;
  for (int i = 0; i < 2; ++i) {
    const char* q = begin[i];
    result.count += end[i] - q;
    if constexpr (std::endian::native == std::endian::little) {
      while (end[i] - q >= 8) {
        if (chunk_digits > 11) {
          result.value.MulAdd(kPowers[chunk_digits], chunk);
          chunk = 0;
          chunk_digits = 0;
        }
        limb_t chars;
        std::memcpy(&chars, q, 8);
        chunk = chunk * 100000000u + EightDigits(chars);
        chunk_digits += 8;
        q += 8;
        // A full chunk takes no more digits.
        if (chunk_digits == 19) {
          result.value.MulAdd(kPowers[19], chunk);
          chunk = 0;
          chunk_digits = 0;
        }
      }
    }
    for (; q != end[i]; ++q) {
      chunk = chunk * 10 + static_cast<limb_t>(*q - '0');
      if (++chunk_digits == 19) {
        result.value.MulAdd(kPowers[19], chunk);
        chunk = 0;
        chunk_digits = 0;
      }
    }
  }
  if (chunk_digits > 0) result.value.MulAdd(kPowers[chunk_digits], chunk);
  return result;
}

// A nonzero integer shifted so that its top bit is set.
inline BigUint Normalized(BigUint value) {
  value.ShiftLeft(static_cast<std::size_t>(countl_zero(value.limbs().back())));
  return value;
}

// Rounds the nonzero value `integer` * 2^exponent2.
template <class F, class Env>
F RoundInteger(bool negative, const BigUint& integer, std::int64_t exponent2,
               bool sticky, const Env& env) {
  const BigUint normalized = Normalized(integer);
  const std::int64_t top = SaturatingAdd(
      static_cast<std::int64_t>(integer.bit_length()) - 1, exponent2);
  return FloatAccess::Round<F>(negative, top, normalized.limbs().data(),
                               normalized.limbs().size(), sticky, env);
}

// The significand of a finite nonzero x as an integer f with x = f * 2^e,
// where e is returned.
template <int P, int E>
BigUint Significand(const Float<P, E>& x, std::int64_t& e) {
  const auto& limbs = x.raw_limbs();
  BigUint f = BigUint::FromLimbs(limbs.data(), limbs.size());
  constexpr int kSpare = static_cast<int>(Float<P, E>::kLimbs) * kLimbBits - P;
  f.ShiftRight(kSpare);
  e = std::int64_t{x.raw_exponent()} - (P - 1);
  return f;
}

// The odd part of the significand of a finite nonzero x, with x = m * 2^e.
template <int P, int E>
BigUint OddSignificand(const Float<P, E>& x, std::int64_t& e) {
  BigUint m = Significand(x, e);
  const std::size_t zeros = m.trailing_zeros();
  m.ShiftRight(zeros);
  e += static_cast<std::int64_t>(zeros);
  return m;
}

inline constexpr double kLog2Of5 = 2.32192809488736234787;

// --- Fast paths ---
//
// Interval arithmetic in a wider format, with each operation rounded down for
// the lower bound and up for the upper one, brackets the exact value of a
// conversion. When the whole interval rounds to the same result, that result
// is correct. When it does not, which happens only very close to a rounding
// boundary, the exact path decides. Exactness, and so the inexact flag, is
// settled by number theory on the operands, which is cheap.

// The working formats have an exponent range that holds every power of ten
// that F's conversions need. The first has at least 16 bits more than F, which
// decides nearly every case at the least cost; the second at least 64 more.
template <class F>
using NarrowWorking =
    Float<64 * ((F::kPrecision + 16 + kLimbBits - 1) / kLimbBits), 29>;
template <class F>
using Working = Float<static_cast<int>(64 * (F::kLimbs + 1)), 29>;
template <class F>
inline constexpr bool kHasNarrowWorking =
    !std::is_same_v<NarrowWorking<F>, Working<F>>;

template <class F>
inline constexpr bool kHasFastPath = F::kExponentBits <= 26;

// The fast paths take powers of ten below 10^(2^24).
inline constexpr std::int64_t kFastPowerLimit = std::int64_t{1} << 24;

// Lower and upper bounds on a positive value.
template <class T>
struct Bounds {
  T lower;
  T upper;
};

// Bounds on powers of ten and their reciprocals. 10^n for |n| < 2048 is the
// product of two table entries, 10^(32 a) and 10^b, which costs one
// multiplication per bound. Larger n add factors of 10^(2^k).
template <class T>
struct PowerOfTenTable {
  static constexpr int kUnits = 32;
  static constexpr int kThirtyTwos = 64;
  static constexpr int kFirstLevel = 11;  // 2^11 = kUnits * kThirtyTwos
  static constexpr int kLevels = 24;

  // [0] for 10^n and [1] for 10^-n.
  Bounds<T> units[2][kUnits];
  Bounds<T> thirty_twos[2][kThirtyTwos];
  Bounds<T> levels[2][kLevels];

  PowerOfTenTable() {
    // The entries come from long chains of products, so they are computed
    // 64 bits wider and rounded once. Each bound is then within an ulp or two.
    using Wide = Float<T::kPrecision + 64, 29>;
    Context down{Rounding::kDown, 0};
    Context up{Rounding::kUp, 0};
    const Wide one(1);
    const Wide ten(10);
    const auto reciprocal = [&](const Bounds<Wide>& x) {
      return Bounds<Wide>{div(one, x.upper, down), div(one, x.lower, up)};
    };
    const auto product = [&](const Bounds<Wide>& a, const Bounds<Wide>& b) {
      return Bounds<Wide>{mul(a.lower, b.lower, down),
                          mul(a.upper, b.upper, up)};
    };
    const auto narrow = [&](const Bounds<Wide>& x) {
      return Bounds<T>{convert<T>(x.lower, down), convert<T>(x.upper, up)};
    };
    const auto fill = [&](Bounds<T>(&target)[2], const Bounds<Wide>& x) {
      target[0] = narrow(x);
      target[1] = narrow(reciprocal(x));
    };
    Bounds<T> pair[2];
    // 10^31 has 72 significant bits, so the units powers are exact.
    Bounds<Wide> power{one, one};
    for (int b = 0; b < kUnits; ++b) {
      fill(pair, power);
      units[0][b] = pair[0];
      units[1][b] = pair[1];
      power = product(power, {ten, ten});
    }
    const Bounds<Wide> step = power;  // 10^32
    power = {one, one};
    for (int a = 0; a < kThirtyTwos; ++a) {
      fill(pair, power);
      thirty_twos[0][a] = pair[0];
      thirty_twos[1][a] = pair[1];
      power = product(power, step);
    }
    power = {ten, ten};
    for (int k = 0; k < kLevels; ++k) {
      fill(pair, power);
      levels[0][k] = pair[0];
      levels[1][k] = pair[1];
      power = product(power, power);
    }
  }

  static const PowerOfTenTable& Get() {
    static const PowerOfTenTable table;
    return table;
  }
};

// Bounds on 10^n, for |n| < kFastPowerLimit.
template <class T>
void PowerOfTenBounds(std::int64_t n, T& lower, T& upper) {
  using Table = PowerOfTenTable<T>;
  const Table& table = Table::Get();
  const int sign = n < 0 ? 1 : 0;
  const std::int64_t m = n < 0 ? -n : n;
  const Bounds<T>& units = table.units[sign][m % Table::kUnits];
  const Bounds<T>& thirty_twos =
      table.thirty_twos[sign][(m / Table::kUnits) % Table::kThirtyTwos];
  Context down{Rounding::kDown, 0};
  Context up{Rounding::kUp, 0};
  if (m % Table::kUnits == 0) {
    lower = thirty_twos.lower;
    upper = thirty_twos.upper;
  } else if (m / Table::kUnits % Table::kThirtyTwos == 0) {
    lower = units.lower;
    upper = units.upper;
  } else {
    lower = mul(thirty_twos.lower, units.lower, down);
    upper = mul(thirty_twos.upper, units.upper, up);
  }
  std::int64_t high = m >> Table::kFirstLevel;
  for (int k = Table::kFirstLevel; high != 0; ++k, high >>= 1) {
    if ((high & 1) != 0) {
      lower = mul(lower, table.levels[sign][k].lower, down);
      upper = mul(upper, table.levels[sign][k].upper, up);
    }
  }
}

// The integer part of a nonnegative finite value below 2^(64 L).
template <int P, int E>
BigUint IntegerPart(const Float<P, E>& x) {
  if (x.is_zero() || x.raw_exponent() < 0) return BigUint();
  const auto& limbs = x.raw_limbs();
  BigUint value = BigUint::FromLimbs(limbs.data(), limbs.size());
  value.ShiftRight(static_cast<std::size_t>(
      static_cast<std::int64_t>(limbs.size()) * kLimbBits - 1 -
      x.raw_exponent()));
  return value;
}

// The exponent of the lowest set bit of a finite nonzero x: x = m * 2^e with
// m odd.
template <int P, int E>
std::int64_t LowestBitExponent(const Float<P, E>& x) {
  const auto& limbs = x.raw_limbs();
  std::size_t i = 0;
  while (limbs[i] == 0) ++i;
  const auto zeros =
      static_cast<std::int64_t>(i * kLimbBits) + countr_zero(limbs[i]);
  return std::int64_t{x.raw_exponent()} -
         (static_cast<std::int64_t>(limbs.size()) * kLimbBits - 1) + zeros;
}

// Whether digits * 10^t equals x, for nonzero digits and a finite nonzero x.
// The powers of two must agree first, which rules out almost every inexact
// case before any big multiplication.
template <int P, int E>
bool EqualsDecimal(const Float<P, E>& x, const BigUint& digits,
                   std::int64_t t) {
  const std::int64_t e = LowestBitExponent(x);
  const auto d_zeros = static_cast<std::int64_t>(digits.trailing_zeros());
  // digits * 10^t = d * 5^t * 2^(d_zeros + t), with d odd.
  if (d_zeros + t != e) return false;
  std::int64_t unused;
  BigUint m = OddSignificand(x, unused);
  BigUint d = digits;
  d.ShiftRight(static_cast<std::size_t>(d_zeros));
  if (t >= 0) {
    // m < 2^P, so it cannot be a multiple of a larger power of five.
    if (static_cast<double>(t) * kLog2Of5 > P) return false;
    return Compare(d * BigUint::Pow5(static_cast<std::uint64_t>(t)), m) == 0;
  }
  // Then d = m * 5^u for u = -t.
  const std::int64_t u = -t;
  if (static_cast<double>(u) * kLog2Of5 > static_cast<double>(d.bit_length())) {
    return false;
  }
  return Compare(m * BigUint::Pow5(static_cast<std::uint64_t>(u)), d) == 0;
}

// Whether |x| / 10^t is an integer, for finite nonzero x.
template <int P, int E>
bool IsMultipleOfPowerOfTen(const Float<P, E>& x, std::int64_t t) {
  const std::int64_t e = LowestBitExponent(x);
  if (t < 0) return e - t >= 0;  // m * 5^-t * 2^(e - t)
  if (e < t) return false;
  std::int64_t unused;
  const BigUint m = OddSignificand(x, unused);
  if (static_cast<double>(t) * kLog2Of5 > P) return false;
  BigUint quotient;
  BigUint remainder;
  BigUint::DivRem(m, BigUint::Pow5(static_cast<std::uint64_t>(t)), quotient,
                  remainder);
  return remainder.is_zero();
}

// Rounds digits * 10^t on the fast path, working in T. Returns false if it
// cannot decide.
template <class T, class F, class Env>
bool RoundDecimalFast(bool negative, const BigUint& digits, std::int64_t t,
                      const Env& env, F& result) {
  if constexpr (!kHasFastPath<F>) {
    return false;
  } else {
    if (t >= kFastPowerLimit || t <= -kFastPowerLimit) return false;
    const BigUint normalized = Normalized(digits);
    const LimbVector& limbs = normalized.limbs();
    const auto top = static_cast<std::int64_t>(digits.bit_length()) - 1;
    const T d_low =
        FloatAccess::Round<T>(false, top, limbs.data(), limbs.size(), false,
                              SilentEnv{Rounding::kDown});
    const T d_high =
        limbs.size() <= T::kLimbs
            ? d_low
            : FloatAccess::Round<T>(false, top, limbs.data(), limbs.size(),
                                    false, SilentEnv{Rounding::kUp});
    T p_low;
    T p_high;
    PowerOfTenBounds(t, p_low, p_high);
    Context down{Rounding::kDown, 0};
    Context up{Rounding::kUp, 0};
    const T low = mul(d_low, p_low, down);
    const T high = mul(d_high, p_high, up);

    const T signed_low = negative ? -high : low;
    const T signed_high = negative ? -low : high;
    const SilentEnv rounding{env.rounding()};

    // Rounded to P bits with a wider exponent range, which decides overflow
    // and tininess after rounding, as IEEE 754 defines them.
    using Unbounded = Float<F::kPrecision, F::kExponentBits + 1>;
    const Unbounded wide =
        FloatAccess::Convert<Unbounded>(signed_low, rounding);
    if (wide.to_bits() !=
        FloatAccess::Convert<Unbounded>(signed_high, rounding).to_bits()) {
      return false;
    }
    const Unbounded magnitude = abs(wide);
    if (magnitude > Unbounded(F::MaxFinite())) {
      // Raises overflow and inexact, and picks infinity or MaxFinite.
      result = FloatAccess::Convert<F>(wide, env);
      return true;
    }
    if (magnitude >= Unbounded(F::MinNormal())) {
      // Normal, so the wide result is exact in F.
      result = FloatAccess::Convert<F>(wide, rounding);
      if (!EqualsDecimal(result, digits, t)) env.raise(kFlagInexact);
      return true;
    }
    // Tiny: rounding to the subnormal grid must be decided on its own, to
    // avoid rounding twice.
    const F a = FloatAccess::Convert<F>(signed_low, rounding);
    if (a.to_bits() !=
        FloatAccess::Convert<F>(signed_high, rounding).to_bits()) {
      return false;
    }
    if (a.is_zero() || !EqualsDecimal(a, digits, t)) {
      env.raise(kFlagInexact | kFlagUnderflow);
    }
    result = a;
    return true;
  }
}

// Rounds digits * 10^exponent10 for nonzero digits with `count` decimal
// digits.
template <class F, class Env>
F RoundDecimal(bool negative, const BigUint& digits, std::int64_t count,
               std::int64_t exponent10, const Env& env) {
  // The value lies in [10^leading, 10^(leading + 1)).
  const std::int64_t leading = SaturatingAdd(count - 1, exponent10);
  const limb_t top_bit = limb_t{1} << (kLimbBits - 1);
  if (static_cast<double>(leading) >
      static_cast<double>(F::kEmax + 1) * kLog10Of2 + 1) {
    // At least 2^(kEmax + 1): every mode overflows.
    return FloatAccess::Round<F>(negative, std::int64_t{F::kEmax} + 1, &top_bit,
                                 1, false, env);
  }
  if (static_cast<double>(leading + 1) <
      static_cast<double>(F::kEmin - F::kPrecision) * kLog10Of2 - 1) {
    // Below half the smallest subnormal: any value there rounds alike.
    return FloatAccess::Round<F>(negative,
                                 std::int64_t{F::kEmin} - F::kPrecision - 4,
                                 &top_bit, 1, true, env);
  }

  F fast;
  if constexpr (kHasNarrowWorking<F>) {
    if (RoundDecimalFast<NarrowWorking<F>>(negative, digits, exponent10, env,
                                           fast)) {
      return fast;
    }
  }
  if (RoundDecimalFast<Working<F>>(negative, digits, exponent10, env, fast)) {
    return fast;
  }

  if (exponent10 >= 0) {
    // digits * 5^e * 2^e is exact.
    const BigUint scaled =
        digits * BigUint::Pow5(static_cast<std::uint64_t>(exponent10));
    return RoundInteger<F>(negative, scaled, exponent10, false, env);
  }

  // digits / 10^u = (digits * 2^s / 5^u) * 2^(-s - u), with s chosen so that
  // the quotient has at least P + 2 bits.
  const std::int64_t u = -exponent10;
  const BigUint divisor = BigUint::Pow5(static_cast<std::uint64_t>(u));
  const std::int64_t shift = std::max<std::int64_t>(
      0, F::kPrecision + 3 + static_cast<std::int64_t>(divisor.bit_length()) -
             static_cast<std::int64_t>(digits.bit_length()));
  BigUint numerator = digits;
  numerator.ShiftLeft(static_cast<std::size_t>(shift));
  BigUint quotient;
  BigUint remainder;
  BigUint::DivRem(numerator, divisor, quotient, remainder);
  return RoundInteger<F>(negative, quotient, -shift - u, !remainder.is_zero(),
                         env);
}

// Parses in the syntax of std::from_chars: an optional minus sign, no leading
// whitespace, and no 0x prefix for hexadecimal.
template <class F, class Env>
std::from_chars_result Parse(const char* first, const char* last, F& value,
                             std::chars_format fmt, const Env& env) {
  const char* p = first;
  const bool negative = p != last && *p == '-';
  if (negative) ++p;

  if (StartsWithIgnoringCase(p, last, "inf")) {
    p += 3;
    if (StartsWithIgnoringCase(p, last, "inity")) p += 5;
    value = F::Infinity(negative);
    return {p, std::errc{}};
  }
  if (StartsWithIgnoringCase(p, last, "nan")) {
    p += 3;
    if (p != last && *p == '(') {
      const char* q = p + 1;
      while (q != last &&
             (DigitValue(*q, false) >= 0 ||
              (ToLower(*q) >= 'a' && ToLower(*q) <= 'z') || *q == '_')) {
        ++q;
      }
      if (q != last && *q == ')') p = q + 1;
    }
    value = F::QuietNaN(negative);
    return {p, std::errc{}};
  }

  const bool hex = fmt == std::chars_format::hex;
  // The integer digits, then the fraction digits after an optional point.
  const char* integer_begin = p;
  while (p != last && DigitValue(*p, hex) >= 0) ++p;
  const char* integer_end = p;
  const char* fraction_begin = p;
  const char* fraction_end = p;
  if (p != last && *p == '.') {
    fraction_begin = p + 1;
    fraction_end = fraction_begin;
    while (fraction_end != last && DigitValue(*fraction_end, hex) >= 0) {
      ++fraction_end;
    }
    if (integer_end != integer_begin || fraction_end != fraction_begin) {
      p = fraction_end;
    }
  }
  if (integer_end == integer_begin && fraction_end == fraction_begin) {
    return {first, std::errc::invalid_argument};
  }

  std::int64_t exponent = 0;
  bool has_exponent = false;
  const bool allows_exponent =
      hex || HasFormat(fmt, std::chars_format::scientific);
  if (allows_exponent && p != last && ToLower(*p) == (hex ? 'p' : 'e')) {
    const char* q = p + 1;
    bool exponent_negative = false;
    if (q != last && (*q == '+' || *q == '-')) {
      exponent_negative = *q == '-';
      ++q;
    }
    if (q != last && DigitValue(*q, false) >= 0) {
      for (; q != last && DigitValue(*q, false) >= 0; ++q) {
        if (exponent < kExponentLimit) {
          exponent = exponent * 10 + DigitValue(*q, false);
        }
      }
      if (exponent_negative) exponent = -exponent;
      has_exponent = true;
      p = q;
    }
  }
  if (!hex && !has_exponent && !HasFormat(fmt, std::chars_format::fixed)) {
    return {first, std::errc::invalid_argument};
  }

  if (hex) {
    DigitAccumulator accumulator(true);
    for (const char* q = integer_begin; q != integer_end; ++q) {
      accumulator.Push(DigitValue(*q, true), false);
    }
    for (const char* q = fraction_begin; q != fraction_end; ++q) {
      accumulator.Push(DigitValue(*q, true), true);
    }
    const BigUint& digits = accumulator.Finish();
    value =
        digits.is_zero()
            ? F::Zero(negative)
            : RoundInteger<F>(negative, digits,
                              SaturatingAdd(4 * accumulator.adjust(), exponent),
                              false, env);
    return {p, std::errc{}};
  }

  const DecimalDigits digits = ReadDecimalDigits(integer_begin, integer_end,
                                                 fraction_begin, fraction_end);
  value = digits.value.is_zero()
              ? F::Zero(negative)
              : RoundDecimal<F>(negative, digits.value, digits.count,
                                SaturatingAdd(digits.exponent, exponent), env);
  return {p, std::errc{}};
}

// Parses a whole string: an optional sign, then a decimal number, or a
// hexadecimal one after 0x. Anything else, trailing characters included, gives
// the default NaN and raises invalid.
template <class F, class Env>
F ParseText(std::string_view text, const Env& env) {
  const char* first = text.data();
  const char* last = first + text.size();
  const char* p = first;
  bool negative = false;
  if (p != last && (*p == '+' || *p == '-')) {
    negative = *p == '-';
    ++p;
  }
  std::chars_format fmt = std::chars_format::general;
  if (last - p >= 2 && p[0] == '0' && ToLower(p[1]) == 'x') {
    fmt = std::chars_format::hex;
    p += 2;
  }
  F value = F::Zero();
  // Parse would accept a second minus sign.
  const bool second_sign = p != last && (*p == '+' || *p == '-');
  const std::from_chars_result result =
      second_sign ? std::from_chars_result{p, std::errc::invalid_argument}
                  : Parse(p, last, value, fmt, env);
  if (result.ec != std::errc{} || result.ptr != last) {
    env.raise(kFlagInvalid);
    return F::QuietNaN();
  }
  return negative ? -value : value;
}

// --- Output ---

// A decimal value d1.d2d3... * 10^exponent, as a digit string.
struct Decimal {
  std::string digits;
  std::int64_t exponent;
};

inline std::string ToDecimalString(BigUint value) {
  if (value.is_zero()) return "0";
  std::string reversed;
  while (!value.is_zero()) {
    limb_t chunk = value.DivBy10To19();
    const bool last = value.is_zero();
    for (int i = 0; i < 19 && (!last || chunk != 0); ++i) {
      reversed.push_back(static_cast<char>('0' + chunk % 10));
      chunk /= 10;
    }
  }
  return std::string(reversed.rbegin(), reversed.rend());
}

inline BigUint PowerOfTen(std::int64_t exponent) {
  BigUint power = BigUint::Pow5(static_cast<std::uint64_t>(exponent));
  power.ShiftLeft(static_cast<std::size_t>(exponent));
  return power;
}

// |x| / 10^t as an integer part and the rounding information of the rest.
struct Scaled {
  BigUint floor;
  bool round_bit;
  bool sticky;
};

// For finite nonzero x.
template <int P, int E>
Scaled ScaleByPowerOfTen(const Float<P, E>& x, std::int64_t t) {
  std::int64_t e;
  BigUint numerator = Significand(x, e);
  const std::int64_t twos = e - t;
  const std::int64_t fives = -t;
  BigUint denominator(1);
  if (fives >= 0) {
    numerator = numerator * BigUint::Pow5(static_cast<std::uint64_t>(fives));
  } else {
    denominator = BigUint::Pow5(static_cast<std::uint64_t>(-fives));
  }

  BigUint quotient;
  bool round_bit;
  bool sticky;
  if (twos >= 0) {
    numerator.ShiftLeft(static_cast<std::size_t>(twos));
  }
  if (twos < 0 && Compare(denominator, BigUint(1)) == 0) {
    // Dividing by a power of two only needs a shift.
    const auto count = static_cast<std::size_t>(-twos);
    round_bit = numerator.bit(count - 1);
    sticky = numerator.any_bit_below(count - 1);
    quotient = std::move(numerator);
    quotient.ShiftRight(count);
  } else {
    if (twos < 0) denominator.ShiftLeft(static_cast<std::size_t>(-twos));
    BigUint remainder;
    BigUint::DivRem(numerator, denominator, quotient, remainder);
    remainder.ShiftLeft(1);
    const int half = Compare(remainder, denominator);
    round_bit = half >= 0;
    sticky = half != 0 && !remainder.is_zero();
  }
  return {std::move(quotient), round_bit, sticky};
}

// Rounds a scaled value to an integer, raising inexact through `env`.
template <class Env>
BigUint RoundScaled(Scaled scaled, bool negative, Rounding mode,
                    const Env& env) {
  if (scaled.round_bit || scaled.sticky) env.raise(kFlagInexact);
  if (ShouldRoundUp(mode, negative, scaled.floor.bit(0), scaled.round_bit,
                    scaled.sticky)) {
    scaled.floor.Add(BigUint(1));
  }
  return std::move(scaled.floor);
}

// A second working format for output with many digits.
template <class F>
using WideWorking = Float<static_cast<int>(64 * (F::kLimbs + 33)), 29>;

// Rounds x to `count` significant digits on the fast path, working in T and
// starting from an estimate of the decimal exponent. Returns false if it
// cannot decide.
// The decimal digits of a nonnegative integer x below 10^count, exactly
// `count` of them, with leading zeros.
template <int P, int E>
std::string IntegerDigits(const Float<P, E>& x, std::int64_t count) {
  std::string digits(static_cast<std::size_t>(count), '0');
  if (x.is_zero()) return digits;
  auto limbs = x.raw_limbs();
  constexpr int kTotal = static_cast<int>(Float<P, E>::kLimbs) * kLimbBits;
  ShiftRight(limbs, kTotal - 1 - x.raw_exponent());
  std::size_t used = limbs.size();
  auto position = static_cast<std::size_t>(count);
  while (position > 0) {
    while (used > 0 && limbs[used - 1] == 0) --used;
    if (used == 0) break;
    limb_t chunk = DivideLimbsBy10To19(limbs.data(), used);
    for (int i = 0; i < 19 && position > 0; ++i) {
      digits[--position] = static_cast<char>('0' + chunk % 10);
      chunk /= 10;
    }
  }
  return digits;
}

template <class T, int P, int E, class Env>
bool ToSignificantDigitsFast(const Float<P, E>& x, std::int64_t count,
                             std::int64_t exponent, Rounding mode,
                             const Env& env, Decimal& result) {
  using F = Float<P, E>;
  if constexpr (!kHasFastPath<F>) {
    return false;
  } else {
    // The scaled value, below 10^count, must leave fraction bits to spare.
    // Then 10^count is exact in T, and so are the integers compared below.
    if (static_cast<double>(count) * 3.33 + 24 > T::kPrecision) return false;
    // The digit count rarely changes from call to call.
    struct Limits {
      std::int64_t count = -1;
      T lower;
      T upper;
    };
    thread_local Limits limits;
    if (limits.count != count) {
      T unused;
      PowerOfTenBounds(count - 1, limits.lower, unused);
      PowerOfTenBounds(count, limits.upper, unused);
      limits.count = count;
    }
    const T lower = limits.lower;
    const T upper = limits.upper;
    const T magnitude = abs(x);
    Context down{Rounding::kDown, 0};
    Context up{Rounding::kUp, 0};
    for (int attempt = 0; attempt < 3; ++attempt) {
      const std::int64_t t = exponent - count + 1;
      if (t >= kFastPowerLimit || t <= -kFastPowerLimit) return false;
      T p_low;
      T p_high;
      PowerOfTenBounds(-t, p_low, p_high);
      const T low = mul(magnitude, p_low, down);
      const T high = mul(magnitude, p_high, up);
      // The exponent is right when the integer part has `count` digits.
      if (high < lower) {
        --exponent;
        continue;
      }
      if (low >= upper) {
        ++exponent;
        continue;
      }
      if (low < lower || high >= upper) return false;
      Context low_context{mode, 0};
      Context high_context{mode, 0};
      const bool negative = x.signbit();
      const T q = abs(round_to_integral(negative ? -high : low, low_context));
      if (q != abs(round_to_integral(negative ? -low : high, high_context))) {
        return false;
      }
      // The scaled value is either the integer q itself or not an integer.
      if (!IsMultipleOfPowerOfTen(x, t)) env.raise(kFlagInexact);
      if (q == upper) {
        result = {"1" + std::string(static_cast<std::size_t>(count - 1), '0'),
                  exponent + 1};
      } else {
        result = {IntegerDigits(q, count), exponent};
      }
      return true;
    }
    return false;
  }
}

// x rounded to `count` >= 1 significant digits.
template <int P, int E, class Env>
Decimal ToSignificantDigits(const Float<P, E>& x, std::int64_t count,
                            Rounding mode, const Env& env) {
  if (x.is_zero()) {
    return {std::string(static_cast<std::size_t>(count), '0'), 0};
  }
  // An estimate of floor(log10 |x|) from the top bits, corrected below.
  const double top = static_cast<double>(x.raw_limbs().back()) * 0x1p-63;
  std::int64_t exponent = static_cast<std::int64_t>(std::floor(
      static_cast<double>(x.raw_exponent()) * kLog10Of2 + std::log10(top)));
  using F = Float<P, E>;
  Decimal fast;
  if constexpr (kHasNarrowWorking<F>) {
    if (ToSignificantDigitsFast<NarrowWorking<F>>(x, count, exponent, mode, env,
                                                  fast)) {
      return fast;
    }
  }
  if (ToSignificantDigitsFast<Working<F>>(x, count, exponent, mode, env,
                                          fast) ||
      ToSignificantDigitsFast<WideWorking<F>>(x, count, exponent, mode, env,
                                              fast)) {
    return fast;
  }
  const BigUint lower = PowerOfTen(count - 1);
  const BigUint upper = PowerOfTen(count);
  for (;;) {
    Scaled scaled = ScaleByPowerOfTen(x, exponent - count + 1);
    // The exponent is right when the integer part has `count` digits.
    if (Compare(scaled.floor, lower) < 0) {
      --exponent;
      continue;
    }
    if (Compare(scaled.floor, upper) >= 0) {
      ++exponent;
      continue;
    }
    BigUint q = RoundScaled(std::move(scaled), x.signbit(), mode, env);
    if (Compare(q, upper) == 0) {
      // Rounded up to the next power of ten.
      return {"1" + std::string(static_cast<std::size_t>(count - 1), '0'),
              exponent + 1};
    }
    return {ToDecimalString(std::move(q)), exponent};
  }
}

// Whether a decimal reads back as x under round-to-nearest-even.
template <int P, int E>
bool ReadsBack(const Float<P, E>& x, const Decimal& decimal) {
  DigitAccumulator accumulator(false);
  for (const char c : decimal.digits) accumulator.Push(c - '0', false);
  const BigUint& digits = accumulator.Finish();
  if (digits.is_zero()) return false;
  const auto length = static_cast<std::int64_t>(decimal.digits.size());
  const Float<P, E> parsed = RoundDecimal<Float<P, E>>(
      x.signbit(), digits, accumulator.count(),
      decimal.exponent - (length - 1) + accumulator.adjust(), NearestEnv{});
  return parsed.to_bits() == x.to_bits();
}

// The shortest decimal that reads back as x, for x far from 1, where the
// exact method below works with huge integers. Whether some n-digit decimal
// reads back is monotonic in n, so a binary search finds the shortest length,
// and only the two n-digit decimals around x are candidates. The nearer one
// is preferred, as in the exact method.
template <int P, int E>
Decimal ToShortestDigitsBySearch(const Float<P, E>& x) {
  const auto fits = [&x](std::int64_t count, Decimal& result) {
    result =
        ToSignificantDigits(x, count, Rounding::kNearestEven, NearestEnv{});
    if (ReadsBack(x, result)) return true;
    const Decimal toward_zero =
        ToSignificantDigits(x, count, Rounding::kTowardZero, NearestEnv{});
    const Rounding away = x.signbit() ? Rounding::kDown : Rounding::kUp;
    result = toward_zero.digits == result.digits &&
                     toward_zero.exponent == result.exponent
                 ? ToSignificantDigits(x, count, away, NearestEnv{})
                 : toward_zero;
    return ReadsBack(x, result);
  };
  // P bits always need at most this many digits.
  std::int64_t low = 1;
  std::int64_t high = static_cast<std::int64_t>(std::ceil(P * kLog10Of2)) + 1;
  Decimal result;
  while (low < high) {
    const std::int64_t middle = (low + high) / 2;
    if (fits(middle, result)) {
      high = middle;
    } else {
      low = middle + 1;
    }
  }
  fits(low, result);
  return result;
}

// Adds one to a decimal digit string, which may grow by a digit.
inline void IncrementDigits(std::string& digits) {
  for (std::size_t i = digits.size(); i-- > 0;) {
    if (digits[i] != '9') {
      ++digits[i];
      return;
    }
    digits[i] = '0';
  }
  digits.insert(digits.begin(), '1');
}

// The shortest decimal that reads back as x on the fast path, working in T.
// Returns false if it cannot decide.
//
// With |x| scaled by 10^-t to an integer part of about M + 2 digits, where M
// digits always suffice, every candidate is an integer at this scale. The
// integers strictly inside the rounding interval of x form a range [A, B]:
// its shortest members are the multiples of the largest power of ten in it,
// found where the digits of A - 1 and B first differ, and the nearest of
// those to x is the answer. Interval bounds on the scaled interval ends and on
// x settle all of this unless one of them lies very close to an integer.
template <class T, int P, int E>
bool ToShortestDigitsFast(const Float<P, E>& x, Decimal& result) {
  using F = Float<P, E>;
  if constexpr (!kHasFastPath<F>) {
    return false;
  } else {
    constexpr std::int64_t kMaxDigits =
        static_cast<std::int64_t>(P * kLog10Of2) + 2;
    // The estimate of the decimal exponent is off by at most one, so the
    // scaled values stay below 10^width.
    constexpr std::int64_t kWidth = kMaxDigits + 4;
    if constexpr (static_cast<double>(kWidth) * 3.33 + 24 > T::kPrecision) {
      return false;
    } else {
      const double top = static_cast<double>(x.raw_limbs().back()) * 0x1p-63;
      const auto estimate = static_cast<std::int64_t>(std::floor(
          static_cast<double>(x.raw_exponent()) * kLog10Of2 + std::log10(top)));
      const std::int64_t t = estimate - kMaxDigits - 1;
      if (t >= kFastPowerLimit || t <= -kFastPowerLimit) return false;

      // The rounding interval of |x|, whose ends are exact in T. The gap
      // below is half the gap above at a power of two, except at the bottom
      // of the normal range.
      const T value = abs(x);
      const auto grid = static_cast<int>(
          std::max<std::int64_t>(x.raw_exponent(), F::kEmin) - (P - 1));
      const auto& limbs = x.raw_limbs();
      bool power_of_two = limbs.back() == limb_t{1} << (kLimbBits - 1);
      for (std::size_t i = 0; i + 1 < limbs.size(); ++i) {
        power_of_two = power_of_two && limbs[i] == 0;
      }
      const bool lower_closer = power_of_two && x.raw_exponent() > F::kEmin;
      const T half_gap = scalbn(T(1), grid - 1);
      const T low_end =
          value - (lower_closer ? scalbn(T(1), grid - 2) : half_gap);
      const T high_end = value + half_gap;

      T p_low;
      T p_high;
      PowerOfTenBounds(-t, p_low, p_high);
      Context down{Rounding::kDown, 0};
      Context up{Rounding::kUp, 0};
      const auto floor_of = [](const T& y) {
        Context context{Rounding::kDown, 0};
        return round_to_integral(y, context);
      };
      const T low_1 = mul(low_end, p_low, down);
      const T low_2 = mul(low_end, p_high, up);
      const T high_1 = mul(high_end, p_low, down);
      const T high_2 = mul(high_end, p_high, up);
      const T x_1 = mul(value, p_low, down);
      const T x_2 = mul(value, p_high, up);
      const T low_floor = floor_of(low_2);
      const T high_floor = floor_of(high_1);
      const T x_floor = floor_of(x_1);
      // No integer may lie within the bounds of an end, nor may an end be an
      // integer: whether it counts depends on the parity of x.
      if (floor_of(low_1) != low_floor || floor_of(high_2) != high_floor ||
          floor_of(x_2) != x_floor || low_floor == low_1 ||
          high_floor == high_1 || high_floor <= low_floor) {
        return false;
      }

      // Candidates are the integers in [low_floor + 1, high_floor].
      const std::string below = IntegerDigits(low_floor, kWidth);
      const std::string above = IntegerDigits(high_floor, kWidth);
      const std::string nearby = IntegerDigits(x_floor, kWidth);
      std::size_t keep = 0;
      while (below[keep] == above[keep]) ++keep;
      ++keep;
      // Multiples of 10^dropped, as prefixes of `keep` digits.
      const auto dropped = static_cast<std::size_t>(kWidth) - keep;
      std::string lowest = below.substr(0, keep);
      IncrementDigits(lowest);
      const std::string highest = above.substr(0, keep);
      std::string floor_prefix = nearby.substr(0, keep);
      std::string ceiling_prefix = floor_prefix;
      IncrementDigits(ceiling_prefix);

      std::string chosen;
      if (floor_prefix < lowest) {
        chosen = lowest;
      } else if (ceiling_prefix.size() > keep || ceiling_prefix > highest) {
        chosen = highest;
      } else {
        // Compare x with the midpoint between the two prefixes.
        int order;
        if (dropped == 0) {
          const T half(0.5);
          if (x_1 >= x_floor + half) {
            order = 1;
          } else if (x_2 < x_floor + half) {
            order = -1;
          } else {
            return false;
          }
        } else {
          const std::string rest = nearby.substr(keep);
          const std::string middle = "5" + std::string(dropped - 1, '0');
          order = rest.compare(middle);
          if (order == 0 && !IsMultipleOfPowerOfTen(x, t)) order = 1;
        }
        if (order == 0) {
          // A tie: the even last digit.
          order = (floor_prefix.back() - '0') % 2 == 0 ? -1 : 1;
        }
        chosen = order < 0 ? floor_prefix : ceiling_prefix;
      }

      const std::size_t first = chosen.find_first_not_of('0');
      const std::size_t last = chosen.find_last_not_of('0');
      result.digits = chosen.substr(first, last - first + 1);
      result.exponent = t + static_cast<std::int64_t>(dropped) +
                        static_cast<std::int64_t>(chosen.size() - first) - 1;
      return true;
    }
  }
}

// The shortest decimal that reads back as x, for finite nonzero x, after
// R. G. Burger and R. K. Dybvig, "Printing floating-point numbers quickly and
// accurately", PLDI 1996, with exact integers throughout.
template <int P, int E>
Decimal ToShortestDigits(const Float<P, E>& x) {
  using F = Float<P, E>;
  Decimal fast;
  if (ToShortestDigitsFast<Working<F>>(x, fast)) return fast;
  // The exact method's integers grow with the exponent.
  if (kHasFastPath<F> &&
      (x.raw_exponent() > 4096 || x.raw_exponent() < -4096)) {
    return ToShortestDigitsBySearch(x);
  }
  std::int64_t e;
  BigUint f = Significand(x, e);
  // Move subnormals onto their grid, where neighbours are one unit apart.
  if (x.raw_exponent() < F::kEmin) {
    const std::int64_t shift = F::kEmin - std::int64_t{x.raw_exponent()};
    f.ShiftRight(static_cast<std::size_t>(shift));
    e += shift;
  }
  // The gap below is half the gap above at a power of two, except at the
  // bottom of the normal range.
  const bool lower_closer = x.raw_exponent() > F::kEmin &&
                            f.bit_length() == static_cast<std::size_t>(P) &&
                            !f.any_bit_below(P - 1);
  // Ties read back to even significands, so the bounds count if f is even.
  const bool inclusive = !f.bit(0);

  // x = r / s, and its rounding interval is [(r - m_minus) / s,
  // (r + m_plus) / s].
  BigUint r = f;
  BigUint s(1);
  BigUint m_plus(1);
  BigUint m_minus(1);
  const std::size_t extra = lower_closer ? 1 : 0;
  if (e >= 0) {
    r.ShiftLeft(static_cast<std::size_t>(e) + 1 + extra);
    s.ShiftLeft(1 + extra);
    m_plus.ShiftLeft(static_cast<std::size_t>(e) + extra);
    m_minus.ShiftLeft(static_cast<std::size_t>(e));
  } else {
    r.ShiftLeft(1 + extra);
    s.ShiftLeft(static_cast<std::size_t>(-e) + 1 + extra);
    m_plus.ShiftLeft(extra);
  }

  // An estimate of ceil(log10 x) that may be too low, but never too high.
  std::int64_t k = static_cast<std::int64_t>(
      std::ceil(static_cast<double>(x.raw_exponent()) * kLog10Of2 - 1e-6));
  if (k >= 0) {
    s = s * PowerOfTen(k);
  } else {
    const BigUint scale = PowerOfTen(-k);
    r = r * scale;
    m_plus = m_plus * scale;
    m_minus = m_minus * scale;
  }
  const auto high_reaches = [&](const BigUint& value) {
    BigUint high = value;
    high.Add(m_plus);
    const int order = Compare(high, s);
    return inclusive ? order >= 0 : order > 0;
  };
  while (high_reaches(r)) {
    s.MulAdd(10, 0);
    ++k;
  }

  std::string digits;
  for (;;) {
    r.MulAdd(10, 0);
    m_plus.MulAdd(10, 0);
    m_minus.MulAdd(10, 0);
    int digit = 0;
    while (Compare(r, s) >= 0) {
      r.Sub(s);
      ++digit;
    }
    const int low_order = Compare(r, m_minus);
    const bool low = inclusive ? low_order <= 0 : low_order < 0;
    const bool high = high_reaches(r);
    if (!low && !high) {
      digits.push_back(static_cast<char>('0' + digit));
      continue;
    }
    if (low && high) {
      // Both ends terminate: take the nearer, or the even digit on a tie.
      BigUint twice = r;
      twice.ShiftLeft(1);
      const int order = Compare(twice, s);
      if (order > 0 || (order == 0 && digit % 2 != 0)) ++digit;
    } else if (high) {
      ++digit;
    }
    digits.push_back(static_cast<char>('0' + digit));
    break;
  }
  // The last digit never carries: while the loop runs, r + m_plus stays below
  // s, or at most s for an exclusive bound, so a digit that is incremented is
  // at most 8.
  return {digits, k - 1};
}

// d1.d2d3... * 10^exponent in plain positional notation.
inline std::string FixedNotation(const Decimal& decimal) {
  const std::string& digits = decimal.digits;
  const auto length = static_cast<std::int64_t>(digits.size());
  if (decimal.exponent < 0) {
    return "0." +
           std::string(static_cast<std::size_t>(-decimal.exponent - 1), '0') +
           digits;
  }
  if (decimal.exponent + 1 >= length) {
    return digits +
           std::string(static_cast<std::size_t>(decimal.exponent + 1 - length),
                       '0');
  }
  const auto point = static_cast<std::size_t>(decimal.exponent + 1);
  return digits.substr(0, point) + "." + digits.substr(point);
}

inline std::string ExponentSuffix(char marker, std::int64_t exponent,
                                  int min_digits) {
  std::string text(1, marker);
  text.push_back(exponent < 0 ? '-' : '+');
  std::string magnitude = std::to_string(exponent < 0 ? -exponent : exponent);
  if (static_cast<int>(magnitude.size()) < min_digits) {
    magnitude.insert(0, static_cast<std::size_t>(min_digits) - magnitude.size(),
                     '0');
  }
  return text + magnitude;
}

inline std::string ScientificNotation(const Decimal& decimal) {
  std::string text(1, decimal.digits[0]);
  if (decimal.digits.size() > 1) {
    text += ".";
    text += decimal.digits.substr(1);
  }
  return text + ExponentSuffix('e', decimal.exponent, 2);
}

inline void StripTrailingZeros(Decimal& decimal) {
  const std::size_t end = decimal.digits.find_last_not_of('0');
  decimal.digits.resize(end == std::string::npos ? 1 : end + 1);
}

// Formats nonfinite values and zeros, which need no digits, or returns false.
template <int P, int E>
bool FormatSpecial(const Float<P, E>& x, std::string& text) {
  if (x.is_nan()) {
    text = x.signbit() ? "-nan" : "nan";
    return true;
  }
  if (x.is_inf()) {
    text = x.signbit() ? "-inf" : "inf";
    return true;
  }
  return false;
}

// Hexadecimal: 1.hhh...p+e, with `precision` digits after the point, or as
// many as needed for an exact result if it is negative.
template <int P, int E, class Env>
std::string FormatHex(const Float<P, E>& x, int precision, Rounding mode,
                      const Env& env) {
  constexpr int kFractionBits = P - 1;
  constexpr int kFullDigits = (kFractionBits + 3) / 4;
  std::int64_t exponent = 0;
  // The fraction bits, left-aligned to whole hexadecimal digits.
  BigUint fraction;
  if (!x.is_zero()) {
    std::int64_t e;
    fraction = Significand(x, e);
    exponent = x.raw_exponent();
    // Drop the leading 1.
    BigUint leading(1);
    leading.ShiftLeft(kFractionBits);
    fraction.Sub(leading);
    fraction.ShiftLeft(
        static_cast<std::size_t>(4 * kFullDigits - kFractionBits));
  }
  int digits = kFullDigits;
  int leading_digit = x.is_zero() ? 0 : 1;
  if (precision >= 0 && precision < kFullDigits) {
    const auto dropped =
        static_cast<std::size_t>(4 * (kFullDigits - precision));
    const bool round_bit = fraction.bit(dropped - 1);
    const bool sticky = fraction.any_bit_below(dropped - 1);
    fraction.ShiftRight(dropped);
    if (round_bit || sticky) env.raise(kFlagInexact);
    const bool is_odd =
        precision == 0 ? leading_digit % 2 != 0 : fraction.bit(0);
    if (ShouldRoundUp(mode, x.signbit(), is_odd, round_bit, sticky)) {
      fraction.Add(BigUint(1));
      if (fraction.bit_length() > static_cast<std::size_t>(4 * precision)) {
        // 1.fff... rounded up to 2: renormalize to 1.000... * 2^(e + 1).
        fraction = BigUint();
        ++exponent;
      }
    }
    digits = precision;
  }
  std::string hex_digits(static_cast<std::size_t>(digits), '0');
  for (int i = digits; i-- > 0;) {
    const auto nibble = static_cast<int>(
        fraction.limbs().empty() ? 0 : fraction.limbs()[0] & 15);
    hex_digits[static_cast<std::size_t>(i)] = "0123456789abcdef"[nibble];
    fraction.ShiftRight(4);
  }
  if (precision < 0) {
    const std::size_t end = hex_digits.find_last_not_of('0');
    hex_digits.resize(end == std::string::npos ? 0 : end + 1);
  } else if (precision > digits) {
    hex_digits.append(static_cast<std::size_t>(precision - digits), '0');
  }
  std::string text(1, static_cast<char>('0' + leading_digit));
  if (!hex_digits.empty()) text += "." + hex_digits;
  return text + ExponentSuffix('p', exponent, 1);
}

template <int P, int E, class Env>
std::string FormatPrecision(const Float<P, E>& x, std::chars_format fmt,
                            int precision, bool alternate, Rounding mode,
                            const Env& env);

// Shortest output. `plain` picks fixed or scientific notation, whichever is
// shorter, as std::to_chars does without a format.
template <int P, int E>
std::string FormatShortest(const Float<P, E>& x, std::chars_format fmt,
                           bool plain) {
  std::string text;
  if (FormatSpecial(x, text)) return text;
  const std::string sign = x.signbit() ? "-" : "";
  if (fmt == std::chars_format::hex) {
    return sign + FormatHex(x, -1, Rounding::kNearestEven, NearestEnv{});
  }
  const Decimal decimal = x.is_zero() ? Decimal{"0", 0} : ToShortestDigits(x);
  // Fixed notation is %f with the fewest decimals that read back, as in
  // std::to_chars. The integer part is then exact rather than padded with
  // zeros, as the shortest digits are where the integer is no shorter.
  const auto fixed = [&] {
    const auto length = static_cast<std::int64_t>(decimal.digits.size());
    const std::int64_t places =
        std::max<std::int64_t>(0, length - 1 - decimal.exponent);
    if (places > 0 || x.is_zero()) return sign + FixedNotation(decimal);
    return FormatPrecision(x, std::chars_format::fixed, 0, false,
                           Rounding::kNearestEven, NearestEnv{});
  };
  if (plain) {
    // Compares lengths without building the fixed text, which can have
    // millions of digits.
    const auto length = static_cast<std::int64_t>(decimal.digits.size());
    const std::int64_t exponent = decimal.exponent;
    std::int64_t fixed_length = exponent + 1;
    if (exponent < 0) {
      fixed_length = 1 - exponent + length;
    } else if (length > exponent + 1) {
      fixed_length = length + 1;
    }
    const std::string scientific = ScientificNotation(decimal);
    if (static_cast<std::int64_t>(scientific.size()) < fixed_length) {
      return sign + scientific;
    }
    return fixed();
  }
  if (fmt == std::chars_format::fixed) return fixed();
  // %g with its default precision of 6 decides the notation, as in
  // std::to_chars.
  if (fmt == std::chars_format::general && decimal.exponent >= -4 &&
      decimal.exponent < 6) {
    return fixed();
  }
  return sign + ScientificNotation(decimal);
}

// Output with a precision, as printf's %.*e, %.*f, %.*g and %.*a do. With
// `alternate`, %g keeps trailing zeros, as with the # flag.
template <int P, int E, class Env>
std::string FormatPrecision(const Float<P, E>& x, std::chars_format fmt,
                            int precision, bool alternate, Rounding mode,
                            const Env& env) {
  std::string text;
  if (FormatSpecial(x, text)) return text;
  const std::string sign = x.signbit() ? "-" : "";
  if (precision < 0) precision = fmt == std::chars_format::hex ? -1 : 6;
  if (fmt == std::chars_format::hex) {
    return sign + FormatHex(x, precision, mode, env);
  }

  if (fmt == std::chars_format::fixed) {
    std::string digits = "0";
    if (!x.is_zero()) {
      digits = ToDecimalString(RoundScaled(ScaleByPowerOfTen(x, -precision),
                                           x.signbit(), mode, env));
    }
    const auto places = static_cast<std::size_t>(precision);
    if (digits.size() <= places) {
      digits.insert(0, places + 1 - digits.size(), '0');
    }
    if (places == 0) return sign + digits;
    return sign + digits.substr(0, digits.size() - places) + "." +
           digits.substr(digits.size() - places);
  }

  if (fmt == std::chars_format::scientific) {
    const Decimal decimal =
        ToSignificantDigits(x, std::int64_t{precision} + 1, mode, env);
    return sign + ScientificNotation(decimal);
  }

  // %g: precision significant digits, in fixed notation for exponents from -4
  // to below the precision.
  const std::int64_t significant = precision == 0 ? 1 : precision;
  Decimal decimal = ToSignificantDigits(x, significant, mode, env);
  if (!alternate) StripTrailingZeros(decimal);
  if (decimal.exponent >= -4 && decimal.exponent < significant) {
    std::string fixed = FixedNotation(decimal);
    if (alternate && fixed.find('.') == std::string::npos) fixed += ".";
    return sign + fixed;
  }
  std::string scientific = ScientificNotation(decimal);
  if (alternate && decimal.digits.size() == 1) {
    scientific.insert(1, ".");
  }
  return sign + scientific;
}

inline std::to_chars_result CopyOut(char* first, char* last,
                                    const std::string& text) {
  if (static_cast<std::size_t>(last - first) < text.size()) {
    return {last, std::errc::value_too_large};
  }
  std::memcpy(first, text.data(), text.size());
  return {first + text.size(), std::errc{}};
}

}  // namespace detail

// --- Output, as std::to_chars ---

// The shortest representation that reads back as x, in fixed or scientific
// notation, whichever is shorter.
template <int P, int E>
std::to_chars_result to_chars(char* first, char* last, const Float<P, E>& x) {
  return detail::CopyOut(
      first, last, detail::FormatShortest(x, std::chars_format::general, true));
}

// The shortest representation in the given format.
template <int P, int E>
std::to_chars_result to_chars(char* first, char* last, const Float<P, E>& x,
                              std::chars_format fmt) {
  return detail::CopyOut(first, last, detail::FormatShortest(x, fmt, false));
}

// With a precision, as printf would format it, rounded to nearest-even.
template <int P, int E>
std::to_chars_result to_chars(char* first, char* last, const Float<P, E>& x,
                              std::chars_format fmt, int precision) {
  return detail::CopyOut(
      first, last,
      detail::FormatPrecision(x, fmt, precision, false, Rounding::kNearestEven,
                              detail::NearestEnv{}));
}

// With a precision, rounded in context.rounding, raising inexact.
template <int P, int E>
std::to_chars_result to_chars(char* first, char* last, const Float<P, E>& x,
                              std::chars_format fmt, int precision,
                              Context& context) {
  return detail::CopyOut(
      first, last,
      detail::FormatPrecision(x, fmt, precision, false, context.rounding,
                              detail::ContextEnv{context}));
}

// The shortest representation, as to_chars without a format.
template <int P, int E>
std::string to_string(const Float<P, E>& x) {
  return detail::FormatShortest(x, std::chars_format::general, true);
}

// --- Input, as std::from_chars ---

// Rounds to nearest-even.
template <int P, int E>
std::from_chars_result from_chars(
    const char* first, const char* last, Float<P, E>& value,
    std::chars_format fmt = std::chars_format::general) {
  return detail::Parse(first, last, value, fmt, detail::NearestEnv{});
}

// Rounds in context.rounding and raises flags.
template <int P, int E>
std::from_chars_result from_chars(const char* first, const char* last,
                                  Float<P, E>& value, std::chars_format fmt,
                                  Context& context) {
  return detail::Parse(first, last, value, fmt, detail::ContextEnv{context});
}

// Parses a whole string with an optional sign, in decimal, or in hexadecimal
// after 0x. Invalid input gives the default NaN.
template <class F>
  requires detail::kIsFloat<F>
F from_string(std::string_view text) {
  return detail::ParseText<F>(text, detail::NearestEnv{});
}

// As above, rounding in context.rounding. Invalid input raises invalid.
template <class F>
  requires detail::kIsFloat<F>
F from_string(std::string_view text, Context& context) {
  return detail::ParseText<F>(text, detail::ContextEnv{context});
}

// --- Streams ---

// Follows the stream's precision and its fixed, scientific, hexfloat,
// showpoint, showpos and uppercase flags. The default format is %g.
template <int P, int E>
std::ostream& operator<<(std::ostream& stream, const Float<P, E>& x) {
  const std::ios_base::fmtflags flags = stream.flags();
  const std::ios_base::fmtflags field = flags & std::ios_base::floatfield;
  const int precision = static_cast<int>(stream.precision());
  std::string text;
  if (field == (std::ios_base::fixed | std::ios_base::scientific)) {
    text = detail::FormatShortest(x, std::chars_format::hex, false);
    if (x.is_finite()) text.insert(x.signbit() ? 1 : 0, "0x");
  } else {
    std::chars_format fmt = std::chars_format::general;
    if (field == std::ios_base::fixed) fmt = std::chars_format::fixed;
    if (field == std::ios_base::scientific) fmt = std::chars_format::scientific;
    text = detail::FormatPrecision(
        x, fmt, precision, (flags & std::ios_base::showpoint) != 0,
        Rounding::kNearestEven, detail::NearestEnv{});
  }
  if ((flags & std::ios_base::showpos) != 0 && text[0] != '-') {
    text.insert(0, "+");
  }
  if ((flags & std::ios_base::uppercase) != 0) {
    for (char& c : text) {
      if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
    }
  }
  return stream << text;
}

// Reads a token of number characters and parses it as from_string does. A
// token that does not parse sets failbit.
template <int P, int E>
std::istream& operator>>(std::istream& stream, Float<P, E>& x) {
  const std::istream::sentry sentry(stream);
  if (!sentry) return stream;
  std::string token;
  for (;;) {
    // At the end, peek sets eofbit.
    const auto next = stream.peek();
    if (next == std::char_traits<char>::eof()) break;
    const char c = static_cast<char>(next);
    const bool part =
        (c >= '0' && c <= '9') ||
        (detail::ToLower(c) >= 'a' && detail::ToLower(c) <= 'z') || c == '.' ||
        c == '+' || c == '-' || c == '(' || c == ')' || c == '_';
    if (!part) break;
    token.push_back(c);
    stream.get();
  }
  Context context;
  const Float<P, E> value = from_string<Float<P, E>>(token, context);
  if ((context.flags & kFlagInvalid) != 0) {
    stream.setstate(std::ios_base::failbit);
  } else {
    x = value;
  }
  return stream;
}

}  // namespace mpfix

#endif  // MPFIX_CHARCONV_H_
