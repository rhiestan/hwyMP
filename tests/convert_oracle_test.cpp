// Copyright 2026 The hwyMP Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Differential tests of the conversions against MPFR: between formats, from
// and to integers, and from and to decimal strings, in all five rounding
// modes and with flags. Set HWYMP_TEST_SCALE to multiply the iteration counts.

#include <bit>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include "hwymp/charconv.h"
#include "oracle.h"

namespace hwymp {
namespace {

using namespace test;  // NOLINT

template <class F>
std::string Describe(const F& x) {
  return Hex(x) + " (" + to_string(x) + ")";
}

// --- Between formats ---

template <class To, class From>
Expected<To> ExpectedConversion(const From& x, Rounding mode) {
  if (x.is_nan()) {
    return {To::QuietNaN(), x.is_signaling() ? unsigned{kFlagInvalid} : 0u};
  }
  if (x.is_inf()) return {To::Infinity(x.signbit()), 0u};
  if (x.is_zero()) return {To::Zero(x.signbit()), 0u};
  const View<From> view(x);
  return RoundExact<To>(
      [&view](mpfr_ptr result, mpfr_rnd_t rnd) {
        return mpfr_set(result, view.get(), rnd);
      },
      mode);
}

template <class To, class From>
bool CheckConversion(const From& x, Rounding mode) {
  const Expected<To> expected = ExpectedConversion<To>(x, mode);
  Context context{mode};
  const To actual = convert<To>(x, context);
  bool ok =
      SameResult(actual, expected.value) && context.flags == expected.flags;
  if (mode == Rounding::kNearestEven) {
    ok = ok && SameResult(To(x), expected.value);
  }
  if constexpr (std::is_same_v<To, Float<53, 11>>) {
    // The built-in double gives the same bits.
    Context native{mode};
    const double d = convert<double>(x, native);
    ok = ok && native.flags == expected.flags &&
         (x.is_nan()
              ? d != d
              : std::bit_cast<std::uint64_t>(d) == expected.value.to_bits()[0]);
  }
  if (!ok) {
    ADD_FAILURE() << "P=" << From::kPrecision << " E=" << From::kExponentBits
                  << " to P=" << To::kPrecision << " E=" << To::kExponentBits
                  << " mode=" << static_cast<int>(mode)
                  << "\n  x        = " << Hex(x)
                  << "\n  actual   = " << Hex(actual) << " flags "
                  << context.flags << "\n  expected = " << Hex(expected.value)
                  << " flags " << expected.flags;
  }
  return ok;
}

template <class To, class From>
bool CheckConversionAllModes(const From& x) {
  for (const Rounding mode : kModes) {
    if (!CheckConversion<To>(x, mode)) return false;
  }
  return true;
}

template <class From, class... To>
void CheckFormatsExhaustively() {
  const limb_t count = limb_t{1} << (From::kPrecision + From::kExponentBits);
  for (limb_t bits = 0; bits < count; ++bits) {
    const From x = From::from_bits({bits});
    if (!(CheckConversionAllModes<To>(x) && ...)) return;
  }
}

TEST(ConvertOracle, ExhaustiveSmallFormats) {
  using A = Float<3, 2>;
  using B = Float<4, 3>;
  using C = Float<5, 4>;
  using D = Float<8, 5>;
  using G = Float<11, 4>;
  CheckFormatsExhaustively<A, B, C, D, G>();
  CheckFormatsExhaustively<B, A, C, D, G>();
  CheckFormatsExhaustively<C, A, B, D, G>();
  CheckFormatsExhaustively<D, A, B, C, G, Float<53, 11>>();
  CheckFormatsExhaustively<G, A, B, C, D>();
}

template <class From, class... To>
void CheckFormatsRandomly(int iterations, limb_t seed) {
  Random random(seed);
  const int total = iterations * TestScale();
  for (int i = 0; i < total; ++i) {
    const From x = RandomValue<From>(random);
    if (!(CheckConversionAllModes<To>(x) && ...)) return;
  }
}

TEST(ConvertOracle, RandomWideFormats) {
  using Double = Float<53, 11>;
  using Extended = Float<64, 15>;
  CheckFormatsRandomly<binary128, Double, Extended, binary256, binary512,
                       Float<24, 8>, Float<8, 5>>(20000, 1);
  CheckFormatsRandomly<binary256, Double, binary128, binary512, Extended,
                       Float<128, 15>>(20000, 2);
  CheckFormatsRandomly<binary512, Double, binary128, binary256, Float<256, 19>,
                       Float<448, 21>>(20000, 3);
  CheckFormatsRandomly<Double, Float<24, 8>, Float<8, 5>, binary128,
                       Float<52, 11>, Float<53, 8>>(20000, 4);
  CheckFormatsRandomly<Float<129, 15>, Float<128, 15>, Float<127, 15>,
                       Float<65, 15>, Float<64, 15>, Float<63, 15>>(20000, 5);
}

// --- Integers ---

template <class To, class I>
bool CheckFromInteger(I value, Rounding mode) {
  Mpfr exact(64);
  if constexpr (std::is_signed_v<I>) {
    mpfr_set_sj(exact.get(), value, MPFR_RNDN);
  } else {
    mpfr_set_uj(exact.get(), value, MPFR_RNDN);
  }
  Expected<To> expected{To::Zero(), 0u};
  if (value != 0) {
    expected = RoundExact<To>(
        [&exact](mpfr_ptr result, mpfr_rnd_t rnd) {
          return mpfr_set(result, exact.get(), rnd);
        },
        mode);
  }
  Context context{mode};
  const To actual = convert<To>(value, context);
  const bool ok =
      SameResult(actual, expected.value) && context.flags == expected.flags &&
      (mode != Rounding::kNearestEven || SameResult(To(value), expected.value));
  if (!ok) {
    ADD_FAILURE() << "P=" << To::kPrecision << " from integer "
                  << std::to_string(value) << " mode=" << static_cast<int>(mode)
                  << "\n  actual   = " << Hex(actual) << " flags "
                  << context.flags << "\n  expected = " << Hex(expected.value)
                  << " flags " << expected.flags;
  }
  return ok;
}

template <class I>
I RandomInteger(Random& random) {
  using U = std::make_unsigned_t<I>;
  // Uniform bits, or a short value, or one next to a power of two.
  switch (random.Below(4)) {
    case 0:
      return static_cast<I>(static_cast<U>(random.Next()));
    case 1:
      return static_cast<I>(random.Below(300)) -
             static_cast<I>(std::is_signed_v<I> ? 150 : 0);
    default: {
      const int bits = static_cast<int>(sizeof(I) * 8);
      const U power =
          static_cast<U>(U{1} << random.Below(static_cast<limb_t>(bits)));
      const U nearby = static_cast<U>(power + random.Below(5) - 2);
      return static_cast<I>(random.Below(2) != 0 ? nearby
                                                 : static_cast<U>(0 - nearby));
    }
  }
}

template <class To, class I>
void CheckFromIntegers(int iterations, limb_t seed) {
  Random random(seed);
  const int total = iterations * TestScale();
  for (int i = 0; i < total; ++i) {
    const I value = RandomInteger<I>(random);
    for (const Rounding mode : kModes) {
      if (!CheckFromInteger<To>(value, mode)) return;
    }
  }
}

TEST(ConvertOracle, FromIntegers) {
  CheckFromIntegers<Float<8, 5>, std::int64_t>(20000, 11);
  CheckFromIntegers<Float<8, 5>, std::uint8_t>(2000, 12);
  CheckFromIntegers<Float<53, 11>, std::int64_t>(20000, 13);
  CheckFromIntegers<Float<53, 11>, std::uint64_t>(20000, 14);
  CheckFromIntegers<Float<63, 8>, std::int64_t>(20000, 15);
  CheckFromIntegers<Float<63, 8>, std::uint64_t>(20000, 16);
  CheckFromIntegers<Float<64, 8>, std::uint64_t>(20000, 17);
  CheckFromIntegers<Float<24, 8>, std::int32_t>(20000, 18);
  CheckFromIntegers<Float<4, 3>, std::int16_t>(20000, 19);
  CheckFromIntegers<binary128, std::int64_t>(5000, 20);
}

template <class I, class F>
bool CheckToInteger(const F& x, Rounding mode) {
  using Limits = std::numeric_limits<I>;
  I expected_value = 0;
  unsigned expected_flags = 0;
  if (x.is_nan()) {
    expected_flags = kFlagInvalid;
  } else if (x.is_inf()) {
    expected_flags = kFlagInvalid;
    expected_value = x.signbit() ? Limits::min() : Limits::max();
  } else if (!x.is_zero()) {
    const View<F> view(x);
    Mpfr rounded(F::kPrecision + 1);
    switch (mode) {
      case Rounding::kNearestEven:
        mpfr_rint(rounded.get(), view.get(), MPFR_RNDN);
        break;
      case Rounding::kNearestAway:
        mpfr_round(rounded.get(), view.get());
        break;
      case Rounding::kTowardZero:
        mpfr_rint(rounded.get(), view.get(), MPFR_RNDZ);
        break;
      case Rounding::kUp:
        mpfr_rint(rounded.get(), view.get(), MPFR_RNDU);
        break;
      case Rounding::kDown:
        mpfr_rint(rounded.get(), view.get(), MPFR_RNDD);
        break;
    }
    Mpfr low(64);
    Mpfr high(64);
    if constexpr (Limits::is_signed) {
      mpfr_set_sj(low.get(), Limits::min(), MPFR_RNDN);
      mpfr_set_sj(high.get(), Limits::max(), MPFR_RNDN);
    } else {
      mpfr_set_uj(low.get(), 0, MPFR_RNDN);
      mpfr_set_uj(high.get(), Limits::max(), MPFR_RNDN);
    }
    const bool in_range = mpfr_cmp(rounded.get(), low.get()) >= 0 &&
                          mpfr_cmp(rounded.get(), high.get()) <= 0;
    if (in_range) {
      if constexpr (Limits::is_signed) {
        expected_value = static_cast<I>(mpfr_get_sj(rounded.get(), MPFR_RNDN));
      } else {
        expected_value = static_cast<I>(mpfr_get_uj(rounded.get(), MPFR_RNDN));
      }
    }
    if (!in_range) {
      expected_flags = kFlagInvalid;
      expected_value = x.signbit() ? Limits::min() : Limits::max();
    }
  }
  const bool inexact =
      x.is_finite() && !x.is_zero() && mpfr_integer_p(View<F>(x).get()) == 0;
  const unsigned exact_flags =
      expected_flags | (expected_flags == 0 && inexact ? kFlagInexact : 0u);

  Context context{mode};
  const I actual = to_integer<I>(x, mode, context);
  Context exact_context{mode};
  const I exact_actual = to_integer_exact<I>(x, mode, exact_context);
  bool ok = actual == expected_value && context.flags == expected_flags &&
            exact_actual == expected_value &&
            exact_context.flags == exact_flags &&
            to_integer<I>(x, mode) == expected_value;
  if (mode == Rounding::kTowardZero) {
    ok = ok && static_cast<I>(x) == expected_value;
  }
  if (!ok) {
    ADD_FAILURE() << "P=" << F::kPrecision << " to integer of " << sizeof(I) * 8
                  << (Limits::is_signed ? " signed" : " unsigned")
                  << " bits, mode=" << static_cast<int>(mode)
                  << "\n  x        = " << Describe(x)
                  << "\n  actual   = " << std::to_string(actual) << " flags "
                  << context.flags << ", exact " << std::to_string(exact_actual)
                  << " flags " << exact_context.flags
                  << "\n  expected = " << std::to_string(expected_value)
                  << " flags " << expected_flags << ", exact " << exact_flags;
  }
  return ok;
}

template <class F>
void CheckToIntegers(int iterations, limb_t seed) {
  Random random(seed);
  const int total = iterations * TestScale();
  for (int i = 0; i < total; ++i) {
    F x = RandomValue<F>(random);
    if (random.Below(4) != 0 && x.is_finite() && !x.is_zero()) {
      // Mostly magnitudes where integers and their limits live.
      x = scalbn(x, static_cast<int>(random.Below(70)) - x.raw_exponent() - 3);
    }
    if (random.Below(8) == 0) x = round(x) + F(0.5);
    for (const Rounding mode : kModes) {
      if (!CheckToInteger<std::int64_t>(x, mode) ||
          !CheckToInteger<std::uint64_t>(x, mode) ||
          !CheckToInteger<std::int32_t>(x, mode) ||
          !CheckToInteger<std::uint32_t>(x, mode) ||
          !CheckToInteger<std::int8_t>(x, mode) ||
          !CheckToInteger<std::uint8_t>(x, mode)) {
        return;
      }
    }
  }
}

TEST(ConvertOracle, ToIntegers) {
  CheckToIntegers<binary128>(10000, 21);
  CheckToIntegers<binary256>(5000, 22);
  CheckToIntegers<Float<64, 8>>(10000, 23);
  CheckToIntegers<Float<53, 11>>(10000, 24);
}

TEST(ConvertOracle, ToIntegersExhaustively) {
  using Small = Float<8, 5>;
  for (limb_t bits = 0; bits < (limb_t{1} << 13); ++bits) {
    const Small x = Small::from_bits({bits});
    for (const Rounding mode : kModes) {
      if (!CheckToInteger<std::int8_t>(x, mode) ||
          !CheckToInteger<std::uint8_t>(x, mode) ||
          !CheckToInteger<std::int16_t>(x, mode) ||
          !CheckToInteger<std::uint64_t>(x, mode)) {
        return;
      }
    }
  }
}

// --- Decimal input ---

template <class F>
Expected<F> ExpectedParse(const std::string& text, int base, Rounding mode) {
  return RoundExact<F>(
      [&text, base](mpfr_ptr result, mpfr_rnd_t rnd) {
        return mpfr_strtofr(result, text.c_str(), nullptr, base, rnd);
      },
      mode);
}

struct ParseCoverage {
  long checks = 0;
  long exact = 0;
  long overflow = 0;
  long underflow = 0;
  long ties = 0;
};

ParseCoverage parse_coverage;

// `text` must be a number in a syntax both parsers share.
template <class F>
bool CheckParse(const std::string& text, Rounding mode,
                std::chars_format fmt = std::chars_format::general) {
  const bool hex = fmt == std::chars_format::hex;
  Expected<F> expected = ExpectedParse<F>(text, hex ? 16 : 10, mode);
  // A string of zeros is a zero of its sign, which MPFR may not keep.
  if (expected.value.is_zero() && expected.flags == 0) {
    expected.value = F::Zero(text[0] == '-');
  }

  ++parse_coverage.checks;
  if (expected.flags == 0) ++parse_coverage.exact;
  if ((expected.flags & kFlagOverflow) != 0) ++parse_coverage.overflow;
  if ((expected.flags & kFlagUnderflow) != 0) ++parse_coverage.underflow;

  Context context{mode};
  F actual;
  const auto result =
      from_chars(text.data(), text.data() + text.size(), actual, fmt, context);
  bool ok =
      result.ec == std::errc{} && result.ptr == text.data() + text.size() &&
      SameResult(actual, expected.value) && context.flags == expected.flags;
  if (mode == Rounding::kNearestEven && !hex) {
    ok = ok && SameResult(from_string<F>(text), expected.value);
  }
  if (!ok) {
    ADD_FAILURE() << "P=" << F::kPrecision << " E=" << F::kExponentBits
                  << " parse mode=" << static_cast<int>(mode)
                  << "\n  text     = " << text.substr(0, 400)
                  << (text.size() > 400 ? "..." : "")
                  << "\n  actual   = " << Hex(actual) << " flags "
                  << context.flags << "\n  expected = " << Hex(expected.value)
                  << " flags " << expected.flags;
  }
  return ok;
}

template <class F>
bool CheckParseAllModes(const std::string& text,
                        std::chars_format fmt = std::chars_format::general) {
  for (const Rounding mode : kModes) {
    if (!CheckParse<F>(text, mode, fmt)) return false;
  }
  return true;
}

// The decimal exponents of F's overflow threshold and smallest subnormal.
template <class F>
long DecimalEmax() {
  return static_cast<long>(F::kEmax * 0.30103) + 1;
}
template <class F>
long DecimalEmin() {
  return static_cast<long>((F::kEmin - F::kPrecision) * 0.30103) - 1;
}

template <class F>
long RandomDecimalExponent(Random& random) {
  const auto near = [&random](long center, long spread) {
    return center - spread +
           static_cast<long>(random.Below(static_cast<limb_t>(2 * spread + 1)));
  };
  switch (random.Below(8)) {
    case 0:
      return near(DecimalEmax<F>(), 3);
    case 1:
      // Across the subnormal range and just below it.
      return DecimalEmin<F>() - 2 +
             static_cast<long>(
                 random.Below(static_cast<limb_t>(F::kPrecision / 3 + 6)));
    case 2: {
      const long range = DecimalEmax<F>() - DecimalEmin<F>() + 8;
      return DecimalEmin<F>() - 4 +
             static_cast<long>(random.Below(static_cast<limb_t>(range)));
    }
    default:
      return near(0, 60);
  }
}

std::string RandomDigits(Random& random, int count) {
  std::string digits;
  for (int i = 0; i < count; ++i) {
    // Runs of zeros and nines make carries and near-ties.
    const limb_t kind = random.Below(8);
    const char digit = kind == 0   ? '0'
                       : kind == 1 ? '9'
                                   : static_cast<char>('0' + random.Below(10));
    digits.push_back(digit);
  }
  return digits;
}

template <class F>
std::string RandomDecimal(Random& random) {
  const int max_digits =
      random.Below(4) == 0 ? F::kPrecision / 2 + 20 : F::kPrecision / 3 + 4;
  const int count =
      1 + static_cast<int>(random.Below(static_cast<limb_t>(max_digits)));
  std::string digits = RandomDigits(random, count);
  if (digits[0] == '0' && random.Below(2) == 0) digits[0] = '1';
  std::string text = random.Below(2) != 0 ? "-" : "";
  const auto point = static_cast<std::size_t>(
      random.Below(static_cast<limb_t>(digits.size() + 1)));
  text += digits.substr(0, point);
  if (point < digits.size() || random.Below(2) == 0) {
    text += "." + digits.substr(point);
  }
  const long exponent =
      RandomDecimalExponent<F>(random) - static_cast<long>(point) + 1;
  text += (random.Below(2) != 0 ? "e" : "E") + std::to_string(exponent);
  return text;
}

template <class F>
void CheckParseRandomly(int iterations, limb_t seed) {
  Random random(seed);
  const int total = iterations * TestScale();
  parse_coverage = {};
  for (int i = 0; i < total; ++i) {
    if (!CheckParseAllModes<F>(RandomDecimal<F>(random))) return;
  }
  std::printf("P=%d parse checks=%ld exact=%ld overflow=%ld underflow=%ld\n",
              F::kPrecision, parse_coverage.checks, parse_coverage.exact,
              parse_coverage.overflow, parse_coverage.underflow);
  EXPECT_GT(parse_coverage.overflow, parse_coverage.checks / 50);
  EXPECT_GT(parse_coverage.underflow, parse_coverage.checks / 50);
}

TEST(ConvertOracle, ParseBinary128) { CheckParseRandomly<binary128>(4000, 31); }
TEST(ConvertOracle, ParseBinary256) { CheckParseRandomly<binary256>(2000, 32); }
TEST(ConvertOracle, ParseBinary512) { CheckParseRandomly<binary512>(600, 33); }
TEST(ConvertOracle, ParseOtherFormats) {
  CheckParseRandomly<Float<53, 11>>(4000, 34);
  CheckParseRandomly<Float<64, 15>>(4000, 35);
  CheckParseRandomly<Float<24, 8>>(4000, 36);
  CheckParseRandomly<Float<8, 5>>(4000, 37);
}

// Decimal strings on, just above and just below the midpoints between
// neighbours, which decide ties and need every digit.
template <class F>
void CheckParseMidpoints(int iterations, limb_t seed) {
  using Wide = Float<F::kPrecision + 1, F::kExponentBits>;
  Random random(seed);
  const int total = iterations * TestScale();
  long exact_midpoints = 0;
  for (int i = 0; i < total; ++i) {
    F x = RandomValue<F>(random);
    if (random.Below(8) != 0 && x.is_finite() && !x.is_zero()) {
      // Mostly keep the exact decimal expansion short enough to print in
      // full. Far from 1, a string this close to a midpoint needs the exact
      // path, whose integers grow with the exponent.
      x = scalbn(x,
                 static_cast<int>(random.Below(400)) - 200 - x.raw_exponent());
    }
    x = abs(x);
    if (!x.is_finite() || x == F::MaxFinite()) continue;
    // The midpoint is exact one bit wider. Adding half the gap cannot
    // overflow, unlike halving the sum.
    const Wide below = x;
    const Wide midpoint = below + (Wide(next_up(x)) - below) * Wide(0.5);
    char buffer[4000];
    const auto result = to_chars(buffer, buffer + sizeof(buffer), midpoint,
                                 std::chars_format::scientific, 1500);
    std::string text(buffer, result.ptr);
    const std::size_t e = text.find('e');
    if (e == std::string::npos) {
      ADD_FAILURE() << "Unexpected output " << text.substr(0, 100) << " for "
                    << Hex(x);
      return;
    }
    std::string mantissa = text.substr(0, e);
    const std::string suffix = text.substr(e);
    const std::size_t end = mantissa.find_last_not_of('0');
    mantissa.resize(end + 1);
    if (mantissa.back() == '.') mantissa.pop_back();
    if (mantissa.size() < 1400) ++exact_midpoints;
    // Dropping the last digit gives a value just below the midpoint.
    std::string truncated = mantissa.substr(0, mantissa.size() - 1);
    if (!truncated.empty() && truncated.back() == '.') truncated.pop_back();
    const std::string sign = random.Below(2) == 0 ? "-" : "";
    std::vector<std::string> variants = {
        sign + mantissa + suffix,
        sign + mantissa +
            (mantissa.find('.') == std::string::npos ? ".0" : "") +
            "000000000000000000001" + suffix};
    if (!truncated.empty()) variants.push_back(sign + truncated + suffix);
    for (const std::string& variant : variants) {
      if (!CheckParseAllModes<F>(variant)) return;
    }
  }
  std::printf("P=%d exact midpoints=%ld of %d\n", F::kPrecision,
              exact_midpoints, total);
  EXPECT_GT(exact_midpoints, total / 2);
}

TEST(ConvertOracle, ParseMidpoints) {
  CheckParseMidpoints<binary128>(1000, 41);
  CheckParseMidpoints<binary256>(300, 42);
  CheckParseMidpoints<binary512>(120, 43);
  CheckParseMidpoints<Float<53, 11>>(1000, 44);
  CheckParseMidpoints<Float<8, 5>>(1000, 45);
}

TEST(ConvertOracle, ParseHex) {
  Random random(46);
  const int total = 4000 * TestScale();
  for (int i = 0; i < total; ++i) {
    std::string text = random.Below(2) != 0 ? "-" : "";
    const int count = 1 + static_cast<int>(random.Below(150));
    for (int d = 0; d < count; ++d) {
      text.push_back("0123456789abcdefABCDEF"[random.Below(22)]);
      if (d == 0 && random.Below(2) == 0) text += ".";
    }
    std::int64_t exponent =
        static_cast<std::int64_t>(random.Below(40000)) - 20000;
    // Now and then far beyond any format's range. MPFR on Windows reads at
    // most 32-bit exponents; charconv_test covers larger ones.
    if (random.Below(16) == 0) exponent *= 10000;
    text += "p" + std::to_string(exponent);
    if (!CheckParseAllModes<binary128>(text, std::chars_format::hex) ||
        !CheckParseAllModes<Float<8, 5>>(text, std::chars_format::hex)) {
      return;
    }
  }
}

// --- Decimal output ---

// The n significant decimal digits of x rounded in `mode`, as MPFR gives
// them, with value 0.digits * 10^exponent. MPFR has no ties-away mode, so ties
// are detected from n + 1 digits.
struct MpfrDigits {
  std::string digits;
  long exponent;
};

MpfrDigits GetDigits(mpfr_srcptr x, long n, mpfr_rnd_t rnd) {
  std::vector<char> buffer(static_cast<std::size_t>(n) + 8);
  mpfr_exp_t exponent;
  mpfr_get_str(buffer.data(), &exponent, 10, static_cast<std::size_t>(n), x,
               rnd);
  std::string digits = buffer.data();
  if (digits[0] == '-') digits.erase(0, 1);
  return {digits, static_cast<long>(exponent)};
}

bool operator==(const MpfrDigits& a, const MpfrDigits& b) {
  return a.digits == b.digits && a.exponent == b.exponent;
}

MpfrDigits ExpectedDigits(mpfr_srcptr x, long n, Rounding mode, bool& inexact) {
  inexact = !(GetDigits(x, n, MPFR_RNDZ) == GetDigits(x, n, MPFR_RNDA));
  switch (mode) {
    case Rounding::kNearestEven:
      return GetDigits(x, n, MPFR_RNDN);
    case Rounding::kTowardZero:
      return GetDigits(x, n, MPFR_RNDZ);
    case Rounding::kUp:
      return GetDigits(x, n, MPFR_RNDU);
    case Rounding::kDown:
      return GetDigits(x, n, MPFR_RNDD);
    case Rounding::kNearestAway:
      break;
  }
  const MpfrDigits wider = GetDigits(x, n + 1, MPFR_RNDZ);
  const bool tie =
      wider == GetDigits(x, n + 1, MPFR_RNDA) && wider.digits.back() == '5';
  return GetDigits(x, n, tie ? MPFR_RNDA : MPFR_RNDN);
}

// Splits our scientific output, d.ddde+XX, into MPFR's form.
MpfrDigits SplitScientific(const std::string& text) {
  std::size_t start = text[0] == '-' ? 1 : 0;
  const std::size_t e = text.find('e');
  std::string digits;
  for (std::size_t i = start; i < e; ++i) {
    if (text[i] != '.') digits.push_back(text[i]);
  }
  return {digits, std::stol(text.substr(e + 1)) + 1};
}

template <class F>
bool CheckDigits(const F& x, int precision, Rounding mode) {
  const View<F> view(x);
  bool inexact;
  const MpfrDigits expected =
      ExpectedDigits(view.get(), precision + 1, mode, inexact);
  std::string buffer(static_cast<std::size_t>(precision) + 64, '\0');
  Context context{mode};
  const auto result =
      to_chars(buffer.data(), buffer.data() + buffer.size(), x,
               std::chars_format::scientific, precision, context);
  const std::string text(buffer.data(), result.ptr);
  const bool ok = result.ec == std::errc{} &&
                  SplitScientific(text) == expected &&
                  (text[0] == '-') == x.signbit() &&
                  context.flags == (inexact ? unsigned{kFlagInexact} : 0u);
  if (!ok) {
    ADD_FAILURE() << "P=" << F::kPrecision << " to_chars precision "
                  << precision << " mode=" << static_cast<int>(mode)
                  << "\n  x        = " << Hex(x) << "\n  actual   = " << text
                  << " flags " << context.flags << "\n  expected = 0."
                  << expected.digits << "e" << expected.exponent
                  << (inexact ? " inexact" : " exact");
  }
  return ok;
}

template <class F>
F RandomFinite(Random& random) {
  for (;;) {
    F x = RandomValue<F>(random);
    if (!x.is_finite() || x.is_zero()) continue;
    if (random.Below(4) != 0) {
      // Mostly moderate magnitudes, where the work is cheap.
      x = scalbn(x,
                 static_cast<int>(random.Below(800)) - 400 - x.raw_exponent());
    }
    // Tiny formats overflow or vanish when scaled.
    if (x.is_finite() && !x.is_zero()) return x;
  }
}

template <class F>
void CheckDigitsRandomly(int iterations, limb_t seed) {
  Random random(seed);
  const int total = iterations * TestScale();
  for (int i = 0; i < total; ++i) {
    F x = RandomFinite<F>(random);
    if (random.Below(8) == 0) {
      // Exactly representable decimals, which give exact and tied output.
      const int digits = 1 + static_cast<int>(random.Below(30));
      std::string text = RandomDigits(random, digits);
      text = "1" + text + "e" +
             std::to_string(static_cast<int>(random.Below(30)) - 30);
      x = from_string<F>(text);
      // Tiny formats may round it to zero or infinity.
      if (x.is_zero() || !x.is_finite()) continue;
    }
    const int precision = static_cast<int>(
        random.Below(random.Below(4) == 0 ? 300 : F::kPrecision / 3 + 3));
    for (const Rounding mode : kModes) {
      if (!CheckDigits(x, precision, mode)) return;
    }
  }
}

TEST(ConvertOracle, DigitsBinary128) {
  CheckDigitsRandomly<binary128>(3000, 51);
}
TEST(ConvertOracle, DigitsBinary256) {
  CheckDigitsRandomly<binary256>(1500, 52);
}
TEST(ConvertOracle, DigitsBinary512) {
  CheckDigitsRandomly<binary512>(800, 53);
}
TEST(ConvertOracle, DigitsOtherFormats) {
  CheckDigitsRandomly<Float<64, 15>>(3000, 54);
  CheckDigitsRandomly<Float<8, 5>>(3000, 55);
}

TEST(ConvertOracle, DigitsExhaustively) {
  using Small = Float<8, 5>;
  for (limb_t bits = 0; bits < (limb_t{1} << 13); ++bits) {
    const Small x = Small::from_bits({bits});
    if (!x.is_finite() || x.is_zero()) continue;
    for (int precision = 0; precision < 5; ++precision) {
      for (const Rounding mode : kModes) {
        if (!CheckDigits(x, precision, mode)) return;
      }
    }
  }
}

// Shortest output reads back exactly, no shorter decimal does, and among the
// decimals of its length it is the nearest.
template <class F>
bool CheckShortest(const F& x) {
  char buffer[4096];
  const auto result = to_chars(buffer, buffer + sizeof(buffer), x,
                               std::chars_format::scientific);
  const std::string text(buffer, result.ptr);
  const auto reads_back = [&x](const MpfrDigits& digits) {
    const std::string decimal = (x.signbit() ? "-0." : "0.") + digits.digits +
                                "e" + std::to_string(digits.exponent);
    return SameResult(
        ExpectedParse<F>(decimal, 10, Rounding::kNearestEven).value, x);
  };
  const MpfrDigits ours = SplitScientific(text);
  const auto length = static_cast<long>(ours.digits.size());
  const View<F> view(x);
  bool ok = reads_back(ours);
  if (length > 1) {
    ok = ok && !reads_back(GetDigits(view.get(), length - 1, MPFR_RNDD)) &&
         !reads_back(GetDigits(view.get(), length - 1, MPFR_RNDU));
  }
  const MpfrDigits nearest = GetDigits(view.get(), length, MPFR_RNDN);
  if (reads_back(nearest)) ok = ok && ours == nearest;
  // The plain form reads back too.
  ok = ok && SameResult(from_string<F>(to_string(x)), x);
  if (!ok) {
    ADD_FAILURE() << "P=" << F::kPrecision << " shortest"
                  << "\n  x        = " << Hex(x) << "\n  actual   = " << text
                  << "\n  nearest  = 0." << nearest.digits << "e"
                  << nearest.exponent;
  }
  return ok;
}

template <class F>
void CheckShortestRandomly(int iterations, limb_t seed) {
  Random random(seed);
  const int total = iterations * TestScale();
  for (int i = 0; i < total; ++i) {
    if (!CheckShortest(RandomFinite<F>(random))) return;
  }
}

TEST(ConvertOracle, ShortestBinary128) {
  CheckShortestRandomly<binary128>(3000, 61);
}
TEST(ConvertOracle, ShortestBinary256) {
  CheckShortestRandomly<binary256>(1500, 62);
}
TEST(ConvertOracle, ShortestBinary512) {
  CheckShortestRandomly<binary512>(600, 63);
}
TEST(ConvertOracle, ShortestExhaustively) {
  using Small = Float<8, 5>;
  for (limb_t bits = 0; bits < (limb_t{1} << 13); ++bits) {
    const Small x = Small::from_bits({bits});
    if (x.is_finite() && !x.is_zero() && !CheckShortest(x)) return;
  }
  using Tiny = Float<3, 3>;
  for (limb_t bits = 0; bits < (limb_t{1} << 6); ++bits) {
    const Tiny x = Tiny::from_bits({bits});
    if (x.is_finite() && !x.is_zero() && !CheckShortest(x)) return;
  }
}

// --- Hexadecimal output ---

// Rounding to `precision` hexadecimal digits is rounding to 4 * precision + 1
// bits. The exponent range is widened, since hexadecimal output keeps
// subnormals normalized.
template <class F, int kDigits>
bool CheckHexDigits(const F& x, Rounding mode) {
  using Rounded = Float<4 * kDigits + 1, F::kExponentBits + 2>;
  Context expected_context{mode};
  const Rounded expected = convert<Rounded>(x, expected_context);
  char buffer[512];
  Context context{mode};
  const auto result = to_chars(buffer, buffer + sizeof(buffer), x,
                               std::chars_format::hex, kDigits, context);
  Rounded actual;
  from_chars(buffer, result.ptr, actual, std::chars_format::hex);
  const bool ok =
      SameResult(actual, expected) && context.flags == expected_context.flags;
  if (!ok) {
    ADD_FAILURE() << "P=" << F::kPrecision << " hex digits " << kDigits
                  << " mode=" << static_cast<int>(mode)
                  << "\n  x      = " << Hex(x)
                  << "\n  actual = " << std::string(buffer, result.ptr);
  }
  return ok;
}

TEST(ConvertOracle, HexDigits) {
  Random random(71);
  const int total = 3000 * TestScale();
  for (int i = 0; i < total; ++i) {
    const binary128 x = RandomFinite<binary128>(random);
    for (const Rounding mode : kModes) {
      if (!CheckHexDigits<binary128, 1>(x, mode) ||
          !CheckHexDigits<binary128, 5>(x, mode) ||
          !CheckHexDigits<binary128, 16>(x, mode) ||
          !CheckHexDigits<binary128, 27>(x, mode)) {
        return;
      }
    }
  }
}

}  // namespace
}  // namespace hwymp
