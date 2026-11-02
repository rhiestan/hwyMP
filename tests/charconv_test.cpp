// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// String conversions. Float<53, 11> has the values of double, so the C and C++
// libraries' correctly rounded double conversions serve as a reference here.
// The differential tests against MPFR, for the wide formats, are in
// convert_oracle_test.cpp.

#include "mpfix/charconv.h"

#include <gtest/gtest.h>

#include <bit>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <string>

#include "mpfix/bigint.h"

namespace mpfix {
namespace {

using detail::BigUint;
using Double = Float<53, 11>;

// --- BigUint ---

BigUint RandomBig(std::mt19937_64& random, std::size_t limbs) {
  std::vector<limb_t> values(limbs);
  for (limb_t& value : values) value = random();
  // Sparse limbs reach the carry paths more often.
  if (random() % 4 == 0) values[random() % limbs] = ~limb_t{0};
  if (random() % 4 == 0) values[random() % limbs] = 0;
  // Never zero, and of the requested size.
  if (values.back() == 0) values.back() = 1;
  return BigUint::FromLimbs(values.data(), values.size());
}

TEST(BigUint, MultiplicationIdentities) {
  std::mt19937_64 random(1);
  for (int i = 0; i < 300; ++i) {
    // Sizes either side of the Karatsuba threshold, balanced and not.
    const std::size_t an = 1 + random() % 150;
    const std::size_t bn = 1 + random() % 150;
    const BigUint a = RandomBig(random, an);
    const BigUint b = RandomBig(random, bn);
    const BigUint product = a * b;
    ASSERT_EQ(Compare(product, b * a), 0);
    // (a + 1) * b = a * b + b.
    BigUint a1 = a;
    a1.Add(BigUint(1));
    BigUint expected = product;
    expected.Add(b);
    ASSERT_EQ(Compare(a1 * b, expected), 0);
    // Division undoes the multiplication, also with a remainder.
    BigUint c = RandomBig(random, std::max<std::size_t>(1, bn / 2));
    if (Compare(c, b) >= 0) c = BigUint();
    BigUint dividend = product;
    dividend.Add(c);
    BigUint quotient;
    BigUint remainder;
    BigUint::DivRem(dividend, b, quotient, remainder);
    ASSERT_EQ(Compare(quotient, a), 0) << an << " " << bn;
    ASSERT_EQ(Compare(remainder, c), 0);
  }
}

TEST(BigUint, KaratsubaMatchesSchoolbook) {
  // Squaring (B^n - 1) gives B^2n - 2 B^n + 1, whose limbs are known.
  for (std::size_t n : {31u, 32u, 33u, 64u, 100u, 257u}) {
    std::vector<limb_t> ones(n, ~limb_t{0});
    const BigUint x = BigUint::FromLimbs(ones.data(), n);
    const BigUint square = x * x;
    ASSERT_EQ(square.limbs().size(), 2 * n);
    EXPECT_EQ(square.limbs()[0], 1u);
    for (std::size_t i = 1; i < n; ++i) ASSERT_EQ(square.limbs()[i], 0u);
    EXPECT_EQ(square.limbs()[n], ~limb_t{1});
    for (std::size_t i = n + 1; i < 2 * n; ++i) {
      ASSERT_EQ(square.limbs()[i], ~limb_t{0});
    }
  }
}

TEST(BigUint, PowersOfFive) {
  BigUint power(1);
  for (std::uint64_t k = 0; k < 400; ++k) {
    ASSERT_EQ(Compare(BigUint::Pow5(k), power), 0) << k;
    power.MulAdd(5, 0);
  }
}

TEST(BigUint, SmallDivisionAndShifts) {
  std::mt19937_64 random(2);
  for (int i = 0; i < 1000; ++i) {
    const BigUint a = RandomBig(random, 1 + random() % 10);
    limb_t divisor = random() >> (random() % 64);
    if (divisor == 0) divisor = 3;
    BigUint q = a;
    const limb_t r = q.DivSmall(divisor);
    ASSERT_LT(r, divisor);
    q.MulAdd(divisor, r);
    ASSERT_EQ(Compare(q, a), 0);

    const std::size_t shift = random() % 300;
    BigUint shifted = a;
    shifted.ShiftLeft(shift);
    ASSERT_EQ(shifted.bit_length(), a.bit_length() + shift);
    ASSERT_EQ(shifted.any_bit_below(shift), false);
    shifted.ShiftRight(shift);
    ASSERT_EQ(Compare(shifted, a), 0);
  }
}

// --- Output ---

std::string ToChars(const auto& x) {
  char buffer[2048];
  const auto result = to_chars(buffer, buffer + sizeof(buffer), x);
  EXPECT_EQ(result.ec, std::errc{});
  return std::string(buffer, result.ptr);
}

std::string ToChars(const auto& x, std::chars_format fmt, int precision) {
  std::string buffer(5000, '\0');
  const auto result =
      to_chars(buffer.data(), buffer.data() + buffer.size(), x, fmt, precision);
  EXPECT_EQ(result.ec, std::errc{});
  return std::string(buffer.data(), result.ptr);
}

std::string Printf(const char* format, int precision, double x) {
  char buffer[2048];
  std::snprintf(buffer, sizeof(buffer), format, precision, x);
  return buffer;
}

double RandomDouble(std::mt19937_64& random) {
  for (;;) {
    const double x = std::bit_cast<double>(random());
    if (std::isfinite(x)) return x;
  }
}

TEST(Charconv, PrecisionMatchesPrintf) {
  std::mt19937_64 random(3);
  for (int i = 0; i < 20000; ++i) {
    const double x = RandomDouble(random);
    const Double y = x;
    const int precision = static_cast<int>(random() % 40);
    ASSERT_EQ(ToChars(y, std::chars_format::scientific, precision),
              Printf("%.*e", precision, x));
    ASSERT_EQ(ToChars(y, std::chars_format::general, precision),
              Printf("%.*g", precision, x));
    // Fixed notation only for moderate values, to keep the strings short.
    const double moderate =
        std::ldexp(x, -std::ilogb(x) + static_cast<int>(random() % 80) - 40);
    ASSERT_EQ(ToChars(Double(moderate), std::chars_format::fixed, precision),
              Printf("%.*f", precision, moderate));
  }
}

TEST(Charconv, PrecisionKnownValues) {
  EXPECT_EQ(ToChars(Double(0.5), std::chars_format::fixed, 0), "0");
  EXPECT_EQ(ToChars(Double(1.5), std::chars_format::fixed, 0), "2");
  EXPECT_EQ(ToChars(Double(2.5), std::chars_format::fixed, 0), "2");
  EXPECT_EQ(ToChars(Double(9.5), std::chars_format::scientific, 0), "1e+01");
  EXPECT_EQ(ToChars(Double(0.0), std::chars_format::scientific, 3),
            "0.000e+00");
  EXPECT_EQ(ToChars(Double(-0.0), std::chars_format::general, 6), "-0");
  EXPECT_EQ(ToChars(Double(1e-5), std::chars_format::general, 6), "1e-05");
  EXPECT_EQ(ToChars(Double(1e-4), std::chars_format::general, 6), "0.0001");
  EXPECT_EQ(ToChars(Double(123456), std::chars_format::general, 6), "123456");
  EXPECT_EQ(ToChars(Double(1234567), std::chars_format::general, 6),
            "1.23457e+06");
  EXPECT_EQ(ToChars(binary128::Infinity(true), std::chars_format::fixed, 2),
            "-inf");
  EXPECT_EQ(ToChars(binary128::QuietNaN(), std::chars_format::fixed, 2), "nan");

  // Directed rounding, with inexact.
  Context context{Rounding::kUp, 0};
  char buffer[64];
  auto result = to_chars(buffer, buffer + 64, binary128(1) / binary128(3),
                         std::chars_format::scientific, 3, context);
  EXPECT_EQ(std::string(buffer, result.ptr), "3.334e-01");
  EXPECT_EQ(context.flags, unsigned{kFlagInexact});
  context = {Rounding::kDown, 0};
  result = to_chars(buffer, buffer + 64, -binary128(2) / binary128(3),
                    std::chars_format::fixed, 2, context);
  EXPECT_EQ(std::string(buffer, result.ptr), "-0.67");
  context = {Rounding::kNearestAway, 0};
  result = to_chars(buffer, buffer + 64, binary128(0.125),
                    std::chars_format::fixed, 2, context);
  EXPECT_EQ(std::string(buffer, result.ptr), "0.13");
  context = {Rounding::kNearestEven, 0};
  result = to_chars(buffer, buffer + 64, binary128(0.125),
                    std::chars_format::fixed, 2, context);
  EXPECT_EQ(std::string(buffer, result.ptr), "0.12");
  context = {};
  result = to_chars(buffer, buffer + 64, binary128(0.125),
                    std::chars_format::fixed, 3, context);
  EXPECT_EQ(std::string(buffer, result.ptr), "0.125");
  EXPECT_EQ(context.flags, 0u);

  // Too small a buffer.
  result = to_chars(buffer, buffer + 3, binary128(0.125));
  EXPECT_EQ(result.ec, std::errc::value_too_large);
  EXPECT_EQ(result.ptr, buffer + 3);
}

TEST(Charconv, Hex) {
  EXPECT_EQ(ToChars(Double(1.0), std::chars_format::hex, -1), "1p+0");
  EXPECT_EQ(ToChars(Double(-3.0), std::chars_format::hex, 2), "-1.80p+1");
  // A tie rounds the odd leading 1 up, and 2 is written as 1p+1.
  EXPECT_EQ(ToChars(Double(1.5), std::chars_format::hex, 0), "1p+1");
  EXPECT_EQ(ToChars(Double(1.25), std::chars_format::hex, 0), "1p+0");
  EXPECT_EQ(ToChars(Double(1.75), std::chars_format::hex, 0), "1p+1");
  EXPECT_EQ(ToChars(Double(0.0), std::chars_format::hex, -1), "0p+0");
  char buffer[256];
  // 1/10 has 28 hexadecimal digits after the point, and the last rounds up.
  auto result = to_chars(buffer, buffer + 256, binary128(1) / binary128(10),
                         std::chars_format::hex);
  EXPECT_EQ(std::string(buffer, result.ptr),
            "1." + std::string(27, '9') + "ap-4");
  result = to_chars(buffer, buffer + 256,
                    Double(std::numeric_limits<double>::denorm_min()),
                    std::chars_format::hex);
  EXPECT_EQ(std::string(buffer, result.ptr), "1p-1074");

  // Hexadecimal output reads back exactly, in every format.
  std::mt19937_64 random(4);
  for (int i = 0; i < 2000; ++i) {
    const binary256 x =
        binary256(RandomDouble(random)) / binary256(RandomDouble(random) + 3.0);
    result = to_chars(buffer, buffer + 256, x, std::chars_format::hex);
    binary256 back;
    const auto parsed =
        from_chars(buffer, result.ptr, back, std::chars_format::hex);
    ASSERT_EQ(parsed.ptr, result.ptr);
    ASSERT_EQ(back.to_bits(), x.to_bits()) << std::string(buffer, result.ptr);
  }
}

TEST(Charconv, ShortestRoundTrips) {
  std::mt19937_64 random(5);
  for (int i = 0; i < 20000; ++i) {
    const double x = RandomDouble(random);
    const std::string text = ToChars(Double(x));
    // strtod is correctly rounded on all the supported platforms.
    ASSERT_EQ(std::strtod(text.c_str(), nullptr), x) << text;
#if defined(__cpp_lib_to_chars)
    char buffer[64];
    const auto result = std::to_chars(buffer, buffer + 64, x);
    ASSERT_EQ(text, std::string(buffer, result.ptr));
    for (const auto fmt :
         {std::chars_format::fixed, std::chars_format::scientific,
          std::chars_format::general}) {
      char ours[2048];
      const auto mine = to_chars(ours, ours + sizeof(ours), Double(x), fmt);
      char reference[2048];
      const auto theirs =
          std::to_chars(reference, reference + sizeof(reference), x, fmt);
      ASSERT_EQ(std::string(ours, mine.ptr),
                std::string(reference, theirs.ptr));
    }
#endif
  }
  EXPECT_EQ(ToChars(Double(0.1)), "0.1");
  EXPECT_EQ(ToChars(Double(1e23)), "1e+23");
  EXPECT_EQ(ToChars(Double(5e-324)), "5e-324");
  EXPECT_EQ(ToChars(Double(-0.0)), "-0");
  EXPECT_EQ(ToChars(Double(123.0)), "123");
  EXPECT_EQ(ToChars(Double(1e21)), "1e+21");
  EXPECT_EQ(ToChars(Double(std::numeric_limits<double>::max())),
            "1.7976931348623157e+308");
  EXPECT_EQ(to_string(binary128(0.1)), "0.1000000000000000055511151231257827");
  EXPECT_EQ(to_string(binary128(1) / binary128(3)),
            "0.3333333333333333333333333333333333");
  EXPECT_EQ(to_string(binary128::MaxFinite()),
            "1.189731495357231765085759326628007e+4932");
  EXPECT_EQ(to_string(binary128::MinSubnormal()), "6e-4966");
}

TEST(Charconv, ExtremesOfBinary512) {
  for (const binary512& x :
       {binary512::MaxFinite(), binary512::MinNormal(),
        binary512::MinSubnormal(), next_down(binary512::MaxFinite()),
        next_up(binary512::MinSubnormal())}) {
    const std::string text = to_string(x);
    const binary512 back = from_string<binary512>(text);
    ASSERT_EQ(back.to_bits(), x.to_bits()) << text;
  }
}

// --- Input ---

TEST(Charconv, ParseMatchesStrtod) {
  std::mt19937_64 random(6);
  for (int i = 0; i < 50000; ++i) {
    // Random digit strings with random exponents, including long ones and
    // ones near the edges of the range.
    std::string text;
    if (random() % 2 == 0) text += "-";
    const int digits =
        1 + static_cast<int>(random() % (random() % 4 == 0 ? 60 : 20));
    const int point = static_cast<int>(random() % (digits + 1));
    for (int d = 0; d < digits; ++d) {
      if (d == point && d != 0) text += ".";
      text += static_cast<char>('0' + random() % 10);
    }
    text += "e" + std::to_string(static_cast<int>(random() % 700) - 350);
    const double expected = std::strtod(text.c_str(), nullptr);
    Double value;
    const auto result =
        from_chars(text.data(), text.data() + text.size(), value);
    ASSERT_EQ(result.ec, std::errc{});
    ASSERT_EQ(result.ptr, text.data() + text.size());
    ASSERT_EQ(std::bit_cast<std::uint64_t>(static_cast<double>(value)),
              std::bit_cast<std::uint64_t>(expected))
        << text;
  }
}

TEST(Charconv, ParseWithThePointAnywhere) {
  // The digits are read in chunks across the point, so every split of the
  // same digits, with the exponent compensating, must give the same value.
  std::mt19937_64 random(8);
  for (int i = 0; i < 300; ++i) {
    const int length = 1 + static_cast<int>(random() % 70);
    std::string digits;
    for (int d = 0; d < length; ++d) {
      digits +=
          static_cast<char>('0' + random() % (random() % 3 == 0 ? 2 : 10));
    }
    const std::string reference_text = digits + "e-" + std::to_string(length);
    const binary256 reference = from_string<binary256>(reference_text);
    for (int point = 0; point <= length; ++point) {
      const std::string text =
          digits.substr(0, static_cast<std::size_t>(point)) + "." +
          digits.substr(static_cast<std::size_t>(point)) + "e-" +
          std::to_string(point);
      ASSERT_EQ(from_string<binary256>(text).to_bits(), reference.to_bits())
          << text;
    }
  }
}

TEST(Charconv, ParseHalfwayCases) {
  // Exact midpoints between adjacent doubles, which strtod rounds to even,
  // and the same plus or minus a tiny amount.
  std::mt19937_64 random(7);
  for (int i = 0; i < 3000; ++i) {
    const double x = std::abs(RandomDouble(random));
    if (x == std::numeric_limits<double>::max()) continue;
    const binary128 midpoint =
        (binary128(x) + binary128(std::nextafter(x, INFINITY))) / binary128(2);
    // binary128 holds the midpoint exactly, and enough digits print it exactly.
    const std::string exact =
        ToChars(midpoint, std::chars_format::scientific, 800);
    Double parsed;
    from_chars(exact.data(), exact.data() + exact.size(), parsed);
    ASSERT_EQ(static_cast<double>(parsed), std::strtod(exact.c_str(), nullptr))
        << exact;
    const std::string above =
        exact.substr(0, exact.find('e')) + "1" + exact.substr(exact.find('e'));
    from_chars(above.data(), above.data() + above.size(), parsed);
    ASSERT_EQ(static_cast<double>(parsed), std::nextafter(x, INFINITY));
  }
}

TEST(Charconv, ParseSyntax) {
  const auto parse = [](std::string_view text,
                        std::chars_format fmt = std::chars_format::general) {
    Double value = 12345.0;
    const auto result =
        from_chars(text.data(), text.data() + text.size(), value, fmt);
    return std::pair{result.ec == std::errc{}
                         ? static_cast<int>(result.ptr - text.data())
                         : -1,
                     static_cast<double>(value)};
  };
  EXPECT_EQ(parse("1.5"), std::pair(3, 1.5));
  EXPECT_EQ(parse("-1.5e3x"), std::pair(6, -1500.0));
  EXPECT_EQ(parse(".5"), std::pair(2, 0.5));
  EXPECT_EQ(parse("5."), std::pair(2, 5.0));
  EXPECT_EQ(parse("."), std::pair(-1, 12345.0));
  EXPECT_EQ(parse("+1"), std::pair(-1, 12345.0));
  EXPECT_EQ(parse(" 1"), std::pair(-1, 12345.0));
  EXPECT_EQ(parse("1e"), std::pair(1, 1.0));
  EXPECT_EQ(parse("1e+"), std::pair(1, 1.0));
  EXPECT_EQ(parse("1e5", std::chars_format::fixed), std::pair(1, 1.0));
  EXPECT_EQ(parse("1", std::chars_format::scientific), std::pair(-1, 12345.0));
  EXPECT_EQ(parse("1e2", std::chars_format::scientific), std::pair(3, 100.0));
  EXPECT_EQ(parse("1.8p1", std::chars_format::hex), std::pair(5, 3.0));
  EXPECT_EQ(parse("ff", std::chars_format::hex), std::pair(2, 255.0));
  EXPECT_EQ(parse("0x1", std::chars_format::hex), std::pair(1, 0.0));
  EXPECT_EQ(parse("000000000000000000000000000000001.25"), std::pair(36, 1.25));
  EXPECT_EQ(parse("0.0000000000000000000000000000000000125e35"),
            std::pair(42, 1.25));
  EXPECT_EQ(parse("1e99999999999999999999999"),
            std::pair(25, std::numeric_limits<double>::infinity()));
  EXPECT_EQ(parse("1e-99999999999999999999999"), std::pair(26, 0.0));
  EXPECT_EQ(parse("0e99999999999999999999999"), std::pair(25, 0.0));
  EXPECT_EQ(parse("InFiNiTy"),
            std::pair(8, std::numeric_limits<double>::infinity()));
  EXPECT_EQ(parse("-inf"),
            std::pair(4, -std::numeric_limits<double>::infinity()));
  EXPECT_EQ(parse("infinit").first, 3);
  EXPECT_EQ(parse("nan(abc_1)").first, 10);
  EXPECT_EQ(parse("nan(abc").first, 3);
  EXPECT_TRUE(std::isnan(parse("NaN").second));

  // A string of many digits, rounded only at the far end.
  std::string many = "1." + std::string(5000, '0') + "1";
  Context context{Rounding::kUp, 0};
  Double value;
  from_chars(many.data(), many.data() + many.size(), value,
             std::chars_format::general, context);
  EXPECT_EQ(static_cast<double>(value), std::nextafter(1.0, 2.0));
  EXPECT_EQ(context.flags, unsigned{kFlagInexact});

  context = {};
  const std::string huge = "1e400";
  from_chars(huge.data(), huge.data() + huge.size(), value,
             std::chars_format::general, context);
  EXPECT_TRUE(value.is_inf());
  EXPECT_EQ(context.flags, kFlagOverflow | kFlagInexact);
  // Binary exponents beyond 32 bits, in both directions. 2^32 would wrap to 0
  // in 32-bit arithmetic.
  for (const char* text : {"1p4294967296", "1p-4294967296"}) {
    context = {Rounding::kTowardZero, 0};
    const std::string_view view(text);
    from_chars(view.data(), view.data() + view.size(), value,
               std::chars_format::hex, context);
    if (text[2] == '-') {
      EXPECT_TRUE(value.is_zero()) << text;
      EXPECT_EQ(context.flags, kFlagUnderflow | kFlagInexact) << text;
    } else {
      EXPECT_EQ(value, Double::MaxFinite()) << text;
      EXPECT_EQ(context.flags, kFlagOverflow | kFlagInexact) << text;
    }
  }

  context = {};
  const std::string tiny = "1e-400";
  from_chars(tiny.data(), tiny.data() + tiny.size(), value,
             std::chars_format::general, context);
  EXPECT_TRUE(value.is_zero());
  EXPECT_EQ(context.flags, kFlagUnderflow | kFlagInexact);
}

TEST(Charconv, FromString) {
  EXPECT_EQ(from_string<binary128>("+2.5"), binary128(2.5));
  EXPECT_EQ(from_string<binary128>("-0x1.8p1"), binary128(-3));
  EXPECT_EQ(from_string<binary128>("0X10"), binary128(16));
  Context context;
  EXPECT_TRUE(from_string<binary128>("1.5x", context).is_nan());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  EXPECT_TRUE(from_string<binary128>("+-1").is_nan());
  EXPECT_TRUE(from_string<binary128>("").is_nan());
  EXPECT_TRUE(from_string<binary128>("-inf").is_inf());
  EXPECT_EQ(from_string<binary256>("0.1"), binary256(1) / binary256(10));
}

// The example in README.md.
TEST(Charconv, ReadmeExample) {
  using mpfix::binary256;

  binary256 x = 1.5;
  binary256 third = binary256(1) / binary256(3);
  double d = static_cast<double>(third);

  mpfix::Context context{mpfix::Rounding::kUp};
  float f = mpfix::convert<float>(third, context);
  long n = mpfix::to_integer<long>(x, mpfix::Rounding::kNearestEven);

  std::string text = mpfix::to_string(third);
  binary256 y = mpfix::from_string<binary256>("0.1");
  std::ostringstream out;
  out << std::setprecision(40) << y << '\n';

  EXPECT_EQ(d, 1.0 / 3.0);
  // The nearest float is above 1/3, so rounding up gives it too.
  EXPECT_EQ(f, 1.0f / 3.0f);
  EXPECT_EQ(n, 2);
  EXPECT_EQ(mpfix::from_string<binary256>(text), third);
  // %g drops trailing zeros.
  EXPECT_EQ(out.str(), "0.1\n");
}

TEST(Charconv, Streams) {
  std::ostringstream out;
  out << binary128(1) / binary128(3);
  EXPECT_EQ(out.str(), "0.333333");
  out.str("");
  out.precision(34);
  out << binary128(1) / binary128(3);
  EXPECT_EQ(out.str(), "0.3333333333333333333333333333333333");
  out.str("");
  out << std::scientific << std::setprecision(3) << std::uppercase
      << binary128(12345);
  EXPECT_EQ(out.str(), "1.234E+04");
  out.str("");
  out << std::fixed << std::nouppercase << std::showpos << binary128(2.5);
  EXPECT_EQ(out.str(), "+2.500");
  out.str("");
  out << std::hexfloat << std::noshowpos << binary128(-3);
  EXPECT_EQ(out.str(), "-0x1.8p+1");
  out.str("");
  out << std::defaultfloat << std::setw(8) << binary128(1.5);
  EXPECT_EQ(out.str(), "     1.5");

  std::istringstream in("1.25 -0x1p3 bad");
  binary128 a;
  binary128 b;
  binary128 c = 7;
  in >> a >> b;
  EXPECT_EQ(a, binary128(1.25));
  EXPECT_EQ(b, binary128(-8));
  in >> c;
  EXPECT_TRUE(in.fail());
  EXPECT_EQ(c, binary128(7));
}

}  // namespace
}  // namespace mpfix
