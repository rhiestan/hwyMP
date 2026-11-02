// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Conversions between formats, built-in floating-point types and integers. The
// differential tests against MPFR are in convert_oracle_test.cpp.

#include <gtest/gtest.h>

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <type_traits>

#include "mpfix/float.h"

namespace mpfix {
namespace {

using Double = Float<53, 11>;
using Single = Float<24, 8>;
using Small = Float<4, 3>;

// Exact conversions are implicit, and rounding ones explicit.
static_assert(std::is_convertible_v<double, binary128>);
static_assert(std::is_convertible_v<float, binary128>);
static_assert(std::is_convertible_v<std::int64_t, binary128>);
static_assert(std::is_convertible_v<std::uint64_t, binary128>);
static_assert(std::is_convertible_v<int, binary128>);
static_assert(std::is_convertible_v<binary128, binary256>);
static_assert(std::is_convertible_v<binary256, binary512>);
static_assert(!std::is_convertible_v<binary256, binary128>);
static_assert(std::is_constructible_v<binary128, binary256>);
static_assert(!std::is_convertible_v<double, Small>);
static_assert(std::is_constructible_v<Small, double>);
static_assert(!std::is_convertible_v<binary128, double>);
static_assert(std::is_constructible_v<double, binary128>);
static_assert(!std::is_convertible_v<bool, binary128>);
static_assert(std::is_trivially_copyable_v<binary128>);
static_assert(std::is_same_v<detail::NativeFormat<double>, Double>);
static_assert(std::is_same_v<detail::NativeFormat<float>, Single>);

template <class T>
T Bits(double x) {
  return std::bit_cast<T>(x);
}

TEST(Convert, DoubleRoundTripsThroughEveryFormat) {
  std::mt19937_64 random(1);
  for (int i = 0; i < 200000; ++i) {
    const std::uint64_t bits = random();
    const double x = std::bit_cast<double>(bits);
    const binary128 wide = x;
    const double back = static_cast<double>(wide);
    if (std::isnan(x)) {
      // Signaling NaNs come back quiet; the payload is otherwise kept.
      EXPECT_EQ(std::bit_cast<std::uint64_t>(back),
                bits | (std::uint64_t{1} << 51));
      continue;
    }
    ASSERT_EQ(std::bit_cast<std::uint64_t>(back), bits) << x;
    ASSERT_EQ(std::bit_cast<std::uint64_t>(
                  static_cast<double>(binary512(binary256(wide)))),
              bits);
  }
}

TEST(Convert, FloatRoundTrips) {
  std::mt19937 random(2);
  for (int i = 0; i < 200000; ++i) {
    const std::uint32_t bits = random();
    const float x = std::bit_cast<float>(bits);
    if (std::isnan(x)) continue;
    const binary128 wide = x;
    ASSERT_EQ(std::bit_cast<std::uint32_t>(static_cast<float>(wide)), bits);
    ASSERT_EQ(static_cast<double>(wide), static_cast<double>(x));
  }
}

TEST(Convert, DoubleToFloatMatchesHardware) {
  // The host rounds to nearest-even, as the conversion operator does.
  std::mt19937_64 random(3);
  for (int i = 0; i < 200000; ++i) {
    const double x = std::bit_cast<double>(random());
    if (std::isnan(x)) continue;
    const float expected = static_cast<float>(x);
    const float actual = static_cast<float>(Double(x));
    ASSERT_EQ(std::bit_cast<std::uint32_t>(actual),
              std::bit_cast<std::uint32_t>(expected))
        << x;
  }
}

TEST(Convert, NarrowingRounds) {
  // 1 + 2^-60 lies just above the double 1.
  const binary128 one = 1.0;
  const binary128 tiny = std::ldexp(1.0, -60);
  const binary128 x = one + tiny;
  EXPECT_EQ(static_cast<double>(x), 1.0);
  Context up{Rounding::kUp, 0};
  EXPECT_EQ(convert<double>(x, up), std::nextafter(1.0, 2.0));
  EXPECT_EQ(up.flags, unsigned{kFlagInexact});

  // A tie rounds to even, and away in ties-away mode.
  const binary128 tie = one + binary128(std::ldexp(1.0, -53));
  EXPECT_EQ(static_cast<double>(tie), 1.0);
  Context away{Rounding::kNearestAway, 0};
  EXPECT_EQ(convert<double>(tie, away), std::nextafter(1.0, 2.0));

  // Overflow and underflow.
  Context context;
  const binary128 huge = binary128(std::numeric_limits<double>::max()) * 2.0;
  EXPECT_EQ(convert<double>(huge, context),
            std::numeric_limits<double>::infinity());
  EXPECT_EQ(context.flags, kFlagOverflow | kFlagInexact);
  context = {Rounding::kTowardZero, 0};
  EXPECT_EQ(convert<double>(-huge, context),
            -std::numeric_limits<double>::max());

  context = {};
  const binary128 small =
      binary128(std::numeric_limits<double>::denorm_min()) * binary128(0.75);
  EXPECT_EQ(convert<double>(small, context),
            std::numeric_limits<double>::denorm_min());
  EXPECT_EQ(context.flags, kFlagUnderflow | kFlagInexact);
  context = {};
  EXPECT_EQ(convert<double>(small * binary128(0.5), context), 0.0);
  EXPECT_EQ(context.flags, kFlagUnderflow | kFlagInexact);

  // Exact subnormals raise nothing.
  context = {};
  const double subnormal = std::numeric_limits<double>::denorm_min() * 12345;
  EXPECT_EQ(convert<double>(binary128(subnormal), context), subnormal);
  EXPECT_EQ(context.flags, 0u);
}

TEST(Convert, NaNs) {
  Context context;
  const binary128 signaling = binary128::SignalingNaN(true);
  const double quieted = convert<double>(signaling, context);
  EXPECT_TRUE(std::isnan(quieted));
  EXPECT_TRUE(std::signbit(quieted));
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});

  // A payload in the top bits survives narrowing and widening.
  binary512 nan = binary512::QuietNaN();
  auto bits = nan.to_bits();
  bits[7] |= limb_t{0x5A5A} << 40;  // Bits 488 down: the top of the payload.
  nan = binary512::from_bits(bits);
  const binary128 narrow(nan);
  EXPECT_TRUE(narrow.is_nan());
  EXPECT_FALSE(narrow.is_signaling());
  const binary512 wide = narrow;
  EXPECT_EQ(wide.to_bits(), nan.to_bits());

  context = {};
  const binary128 quiet = convert<binary128>(nan, context);
  EXPECT_TRUE(quiet.is_nan());
  EXPECT_EQ(context.flags, 0u);
}

TEST(Convert, SpecialValues) {
  EXPECT_TRUE(binary128(-0.0).is_zero());
  EXPECT_TRUE(binary128(-0.0).signbit());
  EXPECT_TRUE(binary128(std::numeric_limits<double>::infinity()).is_inf());
  EXPECT_TRUE(std::signbit(static_cast<double>(binary256::Zero(true))));
  EXPECT_EQ(static_cast<float>(binary256::Infinity(true)),
            -std::numeric_limits<float>::infinity());
  EXPECT_TRUE(binary128(0).is_zero());
  EXPECT_FALSE(binary128(0).signbit());
}

TEST(Convert, Integers) {
  EXPECT_EQ(binary128(1), binary128(1.0));
  EXPECT_EQ(binary128(-7), binary128(-7.0));
  const binary128 most_negative = std::numeric_limits<std::int64_t>::min();
  EXPECT_EQ(most_negative, -binary128(std::ldexp(1.0, 63)));
  const binary128 most_positive = std::numeric_limits<std::uint64_t>::max();
  EXPECT_EQ(most_positive + binary128(1), binary128(std::ldexp(1.0, 64)));
  EXPECT_EQ(static_cast<std::uint64_t>(most_positive),
            std::numeric_limits<std::uint64_t>::max());
  EXPECT_EQ(static_cast<std::int64_t>(most_negative),
            std::numeric_limits<std::int64_t>::min());

  // Rounding to a narrow format, whose largest value is 240.
  using Narrow = Float<4, 4>;
  Context context;
  EXPECT_EQ(convert<Narrow>(17, context), Narrow(16.0));
  EXPECT_EQ(context.flags, unsigned{kFlagInexact});
  EXPECT_EQ(convert<Narrow>(19, context), Narrow(20.0));
  context = {};
  EXPECT_TRUE(convert<Narrow>(1000, context).is_inf());
  EXPECT_EQ(context.flags, kFlagOverflow | kFlagInexact);
  using Narrower = Float<3, 4>;
  EXPECT_EQ(Narrower(static_cast<unsigned char>(9)), Narrower(8.0));
  EXPECT_EQ(Narrower(static_cast<signed char>(-11)), Narrower(-12.0));

  std::mt19937_64 random(4);
  for (int i = 0; i < 100000; ++i) {
    const auto value = static_cast<std::int64_t>(random());
    const binary128 x = value;
    ASSERT_EQ(static_cast<std::int64_t>(x), value);
    ASSERT_EQ(to_integer<std::int64_t>(x, Rounding::kNearestEven), value);
    // A double conversion of the same integer rounds the same way.
    ASSERT_EQ(static_cast<double>(x), static_cast<double>(value));
    ASSERT_EQ(static_cast<double>(Double(value)), static_cast<double>(value));
    const auto narrow = static_cast<std::int32_t>(value);
    ASSERT_EQ(static_cast<std::int32_t>(binary128(narrow)), narrow);
  }
}

TEST(Convert, ToInteger) {
  const binary128 x = 2.5;
  Context context;
  EXPECT_EQ(to_integer<int>(x, Rounding::kNearestEven, context), 2);
  EXPECT_EQ(to_integer<int>(x, Rounding::kNearestAway, context), 3);
  EXPECT_EQ(to_integer<int>(x, Rounding::kTowardZero, context), 2);
  EXPECT_EQ(to_integer<int>(x, Rounding::kUp, context), 3);
  EXPECT_EQ(to_integer<int>(x, Rounding::kDown, context), 2);
  EXPECT_EQ(to_integer<int>(-x, Rounding::kDown, context), -3);
  EXPECT_EQ(context.flags, 0u);
  EXPECT_EQ(to_integer_exact<int>(x, Rounding::kDown, context), 2);
  EXPECT_EQ(context.flags, unsigned{kFlagInexact});
  context = {};
  EXPECT_EQ(to_integer_exact<int>(binary128(-4), Rounding::kDown, context), -4);
  EXPECT_EQ(context.flags, 0u);

  // Out of range: invalid, saturated, and never inexact.
  EXPECT_EQ(to_integer_exact<std::int8_t>(binary128(127.5),
                                          Rounding::kNearestEven, context),
            127);
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = {};
  EXPECT_EQ(to_integer<std::int8_t>(binary128(-128.5), Rounding::kTowardZero,
                                    context),
            -128);
  EXPECT_EQ(context.flags, 0u);
  EXPECT_EQ(
      to_integer<std::int8_t>(binary128(-128.5), Rounding::kDown, context),
      -128);
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = {};
  EXPECT_EQ(
      to_integer<unsigned>(binary128(-0.5), Rounding::kTowardZero, context),
      0u);
  EXPECT_EQ(context.flags, 0u);
  EXPECT_EQ(to_integer<unsigned>(binary128(-1), Rounding::kTowardZero, context),
            0u);
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = {};
  EXPECT_EQ(to_integer<std::uint64_t>(binary128(std::ldexp(1.0, 64)),
                                      Rounding::kTowardZero, context),
            std::numeric_limits<std::uint64_t>::max());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = {};
  EXPECT_EQ(
      to_integer<int>(binary128::QuietNaN(), Rounding::kTowardZero, context),
      0);
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  EXPECT_EQ(static_cast<int>(binary128::Infinity(true)),
            std::numeric_limits<int>::min());
  EXPECT_EQ(static_cast<long long>(binary512(-1e300)),
            std::numeric_limits<long long>::min());
  EXPECT_EQ(static_cast<int>(binary128(-2.75)), -2);
  EXPECT_TRUE(static_cast<bool>(binary128(0.25)));
  EXPECT_FALSE(static_cast<bool>(binary128::Zero(true)));
}

TEST(Convert, LongDouble) {
  std::mt19937_64 random(5);
  for (int i = 0; i < 100000; ++i) {
    const double d = std::bit_cast<double>(random());
    if (std::isnan(d)) continue;
    // A long double holds every double exactly on every platform.
    long double x = d;
    if constexpr (std::numeric_limits<long double>::digits > 53) {
      if (std::isfinite(d)) {
        x += static_cast<long double>(d) *
             std::ldexp(static_cast<long double>(1), -60);
      }
    }
    const binary128 wide = x;
    ASSERT_EQ(static_cast<long double>(wide), x);
    if constexpr (std::numeric_limits<long double>::digits <= 113) {
      ASSERT_EQ(binary128(static_cast<double>(x)),
                binary128(static_cast<double>(wide)));
    }
  }
  const long double denorm = std::numeric_limits<long double>::denorm_min();
  EXPECT_EQ(static_cast<long double>(binary128(denorm)), denorm);
  EXPECT_TRUE(std::isinf(static_cast<long double>(binary128::Infinity())));
}

}  // namespace
}  // namespace mpfix
