// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Known-answer tests that need no oracle: encodings, special values, the NaN
// rules (which MPFR does not model) and comparisons.

#include "mpfix/float.h"

#include <gtest/gtest.h>

#include <cmath>
#include <compare>
#include <cstddef>
#include <cstring>
#include <limits>
#include <tuple>
#include <type_traits>
#include <vector>

namespace mpfix {
namespace {

using Tiny = Float<8, 4>;  // 12 bits: 1 sign, 4 exponent, 7 trailing.

static_assert(sizeof(binary128) == 24);
static_assert(sizeof(binary256) == 40);
static_assert(sizeof(binary512) == 72);
static_assert(alignof(binary256) == 8);
static_assert(std::is_trivially_copyable_v<binary256>);
static_assert(std::is_trivially_default_constructible_v<binary256>);
static_assert(binary128::kBitsLimbs == 2);
static_assert(binary256::kBitsLimbs == 4);
static_assert(binary512::kBitsLimbs == 8);
static_assert(binary128::kEmax == 16383);
static_assert(binary256::kEmax == 262143);
static_assert(binary512::kEmax == 4194303);

binary128 B128(limb_t high, limb_t low) {
  return binary128::from_bits({low, high});
}

constexpr limb_t kOneHigh = 0x3FFF000000000000u;
constexpr limb_t kTwoHigh = 0x4000000000000000u;
constexpr limb_t kThreeHigh = 0x4000800000000000u;
constexpr limb_t kSixHigh = 0x4001800000000000u;
constexpr limb_t kInfinityHigh = 0x7FFF000000000000u;
constexpr limb_t kQuietHigh = 0x7FFF800000000000u;
constexpr limb_t kSignHigh = 0x8000000000000000u;

void ExpectBits(const binary128& x, limb_t high, limb_t low) {
  const binary128::Bits bits = x.to_bits();
  EXPECT_EQ(bits[1], high);
  EXPECT_EQ(bits[0], low);
}

TEST(Float, Binary128KnownArithmetic) {
  const binary128 one = B128(kOneHigh, 0);
  const binary128 two = B128(kTwoHigh, 0);
  const binary128 three = B128(kThreeHigh, 0);
  ExpectBits(one + one, kTwoHigh, 0);
  ExpectBits(one + two, kThreeHigh, 0);
  ExpectBits(three - two, kOneHigh, 0);
  ExpectBits(two - three, kOneHigh | kSignHigh, 0);
  ExpectBits(two * three, kSixHigh, 0);
  ExpectBits(-two * three, kSixHigh | kSignHigh, 0);
}

// Libraries such as Eigen clear arrays of trivially copyable scalars with
// memset, so zero bytes must be +0.
TEST(Float, ZeroBytesArePositiveZero) {
  const auto check = [](auto x) {
    std::memset(static_cast<void*>(&x), 0, sizeof(x));
    EXPECT_TRUE(x.is_zero());
    EXPECT_FALSE(x.signbit());
    EXPECT_EQ(x.to_bits(), decltype(x)::Zero().to_bits());
    EXPECT_EQ(x + decltype(x)::MinNormal(), decltype(x)::MinNormal());
  };
  check(binary128());
  check(binary256());
  check(binary512());
  check(Float<4, 3>());
}

TEST(Float, CompoundAssignment) {
  binary128 x = B128(kOneHigh, 0);
  x += B128(kTwoHigh, 0);
  ExpectBits(x, kThreeHigh, 0);
  x *= B128(kTwoHigh, 0);
  ExpectBits(x, kSixHigh, 0);
  x -= B128(kThreeHigh, 0);
  ExpectBits(x, kThreeHigh, 0);
}

template <class F>
void CheckBitsRoundTrip() {
  for (const F& x : {F::Zero(), F::Zero(true), F::Infinity(), F::Infinity(true),
                     F::QuietNaN(), F::SignalingNaN(), F::MaxFinite(),
                     F::MaxFinite(true), F::MinNormal(), F::MinSubnormal(),
                     F::MinSubnormal(true)}) {
    EXPECT_EQ(F::from_bits(x.to_bits()).to_bits(), x.to_bits());
  }
}

TEST(Float, BitsRoundTrip) {
  CheckBitsRoundTrip<Tiny>();
  CheckBitsRoundTrip<Float<64, 8>>();
  CheckBitsRoundTrip<Float<65, 8>>();
  CheckBitsRoundTrip<binary128>();
  CheckBitsRoundTrip<binary256>();
  CheckBitsRoundTrip<binary512>();
}

TEST(Float, TinyEncodingsRoundTripExhaustively) {
  for (limb_t bits = 0; bits < (limb_t{1} << 12); ++bits) {
    EXPECT_EQ(Tiny::from_bits({bits}).to_bits()[0], bits);
  }
}

TEST(Float, Encodings) {
  ExpectBits(binary128::Infinity(true), kInfinityHigh | kSignHigh, 0);
  ExpectBits(binary128::QuietNaN(), kQuietHigh, 0);
  ExpectBits(binary128::SignalingNaN(), kInfinityHigh, 1);
  ExpectBits(binary128::MinSubnormal(), 0, 1);
  ExpectBits(binary128::MinNormal(), 0x0001000000000000u, 0);
  ExpectBits(binary128::MaxFinite(), 0x7FFEFFFFFFFFFFFFu, ~limb_t{0});
}

TEST(Float, Classification) {
  EXPECT_EQ(binary256::Zero().classify(), FloatClass::kZero);
  EXPECT_EQ(binary256::Infinity().classify(), FloatClass::kInfinite);
  EXPECT_EQ(binary256::QuietNaN().classify(), FloatClass::kNaN);
  EXPECT_EQ(binary256::MinSubnormal().classify(), FloatClass::kSubnormal);
  EXPECT_EQ(binary256::MinNormal().classify(), FloatClass::kNormal);
  EXPECT_TRUE(binary256::Zero().is_finite());
  EXPECT_FALSE(binary256::Zero().is_normal());
  EXPECT_FALSE(binary256::Infinity().is_finite());
  EXPECT_TRUE(binary256::SignalingNaN().is_signaling());
  EXPECT_FALSE(binary256::QuietNaN().is_signaling());
  EXPECT_TRUE(binary256::Zero(true).signbit());
}

TEST(Float, SignOperationsDoNotQuietNaN) {
  const binary128 s = binary128::SignalingNaN();
  EXPECT_TRUE((-s).is_signaling());
  EXPECT_TRUE((-s).signbit());
  EXPECT_FALSE(abs(-s).signbit());
  EXPECT_TRUE(copysign(s, binary128::Zero(true)).signbit());
}

TEST(Float, SignedZeroRules) {
  const binary128 zero = binary128::Zero();
  const binary128 negative_zero = binary128::Zero(true);
  const binary128 one = B128(kOneHigh, 0);

  EXPECT_FALSE((zero + negative_zero).signbit());
  EXPECT_TRUE((negative_zero + negative_zero).signbit());
  EXPECT_FALSE((one - one).signbit());
  EXPECT_TRUE((negative_zero * one).signbit());

  Context down{Rounding::kDown};
  EXPECT_TRUE(sub(one, one, down).signbit());
  EXPECT_TRUE(add(zero, negative_zero, down).signbit());
  EXPECT_EQ(down.flags, 0u);
}

TEST(Float, InvalidOperationsReturnDefaultNaN) {
  const binary128 infinity = binary128::Infinity();
  Context context;
  ExpectBits(sub(infinity, infinity, context), kQuietHigh, 0);
  EXPECT_EQ(context.flags, kFlagInvalid);

  context = {};
  ExpectBits(add(infinity, -infinity, context), kQuietHigh, 0);
  EXPECT_EQ(context.flags, kFlagInvalid);

  context = {};
  ExpectBits(mul(infinity, binary128::Zero(true), context), kQuietHigh, 0);
  EXPECT_EQ(context.flags, kFlagInvalid);
}

TEST(Float, NaNPropagationKeepsFirstOperandQuieted) {
  const binary128 quiet_a = B128(kQuietHigh | 5, 7);
  const binary128 quiet_b = B128(kQuietHigh | kSignHigh | 9, 11);
  const binary128 signaling = B128(kInfinityHigh | 3, 1);
  const binary128 one = B128(kOneHigh, 0);

  Context context;
  ExpectBits(add(quiet_a, quiet_b, context), kQuietHigh | 5, 7);
  ExpectBits(mul(quiet_b, quiet_a, context), kQuietHigh | kSignHigh | 9, 11);
  ExpectBits(sub(one, quiet_a, context), kQuietHigh | 5, 7);
  EXPECT_EQ(context.flags, 0u);

  ExpectBits(add(signaling, one, context), kQuietHigh | 3, 1);
  EXPECT_EQ(context.flags, kFlagInvalid);

  context = {};
  ExpectBits(mul(one, signaling, context), kQuietHigh | 3, 1);
  EXPECT_EQ(context.flags, kFlagInvalid);

  context = {};
  ExpectBits(add(quiet_a, signaling, context), kQuietHigh | 5, 7);
  EXPECT_EQ(context.flags, kFlagInvalid);
}

TEST(Float, Comparisons) {
  const binary128 one = B128(kOneHigh, 0);
  const binary128 two = B128(kTwoHigh, 0);
  const binary128 nan = binary128::QuietNaN();

  EXPECT_LT(one, two);
  EXPECT_LE(one, one);
  EXPECT_GT(-one, -two);
  EXPECT_LT(-two, one);
  EXPECT_LT(binary128::MaxFinite(), binary128::Infinity());
  EXPECT_LT(binary128::Zero(), binary128::MinSubnormal());
  EXPECT_LT(binary128::MinSubnormal(), binary128::MinNormal());
  EXPECT_GT(binary128::Zero(), -binary128::MinSubnormal());
  EXPECT_EQ(binary128::Zero(), binary128::Zero(true));
  EXPECT_LT(-binary128::Infinity(), binary128::Zero());

  EXPECT_FALSE(nan == nan);
  EXPECT_TRUE(nan != nan);
  EXPECT_FALSE(nan < one);
  EXPECT_FALSE(nan > one);
  EXPECT_EQ(nan <=> one, std::partial_ordering::unordered);
}

TEST(Float, ComparisonFlags) {
  const binary128 one = B128(kOneHigh, 0);
  Context context;
  (void)compare_quiet(one, binary128::QuietNaN(), context);
  EXPECT_EQ(context.flags, 0u);
  (void)compare_quiet(one, binary128::SignalingNaN(), context);
  EXPECT_EQ(context.flags, kFlagInvalid);

  context = {};
  (void)compare_signaling(one, one, context);
  EXPECT_EQ(context.flags, 0u);
  (void)compare_signaling(binary128::QuietNaN(), one, context);
  EXPECT_EQ(context.flags, kFlagInvalid);
}

TEST(Float, OverflowDependsOnRoundingMode) {
  const binary128 max = binary128::MaxFinite();
  const binary128 two = B128(kTwoHigh, 0);
  for (const auto& [mode, positive_infinity, negative_infinity] :
       {std::tuple{Rounding::kNearestEven, true, true},
        std::tuple{Rounding::kNearestAway, true, true},
        std::tuple{Rounding::kTowardZero, false, false},
        std::tuple{Rounding::kUp, true, false},
        std::tuple{Rounding::kDown, false, true}}) {
    Context positive{mode};
    EXPECT_EQ(mul(max, two, positive).is_inf(), positive_infinity);
    EXPECT_EQ(positive.flags, kFlagOverflow | kFlagInexact);
    Context negative{mode};
    EXPECT_EQ(mul(-max, two, negative).is_inf(), negative_infinity);
    EXPECT_EQ(negative.flags, kFlagOverflow | kFlagInexact);
  }
}

TEST(Float, SubnormalArithmeticIsExact) {
  const binary128 min = binary128::MinSubnormal();
  Context context;
  ExpectBits(add(min, min, context), 0, 2);
  ExpectBits(sub(add(min, min, context), min, context), 0, 1);
  EXPECT_EQ(context.flags, 0u);
}

TEST(Float, UnderflowToZeroKeepsSign) {
  const binary128 min = binary128::MinSubnormal();
  const binary128 half = B128(0x3FFE000000000000u, 0);
  Context context;
  const binary128 result = mul(-min, half, context);
  EXPECT_TRUE(result.is_zero());
  EXPECT_TRUE(result.signbit());
  EXPECT_EQ(context.flags, kFlagUnderflow | kFlagInexact);

  Context up{Rounding::kUp};
  ExpectBits(mul(min, half, up), 0, 1);
  EXPECT_EQ(up.flags, kFlagUnderflow | kFlagInexact);
}

TEST(Float, ExactSubnormalResultRaisesNoUnderflow) {
  const binary128 two_min = binary128::from_bits({2, 0});
  const binary128 half = B128(0x3FFE000000000000u, 0);
  Context context;
  ExpectBits(mul(two_min, half, context), 0, 1);
  EXPECT_EQ(context.flags, 0u);
}

using Small = Float<4, 3>;

// Every non-NaN value of Small in increasing order, with -0 before +0.
std::vector<Small> SortedValues() {
  std::vector<Small> positive;
  const limb_t count = limb_t{1} << 6;  // Every encoding with the sign clear.
  for (limb_t bits = 0; bits < count; ++bits) {
    const Small x = Small::from_bits({bits});
    if (!x.is_nan()) positive.push_back(x);
  }
  std::vector<Small> sorted;
  for (std::size_t i = positive.size(); i-- > 0;) {
    sorted.push_back(-positive[i]);
  }
  sorted.insert(sorted.end(), positive.begin(), positive.end());
  return sorted;
}

template <int P, int E>
bool SameBits(const Float<P, E>& a, const Float<P, E>& b) {
  return a.to_bits() == b.to_bits();
}

TEST(FloatNeighbours, StepsThroughEveryValueInOrder) {
  const std::vector<Small> sorted = SortedValues();
  ASSERT_GT(sorted.size(), 100u);
  for (std::size_t i = 0; i + 1 < sorted.size(); ++i) {
    const Small& x = sorted[i];
    const Small& next = sorted[i + 1];
    if (x.is_zero() && next.is_zero()) continue;
    if (x.is_zero()) {
      EXPECT_TRUE(SameBits(next_up(x), Small::MinSubnormal()));
      continue;
    }
    if (next.is_zero()) {
      EXPECT_TRUE(SameBits(next_up(x), next)) << i;
      continue;
    }
    EXPECT_TRUE(SameBits(next_up(x), next)) << i;
    EXPECT_TRUE(SameBits(next_down(next), x)) << i;
  }
  EXPECT_TRUE(SameBits(next_up(Small::Infinity()), Small::Infinity()));
  EXPECT_TRUE(SameBits(next_down(Small::Infinity(true)), Small::Infinity(true)));
  EXPECT_TRUE(SameBits(next_down(Small::Zero(true)), Small::MinSubnormal(true)));
  EXPECT_TRUE(SameBits(next_down(Small::Zero()), Small::MinSubnormal(true)));
  EXPECT_TRUE(SameBits(next_up(Small::MaxFinite()), Small::Infinity()));
  EXPECT_TRUE(SameBits(next_up(Small::MinSubnormal(true)), Small::Zero(true)));
}

TEST(FloatNeighbours, RaiseInvalidOnlyForSignalingNaN) {
  Context context;
  EXPECT_TRUE(next_up(Small::QuietNaN(), context).is_nan());
  EXPECT_EQ(context.flags, 0u);
  const Small result = next_down(Small::SignalingNaN(), context);
  EXPECT_TRUE(result.is_nan());
  EXPECT_FALSE(result.is_signaling());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  next_up(Small::MaxFinite(), context = Context{});
  EXPECT_EQ(context.flags, 0u);
}

TEST(FloatMinMax, FollowsTheOrderOfValues) {
  const std::vector<Small> sorted = SortedValues();
  for (std::size_t i = 0; i < sorted.size(); ++i) {
    for (std::size_t j = 0; j < sorted.size(); ++j) {
      const Small& low = sorted[i < j ? i : j];
      const Small& high = sorted[i < j ? j : i];
      EXPECT_TRUE(SameBits(minimum(sorted[i], sorted[j]), low));
      EXPECT_TRUE(SameBits(maximum(sorted[i], sorted[j]), high));
      EXPECT_TRUE(SameBits(minimum_number(sorted[i], sorted[j]), low));
      EXPECT_TRUE(SameBits(maximum_number(sorted[i], sorted[j]), high));
    }
  }
}

TEST(FloatMinMax, HandlesNaN) {
  const binary128 one = binary128::MinNormal();
  const binary128 quiet = binary128::QuietNaN();
  const binary128 signaling = binary128::SignalingNaN();

  Context context;
  EXPECT_TRUE(minimum(one, quiet, context).is_nan());
  EXPECT_TRUE(maximum(quiet, one, context).is_nan());
  EXPECT_EQ(context.flags, 0u);
  EXPECT_FALSE(minimum(one, signaling, context).is_signaling());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});

  context = Context{};
  EXPECT_TRUE(minimum_number(one, quiet, context) == one);
  EXPECT_TRUE(maximum_number(quiet, one, context) == one);
  EXPECT_EQ(context.flags, 0u);
  EXPECT_TRUE(minimum_number(signaling, one, context) == one);
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});

  context = Context{};
  EXPECT_TRUE(maximum_number(quiet, quiet, context).is_nan());
  EXPECT_EQ(context.flags, 0u);
  const binary128 both = maximum_number(signaling, quiet, context);
  EXPECT_TRUE(both.is_nan());
  EXPECT_FALSE(both.is_signaling());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
}

TEST(FloatExponent, Ilogb) {
  EXPECT_EQ(ilogb(binary128::MinNormal()), binary128::kEmin);
  EXPECT_EQ(ilogb(binary128::MaxFinite()), binary128::kEmax);
  EXPECT_EQ(ilogb(binary128::MinSubnormal()), binary128::kEmin - 112);

  Context context;
  EXPECT_EQ(ilogb(binary128::Zero(), context), binary128::kIlogbZero);
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = Context{};
  EXPECT_EQ(ilogb(binary128::Infinity(true), context),
            binary128::kIlogbInfinity);
  EXPECT_EQ(ilogb(binary128::QuietNaN(), context), binary128::kIlogbNaN);
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = Context{};
  EXPECT_EQ(ilogb(binary128::MaxFinite(), context), binary128::kEmax);
  EXPECT_EQ(context.flags, 0u);
}

TEST(FloatExponent, Logb) {
  const binary128 min_normal = binary128::MinNormal();
  const binary128 one = scalbn(min_normal, -binary128::kEmin);
  EXPECT_TRUE(logb(one).is_zero());
  EXPECT_FALSE(logb(one).signbit());
  EXPECT_TRUE(logb(scalbn(one, 1)) == one);
  EXPECT_TRUE(logb(scalbn(one, -1)) == -one);
  const binary128 two = scalbn(one, 1);
  const binary128 emin = -(scalbn(one, 14) - two);  // kEmin is -16382.
  EXPECT_TRUE(logb(min_normal) == emin);

  Context context{Rounding::kUp};
  const binary128 at_zero = logb(binary128::Zero(true), context);
  EXPECT_TRUE(at_zero.is_inf());
  EXPECT_TRUE(at_zero.signbit());
  EXPECT_EQ(context.flags, unsigned{kFlagDivByZero});

  context = Context{};
  EXPECT_TRUE(logb(binary128::Infinity(true), context) ==
              binary128::Infinity());
  EXPECT_TRUE(logb(binary128::QuietNaN(), context).is_nan());
  EXPECT_EQ(context.flags, 0u);
  EXPECT_FALSE(logb(binary128::SignalingNaN(), context).is_signaling());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
}

TEST(FloatScalbn, ExtremeCountsSaturate) {
  const int huge = std::numeric_limits<int>::max();
  const binary512 one = scalbn(binary512::MinNormal(), -binary512::kEmin);
  EXPECT_TRUE(scalbn(one, huge).is_inf());
  EXPECT_TRUE(scalbn(one, -huge).is_zero());
  EXPECT_TRUE(scalbn(binary512::MinSubnormal(), huge).is_inf());
  EXPECT_TRUE(scalbn(binary512::MaxFinite(), std::numeric_limits<int>::min())
                  .is_zero());
  EXPECT_TRUE(scalbn(one, 0) == one);
}

TEST(FloatScalbn, LdexpAndScalblnMatchScalbn) {
  const binary256 one = scalbn(binary256::MinNormal(), -binary256::kEmin);
  for (const int n : {-300000, -262378, -1000, -1, 0, 1, 1000, 262143, 300000}) {
    EXPECT_TRUE(SameBits(ldexp(one, n), scalbn(one, n)));
    EXPECT_TRUE(SameBits(scalbln(one, n), scalbn(one, n)));
  }
  // The count is clamped before it is added, so no long overflows.
  const long huge = std::numeric_limits<long>::max();
  const long tiny = std::numeric_limits<long>::min();
  EXPECT_TRUE(scalbln(binary256::MaxFinite(), huge).is_inf());
  EXPECT_TRUE(scalbln(binary256::MinSubnormal(), tiny).is_zero());
  EXPECT_TRUE(scalbln(binary256::MinSubnormal(), huge).is_inf());
  EXPECT_TRUE(scalbln(binary256::MaxFinite(), tiny).is_zero());

  Context context{Rounding::kTowardZero};
  EXPECT_TRUE(SameBits(ldexp(one, 300000, context), binary256::MaxFinite()));
  EXPECT_EQ(context.flags, unsigned{kFlagOverflow | kFlagInexact});
  context = Context{};
  EXPECT_TRUE(scalbln(binary256::MinNormal(), -1, context).is_subnormal());
  EXPECT_EQ(context.flags, 0u);
}

TEST(FloatFrexp, MatchesTheCLibraryOnEveryValue) {
  for (const Small& x : SortedValues()) {
    int exponent = 12345;
    const Small significand = frexp(x, &exponent);
    int expected_exponent = 12345;
    const double expected =
        std::frexp(static_cast<double>(x), &expected_exponent);
    EXPECT_EQ(static_cast<double>(significand), expected);
    EXPECT_EQ(significand.signbit(), x.signbit());
    if (x.is_inf()) {
      EXPECT_EQ(exponent, 0);  // The C value is unspecified here.
    } else {
      EXPECT_EQ(exponent, expected_exponent);
    }
    if (!x.is_zero() && !x.is_inf()) {
      EXPECT_TRUE(SameBits(ldexp(significand, exponent), x));
    }
  }
}

TEST(FloatFrexp, ExtremesAndNaN) {
  int exponent = 0;
  const binary512 half = frexp(binary512::MaxFinite(), &exponent);
  EXPECT_EQ(exponent, binary512::kEmax + 1);
  EXPECT_TRUE(SameBits(ldexp(half, exponent), binary512::MaxFinite()));
  frexp(binary512::MinSubnormal(true), &exponent);
  EXPECT_EQ(exponent, binary512::kEmin - 488 + 1);
  EXPECT_TRUE(frexp(binary512::MinSubnormal(true), &exponent) ==
              -scalbn(binary512::MinNormal(), -binary512::kEmin - 1));

  Context context;
  exponent = 7;
  EXPECT_TRUE(frexp(binary512::QuietNaN(), &exponent, context).is_nan());
  EXPECT_EQ(exponent, 0);
  EXPECT_EQ(context.flags, 0u);
  EXPECT_FALSE(frexp(binary512::SignalingNaN(), &exponent, context)
                   .is_signaling());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
}

TEST(FloatModf, MatchesTheCLibraryOnEveryValue) {
  for (const Small& x : SortedValues()) {
    Small integral = Small::QuietNaN();
    const Small fraction = modf(x, &integral);
    double expected_integral = 0;
    const double expected =
        std::modf(static_cast<double>(x), &expected_integral);
    EXPECT_EQ(static_cast<double>(fraction), expected);
    EXPECT_EQ(fraction.signbit(), std::signbit(expected));
    EXPECT_EQ(static_cast<double>(integral), expected_integral);
    EXPECT_EQ(integral.signbit(), std::signbit(expected_integral));
  }
}

TEST(FloatModf, LargeValuesAndNaN) {
  const binary128 one = scalbn(binary128::MinNormal(), -binary128::kEmin);
  const binary128 big = scalbn(one, 112) + one;  // An odd integer.
  binary128 integral;
  const binary128 fraction = modf(-big, &integral);
  EXPECT_TRUE(fraction.is_zero());
  EXPECT_TRUE(fraction.signbit());
  EXPECT_TRUE(integral == -big);
  const binary128 x = scalbn(one, 111) + scalbn(one, -1);
  EXPECT_TRUE(modf(x, &integral) == scalbn(one, -1));
  EXPECT_TRUE(integral == scalbn(one, 111));

  Context context;
  EXPECT_FALSE(modf(binary128::SignalingNaN(), &integral, context)
                   .is_signaling());
  EXPECT_TRUE(integral.is_nan());
  EXPECT_FALSE(integral.is_signaling());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = Context{};
  modf(x, &integral, context);
  EXPECT_EQ(context.flags, 0u);
}

TEST(FloatNextAfter, StepsTowardEveryValue) {
  const std::vector<Small> sorted = SortedValues();
  for (std::size_t i = 0; i < sorted.size(); ++i) {
    for (std::size_t j = 0; j < sorted.size(); ++j) {
      const Small& x = sorted[i];
      const Small& y = sorted[j];
      Small expected = y;
      if (!(x == y)) expected = i < j ? next_up(x) : next_down(x);
      Context context;
      const Small result = nextafter(x, y, context);
      EXPECT_TRUE(SameBits(result, expected)) << i << ' ' << j;
      unsigned flags = 0;
      if (!(x == y) && result.is_inf()) flags = kFlagOverflow | kFlagInexact;
      if (!(x == y) && (result.is_zero() || result.is_subnormal())) {
        flags = kFlagUnderflow | kFlagInexact;
      }
      EXPECT_EQ(context.flags, flags) << i << ' ' << j;
      // Every Small is exactly a long double.
      EXPECT_TRUE(SameBits(
          nexttoward(x, static_cast<long double>(static_cast<double>(y))),
          expected));
    }
  }
}

TEST(FloatNextAfter, TowardLongDoublesBetweenValues) {
  const binary128 one = scalbn(binary128::MinNormal(), -binary128::kEmin);
  const long double above = 1.0L + std::ldexp(1.0L, -40);
  EXPECT_TRUE(SameBits(nexttoward(one, above), next_up(one)));
  EXPECT_TRUE(SameBits(nexttoward(-one, -above), next_down(-one)));
  EXPECT_TRUE(SameBits(nexttoward(one, -above), next_down(one)));
  // Small cannot hold `above`, which still compares exactly.
  EXPECT_TRUE(SameBits(nexttoward(Small(1.0), above), next_up(Small(1.0))));
  EXPECT_TRUE(SameBits(nexttoward(Small(1.0), -above), next_down(Small(1.0))));
  EXPECT_TRUE(SameBits(nexttoward(Small(1.0), 1.0L), Small(1.0)));
  EXPECT_TRUE(SameBits(nexttoward(Small::Zero(), -0.0L), Small::Zero(true)));
  EXPECT_TRUE(nexttoward(Small::MaxFinite(),
                         std::numeric_limits<long double>::infinity())
                  .is_inf());
  EXPECT_TRUE(nexttoward(one, std::numeric_limits<long double>::quiet_NaN())
                  .is_nan());

  Context context;
  EXPECT_FALSE(nextafter(binary128::SignalingNaN(), one, context)
                   .is_signaling());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = Context{};
  EXPECT_TRUE(nextafter(one, binary128::QuietNaN(), context).is_nan());
  EXPECT_EQ(context.flags, 0u);
}

TEST(FloatRounding, CStyleFunctions) {
  const binary256 one = scalbn(binary256::MinNormal(), -binary256::kEmin);
  const binary256 half = scalbn(one, -1);
  const binary256 two = one + one;
  const binary256 x = one + half;  // 1.5

  EXPECT_TRUE(floor(x) == one);
  EXPECT_TRUE(ceil(x) == two);
  EXPECT_TRUE(trunc(x) == one);
  EXPECT_TRUE(round(x) == two);
  EXPECT_TRUE(floor(-x) == -two);
  EXPECT_TRUE(ceil(-x) == -one);
  EXPECT_TRUE(trunc(-x) == -one);
  EXPECT_TRUE(round(-x) == -two);

  EXPECT_TRUE(round(half) == one);
  EXPECT_TRUE(trunc(-half).is_zero());
  EXPECT_TRUE(trunc(-half).signbit());
  EXPECT_TRUE(ceil(-half).signbit());
  EXPECT_TRUE(floor(half).is_zero());
  EXPECT_TRUE(floor(binary256::Infinity(true)).is_inf());
}

// The expected limbs come from an exact integer square root. MPFR 4.2.2 returns
// wrong values for this input when rounding down or up, and a wrong ternary value
// when rounding to nearest, which is how it was found.
TEST(Sqrt, Mpfr422Regression) {
  const binary512::Limbs x = {0,
                              0,
                              2,
                              ~limb_t{0},
                              ~limb_t{0},
                              1,
                              0,
                              limb_t{1} << 63};
  const binary512::Limbs expected = {0,
                                     0,
                                     limb_t{1} << 63,
                                     ~limb_t{0},
                                     ~limb_t{0},
                                     0,
                                     0,
                                     limb_t{1} << 63};
  binary512::Limbs above = expected;
  above[0] = limb_t{1} << 23;
  const binary512 input = binary512::FromRaw(false, -2, x);

  Context down{Rounding::kTowardZero};
  const binary512 root = sqrt(input, down);
  EXPECT_EQ(root.raw_exponent(), -1);
  EXPECT_EQ(root.raw_limbs(), expected);
  EXPECT_EQ(down.flags, unsigned{kFlagInexact});

  Context up{Rounding::kUp};
  EXPECT_EQ(sqrt(input, up).raw_limbs(), above);
  EXPECT_EQ(up.flags, unsigned{kFlagInexact});

  // The root is just above the lower neighbour, so nearest goes down.
  EXPECT_EQ(sqrt(input).raw_limbs(), expected);
}

}  // namespace
}  // namespace mpfix
