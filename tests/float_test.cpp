// Copyright 2026 The hwyMP Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Known-answer tests that need no oracle: encodings, special values, the NaN
// rules (which MPFR does not model) and comparisons.

#include "hwymp/float.h"

#include <gtest/gtest.h>

#include <compare>
#include <tuple>
#include <type_traits>

namespace hwymp {
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

}  // namespace
}  // namespace hwymp
