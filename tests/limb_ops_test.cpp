// Copyright 2026 The hwyMP Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

#include "hwymp/limb_ops.h"

#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <vector>

namespace hwymp {
namespace {

constexpr limb_t kMax = ~limb_t{0};

// The primitives must also work in constant expressions.
constexpr limb_t ConstexprAdc() {
  carry_t carry = 1;
  return adc(kMax, 0, carry) + carry;
}
static_assert(ConstexprAdc() == 1);

constexpr limb_t ConstexprMulHigh() {
  limb_t hi = 0;
  (void)mul_wide(kMax, kMax, hi);
  return hi;
}
static_assert(ConstexprMulHigh() == kMax - 1);

constexpr limb_t ConstexprDiv() {
  limb_t remainder = 0;
  return div_2by1(1, 0, 3, remainder) + remainder;
}
static_assert(ConstexprDiv() == 0x5555555555555555u + 1);

class Xoshiro256 {
 public:
  explicit Xoshiro256(limb_t seed) {
    for (limb_t& word : state_) {
      seed += 0x9E3779B97F4A7C15u;
      limb_t z = seed;
      z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9u;
      z = (z ^ (z >> 27)) * 0x94D049BB133111EBu;
      word = z ^ (z >> 31);
    }
  }

  limb_t Next() {
    const limb_t result = Rotl(state_[1] * 5, 7) * 9;
    const limb_t t = state_[1] << 17;
    state_[2] ^= state_[0];
    state_[3] ^= state_[1];
    state_[1] ^= state_[2];
    state_[0] ^= state_[3];
    state_[2] ^= t;
    state_[3] = Rotl(state_[3], 45);
    return result;
  }

  // Biased towards long runs of ones and zeros, where carries are interesting.
  limb_t NextBiased() {
    const limb_t value = Next();
    switch (Next() & 3) {
      case 0:
        return value;
      case 1:
        return value >> (Next() & 63);
      case 2:
        return ~(value >> (Next() & 63));
      default:
        return value << (Next() & 63);
    }
  }

 private:
  static limb_t Rotl(limb_t x, int k) { return (x << k) | (x >> (64 - k)); }
  std::array<limb_t, 4> state_;
};

std::vector<limb_t> EdgeValues() {
  std::vector<limb_t> values = {0, 1, 2, 3, kMax, kMax - 1, kMax - 2};
  for (int bit = 0; bit < 64; ++bit) {
    const limb_t power = limb_t{1} << bit;
    values.push_back(power);
    values.push_back(power - 1);
    values.push_back(power + 1);
    values.push_back(~power);
  }
  return values;
}

constexpr int kRandomIterations = 2000000;

void CheckAddSub(limb_t a, limb_t b) {
  for (carry_t incoming = 0; incoming <= 1; ++incoming) {
    carry_t carry = incoming;
    carry_t expected_carry = incoming;
    const limb_t sum = adc(a, b, carry);
    ASSERT_EQ(sum, detail::portable::adc(a, b, expected_carry));
    ASSERT_EQ(carry, expected_carry);

    carry_t borrow = incoming;
    carry_t expected_borrow = incoming;
    const limb_t difference = sbb(a, b, borrow);
    ASSERT_EQ(difference, detail::portable::sbb(a, b, expected_borrow));
    ASSERT_EQ(borrow, expected_borrow);

    // Subtracting b again must undo the addition.
    carry_t undo = incoming;
    ASSERT_EQ(sbb(sum, b, undo), a);
    ASSERT_EQ(undo, carry);
  }
}

void CheckMul(limb_t a, limb_t b) {
  limb_t hi = 0;
  limb_t expected_hi = 0;
  const limb_t lo = mul_wide(a, b, hi);
  ASSERT_EQ(lo, detail::portable::mul_wide(a, b, expected_hi));
  ASSERT_EQ(hi, expected_hi);
#if defined(__SIZEOF_INT128__)
  // __int128 is allowed in tests only; it gives an independent reference.
  __extension__ typedef unsigned __int128 U128;
  const U128 product = static_cast<U128>(a) * b;
  ASSERT_EQ(lo, static_cast<limb_t>(product));
  ASSERT_EQ(hi, static_cast<limb_t>(product >> 64));
#endif
}

// Checks quotient * d + remainder == hi:lo using only multiply and add.
void CheckDivResult(limb_t hi, limb_t lo, limb_t d, limb_t quotient,
                    limb_t remainder) {
  ASSERT_LT(remainder, d);
  limb_t product_hi = 0;
  const limb_t product_lo = mul_wide(quotient, d, product_hi);
  carry_t carry = 0;
  const limb_t total_lo = adc(product_lo, remainder, carry);
  const limb_t total_hi = adc(product_hi, 0, carry);
  ASSERT_EQ(carry, 0);
  ASSERT_EQ(total_lo, lo);
  ASSERT_EQ(total_hi, hi);
}

void CheckDiv(limb_t hi, limb_t lo, limb_t d) {
  if (d == 0) return;
  hi %= d;
  limb_t remainder = 0;
  const limb_t quotient = div_2by1(hi, lo, d, remainder);
  CheckDivResult(hi, lo, d, quotient, remainder);

  limb_t portable_remainder = 0;
  const limb_t portable_quotient =
      detail::portable::div_2by1(hi, lo, d, portable_remainder);
  ASSERT_EQ(portable_quotient, quotient);
  ASSERT_EQ(portable_remainder, remainder);
}

TEST(LimbOps, AddSubEdgeValues) {
  const std::vector<limb_t> values = EdgeValues();
  for (const limb_t a : values) {
    for (const limb_t b : values) {
      CheckAddSub(a, b);
      if (HasFatalFailure()) return;
    }
  }
}

TEST(LimbOps, AddSubRandom) {
  Xoshiro256 rng(1);
  for (int i = 0; i < kRandomIterations; ++i) {
    CheckAddSub(rng.NextBiased(), rng.NextBiased());
    if (HasFatalFailure()) return;
  }
}

TEST(LimbOps, CarryChainPropagatesAcrossLimbs) {
  std::array<limb_t, 4> a = {kMax, kMax, kMax, kMax};
  const std::array<limb_t, 4> b = {1, 0, 0, 0};
  carry_t carry = 0;
  for (std::size_t i = 0; i < a.size(); ++i) a[i] = adc(a[i], b[i], carry);
  EXPECT_EQ(a, (std::array<limb_t, 4>{0, 0, 0, 0}));
  EXPECT_EQ(carry, 1);

  carry_t borrow = 0;
  for (std::size_t i = 0; i < a.size(); ++i) a[i] = sbb(a[i], b[i], borrow);
  EXPECT_EQ(a, (std::array<limb_t, 4>{kMax, kMax, kMax, kMax}));
  EXPECT_EQ(borrow, 1);
}

TEST(LimbOps, MulEdgeValues) {
  const std::vector<limb_t> values = EdgeValues();
  for (const limb_t a : values) {
    for (const limb_t b : values) {
      CheckMul(a, b);
      if (HasFatalFailure()) return;
    }
  }
}

TEST(LimbOps, MulRandom) {
  Xoshiro256 rng(2);
  for (int i = 0; i < kRandomIterations; ++i) {
    CheckMul(rng.NextBiased(), rng.NextBiased());
    if (HasFatalFailure()) return;
  }
}

TEST(LimbOps, DivEdgeValues) {
  const std::vector<limb_t> values = EdgeValues();
  for (const limb_t hi : values) {
    for (const limb_t lo : values) {
      for (const limb_t d : {limb_t{1}, limb_t{2}, limb_t{3},
                             limb_t{0xFFFFFFFFu}, limb_t{0x100000000u},
                             limb_t{0x100000001u}, limb_t{1} << 63, kMax}) {
        CheckDiv(hi, lo, d);
        if (HasFatalFailure()) return;
      }
    }
  }
}

TEST(LimbOps, DivRandom) {
  Xoshiro256 rng(3);
  for (int i = 0; i < kRandomIterations; ++i) {
    CheckDiv(rng.NextBiased(), rng.NextBiased(), rng.NextBiased());
    if (HasFatalFailure()) return;
  }
}

TEST(LimbOps, CountZeros) {
  EXPECT_EQ(countl_zero(0), 64);
  EXPECT_EQ(countr_zero(0), 64);
  for (int bit = 0; bit < 64; ++bit) {
    const limb_t power = limb_t{1} << bit;
    EXPECT_EQ(countl_zero(power), 63 - bit);
    EXPECT_EQ(countr_zero(power), bit);
    EXPECT_EQ(countl_zero(power | 1), 63 - bit);
    EXPECT_EQ(countr_zero(power | (limb_t{1} << 63)), bit);
  }
}

TEST(LimbOps, DoubleShifts) {
  const limb_t hi = 0x0123456789ABCDEFu;
  const limb_t lo = 0xFEDCBA9876543210u;
  for (int count = 1; count < 64; ++count) {
    EXPECT_EQ(shift_left_2(hi, lo, count),
              (hi << count) | (lo >> (64 - count)));
    EXPECT_EQ(shift_right_2(hi, lo, count),
              (lo >> count) | (hi << (64 - count)));
  }
  EXPECT_EQ(shift_left_2(hi, lo, 4), 0x123456789ABCDEFFu);
  EXPECT_EQ(shift_right_2(hi, lo, 4), 0xFFEDCBA987654321u);
}

}  // namespace
}  // namespace hwymp
