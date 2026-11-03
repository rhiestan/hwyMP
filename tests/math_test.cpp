// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Special values, flags and known results of the exponentials and logarithms.
// math_oracle_test compares them against MPFR.

#include "mpfix/math.h"

#include <gtest/gtest.h>

#include <cmath>
#include <string>

namespace mpfix {
namespace {

template <int P, int E>
bool SameBits(const Float<P, E>& a, const Float<P, E>& b) {
  return a.to_bits() == b.to_bits();
}

TEST(Math, SpecialValues) {
  using F = binary128;
  const F inf = F::Infinity();
  const F one(1);
  EXPECT_TRUE(SameBits(exp(F::Zero(true)), one));
  EXPECT_TRUE(SameBits(exp(inf), inf));
  EXPECT_TRUE(SameBits(exp(-inf), F::Zero()));
  EXPECT_TRUE(SameBits(exp2(-inf), F::Zero()));
  EXPECT_TRUE(SameBits(expm1(-inf), -one));
  EXPECT_TRUE(SameBits(expm1(F::Zero(true)), F::Zero(true)));
  EXPECT_TRUE(SameBits(log1p(F::Zero(true)), F::Zero(true)));
  EXPECT_TRUE(SameBits(log(one), F::Zero()));
  EXPECT_TRUE(SameBits(log2(one), F::Zero()));
  EXPECT_TRUE(SameBits(log10(one), F::Zero()));
  EXPECT_TRUE(SameBits(log(inf), inf));
  EXPECT_TRUE(SameBits(log1p(inf), inf));

  Context context;
  EXPECT_TRUE(SameBits(log(F::Zero(true), context), -inf));
  EXPECT_EQ(context.flags, unsigned{kFlagDivByZero});
  context = Context{};
  EXPECT_TRUE(SameBits(log1p(-one, context), -inf));
  EXPECT_EQ(context.flags, unsigned{kFlagDivByZero});
  context = Context{};
  EXPECT_TRUE(log(-one, context).is_nan());
  EXPECT_TRUE(log2(-inf, context).is_nan());
  EXPECT_TRUE(log1p(F(-2), context).is_nan());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = Context{};
  EXPECT_TRUE(exp(F::QuietNaN(), context).is_nan());
  EXPECT_EQ(context.flags, 0u);
  EXPECT_FALSE(log(F::SignalingNaN(), context).is_signaling());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
}

TEST(Math, ExactResultsRaiseNothing) {
  using F = binary256;
  Context context;
  EXPECT_TRUE(exp2(F(-10), context) == scalbn(F(1), -10));
  EXPECT_TRUE(log2(F(1024), context) == F(10));
  EXPECT_TRUE(log2(F(0.125), context) == F(-3));
  EXPECT_TRUE(log10(F(1000000), context) == F(6));
  EXPECT_TRUE(exp(F::Zero(), context) == F(1));
  EXPECT_EQ(context.flags, 0u);
  EXPECT_TRUE(log10(F(0.001), context) < F(-2.99));
  EXPECT_EQ(context.flags, unsigned{kFlagInexact});
}

TEST(Math, MatchesDoubleClosely) {
  using F = binary128;
  for (const double x : {-700.0, -3.5, -1.0, -1e-5, 1e-300, 0.1, 0.5, 1.0, 2.0,
                         10.0, 123.456, 700.0}) {
    EXPECT_NEAR(static_cast<double>(exp(F(x))), std::exp(x),
                std::abs(std::exp(x)) * 1e-15);
    EXPECT_NEAR(static_cast<double>(exp2(F(x / 2))), std::exp2(x / 2),
                std::abs(std::exp2(x / 2)) * 1e-15);
    EXPECT_NEAR(static_cast<double>(expm1(F(x))), std::expm1(x),
                std::abs(std::expm1(x)) * 1e-15);
    const double y = std::abs(x);
    EXPECT_NEAR(static_cast<double>(log(F(y))), std::log(y),
                std::abs(std::log(y)) * 1e-15 + 1e-300);
    EXPECT_NEAR(static_cast<double>(log2(F(y))), std::log2(y),
                std::abs(std::log2(y)) * 1e-15 + 1e-300);
    EXPECT_NEAR(static_cast<double>(log10(F(y))), std::log10(y),
                std::abs(std::log10(y)) * 1e-15 + 1e-300);
    EXPECT_NEAR(static_cast<double>(log1p(F(y))), std::log1p(y),
                std::abs(std::log1p(y)) * 1e-15);
  }
}

TEST(Math, KnownDigits) {
  using F = binary256;
  EXPECT_EQ(to_string(exp(F(1))),
            to_string(from_string<F>(
                "2.7182818284590452353602874713526624977572470936999595749669"
                "676277240766303535")));
  EXPECT_EQ(to_string(log(F(2))),
            to_string(from_string<F>(
                "0.6931471805599453094172321214581765680755001343602552541206"
                "800094933936219697")));
}

TEST(Math, OverflowAndUnderflow) {
  using F = binary128;
  Context context;
  EXPECT_TRUE(exp(F(12000), context).is_inf());
  EXPECT_EQ(context.flags, unsigned{kFlagOverflow | kFlagInexact});
  context = Context{Rounding::kTowardZero};
  EXPECT_TRUE(SameBits(exp2(F(16384), context), F::MaxFinite()));
  EXPECT_EQ(context.flags, unsigned{kFlagOverflow | kFlagInexact});
  context = Context{};
  EXPECT_TRUE(SameBits(exp(F(-12000), context), F::Zero()));
  EXPECT_EQ(context.flags, unsigned{kFlagUnderflow | kFlagInexact});
  context = Context{Rounding::kUp};
  EXPECT_TRUE(SameBits(exp(F(-12000), context), F::MinSubnormal()));
  context = Context{};
  // 2^-16400 is a subnormal, exactly.
  EXPECT_TRUE(exp2(F(-16400), context).is_subnormal());
  EXPECT_EQ(context.flags, 0u);
  // Tiny arguments: e^x - 1 is just above x.
  const F tiny = F::MinSubnormal();
  context = Context{Rounding::kUp};
  EXPECT_TRUE(SameBits(expm1(tiny, context), next_up(tiny)));
  EXPECT_EQ(context.flags, unsigned{kFlagUnderflow | kFlagInexact});
  context = Context{Rounding::kDown};
  EXPECT_TRUE(SameBits(log1p(tiny, context), F::Zero()));
  EXPECT_TRUE(SameBits(exp(-tiny, context), next_down(F(1))));
}

TEST(Math, PowSpecialValues) {
  using F = binary128;
  const F inf = F::Infinity();
  const F one(1);
  const F nan = F::QuietNaN();
  EXPECT_TRUE(SameBits(pow(nan, F::Zero(true)), one));
  EXPECT_TRUE(SameBits(pow(one, nan), one));
  EXPECT_TRUE(pow(F(2), nan).is_nan());
  EXPECT_TRUE(SameBits(pow(F(-1), inf), one));
  EXPECT_TRUE(SameBits(pow(F(0.5), inf), F::Zero()));
  EXPECT_TRUE(SameBits(pow(F(0.5), -inf), inf));
  EXPECT_TRUE(SameBits(pow(F(-3), inf), inf));
  EXPECT_TRUE(SameBits(pow(-inf, F(3)), -inf));
  EXPECT_TRUE(SameBits(pow(-inf, F(-3)), F::Zero(true)));
  EXPECT_TRUE(SameBits(pow(-inf, F(2.5)), inf));
  EXPECT_TRUE(SameBits(pow(F::Zero(true), F(3)), F::Zero(true)));
  EXPECT_TRUE(SameBits(pow(F::Zero(true), F(0.5)), F::Zero()));
  EXPECT_TRUE(SameBits(pow(F(-1), F(-7)), -one));

  Context context;
  EXPECT_TRUE(SameBits(pow(F::Zero(true), F(-3), context), -inf));
  EXPECT_EQ(context.flags, unsigned{kFlagDivByZero});
  context = Context{};
  EXPECT_TRUE(SameBits(pow(F::Zero(true), -inf, context), inf));
  EXPECT_EQ(context.flags, 0u);
  EXPECT_TRUE(pow(F(-2), F(0.5), context).is_nan());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = Context{};
  EXPECT_TRUE(pow(F::SignalingNaN(), F::Zero(), context).is_nan());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
}

TEST(Math, PowExactResultsRaiseNothing) {
  using F = binary256;
  Context context;
  EXPECT_TRUE(pow(F(9), F(1.5), context) == F(27));
  EXPECT_TRUE(pow(F(-3), F(5), context) == F(-243));
  EXPECT_TRUE(pow(F(2), F(-1000), context) == scalbn(F(1), -1000));
  EXPECT_TRUE(pow(F(16), F(-0.75), context) == F(0.125));
  EXPECT_TRUE(pow(F(0.25), F(0.5), context) == F(0.5));
  EXPECT_TRUE(cbrt(F(-27), context) == F(-3));
  EXPECT_TRUE(cbrt(F(0.125), context) == F(0.5));
  EXPECT_TRUE(hypot(F(-5), F(12), context) == F(13));
  EXPECT_TRUE(hypot(F(7), F::Zero(), context) == F(7));
  EXPECT_EQ(context.flags, 0u);
  EXPECT_TRUE(pow(F(2), F(0.5), context) == sqrt(F(2)));
  EXPECT_EQ(context.flags, unsigned{kFlagInexact});
}

TEST(Math, PowRootsAndHypotMatchDouble) {
  using F = binary128;
  for (const double x : {0.001, 0.3, 0.9, 1.5, 2.0, 7.25, 1e10}) {
    for (const double y : {-50.5, -3.0, -0.5, 0.1, 1.0 / 3, 2.0, 17.75}) {
      EXPECT_NEAR(static_cast<double>(pow(F(x), F(y))), std::pow(x, y),
                  std::pow(x, y) * 1e-15);
      EXPECT_NEAR(static_cast<double>(hypot(F(x), F(y))), std::hypot(x, y),
                  std::hypot(x, y) * 1e-15);
    }
    EXPECT_NEAR(static_cast<double>(cbrt(F(-x))), std::cbrt(-x),
                std::cbrt(x) * 1e-15);
  }
  // No intermediate overflow or underflow.
  EXPECT_TRUE(hypot(F::MaxFinite(), F::MaxFinite()).is_inf());
  EXPECT_TRUE(hypot(scalbn(F(3), 16000), scalbn(F(4), 16000)) ==
              scalbn(F(5), 16000));
  const F small = F::MinSubnormal();
  EXPECT_TRUE(hypot(F(3) * small, F(4) * small) == F(5) * small);
}

}  // namespace
}  // namespace mpfix
