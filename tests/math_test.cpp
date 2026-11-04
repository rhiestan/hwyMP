// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Special values, flags and known results of the elementary functions.
// math_oracle_test compares them against MPFR.

#include "mpfix/math.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numbers>
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

TEST(Math, Hyperbolic) {
  using F = binary128;
  const F inf = F::Infinity();
  const F one(1);
  EXPECT_TRUE(SameBits(sinh(F::Zero(true)), F::Zero(true)));
  EXPECT_TRUE(SameBits(tanh(F::Zero(true)), F::Zero(true)));
  EXPECT_TRUE(SameBits(cosh(F::Zero(true)), one));
  EXPECT_TRUE(SameBits(sinh(-inf), -inf));
  EXPECT_TRUE(SameBits(cosh(-inf), inf));
  EXPECT_TRUE(SameBits(tanh(-inf), -one));
  EXPECT_TRUE(SameBits(tanh(F(1000)), one));

  Context context;
  EXPECT_TRUE(sinh(F::QuietNaN(), context).is_nan());
  EXPECT_TRUE(cosh(F::Zero(), context) == one);
  EXPECT_EQ(context.flags, 0u);
  EXPECT_TRUE(SameBits(sinh(F(-12000), context), -inf));
  EXPECT_EQ(context.flags, unsigned{kFlagOverflow | kFlagInexact});
  context = Context{Rounding::kTowardZero};
  EXPECT_TRUE(SameBits(tanh(F(-1000), context), -next_down(one)));
  EXPECT_TRUE(SameBits(cosh(F(1e-30), context), one));
  EXPECT_EQ(context.flags, unsigned{kFlagInexact});
  context = Context{Rounding::kUp};
  const F tiny = F::MinSubnormal();
  EXPECT_TRUE(SameBits(sinh(tiny, context), next_up(tiny)));
  EXPECT_TRUE(SameBits(tanh(tiny, context), tiny));
  EXPECT_EQ(context.flags, unsigned{kFlagUnderflow | kFlagInexact});

  for (const double x : {-30.0, -2.0, -0.4, -0.3, -1e-5, 1e-20, 0.25, 0.35,
                         0.36, 1.0, 5.0, 40.0, 700.0}) {
    EXPECT_NEAR(static_cast<double>(sinh(F(x))), std::sinh(x),
                std::abs(std::sinh(x)) * 1e-15);
    EXPECT_NEAR(static_cast<double>(cosh(F(x))), std::cosh(x),
                std::cosh(x) * 1e-15);
    EXPECT_NEAR(static_cast<double>(tanh(F(x))), std::tanh(x),
                std::abs(std::tanh(x)) * 1e-15);
  }
}

TEST(Math, Trigonometric) {
  using F = binary128;
  const F inf = F::Infinity();
  const F one(1);
  EXPECT_TRUE(SameBits(sin(F::Zero(true)), F::Zero(true)));
  EXPECT_TRUE(SameBits(tan(F::Zero(true)), F::Zero(true)));
  EXPECT_TRUE(SameBits(cos(F::Zero(true)), one));

  Context context;
  EXPECT_TRUE(sin(F::QuietNaN(), context).is_nan());
  EXPECT_TRUE(cos(F::Zero(), context) == one);
  EXPECT_EQ(context.flags, 0u);
  EXPECT_TRUE(sin(inf, context).is_nan());
  EXPECT_TRUE(cos(-inf, context).is_nan());
  EXPECT_TRUE(tan(inf, context).is_nan());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = Context{Rounding::kUp};
  const F tiny = F::MinSubnormal();
  EXPECT_TRUE(SameBits(sin(tiny, context), tiny));
  EXPECT_TRUE(SameBits(tan(tiny, context), next_up(tiny)));
  EXPECT_EQ(context.flags, unsigned{kFlagUnderflow | kFlagInexact});
  context = Context{Rounding::kUp};
  EXPECT_TRUE(SameBits(cos(F(1e-30), context), one));
  context = Context{Rounding::kDown};
  EXPECT_TRUE(SameBits(cos(F(1e-30), context), next_down(one)));
  EXPECT_EQ(context.flags, unsigned{kFlagInexact});

  for (const double x : {-1e22, -1e6, -30.0, -2.0, -0.5, -1e-3, 1e-20, 0.25,
                         0.75, 0.785, 0.786, 1.0, 1.5707963267948966, 3.0,
                         std::numbers::pi, 5.0, 100.0, 1e10, 1e300}) {
    EXPECT_NEAR(static_cast<double>(sin(F(x))), std::sin(x), 1e-15) << x;
    EXPECT_NEAR(static_cast<double>(cos(F(x))), std::cos(x), 1e-15) << x;
    const double t = std::tan(x);
    EXPECT_NEAR(static_cast<double>(tan(F(x))), t, std::abs(t) * 1e-14) << x;
  }
}

TEST(Math, InverseTrigonometric) {
  using F = binary128;
  const F inf = F::Infinity();
  const F one(1);
  const F pi_f = pi<F>();
  const F half_pi = scalbn(pi_f, -1);
  EXPECT_TRUE(SameBits(asin(F::Zero(true)), F::Zero(true)));
  EXPECT_TRUE(SameBits(atan(F::Zero(true)), F::Zero(true)));
  EXPECT_TRUE(SameBits(acos(one), F::Zero()));
  EXPECT_TRUE(SameBits(acos(-one), pi_f));
  EXPECT_TRUE(SameBits(asin(-one), -half_pi));
  EXPECT_TRUE(SameBits(atan(inf), half_pi));
  EXPECT_TRUE(SameBits(atan2(F::Zero(true), -one), -pi_f));
  EXPECT_TRUE(SameBits(atan2(F::Zero(), F::Zero(true)), pi_f));
  EXPECT_TRUE(SameBits(atan2(F::Zero(true), F::Zero()), F::Zero(true)));
  EXPECT_TRUE(SameBits(atan2(-one, F::Zero()), -half_pi));
  EXPECT_TRUE(SameBits(atan2(inf, inf), scalbn(pi_f, -2)));
  EXPECT_TRUE(SameBits(atan2(one, -inf), pi_f));
  EXPECT_TRUE(SameBits(atan2(-one, inf), F::Zero(true)));

  Context context;
  EXPECT_TRUE(atan2(F::Zero(), one, context) == F::Zero());
  EXPECT_TRUE(acos(one, context) == F::Zero());
  EXPECT_EQ(context.flags, 0u);
  EXPECT_TRUE(asin(F(1.5), context).is_nan());
  EXPECT_TRUE(acos(-inf, context).is_nan());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = Context{Rounding::kUp};
  const F tiny = F::MinSubnormal();
  EXPECT_TRUE(SameBits(asin(tiny, context), next_up(tiny)));
  EXPECT_TRUE(SameBits(atan(tiny, context), tiny));
  EXPECT_EQ(context.flags, unsigned{kFlagUnderflow | kFlagInexact});
  // atan2 of an exact tiny quotient lies just below it.
  context = Context{Rounding::kDown};
  const F q = scalbn(F(3), -400);
  EXPECT_TRUE(SameBits(atan2(q, one, context), next_down(q)));
  context = Context{Rounding::kUp};
  EXPECT_TRUE(SameBits(atan2(q, one, context), q));

  for (const double x : {-1.0, -0.999, -0.5, -1e-3, 1e-20, 0.25, 0.75, 0.99999,
                         1.0}) {
    EXPECT_NEAR(static_cast<double>(asin(F(x))), std::asin(x), 1e-15) << x;
    EXPECT_NEAR(static_cast<double>(acos(F(x))), std::acos(x), 1e-15) << x;
  }
  for (const double x : {-1e300, -30.0, -1.0, -0.3, 1e-20, 0.5, 2.0, 1e10}) {
    EXPECT_NEAR(static_cast<double>(atan(F(x))), std::atan(x), 1e-15) << x;
    for (const double y : {-2.0, -1e-5, 0.5, 7.0}) {
      EXPECT_NEAR(static_cast<double>(atan2(F(y), F(x))), std::atan2(y, x),
                  1e-15)
          << y << " " << x;
    }
  }
}

TEST(Math, C23ExponentialsAndLogarithms) {
  using F = binary128;
  const F inf = F::Infinity();
  Context context;
  EXPECT_TRUE(exp10(F(3), context) == F(1000));
  EXPECT_TRUE(exp2m1(F(10), context) == F(1023));
  EXPECT_TRUE(exp2m1(F(-2), context) == F(-0.75));
  EXPECT_TRUE(exp10m1(F(2), context) == F(99));
  EXPECT_TRUE(log2p1(F(7), context) == F(3));
  EXPECT_TRUE(log2p1(F(-0.75), context) == F(-2));
  EXPECT_TRUE(log10p1(F(999), context) == F(3));
  EXPECT_EQ(context.flags, 0u);
  EXPECT_TRUE(SameBits(exp10(-inf), F::Zero()));
  EXPECT_TRUE(SameBits(exp2m1(-inf), F(-1)));
  EXPECT_TRUE(SameBits(exp10m1(F::Zero(true)), F::Zero(true)));
  EXPECT_TRUE(SameBits(log2p1(F::Zero(true)), F::Zero(true)));
  EXPECT_TRUE(SameBits(log10p1(inf), inf));
  EXPECT_TRUE(SameBits(log2p1(F(-1), context), -inf));
  EXPECT_EQ(context.flags, unsigned{kFlagDivByZero});
  context = Context{};
  EXPECT_TRUE(log10p1(F(-2), context).is_nan());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  // log2(1 + 2^200) is just above 200.
  context = Context{Rounding::kUp};
  EXPECT_TRUE(SameBits(log2p1(scalbn(F(1), 200), context), next_up(F(200))));
  for (const double x : {-30.0, -2.5, -0.3, -1e-5, 1e-20, 0.25, 0.75, 5.5}) {
    EXPECT_NEAR(static_cast<double>(exp10(F(x))), std::pow(10.0, x),
                std::pow(10.0, x) * 1e-15);
    const double expected = std::expm1(x * std::numbers::ln2);
    EXPECT_NEAR(static_cast<double>(exp2m1(F(x))), expected,
                std::abs(expected) * 1e-14);
    if (x > -1) {
      EXPECT_NEAR(static_cast<double>(log2p1(F(x))), std::log1p(x) / std::log(2),
                  std::abs(std::log1p(x)) * 1e-14);
    }
  }
}

TEST(Math, ErrorAndGammaFunctions) {
  using F = binary128;
  const F inf = F::Infinity();
  EXPECT_TRUE(SameBits(erf(F::Zero(true)), F::Zero(true)));
  EXPECT_TRUE(SameBits(erf(-inf), F(-1)));
  EXPECT_TRUE(SameBits(erfc(-inf), F(2)));
  EXPECT_TRUE(SameBits(erfc(inf), F::Zero()));
  EXPECT_TRUE(SameBits(lgamma(F(1)), F::Zero()));
  EXPECT_TRUE(SameBits(lgamma(-inf), inf));
  EXPECT_TRUE(SameBits(tgamma(inf), inf));

  Context context;
  EXPECT_TRUE(tgamma(F(5), context) == F(24));
  EXPECT_TRUE(tgamma(F(21), context) == F(2432902008176640000.0));
  EXPECT_TRUE(lgamma(F(2), context) == F::Zero());
  EXPECT_EQ(context.flags, 0u);
  EXPECT_TRUE(SameBits(tgamma(F::Zero(true), context), -inf));
  EXPECT_TRUE(SameBits(lgamma(F(-3), context), inf));
  EXPECT_EQ(context.flags, unsigned{kFlagDivByZero});
  context = Context{};
  EXPECT_TRUE(tgamma(F(-3), context).is_nan());
  EXPECT_TRUE(tgamma(-inf, context).is_nan());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = Context{};
  EXPECT_TRUE(tgamma(F(2000), context).is_inf());
  EXPECT_EQ(context.flags, unsigned{kFlagOverflow | kFlagInexact});
  // erf saturates just below 1, and erfc underflows.
  context = Context{Rounding::kTowardZero};
  EXPECT_TRUE(SameBits(erf(F(20), context), next_down(F(1))));
  EXPECT_TRUE(SameBits(erfc(F(200), context), F::Zero()));
  EXPECT_EQ(context.flags, unsigned{kFlagUnderflow | kFlagInexact});

  for (const double x : {-0.7, 1e-10, 0.3, 1.0, 2.5, 5.0}) {
    EXPECT_NEAR(static_cast<double>(erf(F(x))), std::erf(x),
                std::abs(std::erf(x)) * 1e-15)
        << x;
    EXPECT_NEAR(static_cast<double>(erfc(F(x))), std::erfc(x),
                std::erfc(x) * 1e-15)
        << x;
  }
  for (const double x : {-3.5, -0.5, 1e-10, 0.25, 1.5, 7.0, 30.5, 170.0}) {
    EXPECT_NEAR(static_cast<double>(tgamma(F(x))), std::tgamma(x),
                std::abs(std::tgamma(x)) * 1e-14)
        << x;
    EXPECT_NEAR(static_cast<double>(lgamma(F(x))), std::lgamma(x),
                std::abs(std::lgamma(x)) * 1e-14 + 1e-16)
        << x;
  }
}

TEST(Math, Pi) {
  EXPECT_EQ(to_string(pi<binary256>()),
            to_string(from_string<binary256>(
                "3.1415926535897932384626433832795028841971693993751058209749"
                "445923078164062862")));
  EXPECT_EQ(static_cast<double>(pi<binary128>()), 3.141592653589793);
  Context context{Rounding::kDown};
  const binary128 low = pi<binary128>(context);
  context = Context{Rounding::kUp};
  const binary128 high = pi<binary128>(context);
  EXPECT_TRUE(next_up(low) == high);
  EXPECT_EQ(context.flags, unsigned{kFlagInexact});
}

TEST(Math, InverseHyperbolic) {
  using F = binary128;
  const F inf = F::Infinity();
  const F one(1);
  EXPECT_TRUE(SameBits(asinh(F::Zero(true)), F::Zero(true)));
  EXPECT_TRUE(SameBits(atanh(F::Zero(true)), F::Zero(true)));
  EXPECT_TRUE(SameBits(acosh(one), F::Zero()));
  EXPECT_TRUE(SameBits(asinh(-inf), -inf));
  EXPECT_TRUE(SameBits(acosh(inf), inf));

  Context context;
  EXPECT_TRUE(asinh(F::QuietNaN(), context).is_nan());
  EXPECT_TRUE(acosh(one, context) == F::Zero());
  EXPECT_EQ(context.flags, 0u);
  EXPECT_TRUE(SameBits(atanh(-one, context), -inf));
  EXPECT_EQ(context.flags, unsigned{kFlagDivByZero});
  context = Context{};
  EXPECT_TRUE(acosh(F(0.5), context).is_nan());
  EXPECT_TRUE(acosh(-inf, context).is_nan());
  EXPECT_TRUE(atanh(F(1.5), context).is_nan());
  EXPECT_TRUE(atanh(inf, context).is_nan());
  EXPECT_EQ(context.flags, unsigned{kFlagInvalid});
  context = Context{Rounding::kUp};
  const F tiny = F::MinSubnormal();
  EXPECT_TRUE(SameBits(asinh(tiny, context), tiny));
  EXPECT_TRUE(SameBits(atanh(tiny, context), next_up(tiny)));
  EXPECT_EQ(context.flags, unsigned{kFlagUnderflow | kFlagInexact});
  // No overflow at the largest finite value.
  EXPECT_NEAR(static_cast<double>(asinh(F::MaxFinite())), 11357.21655, 1e-4);
  EXPECT_NEAR(static_cast<double>(acosh(F::MaxFinite())), 11357.21655, 1e-4);

  for (const double x : {-30.0, -2.0, -0.5, -1e-3, 1e-20, 0.25, 0.75, 0.999,
                         1.0, 1.0 + 1e-12, 1.5, 2.0, 5.0, 1e10, 1e300}) {
    EXPECT_NEAR(static_cast<double>(asinh(F(x))), std::asinh(x),
                std::abs(std::asinh(x)) * 1e-15);
    if (x > 1) {
      EXPECT_NEAR(static_cast<double>(acosh(F(x))), std::acosh(x),
                  std::acosh(x) * 1e-15);
    }
    if (std::abs(x) < 1) {
      EXPECT_NEAR(static_cast<double>(atanh(F(x))), std::atanh(x),
                  std::abs(std::atanh(x)) * 1e-15);
    }
  }
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
