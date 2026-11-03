// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// mpfix::boost_backend: numbers on it give exactly Float's results, convert,
// print and compare as Boost.Multiprecision's types do, and work with
// Boost.Math where that is available.

#include "mpfix/boost.h"

#include <gtest/gtest.h>

#include <boost/multiprecision/cpp_bin_float.hpp>
#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <string>

#ifdef BOOST_MP_MATH_AVAILABLE
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/zeta.hpp>
#endif

namespace mpfix {
namespace {

namespace mp = boost::multiprecision;

using Quad = binary128;
using Num = boost_binary128;

bool Same(const Num& a, const Quad& b) {
  const Quad x = a.backend().value();
  if (x.is_nan() || b.is_nan()) return x.is_nan() && b.is_nan();
  return x.to_bits() == b.to_bits();
}

// Values with all of the significand bits in use, of either sign.
Quad RandomValue(std::mt19937_64& random, int low, int high) {
  std::uniform_real_distribution<double> unit(-1, 1);
  std::uniform_int_distribution<int> exponent(low, high);
  const Quad x = Quad(unit(random)) + scalbn(Quad(unit(random)), -52) +
              scalbn(Quad(unit(random)), -104);
  return scalbn(x, exponent(random));
}

TEST(Boost, ArithmeticIsFloats) {
  std::mt19937_64 random(1);
  for (int i = 0; i < 1000; ++i) {
    const Quad x = RandomValue(random, -40, 40);
    const Quad y = RandomValue(random, -40, 40);
    const Quad z = RandomValue(random, -40, 40);
    const Num a = x, b = y, c = z;
    EXPECT_TRUE(Same(a + b, x + y));
    EXPECT_TRUE(Same(a - b, x - y));
    EXPECT_TRUE(Same(a * b, x * y));
    EXPECT_TRUE(Same(a / b, x / y));
    EXPECT_TRUE(Same(-a, -x));
    EXPECT_TRUE(Same(sqrt(abs(a)), sqrt(abs(x))));
    EXPECT_TRUE(Same(fma(a, b, c), fma(x, y, z)));
    // x - trunc(x / y) y is exact in binary512 for these exponents.
    const binary512 wx(x), wy(y);
    EXPECT_TRUE(Same(fmod(a, b), Quad(wx - trunc(wx / wy) * wy))) << i;
    EXPECT_TRUE(Same(remainder(a, b), remainder(x, y)));
    EXPECT_EQ(a < b, x < y);
    EXPECT_EQ(a == b, x == y);
  }
}

TEST(Boost, FmodIsExact) {
  const Num three = 3;
  EXPECT_TRUE(Same(fmod(Num(7.5), three), Quad(1.5)));
  EXPECT_TRUE(Same(fmod(Num(-7.5), three), Quad(-1.5)));
  EXPECT_TRUE(Same(fmod(Num(5), Num(-3)), Quad(2)));
  EXPECT_TRUE(Same(fmod(Num(-6), three), -Quad::Zero()));
  EXPECT_TRUE(Same(fmod(Num(2), Num(3)), Quad(2)));
  // 2^200 mod 3 is 1 and 2^201 mod 3 is 2, far beyond the precision of the
  // quotients.
  EXPECT_TRUE(Same(fmod(Num(scalbn(Quad(1), 200)), three), Quad(1)));
  EXPECT_TRUE(Same(fmod(Num(scalbn(Quad(1), 201)), three), Quad(2)));
  EXPECT_TRUE(isnan(fmod(Num(1), Num(0))));
  EXPECT_TRUE(Same(fmod(Num(1), std::numeric_limits<Num>::infinity()), Quad(1)));
}

TEST(Boost, FunctionsAreMpfixs) {
  std::mt19937_64 random(2);
  for (int i = 0; i < 100; ++i) {
    const Quad x = abs(RandomValue(random, -3, 3));
    const Quad y = RandomValue(random, -3, 3);
    const Quad u = x / (Quad(1) + x);  // In [0, 1).
    const Num a = x, b = y, c = u;
    EXPECT_TRUE(Same(exp(b), exp(y)));
    EXPECT_TRUE(Same(exp2(b), exp2(y)));
    EXPECT_TRUE(Same(log(a), log(x)));
    EXPECT_TRUE(Same(log2(a), log2(x)));
    EXPECT_TRUE(Same(log10(a), log10(x)));
    EXPECT_TRUE(Same(sin(b), sin(y)));
    EXPECT_TRUE(Same(cos(b), cos(y)));
    EXPECT_TRUE(Same(tan(b), tan(y)));
    EXPECT_TRUE(Same(asin(c), asin(u)));
    EXPECT_TRUE(Same(acos(c), acos(u)));
    EXPECT_TRUE(Same(atan(b), atan(y)));
    EXPECT_TRUE(Same(sinh(b), sinh(y)));
    EXPECT_TRUE(Same(cosh(b), cosh(y)));
    EXPECT_TRUE(Same(tanh(b), tanh(y)));
    EXPECT_TRUE(Same(pow(a, b), pow(x, y)));
    EXPECT_TRUE(Same(atan2(b, a), atan2(y, x)));
    EXPECT_TRUE(Same(hypot(a, b), hypot(x, y)));
    EXPECT_TRUE(Same(mp::asinh(b), asinh(y)));
    EXPECT_TRUE(Same(mp::acosh(a + 1), acosh(x + Quad(1))));
    EXPECT_TRUE(Same(mp::atanh(c), atanh(u)));
    EXPECT_TRUE(Same(mp::cbrt(b), cbrt(y)));
    EXPECT_TRUE(Same(mp::expm1(b), expm1(y)));
    EXPECT_TRUE(Same(mp::log1p(a), log1p(x)));
    EXPECT_TRUE(Same(mp::erf(b), erf(y)));
    EXPECT_TRUE(Same(mp::erfc(b), erfc(y)));
    EXPECT_TRUE(Same(mp::tgamma(b), tgamma(y)));
    EXPECT_TRUE(Same(mp::lgamma(b), lgamma(y)));
  }
}

TEST(Boost, RoundingAndExponents) {
  const Num x = 2.5;
  EXPECT_TRUE(Same(floor(-x), Quad(-3)));
  EXPECT_TRUE(Same(ceil(x), Quad(3)));
  EXPECT_TRUE(Same(trunc(-x), Quad(-2)));
  EXPECT_TRUE(Same(round(x), Quad(3)));
  EXPECT_TRUE(Same(rint(x), Quad(2)));
  EXPECT_TRUE(Same(nearbyint(Num(3.5)), Quad(4)));
  int e = 0;
  EXPECT_TRUE(Same(frexp(Num(12), &e), Quad(0.75)));
  EXPECT_EQ(e, 4);
  EXPECT_TRUE(Same(ldexp(Num(0.75), 4), Quad(12)));
  EXPECT_TRUE(Same(scalbn(Num(1), -16494), Quad::MinSubnormal()));
  EXPECT_EQ(ilogb(Num(12)), 3);
  EXPECT_TRUE(Same(logb(Num(12)), Quad(3)));
  EXPECT_TRUE(Same(mp::nextafter(Num(1), Num(2)), next_up(Quad(1))));
  EXPECT_EQ(fpclassify(Num(0)), FP_ZERO);
  EXPECT_EQ(fpclassify(Num(Quad::MinSubnormal())), FP_SUBNORMAL);
  EXPECT_EQ(fpclassify(std::numeric_limits<Num>::infinity()), FP_INFINITE);
  EXPECT_TRUE(isnan(std::numeric_limits<Num>::quiet_NaN()));
  EXPECT_TRUE(signbit(Num(-Quad::Zero())));
}

TEST(Boost, Conversions) {
  // From built-in types and text, exactly where the value fits.
  EXPECT_TRUE(Same(Num(-7), Quad(-7)));
  EXPECT_TRUE(Same(Num(9007199254740993LL), Quad(9007199254740993LL)));
  EXPECT_TRUE(Same(Num(18446744073709551615ULL), Quad(18446744073709551615ULL)));
  EXPECT_TRUE(Same(Num(0.1), Quad(0.1)));
  EXPECT_TRUE(Same(Num(0.1f), Quad(0.1f)));
  EXPECT_TRUE(Same(Num("0.1"), from_string<Quad>("0.1")));
  EXPECT_TRUE(Same(Num(std::string("-0x1.8p3")), Quad(-12)));
  EXPECT_THROW(Num("zero point one"), std::runtime_error);
  // To built-in types: rounded to nearest, or truncated to integers.
  const Num third = Num(1) / 3;
  EXPECT_EQ(third.convert_to<double>(), 1.0 / 3.0);
  EXPECT_EQ(static_cast<long long>(Num(-7.9)), -7);
  EXPECT_EQ(static_cast<unsigned>(Num(42.5)), 42u);
  EXPECT_EQ(static_cast<float>(third), 1.0f / 3.0f);
  EXPECT_TRUE(static_cast<bool>(third));
  // Between formats: exact when widening.
  const boost_binary256 wide = boost_binary256(third);
  EXPECT_EQ(wide.backend().value(), binary256(third.backend().value()));
  EXPECT_TRUE(Same(Num(wide), third.backend().value()));
  // To Boost's own types, through its generic conversion.
  const mp::cpp_bin_float_quad quad(third);
  EXPECT_EQ(quad, mp::cpp_bin_float_quad(1) / 3);
}

TEST(Boost, StringsAndStreams) {
  const Num third = Num(1) / 3;
  // With no digits, the shortest text that reads back.
  EXPECT_EQ(third.str(), "0.3333333333333333333333333333333333");
  EXPECT_EQ(Num(third.str()), third);
  EXPECT_EQ(Num("1e300").str(0, std::ios_base::scientific), "1e+300");
  EXPECT_EQ(third.str(5), "0.33333");
  EXPECT_EQ(third.str(3, std::ios_base::scientific), "3.333e-01");
  EXPECT_EQ(Num(2.5).str(2, std::ios_base::fixed | std::ios_base::showpos),
            "+2.50");
  std::ostringstream out;
  out << std::setprecision(10) << third << ' ' << std::setw(6) << Num(-2);
  EXPECT_EQ(out.str(), "0.3333333333     -2");
  std::istringstream in("1.25e2 nan");
  Num read, nan;
  in >> read >> nan;
  EXPECT_EQ(read, Num(125));
  EXPECT_TRUE(isnan(nan));
}

TEST(Boost, ComparisonsAndLimits) {
  const Num x = 0.5;
  EXPECT_TRUE(x < 1 && x > 0.25 && x == 0.5 && x != 1u && x >= 0.5L);
  // Exact comparison with integers that do not fit double.
  EXPECT_TRUE(Num(9007199254740993LL) > 9007199254740992LL);
  EXPECT_TRUE(Num(9007199254740993LL) == 9007199254740993LL);
  const Num nan = std::numeric_limits<Num>::quiet_NaN();
  EXPECT_FALSE(nan == nan);
  EXPECT_FALSE(nan < 1 || nan > 1);
  using L = std::numeric_limits<Num>;
  static_assert(L::is_specialized && L::digits == 113 && L::digits10 == 33);
  static_assert(L::max_digits10 == 36 && L::max_exponent == 16384);
  EXPECT_TRUE(Same(L::epsilon(), std::numeric_limits<Quad>::epsilon()));
  EXPECT_TRUE(Same((L::max)(), Quad::MaxFinite()));
  EXPECT_TRUE(Same(L::denorm_min(), Quad::MinSubnormal()));
  EXPECT_EQ(std::hash<Num>()(Num(0)), std::hash<Num>()(-Num(0)));
}

TEST(Boost, ExpressionTemplates) {
  using T = mp::number<boost_backend<237, 19>, mp::et_on>;
  const T a = 3, b = T(1) / 7;
  const T c = a * b + sqrt(a) - b / a;
  const binary256 x = 3, y = binary256(1) / 7;
  EXPECT_EQ(c.backend().value(), x * y + sqrt(x) - y / x);
  EXPECT_EQ(T(exp(a * b)).backend().value(), exp(x * y));
}

// The example in README.md.
TEST(Boost, ReadmeExample) {
  mpfix::boost_binary256 x = 2;
  mpfix::boost_binary256 y = sqrt(x) + boost::multiprecision::erf(x);
  std::ostringstream out;
  out << std::setprecision(70) << y << '\n';
  EXPECT_EQ(y.backend().value(), sqrt(binary256(2)) + erf(binary256(2)));
  // sqrt(2) + erf(2), to 70 digits.
  EXPECT_EQ(out.str().substr(0, 12), "2.4095358273");
  EXPECT_GT(out.str().size(), 60u);
}

#ifdef BOOST_MP_MATH_AVAILABLE
// Boost.Math runs on the backend, at its full precision.
TEST(Boost, BoostMath) {
  const Num pi = boost::math::constants::pi<Num>();
  EXPECT_TRUE(Same(pi, mpfix::pi<Quad>()));
  // zeta(2) = pi^2 / 6, within a few ulps.
  const Num zeta = boost::math::zeta(Num(2));
  const Num expected = pi * pi / 6;
  EXPECT_LT(abs(zeta - expected) / expected, 64 * std::numeric_limits<Num>::epsilon());
  // J0 vanishes at its first zero.
  const Num zero = boost::math::cyl_bessel_j_zero(Num(0), 1);
  EXPECT_LT(abs(boost::math::cyl_bessel_j(0, zero)),
            1024 * std::numeric_limits<Num>::epsilon());
}
#endif

}  // namespace
}  // namespace mpfix
