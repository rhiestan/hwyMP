// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Special values, exact results and basic properties of mpfix::Complex;
// complex_oracle_test compares it against MPC.

#include "mpfix/complex.h"

#include <gtest/gtest.h>

#include <complex>
#include <sstream>

namespace mpfix {
namespace {

bool SameResult(const binary128& a, const binary128& b) {
  if (a.is_nan() || b.is_nan()) return a.is_nan() && b.is_nan();
  return a.to_bits() == b.to_bits();
}

// The special values of Annex G.
TEST(Complex, SpecialValues) {
  using F = binary128;
  using C = Complex<F>;
  const F inf = F::Infinity();
  const F nan = F::QuietNaN();
  const F zero = F::Zero();
  const auto same = [](const C& a, const C& b) {
    return SameResult(a.real(), b.real()) && SameResult(a.imag(), b.imag());
  };
  // Multiplication and division recover infinities from NaNs.
  const C product = C(inf, nan) * C(F(1), F(1));
  EXPECT_TRUE(product.real().is_inf() || product.imag().is_inf());
  const C quotient = C(F(1), F(1)) / C(zero, zero);
  EXPECT_TRUE(quotient.real().is_inf() && quotient.imag().is_inf());
  EXPECT_TRUE(same(C(F(3), F(4)) / C(inf, F(1)), C(zero, zero)));
  // sqrt.
  EXPECT_TRUE(same(sqrt(C(nan, inf)), C(inf, inf)));
  EXPECT_TRUE(same(sqrt(C(-inf, F(2))), C(zero, inf)));
  EXPECT_TRUE(same(sqrt(C(inf, -F(2))), C(inf, -zero)));
  EXPECT_TRUE(same(sqrt(C(F(-4), -zero)), C(zero, F(-2))));
  // exp and log.
  EXPECT_TRUE(same(exp(C(-zero, zero)), C(F(1), zero)));
  EXPECT_TRUE(same(exp(C(-inf, F(1))), C(zero * cos(F(1)), zero * sin(F(1)))));
  EXPECT_TRUE(exp(C(F(1), inf)).real().is_nan());
  EXPECT_TRUE(same(log(C(-zero, zero)), C(-inf, pi<F>())));
  EXPECT_TRUE(same(log(C(inf, nan)), C(inf, nan)));
  EXPECT_TRUE(same(log(C(F(1), zero)), C(zero, zero)));
  // The hyperbolic functions and, through them, the trigonometric ones.
  EXPECT_TRUE(same(sinh(C(zero, zero)), C(zero, zero)));
  EXPECT_TRUE(same(sinh(C(-zero, -zero)), C(-zero, -zero)));
  EXPECT_TRUE(same(cosh(C(zero, zero)), C(F(1), zero)));
  EXPECT_TRUE(same(tanh(C(inf, F(1))), C(F(1), zero)));
  EXPECT_TRUE(same(tanh(C(-inf, F(1))), C(F(-1), zero)));
  EXPECT_TRUE(same(sinh(C(inf, zero)), C(inf, zero)));
  EXPECT_TRUE(same(cos(C(zero, zero)), C(F(1), -zero)));
  // sin(x + i0) = sin x + i cos x sinh(+0), whose zero has the sign of cos x.
  EXPECT_TRUE(same(sin(C(F(2), zero)), C(sin(F(2)), -zero)));
  // proj, conj and polar.
  EXPECT_TRUE(same(proj(C(nan, -inf)), C(inf, -zero)));
  EXPECT_TRUE(same(conj(C(F(1), zero)), C(F(1), -zero)));
  EXPECT_TRUE(same(polar(F(2), zero), C(F(2), zero)));
  // Exact results.
  EXPECT_TRUE(same(C(F(1), F(1)) * C(F(1), F(1)), C(zero, F(2))));
  EXPECT_TRUE(same(C(F(1), F(1)) / C(F(1), F(-1)), C(zero, F(1))));
  EXPECT_TRUE(same(sqrt(C(F(3), F(4))), C(F(2), F(1))));
  EXPECT_TRUE(same(pow(C(F(1), F(1)), C(F(4))), C(F(-4), zero)));
}


TEST(Complex, MatchesDouble) {
  using F = binary128;
  using C = Complex<F>;
  const std::complex<double> values[] = {{0.3, -1.7}, {-2.5, 0.25},
                                         {1e-3, 4.0}, {-0.9, -0.4}};
  const auto close = [](const C& a, std::complex<double> b) {
    const double scale = std::abs(b) * 1e-14 + 1e-300;
    return std::abs(static_cast<double>(a.real()) - b.real()) <= scale &&
           std::abs(static_cast<double>(a.imag()) - b.imag()) <= scale;
  };
  for (const auto& d : values) {
    const C z(F(d.real()), F(d.imag()));
    const C w(F(d.imag()), F(-d.real() / 3));
    const std::complex<double> dw(d.imag(), -d.real() / 3);
    EXPECT_TRUE(close(z * w, d * dw));
    EXPECT_TRUE(close(z / w, d / dw));
    EXPECT_TRUE(close(sqrt(z), std::sqrt(d)));
    EXPECT_TRUE(close(exp(z), std::exp(d)));
    EXPECT_TRUE(close(log(z), std::log(d)));
    EXPECT_TRUE(close(log10(z), std::log10(d)));
    EXPECT_TRUE(close(sin(z), std::sin(d)));
    EXPECT_TRUE(close(cos(z), std::cos(d)));
    EXPECT_TRUE(close(tan(z), std::tan(d)));
    EXPECT_TRUE(close(sinh(z), std::sinh(d)));
    EXPECT_TRUE(close(cosh(z), std::cosh(d)));
    EXPECT_TRUE(close(tanh(z), std::tanh(d)));
    EXPECT_TRUE(close(pow(z, w), std::pow(d, dw)));
    EXPECT_NEAR(static_cast<double>(abs(z)), std::abs(d), 1e-15);
    EXPECT_NEAR(static_cast<double>(arg(z)), std::arg(d), 1e-15);
    EXPECT_NEAR(static_cast<double>(norm(z)), std::norm(d), 1e-14);
  }
}

TEST(Complex, Interface) {
  using F = binary256;
  using C = Complex<F>;
  C z(F(1), F(2));
  z += C(F(3));
  z *= F(2);
  z -= F(1);
  z /= C(F(0), F(1));
  EXPECT_TRUE(z == C(F(4), F(-7)));
  EXPECT_TRUE(real(z) == F(4) && imag(z) == F(-7));
  EXPECT_TRUE(C(F(5)) == F(5));
  EXPECT_TRUE(-z == C(F(-4), F(7)));
  std::ostringstream out;
  out << C(F(1.5), F(-2));
  EXPECT_EQ(out.str(), "(1.5,-2)");
  EXPECT_TRUE(C().real().is_zero() && C().imag().is_zero());
}

}  // namespace
}  // namespace mpfix
