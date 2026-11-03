// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// What std::complex<Float> gives with the standard library at hand, checked
// against mpfix::Complex; docs/COMPLEX.md has the results. std::complex<T> is
// unspecified for T other than float, double and long double. MSVC's STL
// deprecates it for other types (STL4037) and computes some of it in double,
// so the test skips there.

#include <gtest/gtest.h>

#include <complex>

#include "mpfix/complex.h"

#if !defined(_MSVC_STL_VERSION)
#include <Eigen/Eigenvalues>

#include "mpfix/eigen.h"
#endif

namespace mpfix {
namespace {

#if defined(_MSVC_STL_VERSION)

TEST(StdComplex, MsvcStl) {
  GTEST_SKIP() << "MSVC's STL does not support std::complex<Float>";
}

#else

using F = binary128;
using S = std::complex<F>;
using C = Complex<F>;

// Within 2^-100 of |b|: full precision but for a few bits, where double would
// give 53.
::testing::AssertionResult Close(const S& a, const C& b) {
  const F bound = scalbn(abs(b), -100);
  const F error = hypot(a.real() - b.real(), a.imag() - b.imag());
  if (error <= bound) return ::testing::AssertionSuccess();
  return ::testing::AssertionFailure()
         << "(" << a.real() << ", " << a.imag() << ") vs (" << b.real() << ", "
         << b.imag() << "), error " << error;
}

TEST(StdComplex, Arithmetic) {
  const S z(F(0.75), F(-1.25));
  const S w(F(-2.5), F(1) / F(3));
  const C cz(z.real(), z.imag());
  const C cw(w.real(), w.imag());
  EXPECT_TRUE(Close(z + w, cz + cw));
  EXPECT_TRUE(Close(z - w, cz - cw));
  EXPECT_TRUE(Close(z * w, cz * cw));
  EXPECT_TRUE(Close(z / w, cz / cw));
  EXPECT_TRUE(Close(F(2) * z, F(2) * cz));
  EXPECT_TRUE(Close(std::conj(z), conj(cz)));
  EXPECT_TRUE(Close(S(std::abs(z)), C(abs(cz))));
  EXPECT_TRUE(Close(S(std::arg(z)), C(arg(cz))));
  EXPECT_TRUE(Close(S(std::norm(z)), C(norm(cz))));
  EXPECT_TRUE(Close(std::polar(F(2), F(0.5)), polar(F(2), F(0.5))));
}

TEST(StdComplex, Functions) {
  const S z(F(0.75), F(-1.25));
  const C cz(z.real(), z.imag());
  EXPECT_TRUE(Close(std::sqrt(z), sqrt(cz)));
  EXPECT_TRUE(Close(std::exp(z), exp(cz)));
  EXPECT_TRUE(Close(std::log(z), log(cz)));
  EXPECT_TRUE(Close(std::log10(z), log10(cz)));
  EXPECT_TRUE(Close(std::pow(z, S(F(0.5), F(2))), pow(cz, C(F(0.5), F(2)))));
  EXPECT_TRUE(Close(std::sin(z), sin(cz)));
  EXPECT_TRUE(Close(std::cos(z), cos(cz)));
  EXPECT_TRUE(Close(std::tan(z), tan(cz)));
  EXPECT_TRUE(Close(std::sinh(z), sinh(cz)));
  EXPECT_TRUE(Close(std::cosh(z), cosh(cz)));
  EXPECT_TRUE(Close(std::tanh(z), tanh(cz)));
  EXPECT_TRUE(Close(std::asin(z), asin(cz)));
  EXPECT_TRUE(Close(std::acos(z), acos(cz)));
  EXPECT_TRUE(Close(std::atan(z), atan(cz)));
  EXPECT_TRUE(Close(std::asinh(z), asinh(cz)));
  EXPECT_TRUE(Close(std::acosh(z), acosh(cz)));
  EXPECT_TRUE(Close(std::atanh(z), atanh(cz)));
}

// Eigen's eigenvalues of a real matrix are std::complex<F>.
TEST(StdComplex, EigenSolver) {
  using Matrix = Eigen::Matrix<F, 3, 3>;
  Matrix a;
  a << F(1), F(2), F(0), F(-2), F(1), F(0), F(0), F(0), F(3);
  const Eigen::EigenSolver<Matrix> solver(a);
  ASSERT_EQ(solver.info(), Eigen::Success);
  const auto values = solver.eigenvalues();
  for (const C expected : {C(F(1), F(2)), C(F(1), F(-2)), C(F(3))}) {
    bool found = false;
    for (int i = 0; i < 3; ++i) found = found || Close(values(i), expected);
    EXPECT_TRUE(found) << expected;
  }
  const auto vectors = solver.eigenvectors();
  const auto residual = a.cast<S>() * vectors - vectors * values.asDiagonal();
  EXPECT_TRUE(residual.cwiseAbs().maxCoeff() < scalbn(F(1), -100));
}

#endif

}  // namespace
}  // namespace mpfix
