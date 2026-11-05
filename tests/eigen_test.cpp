// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Float as an Eigen scalar: the limits, and the dense decompositions, each
// checked by a residual that must be a small multiple of n * epsilon. Such a
// residual shows that the work was done at the format's full precision.

#include "mpfix/eigen.h"

#include <gtest/gtest.h>

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <bit>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>

namespace mpfix {
namespace {

// --- numeric_limits ---

TEST(Limits, MatchDouble) {
  using L = std::numeric_limits<Float<53, 11>>;
  using D = std::numeric_limits<double>;
  static_assert(L::is_specialized && L::is_signed && !L::is_integer);
  static_assert(L::is_iec559 && L::has_infinity && L::has_quiet_NaN &&
                L::has_signaling_NaN);
  static_assert(L::digits == D::digits);
  static_assert(L::digits10 == D::digits10);
  static_assert(L::max_digits10 == D::max_digits10);
  static_assert(L::min_exponent == D::min_exponent);
  static_assert(L::max_exponent == D::max_exponent);
  static_assert(L::min_exponent10 == D::min_exponent10);
  static_assert(L::max_exponent10 == D::max_exponent10);
  static_assert(L::radix == D::radix);
  const auto same = [](const Float<53, 11>& x, double d) {
    return x.to_bits()[0] == std::bit_cast<std::uint64_t>(d);
  };
  EXPECT_TRUE(same(L::epsilon(), D::epsilon()));
  EXPECT_TRUE(same(L::round_error(), D::round_error()));
  EXPECT_TRUE(same((L::min)(), (D::min)()));
  EXPECT_TRUE(same((L::max)(), (D::max)()));
  EXPECT_TRUE(same(L::lowest(), D::lowest()));
  EXPECT_TRUE(same(L::denorm_min(), D::denorm_min()));
  EXPECT_TRUE(same(L::infinity(), D::infinity()));
  EXPECT_TRUE(L::quiet_NaN().is_nan() && !L::quiet_NaN().is_signaling());
  EXPECT_TRUE(L::signaling_NaN().is_signaling());
}

TEST(Limits, Binary128) {
  using L = std::numeric_limits<binary128>;
  // As for the IEEE binary128 type of other languages.
  static_assert(L::digits == 113);
  static_assert(L::digits10 == 33);
  static_assert(L::max_digits10 == 36);
  static_assert(L::min_exponent == -16381);
  static_assert(L::max_exponent == 16384);
  static_assert(L::min_exponent10 == -4931);
  static_assert(L::max_exponent10 == 4932);
  EXPECT_EQ(L::epsilon(), next_up(binary128(1)) - binary128(1));
  static_assert(std::numeric_limits<binary256>::digits10 == 71);
  static_assert(std::numeric_limits<binary512>::digits10 == 146);
}

// --- Decompositions ---

template <class F>
using Matrix = Eigen::Matrix<F, Eigen::Dynamic, Eigen::Dynamic>;
template <class F>
using Vector = Eigen::Matrix<F, Eigen::Dynamic, 1>;

// A matrix with entries in [-1, 1) that use every significand bit, so that
// nothing is exact by accident.
template <class F>
Matrix<F> RandomMatrix(std::mt19937_64& random, int rows, int cols) {
  std::uniform_real_distribution<double> uniform(-1, 1);
  Matrix<F> m(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      m(i, j) = F(uniform(random)) + scalbn(F(uniform(random)), -50) +
                scalbn(F(uniform(random)), -100) +
                scalbn(F(uniform(random)), -150);
    }
  }
  return m;
}

template <class F>
F Norm(const Matrix<F>& m) {
  return m.cwiseAbs().rowwise().sum().maxCoeff();  // The infinity norm.
}

template <class F>
F Epsilon() {
  return Eigen::NumTraits<F>::epsilon();
}

// Expects |residual| <= factor * n * epsilon * scale.
template <class F>
void ExpectSmall(const Matrix<F>& residual, const F& scale, int n,
                 const char* what) {
  const F bound = F(16) * F(n) * Epsilon<F>() * scale;
  const F error = Norm<F>(residual);
  EXPECT_TRUE(error <= bound)
      << what << ": residual " << std::scientific << std::setprecision(3)
      << error << ", bound " << bound;
  // A double-precision solution would leave a far larger residual.
  EXPECT_TRUE(error < F(1e-20) * scale || Epsilon<F>() > F(1e-25)) << what;
}

template <class F>
class EigenDecompositions : public ::testing::Test {};

using Formats = ::testing::Types<binary128, binary256, binary512>;
TYPED_TEST_SUITE(EigenDecompositions, Formats);

TYPED_TEST(EigenDecompositions, LinearSolvers) {
  using F = TypeParam;
  std::mt19937_64 random(1);
  for (const int n : {1, 5, 12}) {
    const Matrix<F> a = RandomMatrix<F>(random, n, n);
    const Matrix<F> b = RandomMatrix<F>(random, n, 2);
    const auto check = [&](const Matrix<F>& x, const char* what) {
      ExpectSmall<F>(a * x - b, Norm<F>(a) * Norm<F>(x) + Norm<F>(b), n, what);
    };
    check(a.partialPivLu().solve(b), "PartialPivLU");
    check(a.fullPivLu().solve(b), "FullPivLU");
    check(a.householderQr().solve(b), "HouseholderQR");
    check(a.colPivHouseholderQr().solve(b), "ColPivHouseholderQR");
    check(a.fullPivHouseholderQr().solve(b), "FullPivHouseholderQR");
    check(a.completeOrthogonalDecomposition().solve(b),
          "CompleteOrthogonalDecomposition");

    // Inverse and determinant.
    const Matrix<F> inverse = a.inverse();
    ExpectSmall<F>(a * inverse - Matrix<F>::Identity(n, n),
                   Norm<F>(a) * Norm<F>(inverse), n, "inverse");
    const F det = a.determinant();
    const F det_lu = a.fullPivLu().determinant();
    EXPECT_TRUE(abs(det - det_lu) <=
                F(64) * F(n) * Epsilon<F>() * Norm<F>(a) * F(n) * abs(det) +
                    Epsilon<F>())
        << "determinant";
  }
}

TYPED_TEST(EigenDecompositions, Cholesky) {
  using F = TypeParam;
  std::mt19937_64 random(2);
  for (const int n : {1, 6, 12}) {
    const Matrix<F> a = RandomMatrix<F>(random, n, n);
    // Symmetric positive definite.
    const Matrix<F> spd = a.transpose() * a + F(n) * Matrix<F>::Identity(n, n);
    const Matrix<F> b = RandomMatrix<F>(random, n, 1);
    const Eigen::LLT<Matrix<F>> llt(spd);
    ASSERT_EQ(llt.info(), Eigen::Success);
    const Matrix<F> l = llt.matrixL();
    ExpectSmall<F>(l * l.transpose() - spd, Norm<F>(spd), n, "LLT factors");
    const Matrix<F> x = llt.solve(b);
    ExpectSmall<F>(spd * x - b, Norm<F>(spd) * Norm<F>(x), n, "LLT solve");
    const Eigen::LDLT<Matrix<F>> ldlt(spd);
    ASSERT_EQ(ldlt.info(), Eigen::Success);
    ExpectSmall<F>(spd * ldlt.solve(b) - b, Norm<F>(spd) * Norm<F>(x), n,
                   "LDLT solve");
  }
}

TYPED_TEST(EigenDecompositions, SingularValues) {
  using F = TypeParam;
  std::mt19937_64 random(3);
  for (const auto& [rows, cols] : {std::pair{1, 1}, {7, 7}, {10, 6}}) {
    const Matrix<F> a = RandomMatrix<F>(random, rows, cols);
    const auto check = [&](const auto& svd, const char* what) {
      const Matrix<F>& u = svd.matrixU();
      const Matrix<F>& v = svd.matrixV();
      const Vector<F>& s = svd.singularValues();
      const int n = std::max(rows, cols);
      ExpectSmall<F>(u * s.asDiagonal() * v.transpose() - a, Norm<F>(a), n,
                     what);
      ExpectSmall<F>(u.transpose() * u - Matrix<F>::Identity(cols, cols), F(1),
                     n, what);
      for (int i = 1; i < s.size(); ++i) {
        EXPECT_TRUE(s(i - 1) >= s(i)) << what << ": unsorted";
      }
    };
    check(
        Eigen::JacobiSVD<Matrix<F>, Eigen::ComputeThinU | Eigen::ComputeThinV>(
            a),
        "JacobiSVD");
    check(
        Eigen::BDCSVD<Matrix<F>, Eigen::ComputeThinU | Eigen::ComputeThinV>(a),
        "BDCSVD");
  }
}

TYPED_TEST(EigenDecompositions, Eigenvalues) {
  using F = TypeParam;
  std::mt19937_64 random(4);
  for (const int n : {1, 6, 12}) {
    const Matrix<F> a = RandomMatrix<F>(random, n, n);
    const Matrix<F> symmetric = a + a.transpose();
    const Eigen::SelfAdjointEigenSolver<Matrix<F>> solver(symmetric);
    ASSERT_EQ(solver.info(), Eigen::Success);
    const Matrix<F>& v = solver.eigenvectors();
    const Vector<F>& d = solver.eigenvalues();
    ExpectSmall<F>(symmetric * v - v * d.asDiagonal(), Norm<F>(symmetric), n,
                   "SelfAdjointEigenSolver");
    ExpectSmall<F>(v.transpose() * v - Matrix<F>::Identity(n, n), F(1), n,
                   "eigenvectors orthonormal");

    // The real Schur form of a general matrix: A = U T U^T.
    const Eigen::RealSchur<Matrix<F>> schur(a);
    ASSERT_EQ(schur.info(), Eigen::Success);
    const Matrix<F>& u = schur.matrixU();
    const Matrix<F>& t = schur.matrixT();
    ExpectSmall<F>(u * t * u.transpose() - a, Norm<F>(a), n, "RealSchur");
  }
}

TYPED_TEST(EigenDecompositions, MoreAccurateThanDouble) {
  // The same Hilbert system, which is badly conditioned, solved in double
  // and in F: F's answer is closer to the exact one by many digits.
  using F = TypeParam;
  constexpr int n = 8;
  Matrix<F> hilbert(n, n);
  Eigen::MatrixXd hilbert_double(n, n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      hilbert(i, j) = F(1) / F(i + j + 1);
      hilbert_double(i, j) = 1.0 / (i + j + 1);
    }
  }
  // The right-hand side for the solution of all ones.
  const Vector<F> b = hilbert * Vector<F>::Ones(n);
  const Vector<F> x = hilbert.fullPivLu().solve(b);
  const Eigen::VectorXd x_double = hilbert_double.fullPivLu().solve(
      hilbert_double * Eigen::VectorXd::Ones(n));
  const F error = (x - Vector<F>::Ones(n)).cwiseAbs().maxCoeff();
  const double error_double =
      (x_double - Eigen::VectorXd::Ones(n)).cwiseAbs().maxCoeff();
  // The condition number is about 1.5e10.
  EXPECT_GT(error_double, 1e-12);
  EXPECT_TRUE(error < F(1e10) * Epsilon<F>() * F(100)) << error;
}

// --- Complex matrices ---

template <class F>
Matrix<Complex<F>> RandomComplexMatrix(std::mt19937_64& random, int rows,
                                       int cols) {
  const Matrix<F> re = RandomMatrix<F>(random, rows, cols);
  const Matrix<F> im = RandomMatrix<F>(random, rows, cols);
  Matrix<Complex<F>> m(rows, cols);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) m(i, j) = Complex<F>(re(i, j), im(i, j));
  }
  return m;
}

template <class F>
F ComplexNorm(const Matrix<Complex<F>>& m) {
  return m.cwiseAbs().rowwise().sum().maxCoeff();
}

template <class F>
void ExpectSmallComplex(const Matrix<Complex<F>>& residual, const F& scale, int n,
                 const char* what) {
  const F bound = F(16) * F(n) * Epsilon<F>() * scale;
  const F error = ComplexNorm<F>(residual);
  EXPECT_TRUE(error <= bound)
      << what << ": residual " << std::scientific << std::setprecision(3)
      << error << ", bound " << bound;
}

template <class F>
class EigenComplex : public ::testing::Test {};

using ComplexFormats = ::testing::Types<binary128, binary256>;
TYPED_TEST_SUITE(EigenComplex, ComplexFormats);

TYPED_TEST(EigenComplex, LinearSolvers) {
  using F = TypeParam;
  using C = Complex<F>;
  std::mt19937_64 random(5);
  for (const int n : {1, 5, 10}) {
    const Matrix<C> a = RandomComplexMatrix<F>(random, n, n);
    const Matrix<C> b = RandomComplexMatrix<F>(random, n, 2);
    const auto check = [&](const Matrix<C>& x, const char* what) {
      ExpectSmallComplex<F>(a * x - b, ComplexNorm<F>(a) * ComplexNorm<F>(x) + ComplexNorm<F>(b), n, what);
    };
    check(a.partialPivLu().solve(b), "complex PartialPivLU");
    check(a.fullPivLu().solve(b), "complex FullPivLU");
    check(a.householderQr().solve(b), "complex HouseholderQR");
    check(a.colPivHouseholderQr().solve(b), "complex ColPivHouseholderQR");
    const Matrix<C> inverse = a.inverse();
    ExpectSmallComplex<F>(a * inverse - Matrix<C>::Identity(n, n),
                   ComplexNorm<F>(a) * ComplexNorm<F>(inverse), n, "complex inverse");

    // Hermitian positive definite.
    const Matrix<C> hpd = a.adjoint() * a + C(F(n)) * Matrix<C>::Identity(n, n);
    const Eigen::LLT<Matrix<C>> llt(hpd);
    ASSERT_EQ(llt.info(), Eigen::Success);
    const Matrix<C> l = llt.matrixL();
    ExpectSmallComplex<F>(l * l.adjoint() - hpd, ComplexNorm<F>(hpd), n, "complex LLT");
  }
}

TYPED_TEST(EigenComplex, Decompositions) {
  using F = TypeParam;
  using C = Complex<F>;
  std::mt19937_64 random(6);
  for (const int n : {1, 6, 10}) {
    const Matrix<C> a = RandomComplexMatrix<F>(random, n, n);
    const Matrix<C> identity = Matrix<C>::Identity(n, n);

    const Eigen::JacobiSVD<Matrix<C>, Eigen::ComputeThinU | Eigen::ComputeThinV>
        svd(a);
    const Matrix<C>& u = svd.matrixU();
    const Matrix<C>& v = svd.matrixV();
    // The singular values are real, and mix with the complex factors.
    const Vector<F>& s = svd.singularValues();
    ExpectSmallComplex<F>(u * s.template cast<C>().asDiagonal() * v.adjoint() - a,
                   ComplexNorm<F>(a), n, "complex JacobiSVD");
    ExpectSmallComplex<F>(u.adjoint() * u - identity, F(1), n, "complex SVD unitary");

    const Matrix<C> hermitian = a + a.adjoint();
    const Eigen::SelfAdjointEigenSolver<Matrix<C>> self_adjoint(hermitian);
    ASSERT_EQ(self_adjoint.info(), Eigen::Success);
    const Matrix<C>& w = self_adjoint.eigenvectors();
    const Vector<F>& lambda = self_adjoint.eigenvalues();
    ExpectSmallComplex<F>(hermitian * w - w * lambda.template cast<C>().asDiagonal(),
                   ComplexNorm<F>(hermitian), n, "Hermitian eigenvectors");

    // The complex Schur form and the eigenvalues of a general matrix.
    const Eigen::ComplexSchur<Matrix<C>> schur(a);
    ASSERT_EQ(schur.info(), Eigen::Success);
    const Matrix<C>& q = schur.matrixU();
    const Matrix<C>& t = schur.matrixT();
    ExpectSmallComplex<F>(q * t * q.adjoint() - a, ComplexNorm<F>(a), n, "ComplexSchur");
    const Eigen::ComplexEigenSolver<Matrix<C>> solver(a);
    ASSERT_EQ(solver.info(), Eigen::Success);
    const Matrix<C>& vectors = solver.eigenvectors();
    const Vector<C>& values = solver.eigenvalues();
    ExpectSmallComplex<F>(a * vectors - vectors * values.asDiagonal(),
                   ComplexNorm<F>(a) * ComplexNorm<F>(vectors), n, "ComplexEigenSolver");
  }
}

TEST(EigenComplex, MixedExpressions) {
  using F = binary128;
  using C = Complex<F>;
  Eigen::Matrix<C, 2, 2> m;
  m << C(F(1), F(1)), C(F(2)), C(F(0), F(-1)), C(F(3), F(2));
  Eigen::Matrix<F, 2, 2> r;
  r << F(1), F(2), F(3), F(4);
  // Real scalars and matrices times complex ones.
  const Eigen::Matrix<C, 2, 2> scaled = F(2) * m;
  EXPECT_TRUE(scaled(1, 1) == C(F(6), F(4)));
  const Eigen::Matrix<C, 2, 2> product = r * m;
  EXPECT_TRUE(product(0, 0) == C(F(1), F(-1)));
  const Eigen::Matrix<C, 2, 2> sum = m + r.cast<C>();
  EXPECT_TRUE(sum(1, 0) == C(F(3), F(-1)));
  // adjoint, the real and imaginary parts, norms and the trace.
  EXPECT_TRUE(m.adjoint()(0, 1) == C(F(0), F(1)));
  EXPECT_TRUE(m.real()(1, 1) == F(3) && m.imag()(0, 0) == F(1));
  EXPECT_EQ(m.squaredNorm(), F(20));
  EXPECT_EQ(m.norm(), sqrt(F(20)));
  EXPECT_TRUE(m.trace() == C(F(4), F(3)));
  EXPECT_EQ(m.cwiseAbs()(1, 1), sqrt(F(13)));
  EXPECT_TRUE(m.allFinite());
  EXPECT_TRUE(m.isApprox(m * C(F(1), scalbn(F(1), -120))));
  const Eigen::Matrix<C, 2, 1> random = Eigen::Matrix<C, 2, 1>::Random();
  EXPECT_TRUE(abs(random(0).real()) <= F(1) && abs(random(1).imag()) <= F(1));
  std::ostringstream out;
  out << m.row(0);
  EXPECT_EQ(out.str(), "    (1,1)     (2,0)");
}

// The example in README.md.
TEST(Eigen, NumTraits) {
  using T = Eigen::NumTraits<binary128>;
  using C = Eigen::NumTraits<Complex<binary128>>;
  EXPECT_TRUE(T::infinity().is_inf() && !T::infinity().signbit());
  EXPECT_TRUE(T::quiet_NaN().is_nan());
  // Called at run time, as Eigen's printing does.
  volatile int digits = T::max_digits10();
  EXPECT_EQ(digits, 36);
  digits = C::max_digits10();
  EXPECT_EQ(digits, 36);
  EXPECT_EQ(C::epsilon(), T::epsilon());
}

TEST(Eigen, ReadmeExample) {
  using Matrix =
      Eigen::Matrix<mpfix::binary256, Eigen::Dynamic, Eigen::Dynamic>;

  Matrix a = Matrix::Random(10, 10);
  Matrix x = a.partialPivLu().solve(Matrix::Ones(10, 1));
  std::ostringstream out;
  out << x.format(Eigen::IOFormat(Eigen::FullPrecision)) << '\n';

  // Random entries lie in [-1, 1] and use more than double's 53 bits.
  EXPECT_TRUE(a.cwiseAbs().maxCoeff() <= binary256(1));
  int long_entries = 0;
  for (const binary256& entry : a.reshaped()) {
    if (binary256(static_cast<double>(entry)) != entry) ++long_entries;
  }
  EXPECT_GT(long_entries, 90);
  ExpectSmall<binary256>(a * x - Matrix::Ones(10, 1),
                         Norm<binary256>(a) * Norm<binary256>(x), 10,
                         "README solve");
  const std::string text = out.str();
  EXPECT_EQ(std::count(text.begin(), text.end(), '\n'), 10);
}

TEST(Eigen, ExpressionsAndPrinting) {
  Eigen::Matrix<binary128, 2, 2> m;
  m << 1, 2, 3, 4;
  const Eigen::Matrix<binary128, 2, 2> square = m * m;
  EXPECT_EQ(square(0, 0), binary128(7));
  EXPECT_EQ(square(1, 1), binary128(22));
  EXPECT_EQ(m.sum(), binary128(10));
  EXPECT_EQ(m.trace(), binary128(5));
  EXPECT_EQ(m.maxCoeff(), binary128(4));
  EXPECT_TRUE(m.isApprox(m + Eigen::Matrix<binary128, 2, 2>::Constant(
                                 scalbn(binary128(1), -120))));
  EXPECT_FALSE(m.isApprox(m * binary128(1.0001)));
  EXPECT_EQ(m.norm(), sqrt(binary128(30)));
  EXPECT_TRUE(m.allFinite());
  m(0, 1) = binary128::QuietNaN();
  EXPECT_TRUE(m.hasNaN());

  Eigen::Matrix<binary128, 1, 2> row;
  row << binary128(1) / binary128(3), binary128(-2.5);
  std::ostringstream out;
  out << row;
  // At the stream's precision by default.
  EXPECT_EQ(out.str(), "0.333333     -2.5");
  // FullPrecision prints NumTraits::max_digits10 digits, which read back
  // exactly.
  out.str("");
  out << row.format(Eigen::IOFormat(Eigen::FullPrecision));
  std::istringstream in(out.str());
  binary128 first;
  binary128 second;
  in >> first >> second;
  EXPECT_EQ(first, row(0)) << out.str();
  EXPECT_EQ(second, row(1));

  // Casting between scalar types goes through the conversions.
  const Eigen::Matrix<double, 1, 2> narrow = row.cast<double>();
  EXPECT_EQ(narrow(0), 1.0 / 3.0);
  const Eigen::Matrix<binary256, 1, 2> wide = row.cast<binary256>();
  EXPECT_EQ(wide(1), binary256(-2.5));
}

}  // namespace
}  // namespace mpfix
