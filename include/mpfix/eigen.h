// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Eigen support: include this header, then use Float as a matrix scalar,
// for example Eigen::Matrix<mpfix::binary256, Eigen::Dynamic, Eigen::Dynamic>.
// Matrices print through the stream operators of mpfix/charconv.h, at the
// stream's precision, or with Eigen::FullPrecision at max_digits10, which
// reads back exactly.
//
// Eigen finds abs, sqrt, exp, log, sin, atan2 and the other functions of
// mpfix/math.h, isnan, isinf and isfinite by argument-dependent lookup.
//
// mpfix::Complex<F> of mpfix/complex.h is a complex scalar too, which mixes
// with F: Eigen::Matrix<mpfix::complex256, Eigen::Dynamic, Eigen::Dynamic> works
// with the dense decompositions, ComplexSchur and ComplexEigenSolver. Eigen's
// eigenvalues of real matrices, of EigenSolver, RealSchur's eigenvalues() and
// the like, are std::complex<F>, which is unspecified for F and does not
// compile with MSVC's standard library; see docs/COMPLEX.md.

#ifndef MPFIX_EIGEN_H_
#define MPFIX_EIGEN_H_

#include <Eigen/Core>
#include <limits>

#include "mpfix/charconv.h"
#include "mpfix/complex.h"
#include "mpfix/float.h"
#include "mpfix/math.h"

namespace Eigen {

template <int P, int E>
struct NumTraits<mpfix::Float<P, E>> {
  using Real = mpfix::Float<P, E>;
  using NonInteger = Real;
  using Nested = Real;
  using Literal = Real;

  enum {
    IsComplex = 0,
    IsInteger = 0,
    IsSigned = 1,
    // Float is trivially copyable and, like double, may stay uninitialized.
    RequireInitialization = 0,
    ReadCost = static_cast<int>(Real::kLimbs) + 1,
    AddCost = 4 * static_cast<int>(Real::kLimbs),
    MulCost = 4 * static_cast<int>(Real::kLimbs * Real::kLimbs)
  };

  static Real epsilon() { return std::numeric_limits<Real>::epsilon(); }
  // The tolerance of approximate comparisons such as isApprox. As for double,
  // whose 1e-12 is about 4500 epsilons, it leaves about 12 bits.
  static Real dummy_precision() { return scalbn(epsilon(), 12); }
  static Real highest() { return Real::MaxFinite(); }
  static Real lowest() { return Real::MaxFinite(true); }
  static Real infinity() { return Real::Infinity(); }
  static Real quiet_NaN() { return Real::QuietNaN(); }
  static constexpr int digits() { return P; }
  static constexpr int digits10() {
    return std::numeric_limits<Real>::digits10;
  }
  static constexpr int max_digits10() {
    return std::numeric_limits<Real>::max_digits10;
  }
  static constexpr int min_exponent() {
    return std::numeric_limits<Real>::min_exponent;
  }
  static constexpr int max_exponent() {
    return std::numeric_limits<Real>::max_exponent;
  }
};

template <int P, int E>
struct NumTraits<mpfix::Complex<mpfix::Float<P, E>>> {
  using Real = mpfix::Float<P, E>;
  using NonInteger = mpfix::Complex<Real>;
  using Nested = NonInteger;
  // Real and complex operands mix through Real, as double and
  // std::complex<double> do.
  using Literal = Real;

  enum {
    IsComplex = 1,
    IsInteger = 0,
    IsSigned = 1,
    RequireInitialization = 0,
    ReadCost = 2 * NumTraits<Real>::ReadCost,
    AddCost = 2 * NumTraits<Real>::AddCost,
    MulCost = 4 * NumTraits<Real>::MulCost + 2 * NumTraits<Real>::AddCost
  };

  static Real epsilon() { return NumTraits<Real>::epsilon(); }
  static Real dummy_precision() { return NumTraits<Real>::dummy_precision(); }
  static constexpr int digits() { return P; }
  static constexpr int digits10() { return NumTraits<Real>::digits10(); }
  static constexpr int max_digits10() {
    return NumTraits<Real>::max_digits10();
  }
};


namespace internal {

// Random values with all of Float's fraction bits, not just double's, which is
// what Eigen gives custom floating-point types otherwise. Used by Random().
template <int P, int E>
struct random_float_impl<mpfix::Float<P, E>, false> {
  using F = mpfix::Float<P, E>;
  static int mantissaBits() { return P - 1; }
  // A value in [-1, 1) with numRandomBits random fraction bits.
  static F run(int numRandomBits) {
    eigen_assert(numRandomBits >= 0 && numRandomBits <= mantissaBits());
    typename F::Limbs limbs;
    for (auto& limb : limbs) limb = getRandomBits<std::uint64_t>(64);
    // Keep the top bit and the numRandomBits bits below it, which also
    // clears the spare bits below the precision.
    constexpr int kTotal = static_cast<int>(F::kLimbs) * mpfix::kLimbBits;
    const int clear = kTotal - 1 - numRandomBits;
    mpfix::detail::ShiftRight(limbs, clear);
    mpfix::detail::ShiftLeft(limbs, clear);
    limbs[F::kLimbs - 1] |= mpfix::limb_t{1} << (mpfix::kLimbBits - 1);
    // In [2, 4), then shifted to [-1, 1) exactly.
    return F::FromRaw(false, 1, limbs) - F(3);
  }
};

}  // namespace internal
}  // namespace Eigen

#endif  // MPFIX_EIGEN_H_
