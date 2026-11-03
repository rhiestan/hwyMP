// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Tests of mpfix::Complex against MPC, whose results are correctly rounded:
// the operations documented as correctly rounded must match it exactly, and
// the others within an ulp in each component, which they nearly always match
// exactly too. Set MPFIX_TEST_SCALE to multiply the iteration counts.

#include <mpc.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <string>

#include "mpfix/complex.h"
#include "oracle.h"

namespace mpfix {
namespace {

using namespace test;  // NOLINT

// An mpc_t that owns its storage.
class Mpc {
 public:
  explicit Mpc(long precision) { mpc_init2(value_, precision); }
  ~Mpc() { mpc_clear(value_); }
  Mpc(const Mpc&) = delete;
  Mpc& operator=(const Mpc&) = delete;
  mpc_ptr get() { return value_; }

 private:
  mpc_t value_;
};

template <class F>
void SetMpc(mpc_ptr result, const Complex<F>& z) {
  const View<F> re(z.real());
  const View<F> im(z.imag());
  mpc_set_fr_fr(result, re.get(), im.get(), MPC_RNDNN);
}

// An MPFR value of high precision rounded to F, to nearest. The ternary
// value tells on which side of it the exact value lies: one unit of the high
// precision further that way, the value rounds as the exact one, even where
// it was a tie itself.
template <class F>
F Rounded(mpfr_srcptr value, int ternary) {
  Mpfr nudged(static_cast<long>(mpfr_get_prec(value)));
  mpfr_set(nudged.get(), value, MPFR_RNDN);
  if (ternary > 0) mpfr_nextbelow(nudged.get());
  if (ternary < 0) mpfr_nextabove(nudged.get());
  return RoundExact<F>(
             [&](mpfr_ptr result, mpfr_rnd_t rnd) {
               return mpfr_set(result, nudged.get(), rnd);
             },
             Rounding::kNearestEven)
      .value;
}

template <class F>
bool WithinUlp(const F& actual, const F& expected) {
  if (SameResult(actual, expected)) return true;
  if (actual.is_nan() || expected.is_nan()) return false;
  return next_up(actual) == expected || next_down(actual) == expected;
}

enum class Op {
  kAdd, kSub, kMul, kDiv, kSqrt, kExp, kLog, kLog10, kSin, kCos, kTan,
  kSinh, kCosh, kTanh, kPow
};

const char* Name(Op op) {
  constexpr const char* kNames[] = {"add", "sub", "mul", "div", "sqrt",
                                    "exp", "log", "log10", "sin", "cos",
                                    "tan", "sinh", "cosh", "tanh", "pow"};
  return kNames[static_cast<int>(op)];
}

constexpr Op kUnary[] = {Op::kSqrt, Op::kExp, Op::kLog, Op::kLog10,
                         Op::kSin,  Op::kCos, Op::kTan, Op::kSinh,
                         Op::kCosh, Op::kTanh};

template <class F>
Complex<F> ApplyMpfix(Op op, const Complex<F>& z, const Complex<F>& w) {
  switch (op) {
    case Op::kAdd: return z + w;
    case Op::kSub: return z - w;
    case Op::kMul: return z * w;
    case Op::kDiv: return z / w;
    case Op::kSqrt: return sqrt(z);
    case Op::kExp: return exp(z);
    case Op::kLog: return log(z);
    case Op::kLog10: return log10(z);
    case Op::kSin: return sin(z);
    case Op::kCos: return cos(z);
    case Op::kTan: return tan(z);
    case Op::kSinh: return sinh(z);
    case Op::kCosh: return cosh(z);
    case Op::kTanh: return tanh(z);
    case Op::kPow: return pow(z, w);
  }
  return z;
}

int ApplyMpc(Op op, mpc_ptr result, mpc_srcptr z, mpc_srcptr w) {
  switch (op) {
    case Op::kAdd: return mpc_add(result, z, w, MPC_RNDNN);
    case Op::kSub: return mpc_sub(result, z, w, MPC_RNDNN);
    case Op::kMul: return mpc_mul(result, z, w, MPC_RNDNN);
    case Op::kDiv: return mpc_div(result, z, w, MPC_RNDNN);
    case Op::kSqrt: return mpc_sqrt(result, z, MPC_RNDNN);
    case Op::kExp: return mpc_exp(result, z, MPC_RNDNN);
    case Op::kLog: return mpc_log(result, z, MPC_RNDNN);
    case Op::kLog10: return mpc_log10(result, z, MPC_RNDNN);
    case Op::kSin: return mpc_sin(result, z, MPC_RNDNN);
    case Op::kCos: return mpc_cos(result, z, MPC_RNDNN);
    case Op::kTan: return mpc_tan(result, z, MPC_RNDNN);
    case Op::kSinh: return mpc_sinh(result, z, MPC_RNDNN);
    case Op::kCosh: return mpc_cosh(result, z, MPC_RNDNN);
    case Op::kTanh: return mpc_tanh(result, z, MPC_RNDNN);
    case Op::kPow: return mpc_pow(result, z, w, MPC_RNDNN);
  }
  return 0;
}

struct Stats {
  long checks = 0;
  long off_by_one = 0;
};

// Compares one operation; returns false after a failure.
template <class F>
bool Check(Op op, const Complex<F>& z, const Complex<F>& w, Stats& stats) {
  const long precision = F::kPrecision + 256;
  Mpc mz(precision), mw(precision), mr(precision);
  SetMpc(mz.get(), z);
  SetMpc(mw.get(), w);
  const int inexact = ApplyMpc(op, mr.get(), mz.get(), mw.get());
  const F re = Rounded<F>(mpc_realref(mr.get()), MPC_INEX_RE(inexact));
  const F im = Rounded<F>(mpc_imagref(mr.get()), MPC_INEX_IM(inexact));
  const Complex<F> actual = ApplyMpfix(op, z, w);
  const bool exact_op = op == Op::kAdd || op == Op::kSub || op == Op::kMul;
  bool ok;
  if (op == Op::kPow) {
    // Relative to |z^w| as a whole.
    const F scale = hypot(re, im);
    const F tolerance = scale * scalbn(F(1), 8 - F::kPrecision);
    ok = !(abs(actual.real() - re) > tolerance) &&
         !(abs(actual.imag() - im) > tolerance);
  } else if (exact_op) {
    ok = SameResult(actual.real(), re) && SameResult(actual.imag(), im);
  } else {
    ok = WithinUlp(actual.real(), re) && WithinUlp(actual.imag(), im);
    if (ok && !(SameResult(actual.real(), re) && SameResult(actual.imag(), im))) {
      ++stats.off_by_one;
    }
  }
  ++stats.checks;
  if (!ok) {
    ADD_FAILURE() << Name(op) << " P=" << F::kPrecision << "\n  z = ("
                  << Hex(z.real()) << ", " << Hex(z.imag()) << ")\n  w = ("
                  << Hex(w.real()) << ", " << Hex(w.imag())
                  << ")\n  actual   = (" << Hex(actual.real()) << ", "
                  << Hex(actual.imag()) << ")\n  expected = (" << Hex(re)
                  << ", " << Hex(im) << ")";
  }
  return ok;
}

template <class F>
F RandomPart(Random& random, std::int32_t low, std::int32_t high) {
  if (random.Below(16) == 0) return F::Zero(random.Below(2) != 0);
  low = (std::max)(low, F::kEmin - F::kPrecision + 1);
  high = (std::min)(high, F::kEmax);
  const auto span = static_cast<limb_t>(high - low + 1);
  return MakeRegular<F>(random.Below(2) != 0,
                        low + static_cast<std::int32_t>(random.Below(span)),
                        RandomSignificand<F>(random));
}

template <class F>
Complex<F> RandomComplex(Random& random, std::int32_t low, std::int32_t high) {
  return {RandomPart<F>(random, low, high), RandomPart<F>(random, low, high)};
}

// Near the unit circle, where log |z| cancels.
template <class F>
Complex<F> NearUnitCircle(Random& random) {
  const F angle = RandomPart<F>(random, -4, 1);
  Complex<F> z(cos(angle), sin(angle));
  if (random.Below(2) == 0) z.real(next_up(z.real()));
  return z;
}

template <class F>
void CheckRandomly(int iterations, limb_t seed) {
  Random random(seed);
  Stats stats;
  const int total = iterations * TestScale();
  for (int i = 0; i < total; ++i) {
    // Arithmetic over the whole exponent range, with overflow and underflow.
    const std::int32_t low = F::kEmin - F::kPrecision;
    const std::int32_t high = F::kEmax;
    const bool wide = random.Below(4) == 0;
    const Complex<F> a = wide ? RandomComplex<F>(random, low, high)
                              : RandomComplex<F>(random, -30, 30);
    Complex<F> b = wide ? RandomComplex<F>(random, low, high)
                        : RandomComplex<F>(random, -30, 30);
    if (random.Below(8) == 0) b = Complex<F>(-a.imag(), a.real());
    for (const Op op : {Op::kAdd, Op::kSub, Op::kMul, Op::kDiv}) {
      if (op == Op::kDiv && b.real().is_zero() && b.imag().is_zero()) continue;
      if (!Check(op, a, b, stats)) return;
    }
    // The functions on moderate arguments, and log near |z| = 1. A part next
    // to a multiple of pi/2 makes cos or sin nearly vanish, and with it a
    // component of sin, cos, exp or sinh, or brings tan and tanh near a pole.
    Complex<F> z = RandomComplex<F>(random, -12, 5);
    if (random.Below(4) == 0) {
      const F half_pi = scalbn(pi<F>(), -1);
      // Nonzero multiples: next to zero, MPC works in a precision spanning
      // the exponent range, and takes long.
      const double k = static_cast<double>(random.Below(4) + 1);
      F near = half_pi * F(random.Below(2) == 0 ? k : -k);
      if (random.Below(2) == 0) near = next_up(near);
      // With the other part tiny, tan and tanh come right next to a pole,
      // where cosh 2x + cos 2y cancels.
      const F tiny = random.Below(2) == 0
                         ? RandomPart<F>(random, -F::kPrecision,
                                         -F::kPrecision / 4)
                         : F::Zero();
      if (random.Below(2) == 0) {
        z.real(near);
        if (!tiny.is_zero()) z.imag(tiny);
      } else {
        z.imag(near);
        if (!tiny.is_zero()) z.real(tiny);
      }
    }
    for (const Op op : kUnary) {
      const Complex<F> argument =
          (op == Op::kLog || op == Op::kLog10) && random.Below(3) == 0
              ? NearUnitCircle<F>(random)
              : z;
      if ((op == Op::kLog || op == Op::kLog10) && argument.real().is_zero() &&
          argument.imag().is_zero()) {
        continue;
      }
      if (!Check(op, argument, argument, stats)) return;
    }
    const Complex<F> w = RandomComplex<F>(random, -6, 2);
    if (!(z.real().is_zero() && z.imag().is_zero())) {
      if (!Check(Op::kPow, z, w, stats)) return;
    }
  }
  // Off by one ulp should be rare: the functions are within half an ulp plus
  // about 2^-60.
  EXPECT_LT(stats.off_by_one * 1000, stats.checks) << "P=" << F::kPrecision;
  std::printf("P=%d: %ld checks, %ld off by one ulp\n", F::kPrecision,
              stats.checks, stats.off_by_one);
}

TEST(ComplexOracle, RandomBinary128) { CheckRandomly<binary128>(1000, 1); }
TEST(ComplexOracle, RandomBinary256) { CheckRandomly<binary256>(300, 2); }
TEST(ComplexOracle, RandomBinary512) { CheckRandomly<binary512>(100, 3); }
TEST(ComplexOracle, RandomOtherFormats) {
  CheckRandomly<Float<53, 11>>(2000, 4);
  CheckRandomly<Float<24, 8>>(2000, 5);
  CheckRandomly<Float<11, 5>>(2000, 6);
}

// abs, arg and norm are correctly rounded.
TEST(ComplexOracle, AbsArgNorm) {
  using F = binary128;
  Random random(7);
  const long precision = F::kPrecision + 256;
  Mpfr exact(precision);
  for (int i = 0; i < 2000 * TestScale(); ++i) {
    const Complex<F> z = RandomComplex<F>(random, F::kEmin - F::kPrecision,
                                          F::kEmax);
    Mpc mz(precision);
    SetMpc(mz.get(), z);
    int ternary = mpc_abs(exact.get(), mz.get(), MPFR_RNDN);
    EXPECT_TRUE(SameResult(abs(z), Rounded<F>(exact.get(), ternary)));
    ternary = mpc_arg(exact.get(), mz.get(), MPFR_RNDN);
    EXPECT_TRUE(SameResult(arg(z), Rounded<F>(exact.get(), ternary)));
    ternary = mpc_norm(exact.get(), mz.get(), MPFR_RNDN);
    EXPECT_TRUE(SameResult(norm(z), Rounded<F>(exact.get(), ternary)))
        << Hex(z.real()) << " " << Hex(z.imag());
  }
}

}  // namespace
}  // namespace mpfix
