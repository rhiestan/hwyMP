// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Differential tests of the exponentials, logarithms and hyperbolic functions
// against MPFR, in all five rounding modes and with flags, plus direct checks
// of the error bounds the kernels claim, in every working format. Set
// MPFIX_TEST_SCALE to multiply the iteration counts.

#include <bit>
#include <cstdint>
#include <string>
#include <type_traits>

#include "mpfix/math.h"
#include "oracle.h"

namespace mpfix {
namespace {

using namespace test;  // NOLINT

enum class Function {
  kExp,
  kExp2,
  kExpm1,
  kLog,
  kLog2,
  kLog10,
  kLog1p,
  kSinh,
  kCosh,
  kTanh,
  kAsinh,
  kAcosh,
  kAtanh,
  kExp10,
  kExp2m1,
  kExp10m1,
  kLog2p1,
  kLog10p1
};

constexpr Function kFunctions[] = {
    Function::kExp,   Function::kExp2,  Function::kExpm1, Function::kLog,
    Function::kLog2,  Function::kLog10, Function::kLog1p, Function::kSinh,
    Function::kCosh,  Function::kTanh,  Function::kAsinh, Function::kAcosh,
    Function::kAtanh, Function::kExp10, Function::kExp2m1, Function::kExp10m1,
    Function::kLog2p1, Function::kLog10p1};

const char* Name(Function f) {
  switch (f) {
    case Function::kExp:
      return "exp";
    case Function::kExp2:
      return "exp2";
    case Function::kExpm1:
      return "expm1";
    case Function::kLog:
      return "log";
    case Function::kLog2:
      return "log2";
    case Function::kLog10:
      return "log10";
    case Function::kLog1p:
      return "log1p";
    case Function::kSinh:
      return "sinh";
    case Function::kCosh:
      return "cosh";
    case Function::kTanh:
      return "tanh";
    case Function::kAsinh:
      return "asinh";
    case Function::kAcosh:
      return "acosh";
    case Function::kAtanh:
      return "atanh";
    case Function::kExp10:
      return "exp10";
    case Function::kExp2m1:
      return "exp2m1";
    case Function::kExp10m1:
      return "exp10m1";
    case Function::kLog2p1:
      return "log2p1";
    case Function::kLog10p1:
      return "log10p1";
  }
  return "";
}

int ApplyMpfr(Function f, mpfr_ptr result, mpfr_srcptr x, mpfr_rnd_t rnd) {
  switch (f) {
    case Function::kExp:
      return mpfr_exp(result, x, rnd);
    case Function::kExp2:
      return mpfr_exp2(result, x, rnd);
    case Function::kExpm1:
      return mpfr_expm1(result, x, rnd);
    case Function::kLog:
      return mpfr_log(result, x, rnd);
    case Function::kLog2:
      return mpfr_log2(result, x, rnd);
    case Function::kLog10:
      return mpfr_log10(result, x, rnd);
    case Function::kLog1p:
      return mpfr_log1p(result, x, rnd);
    case Function::kSinh:
      return mpfr_sinh(result, x, rnd);
    case Function::kCosh:
      return mpfr_cosh(result, x, rnd);
    case Function::kTanh:
      return mpfr_tanh(result, x, rnd);
    case Function::kAsinh:
      return mpfr_asinh(result, x, rnd);
    case Function::kAcosh:
      return mpfr_acosh(result, x, rnd);
    case Function::kAtanh:
      return mpfr_atanh(result, x, rnd);
    case Function::kExp10:
      return mpfr_exp10(result, x, rnd);
    case Function::kExp2m1:
      return mpfr_exp2m1(result, x, rnd);
    case Function::kExp10m1:
      return mpfr_exp10m1(result, x, rnd);
    case Function::kLog2p1:
      return mpfr_log2p1(result, x, rnd);
    case Function::kLog10p1:
      return mpfr_log10p1(result, x, rnd);
  }
  return 0;
}

// With kLevel > 1, starts in the second or third working format, which random
// arguments almost never reach otherwise.
template <int kLevel = 1, class F>
F ApplyMpfix(Function f, const F& x, Context& context) {
  if constexpr (kLevel == 1) {
    switch (f) {
      case Function::kExp:
        return exp(x, context);
      case Function::kExp2:
        return exp2(x, context);
      case Function::kExpm1:
        return expm1(x, context);
      case Function::kLog:
        return log(x, context);
      case Function::kLog2:
        return log2(x, context);
      case Function::kLog10:
        return log10(x, context);
      case Function::kLog1p:
        return log1p(x, context);
      case Function::kSinh:
        return sinh(x, context);
      case Function::kCosh:
        return cosh(x, context);
      case Function::kTanh:
        return tanh(x, context);
      case Function::kAsinh:
        return asinh(x, context);
      case Function::kAcosh:
        return acosh(x, context);
      case Function::kAtanh:
        return atanh(x, context);
      case Function::kExp10:
        return exp10(x, context);
      case Function::kExp2m1:
        return exp2m1(x, context);
      case Function::kExp10m1:
        return exp10m1(x, context);
      case Function::kLog2p1:
        return log2p1(x, context);
      case Function::kLog10p1:
        return log10p1(x, context);
    }
  } else {
    using Env = detail::ContextEnv;
    const Env env{context};
    switch (f) {
      case Function::kExp:
        return detail::ExpImpl<F, Env, kLevel>(x, env);
      case Function::kExp2:
        return detail::Exp2Impl<F, Env, kLevel>(x, env);
      case Function::kExpm1:
        return detail::Expm1Impl<F, Env, kLevel>(x, env);
      case Function::kLog:
        return detail::LogImpl<F, Env, kLevel>(x, env);
      case Function::kLog2:
        return detail::Log2Impl<F, Env, kLevel>(x, env);
      case Function::kLog10:
        return detail::Log10Impl<F, Env, kLevel>(x, env);
      case Function::kLog1p:
        return detail::Log1pImpl<F, Env, kLevel>(x, env);
      case Function::kSinh:
        return detail::SinhImpl<F, Env, kLevel>(x, env);
      case Function::kCosh:
        return detail::CoshImpl<F, Env, kLevel>(x, env);
      case Function::kTanh:
        return detail::TanhImpl<F, Env, kLevel>(x, env);
      case Function::kAsinh:
        return detail::AsinhImpl<F, Env, kLevel>(x, env);
      case Function::kAcosh:
        return detail::AcoshImpl<F, Env, kLevel>(x, env);
      case Function::kAtanh:
        return detail::AtanhImpl<F, Env, kLevel>(x, env);
      case Function::kExp10:
        return detail::Exp10Impl<F, Env, kLevel>(x, env);
      case Function::kExp2m1:
        return detail::Exp2m1Impl<F, Env, kLevel>(x, env);
      case Function::kExp10m1:
        return detail::Exp10m1Impl<F, Env, kLevel>(x, env);
      case Function::kLog2p1:
        return detail::LogPlus1Impl<2, F, Env, kLevel>(x, env);
      case Function::kLog10p1:
        return detail::LogPlus1Impl<10, F, Env, kLevel>(x, env);
    }
  }
  return x;
}

template <class F>
F ApplyMpfix(Function f, const F& x) {
  switch (f) {
    case Function::kExp:
      return exp(x);
    case Function::kExp2:
      return exp2(x);
    case Function::kExpm1:
      return expm1(x);
    case Function::kLog:
      return log(x);
    case Function::kLog2:
      return log2(x);
    case Function::kLog10:
      return log10(x);
    case Function::kLog1p:
      return log1p(x);
    case Function::kSinh:
      return sinh(x);
    case Function::kCosh:
      return cosh(x);
    case Function::kTanh:
      return tanh(x);
    case Function::kAsinh:
      return asinh(x);
    case Function::kAcosh:
      return acosh(x);
    case Function::kAtanh:
      return atanh(x);
    case Function::kExp10:
      return exp10(x);
    case Function::kExp2m1:
      return exp2m1(x);
    case Function::kExp10m1:
      return exp10m1(x);
    case Function::kLog2p1:
      return log2p1(x);
    case Function::kLog10p1:
      return log10p1(x);
  }
  return x;
}

template <class F>
Expected<F> ExpectedResult(Function f, const F& x, Rounding mode) {
  if (x.is_nan()) {
    return {F::QuietNaN(), x.is_signaling() ? unsigned{kFlagInvalid} : 0u};
  }
  // The poles, which MPFR reports through its own flags.
  const bool is_log = f == Function::kLog || f == Function::kLog2 ||
                      f == Function::kLog10;
  const bool is_log_plus_1 = f == Function::kLog1p ||
                             f == Function::kLog2p1 ||
                             f == Function::kLog10p1;
  if ((is_log && x.is_zero()) || (is_log_plus_1 && x == F(-1))) {
    return {F::Infinity(true), unsigned{kFlagDivByZero}};
  }
  if (f == Function::kAtanh && abs(x) == F(1)) {
    return {F::Infinity(x.signbit()), unsigned{kFlagDivByZero}};
  }
  const View<F> view(x);
  const auto apply = [&](mpfr_ptr result, mpfr_rnd_t rnd) {
    return ApplyMpfr(f, result, view.get(), rnd);
  };
  // A finite x whose result is beyond even MPFR's range overflows. Rounding
  // away from zero gives infinity there in MPFR.
  Mpfr probe(F::kPrecision);
  apply(probe.get(), MPFR_RNDA);
  if (mpfr_inf_p(probe.get()) && !x.is_inf()) {
    return {OverflowResult<F>(mpfr_signbit(probe.get()) != 0, mode),
            unsigned{kFlagOverflow} | kFlagInexact};
  }
  return RoundExact<F>(apply, mode);
}

template <int kLevel = 1, class F>
bool Check(Function f, const F& x, Rounding mode) {
  const Expected<F> expected = ExpectedResult(f, x, mode);
  Context context{mode};
  const F actual = ApplyMpfix<kLevel>(f, x, context);
  bool ok = SameResult(actual, expected.value) && context.flags == expected.flags;
  if constexpr (kLevel == 1) {
    if (mode == Rounding::kNearestEven) {
      ok = ok && SameResult(ApplyMpfix(f, x), expected.value);
    }
  }
  if (!ok) {
    ADD_FAILURE() << Name(f) << " level " << kLevel << " P=" << F::kPrecision
                  << " E=" << F::kExponentBits
                  << " mode=" << static_cast<int>(mode)
                  << "\n  x        = " << Hex(x)
                  << "\n  actual   = " << Hex(actual) << " flags "
                  << context.flags << "\n  expected = " << Hex(expected.value)
                  << " flags " << expected.flags;
  }
  return ok;
}

template <int kLevel = 1, class F>
bool CheckAllModes(Function f, const F& x) {
  for (const Rounding mode : kModes) {
    if (!Check<kLevel>(f, x, mode)) return false;
  }
  return true;
}

template <class F>
void CheckExhaustively() {
  const limb_t count = limb_t{1} << (F::kPrecision + F::kExponentBits);
  for (const Function f : kFunctions) {
    for (limb_t bits = 0; bits < count; ++bits) {
      if (!CheckAllModes(f, F::from_bits({bits}))) return;
    }
  }
}

TEST(MathOracle, ExhaustiveSmallFormats) {
  CheckExhaustively<Float<4, 3>>();
  CheckExhaustively<Float<8, 4>>();
  CheckExhaustively<Float<8, 5>>();
  CheckExhaustively<Float<11, 5>>();  // binary16
}

// Arguments where the functions do something interesting: around zero, where
// they are linear, out to where they overflow, and near one for logarithms.
template <class F>
F RandomRawArgument(Function f, Random& random) {
  const bool negative = random.Below(2) != 0;
  const auto regular = [&](std::int32_t exponent) {
    return MakeRegular<F>(negative, exponent, RandomSignificand<F>(random));
  };
  switch (random.Below(8)) {
    case 0:
      return RandomValue<F>(random);
    case 1: {
      // Next to one, or next to minus one for log1p, log2p1 and log10p1.
      const bool plus_1 = f == Function::kLog1p || f == Function::kLog2p1 ||
                          f == Function::kLog10p1;
      const F one = plus_1 && negative ? F(-1) : F(1);
      const auto offset = static_cast<std::int32_t>(
          random.Below(static_cast<limb_t>(F::kPrecision) + 8));
      return one + MakeRegular<F>(random.Below(2) != 0, -offset,
                                  RandomSignificand<F>(random));
    }
    case 2:
      // Around the thresholds for tiny arguments.
      return regular(-F::kPrecision - 6 +
                     static_cast<std::int32_t>(random.Below(8)));
    case 3:
      // Near overflow and underflow of the exponentials.
      return regular(F::kExponentBits - 3 +
                     static_cast<std::int32_t>(random.Below(4)));
    case 4:
      // Around ln 2 / 2, where the hyperbolic functions change methods, and
      // around P/2, where tanh reaches one.
      return random.Below(2) == 0
                 ? regular(-3 + static_cast<std::int32_t>(random.Below(3)))
                 : regular(std::bit_width(static_cast<unsigned>(F::kPrecision)) -
                           3 + static_cast<std::int32_t>(random.Below(3)));
    case 5: {
      // Small integers and halves, where exp10, exp2m1 and exp10m1 may be
      // exact or tied.
      const auto span = static_cast<limb_t>(4 * F::kPrecision);
      return scalbn(F(static_cast<double>(random.Below(span)) -
                      static_cast<double>(span) / 2),
                    -static_cast<int>(random.Below(2)));
    }
    default: {
      const auto span =
          static_cast<limb_t>(F::kPrecision + F::kExponentBits + 12);
      return regular(static_cast<std::int32_t>(random.Below(span)) -
                     (F::kPrecision + 8));
    }
  }
}

// Mostly within the domains of acosh and atanh, keeping arguments near one
// there.
template <class F>
F RandomArgument(Function f, Random& random) {
  F x = RandomRawArgument<F>(f, random);
  if (random.Below(4) != 0) {
    if (f == Function::kAcosh && x < F(1)) x = F(2) - x;
    if (f == Function::kAtanh && abs(x) > F(1)) x = F(1) / x;
  }
  return x;
}

template <class F, int kLevel = 1>
void CheckRandomly(int iterations, limb_t seed) {
  Random random(seed);
  const int total = iterations * TestScale();
  for (const Function f : kFunctions) {
    for (int i = 0; i < total; ++i) {
      if (!CheckAllModes<kLevel>(f, RandomArgument<F>(f, random))) return;
    }
  }
}

TEST(MathOracle, RandomBinary128) { CheckRandomly<binary128>(3000, 1); }
TEST(MathOracle, RandomBinary256) { CheckRandomly<binary256>(1000, 2); }
TEST(MathOracle, RandomBinary512) { CheckRandomly<binary512>(300, 3); }
TEST(MathOracle, RandomOtherFormats) {
  CheckRandomly<Float<53, 11>>(2000, 4);
  CheckRandomly<Float<64, 15>>(2000, 5);
  CheckRandomly<Float<200, 20>>(500, 6);
}

TEST(MathOracle, LaterLevels) {
  CheckRandomly<Float<11, 5>, 2>(2000, 7);
  CheckRandomly<Float<11, 5>, 4>(2000, 8);
  CheckRandomly<binary128, 2>(300, 9);
  CheckRandomly<binary128, 4>(300, 10);
  CheckRandomly<binary512, 2>(30, 11);
  CheckRandomly<binary512, 4>(10, 12);
}

// Integers, where exp10, exp2m1 and exp10m1 are exact, tied or just beside a
// rounding boundary, and the arguments where log2p1 and log10p1 are exact or
// just above an integer.
template <class F>
void CheckIntegerCases() {
  for (const Rounding mode : kModes) {
    for (int n = -(F::kPrecision + 8); n <= F::kPrecision + 8; ++n) {
      for (const Function f :
           {Function::kExp10, Function::kExp2m1, Function::kExp10m1}) {
        EXPECT_TRUE(Check(f, F(n), mode)) << Name(f) << " " << n;
      }
      // 2^n - 1, and the powers of two around 2^(P + 4).
      const F power = scalbn(F(1), n);
      EXPECT_TRUE(Check(Function::kLog2p1, power - F(1), mode)) << n;
      EXPECT_TRUE(Check(Function::kLog2p1, scalbn(power, F::kPrecision), mode))
          << n;
    }
    F ten(1);
    for (int n = 0; n < F::kPrecision / 3 + 4; ++n) {
      EXPECT_TRUE(Check(Function::kLog10p1, ten - F(1), mode)) << n;
      EXPECT_TRUE(Check(Function::kLog10p1, ten, mode)) << n;
      EXPECT_TRUE(Check(Function::kLog10p1, next_up(ten - F(1)), mode)) << n;
      ten = ten * F(10);
    }
    for (int n = 1; n < 4; ++n) {
      EXPECT_TRUE(Check(Function::kExp10, F(n * 1000), mode)) << n;
      EXPECT_TRUE(Check(Function::kExp10m1, F(n * 1000 + 1), mode)) << n;
    }
  }
}

TEST(MathOracle, IntegerCases) {
  CheckIntegerCases<Float<11, 5>>();
  CheckIntegerCases<Float<53, 11>>();
  CheckIntegerCases<binary128>();
  CheckIntegerCases<binary256>();
}

TEST(MathOracle, PowersAreExact) {
  for (const Rounding mode : kModes) {
    for (int k = -40; k <= 40; ++k) {
      EXPECT_TRUE(Check(Function::kExp2, binary256(k), mode)) << k;
      EXPECT_TRUE(Check(Function::kLog2, scalbn(binary256(1), k), mode)) << k;
    }
    binary256 power(1);
    for (int k = 0; k <= 102; ++k) {  // 10^102 < 2^339 still fits in 237 bits.
      EXPECT_TRUE(Check(Function::kLog10, power, mode)) << k;
      EXPECT_TRUE(Check(Function::kLog10, next_up(power), mode)) << k;
      power = power * binary256(10);
    }
  }
}

// --- The kernels' error bounds ---

template <class W>
void CheckKernelBounds(int iterations, limb_t seed) {
  Random random(seed);
  const long precision = W::kPrecision + 128;
  Mpfr exact(precision);
  Mpfr error(precision);
  Mpfr bound(precision);
  for (int i = 0; i < iterations * TestScale(); ++i) {
    // |r| < 1/4 suits both kernels' ranges.
    const bool negative = random.Below(2) != 0;
    const std::int32_t exponent =
        -3 - static_cast<std::int32_t>(random.Below(
                 random.Below(4) == 0 ? W::kPrecision : 24));
    const W r = MakeRegular<W>(negative, exponent, RandomSignificand<W>(random));
    const View<W> view(r);
    int units = 0;

    const W q = detail::ExpM1Kernel(r, units);
    mpfr_expm1(exact.get(), view.get(), MPFR_RNDN);
    const View<W> q_view(q);
    mpfr_sub(error.get(), q_view.get(), exact.get(), MPFR_RNDN);
    mpfr_abs(error.get(), error.get(), MPFR_RNDN);
    mpfr_mul_2si(bound.get(), exact.get(), -W::kPrecision, MPFR_RNDN);
    mpfr_mul_si(bound.get(), bound.get(), units, MPFR_RNDN);
    mpfr_abs(bound.get(), bound.get(), MPFR_RNDN);
    EXPECT_LE(mpfr_cmp(error.get(), bound.get()), 0)
        << "expm1 P=" << W::kPrecision << " r=" << Hex(r);

    const W l = detail::Log1pKernel(r, units);
    mpfr_log1p(exact.get(), view.get(), MPFR_RNDN);
    const View<W> l_view(l);
    mpfr_sub(error.get(), l_view.get(), exact.get(), MPFR_RNDN);
    mpfr_abs(error.get(), error.get(), MPFR_RNDN);
    mpfr_mul_2si(bound.get(), exact.get(), -W::kPrecision, MPFR_RNDN);
    mpfr_mul_si(bound.get(), bound.get(), units, MPFR_RNDN);
    mpfr_abs(bound.get(), bound.get(), MPFR_RNDN);
    EXPECT_LE(mpfr_cmp(error.get(), bound.get()), 0)
        << "log1p P=" << W::kPrecision << " d=" << Hex(r);
  }
}

// log(m) through the tables of the first format, for m in [0.705, 1.41) with
// at most F's precision.
template <class F>
void CheckTableLogBounds(int iterations, limb_t seed) {
  using W = detail::MathWorking<F, 1>;
  Random random(seed);
  const long precision = W::kPrecision + 128;
  Mpfr exact(precision);
  Mpfr error(precision);
  Mpfr bound(precision);
  for (int i = 0; i < iterations * TestScale(); ++i) {
    // 1 + d, with d spread over magnitudes down to 2^-(P + 8).
    const auto exponent = -2 - static_cast<std::int32_t>(random.Below(
                                   random.Below(2) == 0 ? 20 : F::kPrecision + 8));
    const F d = MakeRegular<F>(random.Below(2) != 0, exponent,
                               RandomSignificand<F>(random));
    W m = detail::FloatAccess::Convert<W>(F(1) + d, detail::NearestEnv{});
    if (m >= W(1.41)) m = scalbn(m, -1);
    if (m < W(0.705)) m = scalbn(m, 1);
    const View<W> view(m);
    int units = 0;
    const W l = detail::LogOfReduced<W, true>(m, units);
    mpfr_log(exact.get(), view.get(), MPFR_RNDN);
    const View<W> l_view(l);
    mpfr_sub(error.get(), l_view.get(), exact.get(), MPFR_RNDN);
    mpfr_abs(error.get(), error.get(), MPFR_RNDN);
    mpfr_mul_2si(bound.get(), exact.get(), -W::kPrecision, MPFR_RNDN);
    mpfr_mul_si(bound.get(), bound.get(), units, MPFR_RNDN);
    mpfr_abs(bound.get(), bound.get(), MPFR_RNDN);
    EXPECT_LE(mpfr_cmp(error.get(), bound.get()), 0)
        << "table log P=" << W::kPrecision << " m=" << Hex(m);
  }
}

// e^r - 1 through the tables of the first format, for |r| <= 0.35.
template <class F>
void CheckTableExpBounds(int iterations, limb_t seed) {
  using W = detail::MathWorking<F, 1>;
  Random random(seed);
  const long precision = W::kPrecision + 128;
  Mpfr exact(precision);
  Mpfr error(precision);
  Mpfr bound(precision);
  for (int i = 0; i < iterations * TestScale(); ++i) {
    const auto exponent = -3 - static_cast<std::int32_t>(random.Below(
                                   random.Below(2) == 0 ? 20 : W::kPrecision));
    W r = MakeRegular<W>(random.Below(2) != 0, exponent,
                         RandomSignificand<W>(random));
    if (random.Below(2) == 0) r = scalbn(r, 1);  // Up to 0.25 * 1.4.
    if (abs(r) > W(0.35)) r = scalbn(r, -1);
    const View<W> view(r);
    int units = 0;
    const W q = detail::ExpM1Reduced<W, true>(r, units);
    mpfr_expm1(exact.get(), view.get(), MPFR_RNDN);
    const View<W> q_view(q);
    mpfr_sub(error.get(), q_view.get(), exact.get(), MPFR_RNDN);
    mpfr_abs(error.get(), error.get(), MPFR_RNDN);
    mpfr_mul_2si(bound.get(), exact.get(), -W::kPrecision, MPFR_RNDN);
    mpfr_mul_si(bound.get(), bound.get(), units, MPFR_RNDN);
    mpfr_abs(bound.get(), bound.get(), MPFR_RNDN);
    EXPECT_LE(mpfr_cmp(error.get(), bound.get()), 0)
        << "table expm1 P=" << W::kPrecision << " r=" << Hex(r);
  }
}

// The hyperbolic functions and their inverses in the working format W of F,
// with or without tables, against their claimed bounds.
template <class F, class W>
void CheckHyperbolicBounds(int iterations, limb_t seed) {
  constexpr bool kTables = std::is_same_v<W, detail::MathWorking<F, 1>>;
  Random random(seed);
  const long precision = W::kPrecision + 128;
  Mpfr exact(precision);
  Mpfr error(precision);
  Mpfr bound(precision);
  const auto check = [&](const char* name, const W& a, const W& y,
                         std::int64_t units,
                         int (*reference)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t)) {
    const View<W> view(a);
    reference(exact.get(), view.get(), MPFR_RNDN);
    const View<W> y_view(y);
    mpfr_sub(error.get(), y_view.get(), exact.get(), MPFR_RNDN);
    mpfr_abs(error.get(), error.get(), MPFR_RNDN);
    mpfr_mul_2si(bound.get(), exact.get(), -W::kPrecision, MPFR_RNDN);
    mpfr_mul_d(bound.get(), bound.get(), static_cast<double>(units),
               MPFR_RNDN);
    mpfr_abs(bound.get(), bound.get(), MPFR_RNDN);
    EXPECT_LE(mpfr_cmp(error.get(), bound.get()), 0)
        << name << " P=" << W::kPrecision << " tables=" << kTables
        << " a=" << Hex(a);
  };
  for (int i = 0; i < iterations * TestScale(); ++i) {
    // From 2^-(P/2) to 2^14, mostly around the change of method at 0.35,
    // with F's precision.
    const auto exponent =
        random.Below(2) == 0
            ? -3 + static_cast<std::int32_t>(random.Below(3))
            : static_cast<std::int32_t>(random.Below(
                  static_cast<limb_t>(F::kPrecision / 2 + 14))) -
                  F::kPrecision / 2;
    const F value = MakeRegular<F>(false, (std::min)(exponent, F::kEmax - 1),
                                   RandomSignificand<F>(random));
    const auto wide = [](const F& v) {
      return detail::FloatAccess::Convert<W>(v, detail::NearestEnv{});
    };
    const W a = wide(value);
    int units = 0;
    W y = detail::SinhKernel<W, kTables>(a, units);
    check("sinh", a, y, units, mpfr_sinh);
    y = detail::CoshKernel<W, kTables>(a, units);
    check("cosh", a, y, units, mpfr_cosh);
    if (a < W(F::kPrecision / 2 + 2)) {
      y = detail::TanhKernel<W, kTables>(a, units);
      check("tanh", a, y, units, mpfr_tanh);
    }

    std::int64_t wide_units = 0;
    y = detail::AsinhKernel<W, kTables>(a, wide_units);
    check("asinh", a, y, wide_units, mpfr_asinh);
    // Above one for acosh, below for atanh, often next to one.
    const F above = value < F(1) ? F(1) + value : value;
    if (above != F(1)) {
      y = detail::AcoshKernel<W, kTables>(wide(above), wide_units);
      check("acosh", wide(above), y, wide_units, mpfr_acosh);
    }
    const F below = value < F(1) ? value : F(1) - F(1) / value;
    if (below != F(0) && below < F(1)) {
      y = detail::AtanhKernel<W, kTables>(wide(below), wide_units);
      check("atanh", wide(below), y, wide_units, mpfr_atanh);
    }
  }
}

template <class F>
void CheckKernelBoundsForAllLevels(int iterations, limb_t seed) {
  CheckTableLogBounds<F>(iterations, seed + 3);
  CheckTableExpBounds<F>(iterations, seed + 4);
  CheckHyperbolicBounds<F, detail::MathWorking<F, 1>>(iterations, seed + 5);
  CheckHyperbolicBounds<F, detail::MathWorking<F, 2>>(iterations, seed + 6);
  CheckHyperbolicBounds<F, detail::MathWorking<F, 4>>(iterations / 4 + 1,
                                                       seed + 7);
  CheckKernelBounds<detail::MathWorking<F, 1>>(iterations, seed);
  CheckKernelBounds<detail::MathWorking<F, 2>>(iterations, seed + 1);
  CheckKernelBounds<detail::MathWorking<F, 4>>(iterations / 4 + 1, seed + 2);
}

TEST(MathKernels, StayWithinTheirErrorBounds) {
  CheckKernelBoundsForAllLevels<Float<11, 5>>(2000, 10);
  CheckKernelBoundsForAllLevels<binary128>(2000, 20);
  CheckKernelBoundsForAllLevels<binary256>(500, 30);
  CheckKernelBoundsForAllLevels<binary512>(100, 40);
}

}  // namespace
}  // namespace mpfix
