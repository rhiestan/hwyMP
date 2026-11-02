// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Differential tests of the exponentials and logarithms against MPFR, in all
// five rounding modes and with flags, plus direct checks of the error bounds
// the kernels claim, in every working format. Set MPFIX_TEST_SCALE to
// multiply the iteration counts.

#include <cstdint>
#include <string>

#include "mpfix/math.h"
#include "oracle.h"

namespace mpfix {
namespace {

using namespace test;  // NOLINT

enum class Function { kExp, kExp2, kExpm1, kLog, kLog2, kLog10, kLog1p };

constexpr Function kFunctions[] = {Function::kExp,   Function::kExp2,
                                   Function::kExpm1, Function::kLog,
                                   Function::kLog2,  Function::kLog10,
                                   Function::kLog1p};

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
  if ((is_log && x.is_zero()) ||
      (f == Function::kLog1p && x == F(-1))) {
    return {F::Infinity(true), unsigned{kFlagDivByZero}};
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
    return {OverflowResult<F>(false, mode),
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
F RandomArgument(Function f, Random& random) {
  const bool negative = random.Below(2) != 0;
  const auto regular = [&](std::int32_t exponent) {
    return MakeRegular<F>(negative, exponent, RandomSignificand<F>(random));
  };
  switch (random.Below(8)) {
    case 0:
      return RandomValue<F>(random);
    case 1: {
      // Next to one, or next to minus one for log1p.
      const F one = f == Function::kLog1p && negative ? F(-1) : F(1);
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
    default: {
      const auto span =
          static_cast<limb_t>(F::kPrecision + F::kExponentBits + 12);
      return regular(static_cast<std::int32_t>(random.Below(span)) -
                     (F::kPrecision + 8));
    }
  }
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

template <class F>
void CheckKernelBoundsForAllLevels(int iterations, limb_t seed) {
  CheckTableLogBounds<F>(iterations, seed + 3);
  CheckTableExpBounds<F>(iterations, seed + 4);
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
