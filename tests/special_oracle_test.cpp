// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Differential tests of erf, erfc, tgamma and lgamma against MPFR, in all five
// rounding modes and with flags. Set MPFIX_TEST_SCALE to multiply the
// iteration counts.

#include <algorithm>
#include <cstdint>
#include <string>
#include <type_traits>

#include "mpfix/math.h"
#include "oracle.h"

namespace mpfix {
namespace {

using namespace test;  // NOLINT

enum class Function { kErf, kErfc, kTgamma, kLgamma };

constexpr Function kFunctions[] = {Function::kErf, Function::kErfc,
                                   Function::kTgamma, Function::kLgamma};

const char* Name(Function f) {
  switch (f) {
    case Function::kErf:
      return "erf";
    case Function::kErfc:
      return "erfc";
    case Function::kTgamma:
      return "tgamma";
    case Function::kLgamma:
      return "lgamma";
  }
  return "";
}

template <int kLevel = 1, class F>
F ApplyMpfix(Function f, const F& x, Context& context) {
  if constexpr (kLevel == 1) {
    switch (f) {
      case Function::kErf:
        return erf(x, context);
      case Function::kErfc:
        return erfc(x, context);
      case Function::kTgamma:
        return tgamma(x, context);
      case Function::kLgamma:
        return lgamma(x, context);
    }
  } else {
    using Env = detail::ContextEnv;
    const Env env{context};
    switch (f) {
      case Function::kErf:
        return detail::ErfImpl<F, Env, kLevel>(x, env);
      case Function::kErfc:
        return detail::ErfcImpl<F, Env, kLevel>(x, env);
      case Function::kTgamma:
        return detail::TgammaImpl<F, Env, kLevel>(x, env);
      case Function::kLgamma:
        return detail::LgammaImpl<F, Env, kLevel>(x, env);
    }
  }
  return x;
}

template <class F>
F ApplyMpfix(Function f, const F& x) {
  switch (f) {
    case Function::kErf:
      return erf(x);
    case Function::kErfc:
      return erfc(x);
    case Function::kTgamma:
      return tgamma(x);
    case Function::kLgamma:
      return lgamma(x);
  }
  return x;
}

template <class F>
Expected<F> ExpectedResult(Function f, const F& x, Rounding mode) {
  if (x.is_nan()) {
    return {F::QuietNaN(), x.is_signaling() ? unsigned{kFlagInvalid} : 0u};
  }
  // The poles, which MPFR reports through its own flags.
  const bool pole = x.is_zero() || (x.signbit() && !x.is_inf() &&
                                    detail::IsInteger(x));
  if (f == Function::kTgamma && x.is_zero()) {
    return {F::Infinity(x.signbit()), unsigned{kFlagDivByZero}};
  }
  if (f == Function::kLgamma && pole) {
    return {F::Infinity(), unsigned{kFlagDivByZero}};
  }
  const View<F> view(x);
  const auto apply = [&](mpfr_ptr result, mpfr_rnd_t rnd) {
    switch (f) {
      case Function::kErf:
        return mpfr_erf(result, view.get(), rnd);
      case Function::kErfc:
        return mpfr_erfc(result, view.get(), rnd);
      case Function::kTgamma:
        return mpfr_gamma(result, view.get(), rnd);
      case Function::kLgamma: {
        int sign;
        return mpfr_lgamma(result, &sign, view.get(), rnd);
      }
    }
    return 0;
  };
  // A finite x whose result is beyond even MPFR's range overflows.
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
  bool ok =
      SameResult(actual, expected.value) && context.flags == expected.flags;
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

TEST(SpecialOracle, ExhaustiveSmallFormats) {
  CheckExhaustively<Float<4, 3>>();
  CheckExhaustively<Float<8, 4>>();
  CheckExhaustively<Float<8, 5>>();
  CheckExhaustively<Float<11, 5>>();  // binary16
}

template <class F>
F RandomRegular(Random& random, std::int32_t exponent, bool negative) {
  return MakeRegular<F>(negative,
                        (std::max)((std::min)(exponent, F::kEmax),
                                   F::kEmin - F::kPrecision + 1),
                        RandomSignificand<F>(random));
}

// x^2 = target ln 2, rounded to F, for the thresholds of erf and erfc.
template <class F>
F SquareRootOfLog2(double target) {
  return F(std::sqrt(target * 0.69314718055994531));
}

// Arguments where the functions do something interesting: around zero, near
// the thresholds of erf and erfc, near integers and half-integers, near the
// zeros of lgamma, and up to overflow.
template <class F>
F RandomArgument(Function f, Random& random) {
  const bool negative = random.Below(2) != 0;
  const auto near = [&](const F& v) {
    // v, or a value next to it in relative distance 2^-(1 .. P + 8); zero has
    // no such neighbours.
    if (random.Below(4) == 0 || v.is_zero()) return v;
    const auto offset = static_cast<std::int32_t>(
        random.Below(static_cast<limb_t>(F::kPrecision) + 8));
    return v + RandomRegular<F>(random, v.raw_exponent() - 1 - offset,
                                random.Below(2) != 0);
  };
  const bool error_function = f == Function::kErf || f == Function::kErfc;
  switch (random.Below(8)) {
    case 0: {
      // Random values, with exponents below 16.
      F x = RandomValue<F>(random);
      if (x.is_finite() && x.raw_exponent() > 16) {
        x = scalbn(x, 16 - x.raw_exponent());
      }
      return x;
    }
    case 1:
      // Tiny values.
      return RandomRegular<F>(
          random,
          -static_cast<std::int32_t>(random.Below(
              static_cast<limb_t>(3 * F::kPrecision) + 20)),
          negative);
    case 2:
      if (error_function) {
        // Near where erf reaches 1, erfc switches to its asymptotic series,
        // and erfc underflows.
        const double targets[] = {
            F::kPrecision + 2.0, F::kPrecision + 72.0, 2.0 * F::kPrecision,
            F::kPrecision + 3.0 - F::kEmin};
        const double target = targets[random.Below(4)];
        const F x = near(SquareRootOfLog2<F>(target));
        return negative ? -x : x;
      }
      // Near integers and half-integers, positive and negative.
      return near(scalbn(F(static_cast<double>(random.Below(400)) - 200),
                         -static_cast<int>(random.Below(2))));
    case 3:
      if (error_function) {
        return RandomRegular<F>(
            random, static_cast<std::int32_t>(random.Below(6)), negative);
      }
      // Near the zeros of lgamma at 1, 2, -2.457 and -2.748, and the
      // overflow of tgamma.
      switch (random.Below(5)) {
        case 0:
          return near(F(1));
        case 1:
          return near(F(2));
        case 2:
          return near(F(-2.4570247382208005));
        case 3:
          return near(F(-2.7476826467274127));
        default:
          return near(F(detail::GammaOverflowThreshold<F>()));
      }
    default: {
      const std::int32_t top = error_function ? 6 : 10;
      return RandomRegular<F>(
          random,
          static_cast<std::int32_t>(random.Below(
              static_cast<limb_t>(top + 12))) - 12,
          negative);
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

TEST(SpecialOracle, RandomBinary128) { CheckRandomly<binary128>(1500, 1); }
TEST(SpecialOracle, RandomBinary256) { CheckRandomly<binary256>(400, 2); }
TEST(SpecialOracle, RandomBinary512) { CheckRandomly<binary512>(100, 3); }
TEST(SpecialOracle, RandomOtherFormats) {
  CheckRandomly<Float<53, 11>>(2000, 4);
  CheckRandomly<Float<64, 15>>(1000, 5);
  CheckRandomly<Float<200, 20>>(200, 6);
}

TEST(SpecialOracle, LaterLevels) {
  CheckRandomly<Float<11, 5>, 2>(1000, 7);
  CheckRandomly<Float<11, 5>, 4>(1000, 8);
  CheckRandomly<binary128, 2>(100, 9);
  CheckRandomly<binary128, 4>(50, 10);
  CheckRandomly<binary512, 2>(10, 11);
}

TEST(SpecialOracle, SpecialValues) {
  using F = binary128;
  for (const Rounding mode : kModes) {
    for (const F& v : {F::Zero(), F::Zero(true), F::Infinity(),
                       F::Infinity(true), F(1), F(2), F(-1), F(-2), F(0.5),
                       F(-0.5), F(3000), F(3001), F(1755), F(1756)}) {
      for (const Function f : kFunctions) EXPECT_TRUE(Check(f, v, mode));
    }
    // Gamma of tiny arguments, just beside 1/x.
    for (int e : {-5000, -1000, -300, -240, -230, -229, -120}) {
      for (const F& v : {scalbn(F(1), e), scalbn(F(3), e), scalbn(F(-1), e),
                         scalbn(F(-3), e)}) {
        EXPECT_TRUE(Check(Function::kTgamma, v, mode)) << e;
        EXPECT_TRUE(Check(Function::kLgamma, v, mode)) << e;
      }
    }
  }
}

// --- The kernels against their error bounds ---

template <class W>
void CheckBound(const char* name, const detail::Tracked<W>& t, mpfr_ptr exact,
                const W& argument) {
  Mpfr error(W::kPrecision + 128);
  Mpfr bound(W::kPrecision + 128);
  const View<W> view(t.value);
  mpfr_sub(error.get(), view.get(), exact, MPFR_RNDN);
  mpfr_abs(error.get(), error.get(), MPFR_RNDN);
  mpfr_mul_2si(bound.get(), exact, -W::kPrecision, MPFR_RNDN);
  mpfr_mul_d(bound.get(), bound.get(), t.units, MPFR_RNDN);
  mpfr_abs(bound.get(), bound.get(), MPFR_RNDN);
  if (mpfr_cmp(error.get(), bound.get()) > 0) {
    mpfr_div(error.get(), error.get(), bound.get(), MPFR_RNDN);
    ADD_FAILURE() << name << " P=" << W::kPrecision << " x=" << Hex(argument)
                  << " error/bound=" << mpfr_get_d(error.get(), MPFR_RNDN);
  }
}

template <class F, class W>
void CheckSpecialBounds(int iterations, limb_t seed) {
  constexpr bool kFirst = std::is_same_v<W, detail::MathWorking<F, 1>>;
  Random random(seed);
  Mpfr exact(W::kPrecision + 128);
  for (int i = 0; i < iterations * TestScale(); ++i) {
    // erf and erfc of x in (0, 40).
    const W x = MakeRegular<W>(
        false, static_cast<std::int32_t>(random.Below(12)) - 7,
        RandomSignificand<W>(random));
    const View<W> view(x);
    mpfr_erf(exact.get(), view.get(), MPFR_RNDN);
    CheckBound("erf series", detail::ErfPositive<W, kFirst>(x), exact.get(), x);
    mpfr_erfc(exact.get(), view.get(), MPFR_RNDN);
    if (detail::ErfcIsAsymptotic(x)) {
      CheckBound("erfc asymptotic", detail::ErfcAsymptotic<W, kFirst>(x),
                 exact.get(), x);
    } else if constexpr (kFirst) {
      CheckBound("erfc table", detail::ErfcFromTable(x), exact.get(), x);
    }

    // lgamma and gamma of x in (-200, 200), often next to 1, 2 and the
    // negative integers.
    W g = MakeRegular<W>(random.Below(3) == 0,
                         static_cast<std::int32_t>(random.Below(14)) - 6,
                         RandomSignificand<W>(random));
    if (random.Below(3) == 0) {
      g = W(static_cast<double>(random.Below(40)) - 20) +
          scalbn(g, -static_cast<int>(random.Below(60)));
    }
    if (detail::IsInteger(g)) continue;
    const View<W> g_view(g);
    int sign;
    mpfr_lgamma(exact.get(), &sign, g_view.get(), MPFR_RNDN);
    if (!mpfr_zero_p(exact.get())) {
      CheckBound("lgamma", detail::LgammaAny<W, kFirst>(g), exact.get(), g);
    }
    if (!g.signbit() && g < W(1000)) {
      mpfr_gamma(exact.get(), g_view.get(), MPFR_RNDN);
      CheckBound("gamma", detail::GammaPositive<W, kFirst>(
                              detail::Tracked<W>{g, 0}),
                 exact.get(), g);
    }
  }
}

TEST(SpecialKernels, StayWithinTheirErrorBounds) {
  CheckSpecialBounds<Float<11, 5>, detail::MathWorking<Float<11, 5>, 1>>(500,
                                                                         20);
  CheckSpecialBounds<binary128, detail::MathWorking<binary128, 1>>(300, 21);
  CheckSpecialBounds<binary128, detail::MathWorking<binary128, 2>>(100, 22);
  CheckSpecialBounds<binary512, detail::MathWorking<binary512, 1>>(50, 23);
}

}  // namespace
}  // namespace mpfix
