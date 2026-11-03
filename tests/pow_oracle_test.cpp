// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Differential tests of pow, cbrt and hypot against MPFR, in all five rounding
// modes and with flags. Set MPFIX_TEST_SCALE to multiply the iteration counts.

#include <cstdint>

#include "mpfix/math.h"
#include "oracle.h"

namespace mpfix {
namespace {

using namespace test;  // NOLINT

enum class Function { kPow, kHypot, kCbrt };

const char* Name(Function f) {
  switch (f) {
    case Function::kPow:
      return "pow";
    case Function::kHypot:
      return "hypot";
    case Function::kCbrt:
      return "cbrt";
  }
  return "";
}

template <int kLevel = 1, class F>
F ApplyMpfix(Function f, const F& x, const F& y, Context& context) {
  using Env = detail::ContextEnv;
  const Env env{context};
  switch (f) {
    case Function::kPow:
      if constexpr (kLevel == 1) {
        return pow(x, y, context);
      } else {
        return detail::PowImpl<F, Env, kLevel>(x, y, env);
      }
    case Function::kHypot:
      return hypot(x, y, context);
    case Function::kCbrt:
      if constexpr (kLevel == 1) {
        return cbrt(x, context);
      } else {
        return detail::CbrtImpl<F, Env, kLevel>(x, env);
      }
  }
  return x;
}

template <class F>
F ApplyMpfix(Function f, const F& x, const F& y) {
  switch (f) {
    case Function::kPow:
      return pow(x, y);
    case Function::kHypot:
      return hypot(x, y);
    case Function::kCbrt:
      return cbrt(x);
  }
  return x;
}

template <class F>
Expected<F> ExpectedResult(Function f, const F& x, const F& y, Rounding mode) {
  const bool unary = f == Function::kCbrt;
  if (x.is_signaling() || (!unary && y.is_signaling())) {
    return {F::QuietNaN(), unsigned{kFlagInvalid}};
  }
  if (f == Function::kPow) {
    if (y.is_zero() || x == F(1)) return {F(1), 0u};
    if (x.is_zero() && y.signbit() && !y.is_nan()) {
      // The pole, which MPFR reports through its own flag. pow(0, -inf) is
      // +inf without it.
      const bool odd = !y.is_inf() && detail::LowestBitExponent(y) == 0;
      return {F::Infinity(odd && x.signbit()),
              y.is_inf() ? 0u : unsigned{kFlagDivByZero}};
    }
  }
  if (f == Function::kHypot && (x.is_inf() || y.is_inf())) {
    return {F::Infinity(), 0u};
  }
  if (x.is_nan() || (!unary && y.is_nan())) return {F::QuietNaN(), 0u};
  const View<F> vx(x);
  const View<F> vy(y);
  const auto apply = [&](mpfr_ptr result, mpfr_rnd_t rnd) {
    switch (f) {
      case Function::kPow:
        return mpfr_pow(result, vx.get(), vy.get(), rnd);
      case Function::kHypot:
        return mpfr_hypot(result, vx.get(), vy.get(), rnd);
      case Function::kCbrt:
        return mpfr_cbrt(result, vx.get(), rnd);
    }
    return 0;
  };
  // Finite operands whose result is beyond even MPFR's range overflow.
  Mpfr probe(F::kPrecision);
  apply(probe.get(), MPFR_RNDA);
  if (mpfr_inf_p(probe.get()) && !x.is_inf() && (unary || !y.is_inf())) {
    return {OverflowResult<F>(mpfr_signbit(probe.get()) != 0, mode),
            unsigned{kFlagOverflow} | kFlagInexact};
  }
  return RoundExact<F>(apply, mode);
}

template <int kLevel = 1, class F>
bool Check(Function f, const F& x, const F& y, Rounding mode) {
  const Expected<F> expected = ExpectedResult(f, x, y, mode);
  Context context{mode};
  const F actual = ApplyMpfix<kLevel>(f, x, y, context);
  bool ok =
      SameResult(actual, expected.value) && context.flags == expected.flags;
  if constexpr (kLevel == 1) {
    if (mode == Rounding::kNearestEven) {
      ok = ok && SameResult(ApplyMpfix(f, x, y), expected.value);
    }
  }
  if (!ok) {
    ADD_FAILURE() << Name(f) << " level " << kLevel << " P=" << F::kPrecision
                  << " E=" << F::kExponentBits
                  << " mode=" << static_cast<int>(mode)
                  << "\n  x        = " << Hex(x) << "\n  y        = " << Hex(y)
                  << "\n  actual   = " << Hex(actual) << " flags "
                  << context.flags << "\n  expected = " << Hex(expected.value)
                  << " flags " << expected.flags;
  }
  return ok;
}

template <int kLevel = 1, class F>
bool CheckAllModes(Function f, const F& x, const F& y) {
  for (const Rounding mode : kModes) {
    if (!Check<kLevel>(f, x, y, mode)) return false;
  }
  return true;
}

template <class F>
void CheckPairsExhaustively() {
  const limb_t count = limb_t{1} << (F::kPrecision + F::kExponentBits);
  for (const Function f : {Function::kPow, Function::kHypot}) {
    for (limb_t a = 0; a < count; ++a) {
      for (limb_t b = 0; b < count; ++b) {
        if (!CheckAllModes(f, F::from_bits({a}), F::from_bits({b}))) return;
      }
    }
  }
}

template <class F>
void CheckCbrtExhaustively() {
  const limb_t count = limb_t{1} << (F::kPrecision + F::kExponentBits);
  for (limb_t a = 0; a < count; ++a) {
    if (!CheckAllModes(Function::kCbrt, F::from_bits({a}), F::Zero())) return;
  }
}

TEST(PowOracle, ExhaustiveSmallFormats) {
  CheckPairsExhaustively<Float<4, 3>>();
  CheckPairsExhaustively<Float<5, 3>>();
  CheckPairsExhaustively<Float<4, 4>>();
  CheckCbrtExhaustively<Float<4, 3>>();
  CheckCbrtExhaustively<Float<8, 4>>();
  CheckCbrtExhaustively<Float<11, 5>>();
  CheckCbrtExhaustively<Float<16, 6>>();
}

template <class F>
F RandomRegular(Random& random, std::int32_t exponent, bool negative = false) {
  return MakeRegular<F>(negative, (std::min)(exponent, F::kEmax),
                        RandomSignificand<F>(random));
}

// A small integer times a power of two.
template <class F>
F Scaled(limb_t integer, std::int32_t exponent) {
  return scalbn(F(static_cast<double>(integer)), exponent);
}

// Operands where pow does something interesting.
template <class F>
void RandomPowOperands(Random& random, F& x, F& y) {
  const bool negative = random.Below(4) == 0;
  switch (random.Below(10)) {
    case 0:
      x = RandomValue<F>(random);
      y = RandomValue<F>(random);
      return;
    case 1: {
      // Small integer powers, and their exact and tied cases.
      x = random.Below(2) == 0
              ? RandomValue<F>(random)
              : Scaled<F>(random.Below(64) + 1,
                          static_cast<std::int32_t>(random.Below(40)) - 20);
      y = F(static_cast<double>(static_cast<std::int64_t>(random.Below(81)) -
                                40));
      break;
    }
    case 2: {
      // x = z^(2^s) and y = n / 2^s, exact more often than not.
      const int s = 1 + static_cast<int>(random.Below(3));
      F z = Scaled<F>(random.Below(16) + 1,
                      static_cast<std::int32_t>(random.Below(8)) - 4);
      x = z;
      for (int i = 0; i < s; ++i) x = x * x;
      y = scalbn(F(static_cast<double>(random.Below(15) + 1)), -s);
      if (random.Below(4) == 0) x = next_up(x);
      if (random.Below(4) == 0) y = -y;
      break;
    }
    case 3: {
      // Next to one, with large exponents.
      const auto offset = static_cast<std::int32_t>(
          random.Below(static_cast<limb_t>(F::kPrecision) + 8));
      x = F(1) + RandomRegular<F>(random, -offset, random.Below(2) != 0);
      y = RandomRegular<F>(random, static_cast<std::int32_t>(random.Below(
                                       static_cast<limb_t>(F::kPrecision) +
                                       F::kExponentBits)),
                           random.Below(2) != 0);
      break;
    }
    case 4: {
      // Results near the overflow and underflow thresholds.
      x = RandomRegular<F>(random, static_cast<std::int32_t>(random.Below(40)) -
                                       20);
      if (x == F(1)) x = F(3);
      const double target =
          random.Below(2) == 0
              ? F::kEmax + 1.0 - static_cast<double>(random.Below(3))
              : F::kEmin - F::kPrecision + static_cast<double>(random.Below(6)) -
                    2.0;
      y = F(target) / log2(x);
      if (random.Below(2) == 0) y = round(y);
      break;
    }
    case 5:
      // Tiny exponents, where x^y is next to one.
      x = RandomRegular<F>(random, static_cast<std::int32_t>(random.Below(40)) -
                                       20);
      y = RandomRegular<F>(random,
                           -F::kPrecision - 8 +
                               static_cast<std::int32_t>(random.Below(12)),
                           random.Below(2) != 0);
      break;
    default: {
      const auto span = static_cast<limb_t>(F::kExponentBits + 24);
      x = RandomRegular<F>(
          random, static_cast<std::int32_t>(random.Below(2 * span)) -
                      static_cast<std::int32_t>(span));
      y = RandomRegular<F>(random,
                           static_cast<std::int32_t>(random.Below(
                               static_cast<limb_t>(F::kExponentBits) + 8)) -
                               6,
                           random.Below(2) != 0);
      if (random.Below(3) == 0) y = round(y);
      break;
    }
  }
  if (negative) x = -x;
}

template <class F>
void RandomHypotOperands(Random& random, F& x, F& y) {
  switch (random.Below(6)) {
    case 0:
      x = RandomValue<F>(random);
      y = RandomValue<F>(random);
      break;
    case 1: {
      // Pythagorean triples, exact, and their neighbours.
      const limb_t m = random.Below(200) + 2;
      const limb_t n = random.Below(m - 1) + 1;
      const auto e = static_cast<std::int32_t>(
          random.Below(static_cast<limb_t>(F::kEmax) * 2) - F::kEmax);
      x = Scaled<F>(m * m - n * n, e);
      y = Scaled<F>(2 * m * n, e);
      if (random.Below(3) == 0) y = next_up(y);
      break;
    }
    case 2: {
      // Exponents apart by about P.
      const std::int32_t e = RandomExponent<F>(random);
      x = RandomRegular<F>(random, e);
      const std::int32_t gap =
          F::kPrecision - 4 + static_cast<std::int32_t>(random.Below(10));
      y = RandomRegular<F>(random, (std::max)(e - gap, F::kEmin - F::kPrecision));
      break;
    }
    default: {
      const std::int32_t e = RandomExponent<F>(random);
      x = RandomRegular<F>(random, e);
      const auto gap = static_cast<std::int32_t>(random.Below(8));
      y = RandomRegular<F>(random, (std::max)(e - gap, F::kEmin - F::kPrecision));
      break;
    }
  }
  if (random.Below(2) == 0) x = -x;
  if (random.Below(2) == 0) y = -y;
  if (random.Below(2) == 0) std::swap(x, y);
}

template <class F>
F RandomCbrtOperand(Random& random) {
  switch (random.Below(4)) {
    case 0:
      return RandomValue<F>(random);
    case 1: {
      // Perfect cubes and their neighbours.
      const limb_t root = random.Below(1u << 16) + 1;
      F x = Scaled<F>(root, 0);
      x = x * x * x;
      x = scalbn(x, 3 * (static_cast<int>(random.Below(60)) - 30));
      if (random.Below(3) == 0) x = next_down(x);
      return random.Below(2) == 0 ? -x : x;
    }
    default:
      return MakeRegular<F>(random.Below(2) != 0, RandomExponent<F>(random),
                            RandomSignificand<F>(random));
  }
}

template <class F, int kLevel = 1>
void CheckRandomly(int iterations, limb_t seed) {
  Random random(seed);
  const int total = iterations * TestScale();
  for (int i = 0; i < total; ++i) {
    F x;
    F y;
    RandomPowOperands(random, x, y);
    if (!CheckAllModes<kLevel>(Function::kPow, x, y)) return;
    if (!CheckAllModes<kLevel>(Function::kCbrt, RandomCbrtOperand<F>(random),
                               F::Zero())) {
      return;
    }
    if constexpr (kLevel == 1) {
      RandomHypotOperands(random, x, y);
      if (!CheckAllModes(Function::kHypot, x, y)) return;
    }
  }
}

TEST(PowOracle, RandomSmallFormats) {
  CheckRandomly<Float<8, 4>>(20000, 1);
  CheckRandomly<Float<11, 5>>(20000, 2);
}
TEST(PowOracle, RandomBinary128) { CheckRandomly<binary128>(3000, 3); }
TEST(PowOracle, RandomBinary256) { CheckRandomly<binary256>(1000, 4); }
TEST(PowOracle, RandomBinary512) { CheckRandomly<binary512>(300, 5); }
TEST(PowOracle, RandomOtherFormats) {
  CheckRandomly<Float<53, 11>>(3000, 6);
  CheckRandomly<Float<64, 15>>(3000, 7);
  CheckRandomly<Float<128, 20>>(1000, 8);
  CheckRandomly<Float<200, 20>>(500, 9);
}

TEST(PowOracle, LaterLevels) {
  CheckRandomly<Float<11, 5>, 2>(3000, 10);
  CheckRandomly<Float<11, 5>, 4>(3000, 11);
  CheckRandomly<binary128, 2>(300, 12);
  CheckRandomly<binary128, 4>(300, 13);
  CheckRandomly<binary512, 2>(30, 14);
  CheckRandomly<binary512, 4>(10, 15);
}

TEST(PowOracle, ExactCases) {
  using F = binary256;
  for (const Rounding mode : kModes) {
    for (int a = 1; a <= 40; ++a) {
      for (int n = -12; n <= 12; ++n) {
        EXPECT_TRUE(Check(Function::kPow, F(a), F(n), mode)) << a << "^" << n;
        EXPECT_TRUE(Check(Function::kPow, F(-a), F(n), mode))
            << -a << "^" << n;
        // a^2 and a^4 to the powers n/2 and n/4.
        EXPECT_TRUE(
            Check(Function::kPow, F(a * a), scalbn(F(n), -1), mode));
        EXPECT_TRUE(Check(Function::kPow, F(a * a) * F(a * a),
                          scalbn(F(n), -2), mode));
      }
      EXPECT_TRUE(Check(Function::kCbrt, F(a) * F(a) * F(a), F(), mode));
      EXPECT_TRUE(Check(Function::kHypot, F(3 * a), F(4 * a), mode));
    }
    // Exact results of P + 1 = 238 bits, which are ties:
    // (1.5 + 2^-118)^2 = 2.25 + 3 2^-118 + 2^-236 and
    // (1 + 2^-79)^3 = 1 + 3 2^-79 + 3 2^-158 + 2^-237.
    const F square_base = F(1.5) + scalbn(F(1), -118);
    const F cube_base = F(1) + scalbn(F(1), -79);
    EXPECT_TRUE(Check(Function::kPow, square_base, F(2), mode));
    EXPECT_TRUE(Check(Function::kPow, -cube_base, F(3), mode));
    EXPECT_TRUE(Check(Function::kPow, next_up(cube_base), F(3), mode));
  }
}

}  // namespace
}  // namespace mpfix
