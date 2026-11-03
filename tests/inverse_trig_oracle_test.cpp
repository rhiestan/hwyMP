// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Differential tests of asin, acos, atan and atan2 against MPFR, in all five
// rounding modes and with flags, plus checks of the error bounds of the
// arctangent's kernels. Set MPFIX_TEST_SCALE to multiply the iteration counts.

#include <algorithm>
#include <cstdint>
#include <string>
#include <type_traits>

#include "mpfix/math.h"
#include "oracle.h"

namespace mpfix {
namespace {

using namespace test;  // NOLINT

enum class Function { kAsin, kAcos, kAtan, kAtan2 };

constexpr Function kUnary[] = {Function::kAsin, Function::kAcos,
                               Function::kAtan};

const char* Name(Function f) {
  switch (f) {
    case Function::kAsin:
      return "asin";
    case Function::kAcos:
      return "acos";
    case Function::kAtan:
      return "atan";
    case Function::kAtan2:
      return "atan2";
  }
  return "";
}

// f(x), or atan2(x, y).
template <int kLevel = 1, class F>
F ApplyMpfix(Function f, const F& x, const F& y, Context& context) {
  if constexpr (kLevel == 1) {
    switch (f) {
      case Function::kAsin:
        return asin(x, context);
      case Function::kAcos:
        return acos(x, context);
      case Function::kAtan:
        return atan(x, context);
      case Function::kAtan2:
        return atan2(x, y, context);
    }
  } else {
    using Env = detail::ContextEnv;
    const Env env{context};
    switch (f) {
      case Function::kAsin:
        return detail::AsinImpl<F, Env, kLevel>(x, env);
      case Function::kAcos:
        return detail::AcosImpl<F, Env, kLevel>(x, env);
      case Function::kAtan:
        return detail::AtanImpl<F, Env, kLevel>(x, env);
      case Function::kAtan2:
        return detail::Atan2Impl<F, Env, kLevel>(x, y, env);
    }
  }
  return x;
}

template <class F>
F ApplyMpfix(Function f, const F& x, const F& y) {
  switch (f) {
    case Function::kAsin:
      return asin(x);
    case Function::kAcos:
      return acos(x);
    case Function::kAtan:
      return atan(x);
    case Function::kAtan2:
      return atan2(x, y);
  }
  return x;
}

template <class F>
Expected<F> ExpectedResult(Function f, const F& x, const F& y, Rounding mode) {
  const bool binary = f == Function::kAtan2;
  if (x.is_nan() || (binary && y.is_nan())) {
    const bool signaling = x.is_signaling() || (binary && y.is_signaling());
    return {F::QuietNaN(), signaling ? unsigned{kFlagInvalid} : 0u};
  }
  const View<F> vx(x);
  const View<F> vy(y);
  return RoundExact<F>(
      [&](mpfr_ptr result, mpfr_rnd_t rnd) {
        switch (f) {
          case Function::kAsin:
            return mpfr_asin(result, vx.get(), rnd);
          case Function::kAcos:
            return mpfr_acos(result, vx.get(), rnd);
          case Function::kAtan:
            return mpfr_atan(result, vx.get(), rnd);
          case Function::kAtan2:
            return mpfr_atan2(result, vx.get(), vy.get(), rnd);
        }
        return 0;
      },
      mode);
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
void CheckUnaryExhaustively() {
  const limb_t count = limb_t{1} << (F::kPrecision + F::kExponentBits);
  for (const Function f : kUnary) {
    for (limb_t bits = 0; bits < count; ++bits) {
      if (!CheckAllModes(f, F::from_bits({bits}), F::Zero())) return;
    }
  }
}

template <class F>
void CheckAtan2Exhaustively() {
  const limb_t count = limb_t{1} << (F::kPrecision + F::kExponentBits);
  for (limb_t a = 0; a < count; ++a) {
    for (limb_t b = 0; b < count; ++b) {
      if (!CheckAllModes(Function::kAtan2, F::from_bits({a}),
                         F::from_bits({b}))) {
        return;
      }
    }
  }
}

TEST(InverseTrigOracle, ExhaustiveSmallFormats) {
  CheckUnaryExhaustively<Float<4, 3>>();
  CheckUnaryExhaustively<Float<8, 4>>();
  CheckUnaryExhaustively<Float<8, 5>>();
  CheckUnaryExhaustively<Float<11, 5>>();  // binary16
  CheckAtan2Exhaustively<Float<4, 3>>();
  CheckAtan2Exhaustively<Float<5, 3>>();
  CheckAtan2Exhaustively<Float<4, 4>>();
}

template <class F>
F RandomRegular(Random& random, std::int32_t exponent, bool negative) {
  return MakeRegular<F>(negative,
                        (std::max)((std::min)(exponent, F::kEmax),
                                   F::kEmin - F::kPrecision + 1),
                        RandomSignificand<F>(random));
}

// Arguments of the unary functions: mostly within [-1, 1] for asin and acos,
// often next to +-1 and 0, and of any size for atan.
template <class F>
F RandomUnary(Function f, Random& random) {
  const bool negative = random.Below(2) != 0;
  F x;
  switch (random.Below(6)) {
    case 0:
      x = RandomValue<F>(random);
      break;
    case 1: {
      // Next to +-1.
      const auto offset = static_cast<std::int32_t>(
          random.Below(static_cast<limb_t>(F::kPrecision) + 4));
      x = F(1) - RandomRegular<F>(random, -1 - offset, false);
      if (negative) x = -x;
      break;
    }
    case 2:
      // Around the thresholds for tiny arguments.
      x = RandomRegular<F>(random,
                           -F::kPrecision / 2 - 4 +
                               static_cast<std::int32_t>(random.Below(6)),
                           negative);
      break;
    default:
      x = RandomRegular<F>(random,
                           static_cast<std::int32_t>(random.Below(
                               static_cast<limb_t>(F::kPrecision) + 20)) -
                               F::kPrecision - 8,
                           negative);
      if (f == Function::kAtan && random.Below(2) == 0) x = F(1) / x;
      break;
  }
  if (f != Function::kAtan && random.Below(4) != 0 && abs(x) > F(1)) {
    x = F(1) / x;
  }
  return x;
}

// Pairs (y, x) for atan2: random, with exponents far apart and close, with
// exactly representable tiny quotients, and with zeros and infinities.
template <class F>
void RandomAtan2Operands(Random& random, F& y, F& x) {
  const bool negative_y = random.Below(2) != 0;
  const bool negative_x = random.Below(2) != 0;
  switch (random.Below(6)) {
    case 0:
      y = RandomValue<F>(random);
      x = RandomValue<F>(random);
      return;
    case 1: {
      // Exponents apart by about 1.5 P, where atan2 becomes the quotient.
      const std::int32_t e = static_cast<std::int32_t>(random.Below(40)) - 20;
      x = RandomRegular<F>(random, e, negative_x);
      y = RandomRegular<F>(
          random,
          e - 3 * F::kPrecision / 2 - 14 + static_cast<std::int32_t>(random.Below(10)),
          negative_y);
      return;
    }
    case 2: {
      // y = x m 2^-k for a small integer m: quotients of few bits, exact,
      // down to where even the last working format cannot tell atan2 from
      // the quotient.
      x = RandomRegular<F>(random, static_cast<std::int32_t>(random.Below(20)) - 10,
                           negative_x);
      const auto m = static_cast<double>(random.Below(15) + 1);
      const int k = static_cast<int>(
          random.Below(static_cast<limb_t>(4 * F::kPrecision) + 300));
      y = scalbn(x * F(m), -(std::min)(k, F::kEmax - F::kEmin));
      if (negative_y != y.signbit()) y = -y;
      if (random.Below(4) == 0) y = next_up(y);
      return;
    }
    case 3: {
      // Zeros and infinities against anything.
      const F special = random.Below(2) == 0 ? F::Zero(negative_y)
                                             : F::Infinity(negative_y);
      const F other = RandomRegular<F>(
          random, static_cast<std::int32_t>(random.Below(40)) - 20, negative_x);
      if (random.Below(2) == 0) {
        y = special;
        x = random.Below(3) == 0 ? F::Zero(negative_x) : other;
      } else {
        y = other;
        x = random.Below(2) == 0 ? F::Infinity(negative_x) : F::Zero(negative_x);
      }
      return;
    }
    default: {
      const std::int32_t e = static_cast<std::int32_t>(random.Below(
                                 static_cast<limb_t>(F::kExponentBits) + 20)) -
                             10;
      x = RandomRegular<F>(random, e, negative_x);
      y = RandomRegular<F>(
          random, e + static_cast<std::int32_t>(random.Below(16)) - 8, negative_y);
      return;
    }
  }
}

template <class F, int kLevel = 1>
void CheckRandomly(int iterations, limb_t seed) {
  Random random(seed);
  const int total = iterations * TestScale();
  for (int i = 0; i < total; ++i) {
    for (const Function f : kUnary) {
      if (!CheckAllModes<kLevel>(f, RandomUnary<F>(f, random), F::Zero())) {
        return;
      }
    }
    F y;
    F x;
    RandomAtan2Operands(random, y, x);
    if (!CheckAllModes<kLevel>(Function::kAtan2, y, x)) return;
  }
}

TEST(InverseTrigOracle, RandomSmallFormats) {
  CheckRandomly<Float<8, 4>>(20000, 1);
  CheckRandomly<Float<11, 5>>(20000, 2);
}
TEST(InverseTrigOracle, RandomBinary128) { CheckRandomly<binary128>(3000, 3); }
TEST(InverseTrigOracle, RandomBinary256) { CheckRandomly<binary256>(1000, 4); }
TEST(InverseTrigOracle, RandomBinary512) { CheckRandomly<binary512>(300, 5); }
TEST(InverseTrigOracle, RandomOtherFormats) {
  CheckRandomly<Float<53, 11>>(3000, 6);
  CheckRandomly<Float<64, 15>>(2000, 7);
  CheckRandomly<Float<128, 20>>(500, 8);
  CheckRandomly<Float<200, 20>>(500, 9);
}

TEST(InverseTrigOracle, LaterLevels) {
  CheckRandomly<Float<11, 5>, 2>(2000, 10);
  CheckRandomly<Float<11, 5>, 4>(2000, 11);
  CheckRandomly<binary128, 2>(300, 12);
  CheckRandomly<binary128, 4>(300, 13);
  CheckRandomly<binary512, 2>(30, 14);
  CheckRandomly<binary512, 4>(10, 15);
}

TEST(InverseTrigOracle, SpecialValues) {
  using F = binary256;
  const F inf = F::Infinity();
  const F zero = F::Zero();
  for (const Rounding mode : kModes) {
    for (const F& v : {F(1), F(-1), zero, -zero, inf, -inf, F(0.5), F(-0.5)}) {
      for (const Function f : kUnary) EXPECT_TRUE(Check(f, v, zero, mode));
      for (const F& w :
           {F(1), F(-1), zero, -zero, inf, -inf, F(2), F(-3), F(1e-300)}) {
        EXPECT_TRUE(Check(Function::kAtan2, v, w, mode));
      }
    }
    // The largest and smallest quotients.
    EXPECT_TRUE(Check(Function::kAtan2, F::MinSubnormal(), F::MaxFinite(), mode));
    EXPECT_TRUE(Check(Function::kAtan2, F::MaxFinite(), F::MinSubnormal(), mode));
    EXPECT_TRUE(
        Check(Function::kAtan2, -F::MinSubnormal(), -F::MaxFinite(), mode));
    EXPECT_TRUE(Check(Function::kAtan, F::MaxFinite(true), zero, mode));
  }
}

// --- The kernels against their error bounds ---

template <class F, class W>
void CheckAtanBounds(int iterations, limb_t seed) {
  constexpr bool kTables = std::is_same_v<W, detail::MathWorking<F, 1>>;
  Random random(seed);
  const long precision = W::kPrecision + 128;
  Mpfr exact(precision);
  Mpfr error(precision);
  Mpfr bound(precision);
  for (int i = 0; i < iterations * TestScale(); ++i) {
    // From 2^-(P/2) to 2^20, often next to the table points.
    W a = MakeRegular<W>(false,
                         static_cast<std::int32_t>(random.Below(
                             static_cast<limb_t>(W::kPrecision / 2 + 20))) -
                             W::kPrecision / 2,
                         RandomSignificand<W>(random));
    if (random.Below(3) == 0) {
      a = scalbn(W(static_cast<double>(random.Below(129))), -7) +
          scalbn(W(static_cast<double>(random.Below(129)) - 64), -14) +
          scalbn(a, -14 - static_cast<int>(random.Below(30)));
      if (a.signbit() || a.is_zero()) continue;
    }
    const detail::Tracked<W> t = detail::AtanNonnegative<W, kTables>(a, 0);
    const View<W> view(a);
    mpfr_atan(exact.get(), view.get(), MPFR_RNDN);
    const View<W> t_view(t.value);
    mpfr_sub(error.get(), t_view.get(), exact.get(), MPFR_RNDN);
    mpfr_abs(error.get(), error.get(), MPFR_RNDN);
    mpfr_mul_2si(bound.get(), exact.get(), -W::kPrecision, MPFR_RNDN);
    mpfr_mul_d(bound.get(), bound.get(), t.units, MPFR_RNDN);
    if (mpfr_cmp(error.get(), bound.get()) > 0) {
      mpfr_div(error.get(), error.get(), bound.get(), MPFR_RNDN);
      ADD_FAILURE() << "atan P=" << W::kPrecision << " tables=" << kTables
                    << " a=" << Hex(a)
                    << " error/bound=" << mpfr_get_d(error.get(), MPFR_RNDN);
    }
  }
}

template <class F>
void CheckAtanBoundsForAllLevels(int iterations, limb_t seed) {
  CheckAtanBounds<F, detail::MathWorking<F, 1>>(iterations, seed);
  CheckAtanBounds<F, detail::MathWorking<F, 2>>(iterations, seed + 1);
  CheckAtanBounds<F, detail::MathWorking<F, 4>>(iterations / 4 + 1, seed + 2);
}

TEST(InverseTrigKernels, StayWithinTheirErrorBounds) {
  CheckAtanBoundsForAllLevels<Float<11, 5>>(2000, 20);
  CheckAtanBoundsForAllLevels<binary128>(2000, 30);
  CheckAtanBoundsForAllLevels<binary256>(500, 40);
  CheckAtanBoundsForAllLevels<binary512>(200, 50);
}

}  // namespace
}  // namespace mpfix
