// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Differential tests of pi, sin, cos and tan against MPFR, in all five
// rounding modes and with flags, plus checks of the fixed-point pi and 2/pi,
// the big-integer quotients and square roots behind them, and the error
// bounds of the argument reduction and the kernel. Set MPFIX_TEST_SCALE to
// multiply the iteration counts.

#include <cstdint>
#include <string>
#include <type_traits>

#include "mpfix/math.h"
#include "oracle.h"

namespace mpfix {
namespace {

using namespace test;  // NOLINT
using detail::BigUint;

// --- Fixed-point pi and the big-integer helpers ---

// v as an MPFR integer.
void SetMpfr(mpfr_ptr result, const BigUint& v) {
  mpfr_set_prec(result, static_cast<mpfr_prec_t>(v.bit_length() + 64));
  mpfr_set_ui(result, 0, MPFR_RNDN);
  const auto& limbs = v.limbs();
  for (std::size_t i = limbs.size(); i-- > 0;) {
    mpfr_mul_2ui(result, result, 64, MPFR_RNDN);
    mpfr_add_ui(result, result, 0, MPFR_RNDN);
    // Two halves, since unsigned long may have 32 bits.
    Mpfr part(64);
    mpfr_set_ui(part.get(), static_cast<unsigned long>(limbs[i] >> 32),
                MPFR_RNDN);
    mpfr_mul_2ui(part.get(), part.get(), 32, MPFR_RNDN);
    mpfr_add_ui(part.get(), part.get(),
                static_cast<unsigned long>(limbs[i] & 0xFFFFFFFFu), MPFR_RNDN);
    mpfr_add(result, result, part.get(), MPFR_RNDN);
  }
}

// Whether |v - c 2^bits| < 1 for the constant c that `set` gives.
bool WithinOne(const BigUint& v, std::size_t bits,
               int (*set)(mpfr_ptr, mpfr_rnd_t), bool inverse_times_two) {
  const auto precision = static_cast<mpfr_prec_t>(bits + 128);
  Mpfr exact(precision);
  set(exact.get(), MPFR_RNDN);
  if (inverse_times_two) mpfr_ui_div(exact.get(), 2, exact.get(), MPFR_RNDN);
  mpfr_mul_2ui(exact.get(), exact.get(), static_cast<unsigned long>(bits),
               MPFR_RNDN);
  Mpfr value(64);
  SetMpfr(value.get(), v);
  Mpfr difference(precision);
  mpfr_sub(difference.get(), value.get(), exact.get(), MPFR_RNDN);
  return mpfr_cmpabs_ui(difference.get(), 1) < 0;
}

TEST(PiFixedPoint, MatchesMpfr) {
  // The last sizes divide by Newton's method.
  for (const std::size_t bits : {std::size_t{64}, std::size_t{300},
                                 std::size_t{4000}, std::size_t{9000},
                                 std::size_t{40000}}) {
    const detail::PiFixedPoint p = detail::ComputePi(bits);
    EXPECT_TRUE(WithinOne(p.pi, bits, mpfr_const_pi, false)) << bits;
    EXPECT_TRUE(WithinOne(p.two_over_pi, bits, mpfr_const_pi, true)) << bits;
  }
  const auto large = detail::PiAtLeast(detail::kPiBaseBits + 1);
  EXPECT_GE(large->bits, detail::kPiBaseBits + 1);
  EXPECT_TRUE(WithinOne(large->two_over_pi, large->bits, mpfr_const_pi, true));
}

BigUint RandomBig(Random& random, std::size_t limbs) {
  std::vector<limb_t> values(limbs);
  for (limb_t& v : values) v = random.Next();
  values.back() |= limb_t{1} << 63;
  return BigUint::FromLimbs(values.data(), limbs);
}

TEST(PiFixedPoint, QuotientsAndRoots) {
  Random random(1);
  for (int i = 0; i < 40 * TestScale(); ++i) {
    const auto d_limbs = static_cast<std::size_t>(random.Below(150) + 1);
    const auto a_limbs =
        d_limbs + static_cast<std::size_t>(random.Below(150));
    const BigUint a = RandomBig(random, a_limbs);
    BigUint d = RandomBig(random, d_limbs);
    if (random.Below(2) == 0) d.ShiftRight(random.Below(60));
    BigUint expected;
    BigUint remainder;
    BigUint::DivRem(a, d, expected, remainder);
    // Through the reciprocal from the smallest size on.
    EXPECT_EQ(Compare(detail::Quotient(a, d, 1), expected), 0) << i;
    EXPECT_EQ(Compare(detail::Quotient(a, d), expected), 0) << i;

    const BigUint root = detail::SqrtFloor(a);
    BigUint above = root;
    above.Add(BigUint(1));
    EXPECT_LE(Compare(root * root, a), 0) << i;
    EXPECT_GT(Compare(above * above, a), 0) << i;
  }
}

// --- pi in F ---

template <class F>
void CheckPi() {
  for (const Rounding mode : kModes) {
    const Expected<F> expected = RoundExact<F>(
        [](mpfr_ptr result, mpfr_rnd_t rnd) {
          return mpfr_const_pi(result, rnd);
        },
        mode);
    Context context{mode};
    const F actual = pi<F>(context);
    EXPECT_TRUE(SameResult(actual, expected.value) &&
                context.flags == expected.flags)
        << "pi P=" << F::kPrecision << " mode=" << static_cast<int>(mode)
        << "\n  actual   = " << Hex(actual) << " flags " << context.flags
        << "\n  expected = " << Hex(expected.value) << " flags "
        << expected.flags;
  }
}

TEST(TrigOracle, Pi) {
  CheckPi<Float<4, 3>>();
  CheckPi<Float<11, 5>>();
  CheckPi<Float<53, 11>>();
  CheckPi<binary128>();
  CheckPi<binary256>();
  CheckPi<binary512>();
  CheckPi<Float<1000, 20>>();
  CheckPi<Float<3000, 20>>();
}

// --- sin, cos and tan ---

enum class Function { kSin, kCos, kTan };

constexpr Function kFunctions[] = {Function::kSin, Function::kCos,
                                   Function::kTan};

const char* Name(Function f) {
  switch (f) {
    case Function::kSin:
      return "sin";
    case Function::kCos:
      return "cos";
    case Function::kTan:
      return "tan";
  }
  return "";
}

template <int kLevel = 1, class F>
F ApplyMpfix(Function f, const F& x, Context& context) {
  if constexpr (kLevel == 1) {
    switch (f) {
      case Function::kSin:
        return sin(x, context);
      case Function::kCos:
        return cos(x, context);
      case Function::kTan:
        return tan(x, context);
    }
  } else {
    using Env = detail::ContextEnv;
    const Env env{context};
    switch (f) {
      case Function::kSin:
        return detail::SinImpl<F, Env, kLevel>(x, env);
      case Function::kCos:
        return detail::CosImpl<F, Env, kLevel>(x, env);
      case Function::kTan:
        return detail::TanImpl<F, Env, kLevel>(x, env);
    }
  }
  return x;
}

template <class F>
F ApplyMpfix(Function f, const F& x) {
  switch (f) {
    case Function::kSin:
      return sin(x);
    case Function::kCos:
      return cos(x);
    case Function::kTan:
      return tan(x);
  }
  return x;
}

template <class F>
Expected<F> ExpectedResult(Function f, const F& x, Rounding mode) {
  if (x.is_nan()) {
    return {F::QuietNaN(), x.is_signaling() ? unsigned{kFlagInvalid} : 0u};
  }
  const View<F> view(x);
  return RoundExact<F>(
      [&](mpfr_ptr result, mpfr_rnd_t rnd) {
        switch (f) {
          case Function::kSin:
            return mpfr_sin(result, view.get(), rnd);
          case Function::kCos:
            return mpfr_cos(result, view.get(), rnd);
          case Function::kTan:
            return mpfr_tan(result, view.get(), rnd);
        }
        return 0;
      },
      mode);
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

TEST(TrigOracle, ExhaustiveSmallFormats) {
  CheckExhaustively<Float<4, 3>>();
  CheckExhaustively<Float<8, 4>>();
  CheckExhaustively<Float<8, 5>>();
  CheckExhaustively<Float<11, 5>>();  // binary16
}

// The F nearest to n pi/2, with n up to 2^bits, and its neighbours: these
// reduce to tiny arguments.
template <class F>
F NearMultipleOfHalfPi(Random& random, int bits) {
  const auto n = static_cast<double>(random.Below(limb_t{1} << bits) + 1);
  using W = detail::MathWorking<F, 4>;
  const W product = W(n) * detail::GetTrigConstants<W>().half_pi;
  F x = detail::FloatAccess::Convert<F>(product, detail::NearestEnv{});
  switch (random.Below(3)) {
    case 0:
      return next_up(x);
    case 1:
      return next_down(x);
    default:
      return x;
  }
}

// The largest exponent of random arguments: 2^40000, or 2^10000 for wide
// formats, beyond which 2/pi takes long to compute for the first time and
// MPFR takes long to reduce them. KnownHardArguments goes further.
template <class F>
constexpr std::int32_t kMaxRandomExponent =
    (std::min)(F::kEmax, F::kPrecision > 300 ? 10000 : 40000);

// Arguments where the functions do something interesting: around zero, near
// multiples of pi/2, and up to large exponents.
template <class F>
F RandomArgument(Random& random) {
  const bool negative = random.Below(2) != 0;
  const auto regular = [&](std::int32_t exponent) {
    return MakeRegular<F>(negative, (std::min)(exponent, F::kEmax),
                          RandomSignificand<F>(random));
  };
  switch (random.Below(8)) {
    case 0: {
      const F x = RandomValue<F>(random);
      if (!x.is_finite() || x.raw_exponent() <= kMaxRandomExponent<F>) return x;
      return scalbn(x, kMaxRandomExponent<F> - x.raw_exponent());
    }
    case 1:
      // Around the thresholds for tiny arguments.
      return regular(-F::kPrecision / 2 - 4 +
                     static_cast<std::int32_t>(random.Below(6)));
    case 2: {
      // Near multiples of pi/2.
      const F x = NearMultipleOfHalfPi<F>(
          random, 1 + static_cast<int>(random.Below(
                          static_cast<limb_t>((std::min)(F::kPrecision, 50)))));
      return negative ? -x : x;
    }
    case 3: {
      // Huge arguments.
      const std::int32_t top = kMaxRandomExponent<F>;
      return regular(top - static_cast<std::int32_t>(random.Below(
                               static_cast<limb_t>(top) / 4 + 1)));
    }
    case 4:
      // Around 1/2 and pi/4, where the reduction starts.
      return regular(-2 + static_cast<std::int32_t>(random.Below(3)));
    default:
      return regular(static_cast<std::int32_t>(random.Below(
                         static_cast<limb_t>(F::kPrecision) + 40)) -
                     F::kPrecision / 2 - 8);
  }
}

template <class F, int kLevel = 1>
void CheckRandomly(int iterations, limb_t seed) {
  Random random(seed);
  const int total = iterations * TestScale();
  for (const Function f : kFunctions) {
    for (int i = 0; i < total; ++i) {
      if (!CheckAllModes<kLevel>(f, RandomArgument<F>(random))) return;
    }
  }
}

TEST(TrigOracle, RandomBinary128) { CheckRandomly<binary128>(3000, 1); }
TEST(TrigOracle, RandomBinary256) { CheckRandomly<binary256>(1000, 2); }
TEST(TrigOracle, RandomBinary512) { CheckRandomly<binary512>(200, 3); }
TEST(TrigOracle, RandomOtherFormats) {
  CheckRandomly<Float<53, 11>>(3000, 4);
  CheckRandomly<Float<64, 15>>(2000, 5);
  CheckRandomly<Float<128, 20>>(500, 6);
  CheckRandomly<Float<200, 20>>(500, 7);
}

TEST(TrigOracle, LaterLevels) {
  CheckRandomly<Float<11, 5>, 2>(2000, 8);
  CheckRandomly<Float<11, 5>, 4>(2000, 9);
  CheckRandomly<binary128, 2>(300, 10);
  CheckRandomly<binary128, 4>(300, 11);
  CheckRandomly<binary512, 2>(30, 12);
  CheckRandomly<binary512, 4>(10, 13);
}

// The double nearest to a multiple of pi/2, relative to its size: 6381956970095103
// 2^797, whose reduced argument is about 4.7e-19 (Kahan and McDonald).
TEST(TrigOracle, KnownHardArguments) {
  for (const Rounding mode : kModes) {
    const Float<53, 11> hard(std::ldexp(6381956970095103.0, 797));
    for (const Function f : kFunctions) {
      EXPECT_TRUE(Check(f, hard, mode));
      EXPECT_TRUE(Check(f, binary128(std::ldexp(6381956970095103.0, 797)),
                        mode));
    }
    // The largest finite values, and binary512 at 2^1000000, which needs 2/pi
    // to about 1.1 million bits.
    EXPECT_TRUE(Check(Function::kSin, binary128::MaxFinite(), mode));
    EXPECT_TRUE(Check(Function::kTan, binary256::MaxFinite(), mode));
    EXPECT_TRUE(Check(Function::kCos,
                      scalbn(binary512::MaxFinite(true), -(binary512::kEmax -
                                                           1000000)),
                      mode));
  }
}

// Floats very close to k pi/2, found by searching k: below 2^27, within
// 2^-(kPrecision - 68) of the first working format, where Cody and Waite's
// reduction gives way to the exact one; above, within 2^-60, where the exact
// reduction needs more bits of 2/pi than its first attempt takes.
template <class F>
void CheckVeryCloseToHalfPiMultiples(int wanted, limb_t seed) {
  using W = detail::MathWorking<F, 4>;
  constexpr int kSmallThreshold =
      68 - detail::MathWorking<F, 1>::kPrecision;
  const W half_pi = detail::GetTrigConstants<W>().half_pi;
  Random random(seed);
  int found = 0;
  for (int attempt = 0; attempt < 4000000 && found < wanted; ++attempt) {
    const bool small = (attempt & 1) == 0;
    const auto k = static_cast<double>(
        random.Below(limb_t{1} << (small ? 26 : 52)) + 1);
    const W multiple = W(k) * half_pi;
    const F x = detail::FloatAccess::Convert<F>(multiple, detail::NearestEnv{});
    const W gap = detail::FloatAccess::Convert<W>(x, detail::NearestEnv{}) -
                  multiple;
    if (gap.is_zero() ||
        gap.raw_exponent() >= (small ? kSmallThreshold : -60)) {
      continue;
    }
    ++found;
    for (const Rounding mode : kModes) {
      for (const Function f : kFunctions) {
        EXPECT_TRUE(Check(f, x, mode)) << "k=" << k;
        EXPECT_TRUE(Check(f, -x, mode)) << "k=" << k;
      }
    }
  }
  EXPECT_EQ(found, wanted) << "P=" << F::kPrecision;
}

TEST(TrigOracle, VeryCloseToHalfPiMultiples) {
  CheckVeryCloseToHalfPiMultiples<binary128>(20, 14);
  CheckVeryCloseToHalfPiMultiples<binary256>(20, 15);
  CheckVeryCloseToHalfPiMultiples<Float<200, 20>>(20, 16);
}

// --- The reduction and the kernel against their error bounds ---

template <class F, class W>
void CheckTrigBounds(int iterations, limb_t seed) {
  Random random(seed);
  const long precision = W::kPrecision + 128;
  Mpfr exact(precision);
  Mpfr error(precision);
  Mpfr bound(precision);
  const auto check = [&](const char* name, const W& value, double units,
                         const std::string& context) {
    const View<W> value_view(value);
    mpfr_sub(error.get(), value_view.get(), exact.get(), MPFR_RNDN);
    mpfr_abs(error.get(), error.get(), MPFR_RNDN);
    mpfr_mul_2si(bound.get(), exact.get(), -W::kPrecision, MPFR_RNDN);
    mpfr_mul_d(bound.get(), bound.get(), units, MPFR_RNDN);
    mpfr_abs(bound.get(), bound.get(), MPFR_RNDN);
    if (mpfr_cmp(error.get(), bound.get()) > 0) {
      mpfr_div(error.get(), error.get(), bound.get(), MPFR_RNDN);
      ADD_FAILURE() << name << " P=" << W::kPrecision << " " << context
                    << " error/bound=" << mpfr_get_d(error.get(), MPFR_RNDN);
    }
  };
  for (int i = 0; i < iterations * TestScale(); ++i) {
    // The reduction, against MPFR's remainder by pi/2 in enough precision.
    const F x = abs(RandomArgument<F>(random));
    if (!x.is_finite() || x.raw_exponent() < -1) continue;
    int quadrant = -1;
    W r = detail::ReduceByHalfPi<W>(x, quadrant);
    double reduction_units = 3.1;
    // Starting with too few bits, so that the later attempts must decide.
    if (random.Below(4) == 0) {
      int retried_quadrant = -1;
      r = detail::ReduceByHalfPi<W, -64>(x, retried_quadrant);
      EXPECT_EQ(retried_quadrant, quadrant) << "x=" << Hex(x);
    }
    // Below 2^27, through Cody and Waite's method half of the time.
    if (x.raw_exponent() < 27 && random.Below(2) == 0) {
      W small;
      if (detail::ReduceSmallByHalfPi<W>(x, quadrant, small)) {
        r = small;
        reduction_units = 1.5;
      }
    }
    {
      const long wide = static_cast<long>(x.raw_exponent()) + 3 * precision;
      Mpfr half_pi(wide);
      mpfr_const_pi(half_pi.get(), MPFR_RNDN);
      mpfr_div_2ui(half_pi.get(), half_pi.get(), 1, MPFR_RNDN);
      Mpfr quotient(wide);
      const View<F> view(x);
      mpfr_div(quotient.get(), view.get(), half_pi.get(), MPFR_RNDN);
      mpfr_rint(quotient.get(), quotient.get(), MPFR_RNDN);
      Mpfr reduced(wide);
      mpfr_fms(reduced.get(), quotient.get(), half_pi.get(), view.get(),
               MPFR_RNDN);
      mpfr_neg(reduced.get(), reduced.get(), MPFR_RNDN);
      mpfr_set(exact.get(), reduced.get(), MPFR_RNDN);
      Mpfr four(64);
      mpfr_fmod_ui(quotient.get(), quotient.get(), 4, MPFR_RNDN);
      if (mpfr_sgn(quotient.get()) < 0) {
        mpfr_add_ui(quotient.get(), quotient.get(), 4, MPFR_RNDN);
      }
      EXPECT_EQ(mpfr_get_si(quotient.get(), MPFR_RNDN), quadrant)
          << "quadrant x=" << Hex(x);
      check("reduction", r, reduction_units, "x=" + Hex(x));
    }

    // The kernel, on |r| < 0.8 with W's precision.
    const W a = MakeRegular<W>(
        random.Below(2) != 0,
        -1 - static_cast<std::int32_t>(random.Below(
                 random.Below(2) == 0 ? 8 : W::kPrecision / 2)),
        RandomSignificand<W>(random));
    const W kernel_r = abs(a) >= W(0.78) ? scalbn(a, -1) : a;
    W s;
    W c;
    double s_units;
    double c_units;
    const View<W> r_view(kernel_r);
    detail::SinCosByHalving(kernel_r, s, c, s_units, c_units);
    mpfr_sin(exact.get(), r_view.get(), MPFR_RNDN);
    check("sin by halving", s, s_units, "r=" + Hex(kernel_r));
    // cos r - 1 = -2 sin^2(r/2), without cancellation.
    mpfr_div_2ui(exact.get(), r_view.get(), 1, MPFR_RNDN);
    mpfr_sin(exact.get(), exact.get(), MPFR_RNDN);
    mpfr_sqr(exact.get(), exact.get(), MPFR_RNDN);
    mpfr_mul_si(exact.get(), exact.get(), -2, MPFR_RNDN);
    check("cos - 1 by halving", c, c_units, "r=" + Hex(kernel_r));

    // The tables of the first format, often next to their points.
    if constexpr (std::is_same_v<W, detail::MathWorking<F, 1>>) {
      W t = kernel_r;
      if (random.Below(2) == 0) {
        const W point = scalbn(W(static_cast<double>(random.Below(103))), -7) +
                        scalbn(W(static_cast<double>(random.Below(129)) - 64),
                               -14);
        t = point + scalbn(kernel_r, -14 - static_cast<int>(random.Below(20)));
        if (abs(t) >= W(0.8) || t.is_zero()) t = kernel_r;
      }
      const View<W> t_view(t);
      detail::SinCosReduced<W, true>(t, s, c, s_units, c_units);
      mpfr_sin(exact.get(), t_view.get(), MPFR_RNDN);
      check("sin with tables", s, s_units, "r=" + Hex(t));
      mpfr_cos(exact.get(), t_view.get(), MPFR_RNDN);
      check("cos with tables", c, c_units, "r=" + Hex(t));
    }
  }
}

template <class F>
void CheckTrigBoundsForAllLevels(int iterations, limb_t seed) {
  CheckTrigBounds<F, detail::MathWorking<F, 1>>(iterations, seed);
  CheckTrigBounds<F, detail::MathWorking<F, 2>>(iterations, seed + 1);
  CheckTrigBounds<F, detail::MathWorking<F, 4>>(iterations / 4 + 1, seed + 2);
}

TEST(TrigKernels, StayWithinTheirErrorBounds) {
  CheckTrigBoundsForAllLevels<Float<11, 5>>(1000, 20);
  CheckTrigBoundsForAllLevels<binary128>(1000, 30);
  CheckTrigBoundsForAllLevels<binary256>(300, 40);
  CheckTrigBoundsForAllLevels<binary512>(100, 50);
}

}  // namespace
}  // namespace mpfix
