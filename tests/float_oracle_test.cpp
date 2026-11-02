// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Differential tests of Float against MPFR, covering results and flags in all
// five rounding modes. Set MPFIX_TEST_SCALE to multiply the iteration counts
// and, when it is at least 10, to enable the larger exhaustive formats.

#include <cstdint>
#include <cstdio>
#include <string>

#include "oracle.h"

namespace mpfix {
namespace {

using namespace test;  // NOLINT

enum class Op { kAdd, kSub, kMul, kDiv, kRemainder, kFma, kSqrt, kScalbn };

constexpr Op kBinaryOps[] = {Op::kAdd, Op::kSub, Op::kMul, Op::kDiv,
                             Op::kRemainder};

const char* Name(Op op) {
  switch (op) {
    case Op::kAdd:
      return "add";
    case Op::kSub:
      return "sub";
    case Op::kMul:
      return "mul";
    case Op::kDiv:
      return "div";
    case Op::kRemainder:
      return "remainder";
    case Op::kFma:
      return "fma";
    case Op::kSqrt:
      return "sqrt";
    default:
      return "scalbn";
  }
}

// mpfr_sqrt in MPFR 4.2.2 misrounds some radicands whose root has a long run of
// zeros below the rounding position, at the precisions of binary512 at least:
// see the Sqrt.Mpfr422Regression test. So the root is taken at much higher
// precision and proved correct by exact squaring, which needs no trust in
// mpfr_sqrt, before it is rounded to the precision of `result`.
int CertifiedSqrt(mpfr_ptr result, mpfr_srcptr a, mpfr_rnd_t rnd) {
  if (!mpfr_regular_p(a) || mpfr_signbit(a) != 0) {
    return mpfr_sqrt(result, a, rnd);
  }
  const long wide = 2 * static_cast<long>(mpfr_get_prec(result)) +
                    2 * static_cast<long>(mpfr_get_prec(a)) + 64;
  Mpfr root(wide);
  mpfr_sqrt(root.get(), a, MPFR_RNDZ);
  Mpfr above(wide);
  mpfr_set(above.get(), root.get(), MPFR_RNDZ);
  mpfr_nextabove(above.get());

  Mpfr square(2 * wide);
  mpfr_mul(square.get(), root.get(), root.get(), MPFR_RNDZ);
  Mpfr square_above(2 * wide);
  mpfr_mul(square_above.get(), above.get(), above.get(), MPFR_RNDZ);
  const int order = mpfr_cmp(square.get(), a);
  if (order > 0 || mpfr_cmp(square_above.get(), a) <= 0) {
    ADD_FAILURE() << "The wide square root is not the floor of the true root";
  }
  if (order == 0) return mpfr_set(result, root.get(), rnd);

  // Strictly between the two neighbours, as the true root is, so it rounds the
  // same way in every mode.
  Mpfr between(wide + 1);
  mpfr_set(between.get(), root.get(), MPFR_RNDZ);
  mpfr_nextabove(between.get());
  return mpfr_set(result, between.get(), rnd);
}

struct MpfrArgs {
  mpfr_srcptr a;
  mpfr_srcptr b;
  mpfr_srcptr c;
  long n;
};

// Applies `op` to the operands, rounded to the precision of `result`. Returns
// MPFR's ternary value.
int Apply(Op op, mpfr_ptr result, const MpfrArgs& x, mpfr_rnd_t rnd) {
  switch (op) {
    case Op::kAdd:
      return mpfr_add(result, x.a, x.b, rnd);
    case Op::kSub:
      return mpfr_sub(result, x.a, x.b, rnd);
    case Op::kMul:
      return mpfr_mul(result, x.a, x.b, rnd);
    case Op::kDiv:
      return mpfr_div(result, x.a, x.b, rnd);
    case Op::kRemainder:
      return mpfr_remainder(result, x.a, x.b, rnd);
    case Op::kFma:
      return mpfr_fma(result, x.a, x.b, x.c, rnd);
    case Op::kSqrt:
      return CertifiedSqrt(result, x.a, rnd);
    default:
      return mpfr_mul_2si(result, x.a, x.n, rnd);
  }
}

// The operands of one operation. Operations ignore those they do not use.
template <class F>
struct Args {
  F a;
  F b;
  F c;
  int n;
};

template <class F>
Args<F> Binary(const F& a, const F& b) {
  return {a, b, F::Zero(), 0};
}

template <class F>
Expected<F> Compute(Op op, const Args<F>& args, Rounding mode) {
  const F& a = args.a;
  const F& b = args.b;
  const F& c = args.c;
  const bool unary = op == Op::kSqrt || op == Op::kScalbn;
  const bool uses_b = !unary;
  const bool uses_c = op == Op::kFma;
  if (a.is_nan() || (uses_b && b.is_nan()) || (uses_c && c.is_nan())) {
    const bool signaling = a.is_signaling() || (uses_b && b.is_signaling()) ||
                           (uses_c && c.is_signaling());
    return {F::QuietNaN(), signaling ? unsigned{kFlagInvalid} : 0u};
  }
  if (op == Op::kDiv && b.is_zero() && a.is_finite() && !a.is_zero()) {
    return {F::Infinity(a.signbit() != b.signbit()), unsigned{kFlagDivByZero}};
  }

  const View<F> view_a(a);
  const View<F> view_b(b);
  const View<F> view_c(c);
  const MpfrArgs x{view_a.get(), view_b.get(), view_c.get(), args.n};
  const auto apply = [op, &x](mpfr_ptr result, mpfr_rnd_t rnd) {
    return Apply(op, result, x, rnd);
  };
  return RoundExact<F>(apply, mode);
}

template <class F>
F ApplyFloat(Op op, const Args<F>& x, Context& context) {
  switch (op) {
    case Op::kAdd:
      return add(x.a, x.b, context);
    case Op::kSub:
      return sub(x.a, x.b, context);
    case Op::kMul:
      return mul(x.a, x.b, context);
    case Op::kDiv:
      return div(x.a, x.b, context);
    case Op::kRemainder:
      return remainder(x.a, x.b, context);
    case Op::kFma:
      return fma(x.a, x.b, x.c, context);
    case Op::kSqrt:
      return sqrt(x.a, context);
    default:
      return scalbn(x.a, x.n, context);
  }
}

template <class F>
F ApplyOperator(Op op, const Args<F>& x) {
  switch (op) {
    case Op::kAdd:
      return x.a + x.b;
    case Op::kSub:
      return x.a - x.b;
    case Op::kMul:
      return x.a * x.b;
    case Op::kDiv:
      return x.a / x.b;
    case Op::kRemainder:
      return remainder(x.a, x.b);
    case Op::kFma:
      return fma(x.a, x.b, x.c);
    case Op::kSqrt:
      return sqrt(x.a);
    default:
      return scalbn(x.a, x.n);
  }
}

// How many checks reached each kind of outcome, so that a generator which
// stops producing a hard case fails loudly instead of silently testing less.
struct Coverage {
  long checks = 0;
  long invalid = 0;
  long overflow = 0;
  long underflow = 0;
  long inexact = 0;
  long exact_regular = 0;
  long exact_subnormal = 0;
  long subnormal_results = 0;
  long exact_cancellation = 0;
  long exact_division = 0;
  long exact_root = 0;
};

Coverage coverage;

template <class F>
void Record(Op op, const F& a, const F& b, const Expected<F>& expected) {
  ++coverage.checks;
  const unsigned flags = expected.flags;
  if (flags == 0 && expected.value.is_finite() && !expected.value.is_zero()) {
    if (op == Op::kDiv) ++coverage.exact_division;
    if (op == Op::kSqrt) ++coverage.exact_root;
  }
  if ((flags & kFlagInvalid) != 0) ++coverage.invalid;
  if ((flags & kFlagOverflow) != 0) ++coverage.overflow;
  if ((flags & kFlagUnderflow) != 0) ++coverage.underflow;
  if ((flags & kFlagInexact) != 0) ++coverage.inexact;
  const bool regular_inputs =
      a.is_finite() && b.is_finite() && !a.is_zero() && !b.is_zero();
  if (!regular_inputs) return;
  if (expected.value.is_subnormal()) {
    ++coverage.subnormal_results;
    if (flags == 0) ++coverage.exact_subnormal;
  } else if (flags == 0 && expected.value.is_zero()) {
    ++coverage.exact_cancellation;
  } else if (flags == 0) {
    ++coverage.exact_regular;
  }
}

// Returns false, after reporting, if any result or flag differs.
template <class F>
bool Check(Op op, const Args<F>& args, Rounding mode) {
  Context context{mode};
  const F actual = ApplyFloat(op, args, context);
  const Expected<F> expected = Compute(op, args, mode);
  Record(op, args.a, args.b, expected);
  bool ok = SameResult(actual, expected.value) &&
            context.flags == expected.flags;
  if (mode == Rounding::kNearestEven) {
    ok = ok && SameResult(ApplyOperator(op, args), expected.value);
  }
  if (!ok) {
    ADD_FAILURE() << "P=" << F::kPrecision << " E=" << F::kExponentBits << ' '
                  << Name(op) << " mode=" << static_cast<int>(mode)
                  << "\n  a        = " << Hex(args.a)
                  << "\n  b        = " << Hex(args.b)
                  << "\n  c        = " << Hex(args.c) << "\n  n        = "
                  << args.n << "\n  actual   = " << Hex(actual) << " flags "
                  << context.flags << "\n  expected = " << Hex(expected.value)
                  << " flags " << expected.flags;
  }
  return ok;
}

// Checks one operation in every rounding mode.
template <class F>
bool CheckAllModes(Op op, const Args<F>& args) {
  for (const Rounding mode : kModes) {
    if (!Check(op, args, mode)) return false;
  }
  return true;
}

// roundToIntegral in one mode, in both its silent and its exact form.
template <class F>
bool CheckRoundToIntegral(const F& x, Rounding mode) {
  Expected<F> expected{x, 0u};
  bool inexact = false;
  if (x.is_nan()) {
    expected = {F::QuietNaN(),
                x.is_signaling() ? unsigned{kFlagInvalid} : 0u};
  } else if (x.is_finite() && !x.is_zero()) {
    const View<F> view(x);
    inexact = mpfr_integer_p(view.get()) == 0;
    Mpfr rounded(F::kPrecision + 1);
    switch (mode) {
      case Rounding::kNearestEven:
        mpfr_rint(rounded.get(), view.get(), MPFR_RNDN);
        break;
      case Rounding::kNearestAway:
        mpfr_round(rounded.get(), view.get());
        break;
      case Rounding::kTowardZero:
        mpfr_rint(rounded.get(), view.get(), MPFR_RNDZ);
        break;
      case Rounding::kUp:
        mpfr_rint(rounded.get(), view.get(), MPFR_RNDU);
        break;
      case Rounding::kDown:
        mpfr_rint(rounded.get(), view.get(), MPFR_RNDD);
        break;
    }
    if (!mpfr_zero_p(rounded.get()) &&
        mpfr_get_exp(rounded.get()) - 1 > F::kEmax) {
      // Possible only in tiny formats, whose largest values are not integers.
      expected = {OverflowResult<F>(x.signbit(), mode),
                  unsigned{kFlagOverflow}};
    } else {
      expected = {FromMpfr<F>(rounded.get()), 0u};
    }
  }

  Context silent{mode};
  const F silent_result = round_to_integral(x, silent);
  Context exact{mode};
  const F exact_result = round_to_integral_exact(x, exact);
  const unsigned exact_flags =
      expected.flags | (inexact ? unsigned{kFlagInexact} : 0u);
  const bool ok = SameResult(silent_result, expected.value) &&
                  silent.flags == expected.flags &&
                  SameResult(exact_result, expected.value) &&
                  exact.flags == exact_flags;
  if (!ok) {
    ADD_FAILURE() << "P=" << F::kPrecision << " E=" << F::kExponentBits
                  << " roundToIntegral mode=" << static_cast<int>(mode)
                  << "\n  x        = " << Hex(x)
                  << "\n  silent   = " << Hex(silent_result) << " flags "
                  << silent.flags << "\n  exact    = " << Hex(exact_result)
                  << " flags " << exact.flags
                  << "\n  expected = " << Hex(expected.value) << " flags "
                  << exact_flags;
  }
  return ok;
}

// Every operand pair, operation and rounding mode, and every operand of the
// operations that take one.
template <class F>
void CheckExhaustively() {
  const limb_t count = limb_t{1} << (F::kPrecision + F::kExponentBits);
  const int shift_limit = 2 * F::kEmax + 2 * F::kPrecision + 2;
  for (limb_t a_bits = 0; a_bits < count; ++a_bits) {
    const F a = F::from_bits({a_bits});
    for (limb_t b_bits = 0; b_bits < count; ++b_bits) {
      const F b = F::from_bits({b_bits});
      for (const Op op : kBinaryOps) {
        if (!CheckAllModes(op, Binary(a, b))) return;
      }
    }
    if (!CheckAllModes(Op::kSqrt, Binary(a, a))) return;
    for (int n = -shift_limit; n <= shift_limit; ++n) {
      if (!CheckAllModes(Op::kScalbn, Args<F>{a, a, a, n})) return;
    }
    for (const Rounding mode : kModes) {
      if (!CheckRoundToIntegral(a, mode)) return;
    }
  }
}

template <class F>
void CheckFmaExhaustively() {
  const limb_t count = limb_t{1} << (F::kPrecision + F::kExponentBits);
  for (limb_t a_bits = 0; a_bits < count; ++a_bits) {
    const F a = F::from_bits({a_bits});
    for (limb_t b_bits = 0; b_bits < count; ++b_bits) {
      const F b = F::from_bits({b_bits});
      for (limb_t c_bits = 0; c_bits < count; ++c_bits) {
        const F c = F::from_bits({c_bits});
        if (!CheckAllModes(Op::kFma, Args<F>{a, b, c, 0})) return;
      }
    }
  }
}

TEST(FloatOracle, ExhaustiveP3E2) { CheckExhaustively<Float<3, 2>>(); }
TEST(FloatOracle, ExhaustiveP3E3) { CheckExhaustively<Float<3, 3>>(); }
TEST(FloatOracle, ExhaustiveP4E3) { CheckExhaustively<Float<4, 3>>(); }
TEST(FloatOracle, ExhaustiveP5E2) { CheckExhaustively<Float<5, 2>>(); }

TEST(FloatOracle, ExhaustiveFmaP3E2) { CheckFmaExhaustively<Float<3, 2>>(); }
TEST(FloatOracle, ExhaustiveFmaP3E3) { CheckFmaExhaustively<Float<3, 3>>(); }

TEST(FloatOracle, ExhaustiveLargerFormats) {
  if (TestScale() < 10) GTEST_SKIP() << "Needs MPFIX_TEST_SCALE >= 10";
  CheckExhaustively<Float<6, 3>>();
  CheckExhaustively<Float<5, 4>>();
  CheckFmaExhaustively<Float<4, 3>>();
}

// A second operand with a related exponent, so that effective subtraction
// cancels and alignment shifts land on every word and bit boundary.
template <class F>
F RelatedValue(Random& random, const F& a) {
  if (!a.is_finite() || a.is_zero() || random.Below(2) == 0) {
    return RandomValue<F>(random);
  }
  const bool negative = random.Below(2) != 0;
  const limb_t spread = random.Below(2) != 0
                            ? 5
                            : static_cast<limb_t>(F::kPrecision) * 2 + 8;
  const std::int32_t exponent =
      a.raw_exponent() + static_cast<std::int32_t>(spread / 2) -
      static_cast<std::int32_t>(random.Below(spread));
  if (exponent > F::kEmax) return RandomValue<F>(random);

  typename F::Limbs limbs = RandomSignificand<F>(random);
  if (random.Below(2) == 0) {
    // Share a's high bits so the difference cancels deeply.
    limbs = a.raw_limbs();
    const int keep = static_cast<int>(random.Below(
        static_cast<limb_t>(F::kPrecision)));
    const int position = static_cast<int>(F::kLimbs) * kLimbBits - 1 - keep;
    limbs[static_cast<std::size_t>(position) / kLimbBits] ^=
        limb_t{1} << (position % kLimbBits);
    limbs[F::kLimbs - 1] |= limb_t{1} << (kLimbBits - 1);
  }
  return MakeRegular<F>(negative, exponent, limbs);
}

// A scale that sends values to any part of the range, or a small one.
template <class F>
int RandomShift(Random& random) {
  if (random.Below(2) != 0) {
    return static_cast<int>(random.Below(2 * F::kPrecision + 9)) -
           (F::kPrecision + 4);
  }
  const limb_t range = limb_t{4} * F::kEmax + 4 * F::kPrecision;
  return static_cast<int>(random.Below(range)) - static_cast<int>(range / 2);
}

template <class F>
void CheckRandomly(int iterations, limb_t seed) {
  Random random(seed);
  const int total = iterations * TestScale();
  coverage = {};
  for (int i = 0; i < total; ++i) {
    const F a = RandomValue<F>(random);
    const F b = RelatedValue<F>(random, a);
    const F product = a * b;
    const F c = product.is_nan() ? RandomValue<F>(random)
                                 : RelatedValue<F>(random, product);
    const Args<F> args{a, b, c, RandomShift<F>(random)};
    for (const Op op : {Op::kAdd, Op::kSub, Op::kMul, Op::kDiv}) {
      if (!CheckAllModes(op, args)) return;
    }
    if (!CheckAllModes(Op::kFma, args)) return;
    if (!CheckAllModes(Op::kScalbn, args)) return;
    if (!CheckAllModes(Op::kSqrt, Binary(abs(a), b))) return;

    // Remainder is exact, so one mode is enough. Its time grows with the
    // exponent difference, so far-apart operands run only occasionally.
    const Rounding mode = kModes[static_cast<unsigned>(i) % 5];
    const std::int64_t distance =
        std::int64_t{a.raw_exponent()} - b.raw_exponent();
    if ((distance < 4096 && distance > -4096) || random.Below(512) == 0) {
      if (!Check(Op::kRemainder, args, mode)) return;
    }
    if (!CheckRoundToIntegral(a, mode)) return;
    const F near_integer =
        scalbn(a, F::kPrecision / 2 + 1 - a.raw_exponent() % 64);
    if (!CheckRoundToIntegral(near_integer, mode)) return;
  }

  // Each kind of outcome must make up a real share of the checks.
  const long minimum = coverage.checks / 2000;
  std::printf(
      "P=%d checks=%ld invalid=%ld overflow=%ld underflow=%ld inexact=%ld "
      "exact=%ld exact_subnormal=%ld subnormal=%ld cancelled=%ld\n",
      F::kPrecision, coverage.checks, coverage.invalid, coverage.overflow,
      coverage.underflow, coverage.inexact, coverage.exact_regular,
      coverage.exact_subnormal, coverage.subnormal_results,
      coverage.exact_cancellation);
  EXPECT_GE(coverage.checks, 35L * total);
  EXPECT_GT(coverage.invalid, minimum);
  EXPECT_GT(coverage.overflow, minimum);
  EXPECT_GT(coverage.underflow, minimum);
  EXPECT_GT(coverage.inexact, minimum);
  EXPECT_GT(coverage.exact_regular, minimum);
  EXPECT_GT(coverage.exact_subnormal, minimum);
  EXPECT_GT(coverage.subnormal_results, minimum);
  EXPECT_GT(coverage.exact_cancellation, minimum);
}

// A value near one whose low half of the significand is zero, so that products
// of two of them are exact.
template <class F>
F ShortValue(Random& random) {
  typename F::Limbs limbs = RandomSignificand<F>(random);
  constexpr int kTotal = static_cast<int>(F::kLimbs) * kLimbBits;
  const int clear = kTotal - F::kPrecision / 2;
  detail::ShiftRight(limbs, clear);
  detail::ShiftLeft(limbs, clear);
  return F::FromRaw(false, static_cast<std::int32_t>(random.Below(9)) - 4,
                    limbs);
}

template <class F>
F FullValue(Random& random) {
  return F::FromRaw(random.Below(2) != 0,
                    static_cast<std::int32_t>(random.Below(9)) - 4,
                    RandomSignificand<F>(random));
}

// Products of the form 2^e (1 + tiny), for fma. Their tail falls outside the
// double-width window when aligned to a large addend, and only the sticky bit
// remembers it. These come from 2^(4t+2) + 1 = (2^(2t+1) - 2^(t+1) + 1)
// (2^(2t+1) + 2^(t+1) + 1), since no other operands give so long a zero run.
template <class F>
bool CheckTinyTailFma(Random& random) {
  const int t = (F::kPrecision - 2) / 2;
  const F one = scalbn(F::MinNormal(), -F::kEmin);
  const F high = one + scalbn(one, -(2 * t + 1));
  const F step = scalbn(one, -t);
  const F a = high - step;
  const F b = high + step;
  const int scale = static_cast<int>(random.Below(4));
  const F scaled = scalbn(b, scale);
  const int p = F::kPrecision;

  // Addends whose unit in the last place puts the product exactly on a tie or
  // on a carry out of the top, if not for the tail.
  const F tie = scalbn(one, p + scale) + scalbn(one, scale + 1);
  const F carry = next_down(scalbn(one, p + scale));
  const F addends[] = {tie, -tie, next_up(tie), carry, -carry,
                       scalbn(carry, 1)};
  for (const F& addend : addends) {
    if (!CheckAllModes(Op::kFma, Args<F>{a, scaled, addend, 0})) return false;
    if (!CheckAllModes(Op::kFma, Args<F>{-a, scaled, addend, 0})) return false;
  }
  return true;
}

// Operands built to land on, or right next to, the exact and halfway cases that
// random operands almost never reach.
template <class F>
void CheckHardCases(int iterations, limb_t seed) {
  Random random(seed);
  const int total = iterations * TestScale();
  coverage = {};
  for (int i = 0; i < total; ++i) {
    const F s = ShortValue<F>(random);
    const F t = ShortValue<F>(random);
    const F y = abs(FullValue<F>(random));
    const F z = next_up(y);
    const F b = FullValue<F>(random);

    const F exact = s * t;
    const F near_quotient = b * y;
    const F near_midpoint = scalbn(b * y + b * z, -1);
    const F dividends[] = {exact, next_up(exact), next_down(exact),
                           near_quotient, near_midpoint};
    for (const F& dividend : dividends) {
      if (!CheckAllModes(Op::kDiv, Binary(dividend, t))) return;
      if (!CheckAllModes(Op::kDiv, Binary(dividend, b))) return;
    }

    const F square = s * s;
    const F radicands[] = {square, next_up(square), next_down(square), y * y,
                           y * z, scalbn(y * z, 1)};
    for (const F& radicand : radicands) {
      if (!CheckAllModes(Op::kSqrt, Binary(radicand, b))) return;
    }

    const F product = y * b;
    const int scale = static_cast<int>(
        random.Below(static_cast<limb_t>(F::kLimbs) * 2 * kLimbBits + 8));
    const F addends[] = {-product,
                         next_up(-product),
                         next_down(-product),
                         scalbn(product, -scale),
                         -scalbn(product, -scale),
                         scalbn(product, scale),
                         -scalbn(next_up(product), scale)};
    for (const F& addend : addends) {
      if (!CheckAllModes(Op::kFma, Args<F>{y, b, addend, 0})) return;
    }
    if (!CheckTinyTailFma<F>(random)) return;
  }
  std::printf("P=%d checks=%ld exact_division=%ld exact_root=%ld\n",
              F::kPrecision, coverage.checks, coverage.exact_division,
              coverage.exact_root);
  EXPECT_GT(coverage.exact_division, total / 2);
  EXPECT_GT(coverage.exact_root, total / 2);
}

constexpr int kIterations = 40000;

TEST(FloatOracle, RandomBinary128) { CheckRandomly<binary128>(kIterations, 1); }
TEST(FloatOracle, RandomBinary256) { CheckRandomly<binary256>(kIterations, 2); }
TEST(FloatOracle, RandomBinary512) { CheckRandomly<binary512>(kIterations, 3); }

TEST(FloatOracle, HardCasesBinary128) {
  CheckHardCases<binary128>(kIterations / 4, 21);
}
TEST(FloatOracle, HardCasesBinary256) {
  CheckHardCases<binary256>(kIterations / 4, 22);
}
TEST(FloatOracle, HardCasesBinary512) {
  CheckHardCases<binary512>(kIterations / 4, 23);
}
TEST(FloatOracle, HardCasesLimbBoundaries) {
  CheckHardCases<Float<64, 8>>(kIterations / 4, 24);
  CheckHardCases<Float<127, 15>>(kIterations / 4, 25);
  CheckHardCases<Float<128, 15>>(kIterations / 4, 26);
  CheckHardCases<Float<192, 15>>(kIterations / 4, 27);
  CheckHardCases<Float<256, 19>>(kIterations / 4, 28);
}

// Spare-bit counts of 0, 1 and 2 stress the guard logic, and they straddle limb
// boundaries, which the full-limb formats depend on.
TEST(FloatOracle, RandomLimbBoundaries) {
  CheckRandomly<Float<53, 11>>(kIterations, 4);
  CheckRandomly<Float<62, 8>>(kIterations, 5);
  CheckRandomly<Float<63, 8>>(kIterations, 6);
  CheckRandomly<Float<64, 8>>(kIterations, 7);
  CheckRandomly<Float<65, 8>>(kIterations, 8);
  CheckRandomly<Float<127, 15>>(kIterations, 9);
  CheckRandomly<Float<128, 15>>(kIterations, 10);
  CheckRandomly<Float<129, 15>>(kIterations, 11);
  CheckRandomly<Float<192, 15>>(kIterations, 12);
  CheckRandomly<Float<256, 19>>(kIterations, 13);
}

}  // namespace
}  // namespace mpfix
