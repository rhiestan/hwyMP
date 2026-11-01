// Copyright 2026 The hwyMP Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Differential tests of Float against MPFR, covering results and flags in all
// five rounding modes. Set HWYMP_TEST_SCALE to multiply the iteration counts
// and, when it is at least 10, to enable the larger exhaustive formats.
//
// MPFR runs with its wide default exponent range, and the expected IEEE result
// is derived from it here: rounding directly to the precision available at the
// result's exponent avoids any double rounding, and overflow and subnormals are
// resolved explicitly. MPFR has no ties-away mode, so ties are detected by
// rounding towards zero with one extra bit.

#include <gmp.h>
#include <gtest/gtest.h>
#include <mpfr.h>

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "hwymp/float.h"

namespace hwymp {
namespace {

enum class Op { kAdd, kSub, kMul };

constexpr Op kOps[] = {Op::kAdd, Op::kSub, Op::kMul};
constexpr Rounding kModes[] = {Rounding::kNearestEven, Rounding::kNearestAway,
                               Rounding::kTowardZero, Rounding::kUp,
                               Rounding::kDown};

const char* Name(Op op) {
  switch (op) {
    case Op::kAdd:
      return "add";
    case Op::kSub:
      return "sub";
    default:
      return "mul";
  }
}

// An mpfr_t that owns its storage.
class Mpfr {
 public:
  explicit Mpfr(long precision) { mpfr_init2(value_, precision); }
  ~Mpfr() { mpfr_clear(value_); }
  Mpfr(const Mpfr&) = delete;
  Mpfr& operator=(const Mpfr&) = delete;
  mpfr_ptr get() { return value_; }

 private:
  mpfr_t value_;
};

// Applies `op` to the operands, rounded to the precision of `result`. Returns
// MPFR's ternary value.
int Apply(Op op, mpfr_ptr result, mpfr_srcptr a, mpfr_srcptr b,
          mpfr_rnd_t rnd) {
  switch (op) {
    case Op::kAdd:
      return mpfr_add(result, a, b, rnd);
    case Op::kSub:
      return mpfr_sub(result, a, b, rnd);
    default:
      return mpfr_mul(result, a, b, rnd);
  }
}

// Rounds op(a, b) to the precision of `result` in any of the five modes.
int ApplyRounded(Op op, mpfr_ptr result, mpfr_srcptr a, mpfr_srcptr b,
                 Rounding mode) {
  switch (mode) {
    case Rounding::kNearestEven:
      return Apply(op, result, a, b, MPFR_RNDN);
    case Rounding::kTowardZero:
      return Apply(op, result, a, b, MPFR_RNDZ);
    case Rounding::kUp:
      return Apply(op, result, a, b, MPFR_RNDU);
    case Rounding::kDown:
      return Apply(op, result, a, b, MPFR_RNDD);
    case Rounding::kNearestAway:
      break;
  }
  const long precision = mpfr_get_prec(result);
  Mpfr probe(precision + 1);
  const int ternary = Apply(op, probe.get(), a, b, MPFR_RNDZ);
  const bool is_tie = ternary == 0 && mpfr_regular_p(probe.get()) &&
                      static_cast<long>(mpfr_min_prec(probe.get())) ==
                          precision + 1;
  return Apply(op, result, a, b, is_tie ? MPFR_RNDA : MPFR_RNDN);
}

// Exposes a Float's significand to MPFR without copying through strings.
template <class F>
class View {
 public:
  explicit View(const F& x) : limbs_(x.raw_limbs()) {
    const int sign = x.signbit() ? -1 : 1;
    int kind = MPFR_REGULAR_KIND;
    if (x.is_nan()) {
      kind = MPFR_NAN_KIND;
    } else if (x.is_inf()) {
      kind = MPFR_INF_KIND;
    } else if (x.is_zero()) {
      kind = MPFR_ZERO_KIND;
    }
    const mpfr_exp_t exponent =
        kind == MPFR_REGULAR_KIND ? x.raw_exponent() + mpfr_exp_t{1} : 0;
    mpfr_custom_init_set(value_, sign * kind, exponent, F::kPrecision,
                         limbs_.data());
  }
  mpfr_srcptr get() const { return value_; }

 private:
  typename F::Limbs limbs_;
  mpfr_t value_;
};

template <class F>
struct Expected {
  F value;
  unsigned flags;
};

// Converts a regular or zero MPFR value of at most F's precision.
template <class F>
F FromMpfr(mpfr_ptr x) {
  if (mpfr_zero_p(x)) return F::Zero(mpfr_signbit(x) != 0);
  Mpfr wide(F::kPrecision);
  mpfr_set(wide.get(), x, MPFR_RNDN);
  typename F::Limbs limbs;
  const auto* source =
      static_cast<const mp_limb_t*>(mpfr_custom_get_significand(wide.get()));
  for (std::size_t i = 0; i < F::kLimbs; ++i) limbs[i] = source[i];
  return F::FromRaw(mpfr_signbit(x) != 0,
                    static_cast<std::int32_t>(mpfr_get_exp(x) - 1), limbs);
}

template <class F>
F OverflowResult(bool negative, Rounding mode) {
  const bool to_infinity = mode == Rounding::kNearestEven ||
                           mode == Rounding::kNearestAway ||
                           (mode == Rounding::kUp && !negative) ||
                           (mode == Rounding::kDown && negative);
  return to_infinity ? F::Infinity(negative) : F::MaxFinite(negative);
}

template <class F>
Expected<F> Compute(Op op, const F& a, const F& b, Rounding mode) {
  if (a.is_nan() || b.is_nan()) {
    const bool signaling = a.is_signaling() || b.is_signaling();
    return {F::QuietNaN(), signaling ? unsigned{kFlagInvalid} : 0u};
  }

  const View<F> view_a(a);
  const View<F> view_b(b);
  const long precision = F::kPrecision;

  // Rounding towards zero gives the true exponent of the exact result.
  Mpfr truncated(precision);
  const int truncated_ternary =
      Apply(op, truncated.get(), view_a.get(), view_b.get(), MPFR_RNDZ);
  if (mpfr_nan_p(truncated.get())) {
    return {F::QuietNaN(), unsigned{kFlagInvalid}};
  }
  if (mpfr_inf_p(truncated.get())) {
    return {F::Infinity(mpfr_signbit(truncated.get()) != 0), 0u};
  }
  if (mpfr_zero_p(truncated.get())) {
    Mpfr exact(precision);
    ApplyRounded(op, exact.get(), view_a.get(), view_b.get(), mode);
    return {FromMpfr<F>(exact.get()), 0u};
  }

  const bool negative = mpfr_signbit(truncated.get()) != 0;
  const long exponent = static_cast<long>(mpfr_get_exp(truncated.get())) - 1;

  // Rounded with an unbounded exponent, which decides overflow and tininess.
  Mpfr unbounded(precision);
  const int unbounded_ternary =
      ApplyRounded(op, unbounded.get(), view_a.get(), view_b.get(), mode);
  const long unbounded_exponent =
      static_cast<long>(mpfr_get_exp(unbounded.get())) - 1;
  if (unbounded_exponent > F::kEmax) {
    return {OverflowResult<F>(negative, mode),
            unsigned{kFlagOverflow} | kFlagInexact};
  }

  if (exponent >= F::kEmin) {
    return {FromMpfr<F>(unbounded.get()),
            unbounded_ternary != 0 ? unsigned{kFlagInexact} : 0u};
  }

  const bool tiny = unbounded_exponent < F::kEmin;
  const unsigned inexact_flags =
      kFlagInexact | (tiny ? unsigned{kFlagUnderflow} : 0u);
  const long available = precision - (F::kEmin - exponent);

  if (available >= 1) {
    Mpfr rounded(available);
    const int ternary =
        ApplyRounded(op, rounded.get(), view_a.get(), view_b.get(), mode);
    return {FromMpfr<F>(rounded.get()), ternary != 0 ? inexact_flags : 0u};
  }

  // The magnitude is below the smallest subnormal, so the result is inexact
  // and is either zero or the smallest subnormal.
  const bool at_or_above_half = available == 0;
  const bool exactly_half = at_or_above_half && truncated_ternary == 0 &&
                            mpfr_min_prec(truncated.get()) == 1;
  bool round_up = false;
  switch (mode) {
    case Rounding::kNearestEven:
      round_up = at_or_above_half && !exactly_half;
      break;
    case Rounding::kNearestAway:
      round_up = at_or_above_half;
      break;
    case Rounding::kTowardZero:
      break;
    case Rounding::kUp:
      round_up = !negative;
      break;
    case Rounding::kDown:
      round_up = negative;
      break;
  }
  return {round_up ? F::MinSubnormal(negative) : F::Zero(negative),
          inexact_flags};
}

template <class F>
F ApplyFloat(Op op, const F& a, const F& b, Context& context) {
  switch (op) {
    case Op::kAdd:
      return add(a, b, context);
    case Op::kSub:
      return sub(a, b, context);
    default:
      return mul(a, b, context);
  }
}

template <class F>
F ApplyOperator(Op op, const F& a, const F& b) {
  switch (op) {
    case Op::kAdd:
      return a + b;
    case Op::kSub:
      return a - b;
    default:
      return a * b;
  }
}

template <class F>
std::string Hex(const F& x) {
  std::string text = x.is_nan() ? "nan:" : "";
  const typename F::Bits bits = x.to_bits();
  char buffer[17];
  for (std::size_t i = bits.size(); i-- > 0;) {
    std::snprintf(buffer, sizeof(buffer), "%016llx",
                  static_cast<unsigned long long>(bits[i]));
    text += buffer;
  }
  return text;
}

template <class F>
bool SameResult(const F& x, const F& y) {
  if (x.is_nan() || y.is_nan()) {
    return x.is_nan() && y.is_nan() && !x.is_signaling() && !y.is_signaling();
  }
  return x.to_bits() == y.to_bits();
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
};

Coverage coverage;

template <class F>
void Record(const F& a, const F& b, const Expected<F>& expected) {
  ++coverage.checks;
  const unsigned flags = expected.flags;
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
bool Check(Op op, const F& a, const F& b, Rounding mode) {
  Context context{mode};
  const F actual = ApplyFloat(op, a, b, context);
  const Expected<F> expected = Compute(op, a, b, mode);
  Record(a, b, expected);
  bool ok = SameResult(actual, expected.value) &&
            context.flags == expected.flags;
  if (mode == Rounding::kNearestEven) {
    ok = ok && SameResult(ApplyOperator(op, a, b), expected.value);
  }
  if (!ok) {
    ADD_FAILURE() << "P=" << F::kPrecision << " E=" << F::kExponentBits << ' '
                  << Name(op) << " mode=" << static_cast<int>(mode)
                  << "\n  a        = " << Hex(a) << "\n  b        = " << Hex(b)
                  << "\n  actual   = " << Hex(actual) << " flags "
                  << context.flags << "\n  expected = " << Hex(expected.value)
                  << " flags " << expected.flags;
  }
  return ok;
}

int TestScale() {
  const char* text = std::getenv("HWYMP_TEST_SCALE");
  const int scale = text != nullptr ? std::atoi(text) : 1;
  return scale > 0 ? scale : 1;
}

// Every operand pair, operation and rounding mode.
template <class F>
void CheckExhaustively() {
  const limb_t count = limb_t{1} << (F::kPrecision + F::kExponentBits);
  for (limb_t a_bits = 0; a_bits < count; ++a_bits) {
    const F a = F::from_bits({a_bits});
    for (limb_t b_bits = 0; b_bits < count; ++b_bits) {
      const F b = F::from_bits({b_bits});
      for (const Op op : kOps) {
        for (const Rounding mode : kModes) {
          if (!Check(op, a, b, mode)) return;
        }
      }
    }
  }
}

TEST(FloatOracle, ExhaustiveP3E2) { CheckExhaustively<Float<3, 2>>(); }
TEST(FloatOracle, ExhaustiveP3E3) { CheckExhaustively<Float<3, 3>>(); }
TEST(FloatOracle, ExhaustiveP4E3) { CheckExhaustively<Float<4, 3>>(); }
TEST(FloatOracle, ExhaustiveP5E2) { CheckExhaustively<Float<5, 2>>(); }

TEST(FloatOracle, ExhaustiveLargerFormats) {
  if (TestScale() < 10) GTEST_SKIP() << "Needs HWYMP_TEST_SCALE >= 10";
  CheckExhaustively<Float<6, 3>>();
  CheckExhaustively<Float<5, 4>>();
}

class Random {
 public:
  explicit Random(limb_t seed) : state_(seed) {}

  limb_t Next() {
    state_ += 0x9E3779B97F4A7C15u;
    limb_t z = state_;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9u;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBu;
    return z ^ (z >> 31);
  }

  limb_t Below(limb_t bound) { return Next() % bound; }

 private:
  limb_t state_;
};

// A significand limb array of random runs of ones and zeros, in the style of
// Berkeley TestFloat: such patterns exercise carries and rounding boundaries
// far better than uniform bits. The top bit is always set.
template <class F>
typename F::Limbs RandomSignificand(Random& random) {
  typename F::Limbs limbs{};
  constexpr int kTotal = static_cast<int>(F::kLimbs) * kLimbBits;
  constexpr int spare = kTotal - F::kPrecision;
  const auto set_bit = [&limbs](int position) {
    limbs[static_cast<std::size_t>(position) / kLimbBits] |=
        limb_t{1} << (position % kLimbBits);
  };
  if (random.Below(4) == 0) {
    for (limb_t& limb : limbs) limb = random.Next();
  } else {
    bool ones = random.Below(2) != 0;
    int position = kTotal - 1;
    while (position >= spare) {
      const int run = 1 + static_cast<int>(random.Below(
                              random.Below(2) != 0 ? 8 : F::kPrecision));
      for (int i = 0; i < run && position >= spare; ++i, --position) {
        if (ones) set_bit(position);
      }
      ones = !ones;
    }
  }
  if constexpr (spare != 0) limbs[0] &= ~((limb_t{1} << spare) - 1);
  limbs[F::kLimbs - 1] |= limb_t{1} << (kLimbBits - 1);
  return limbs;
}

// Rounds an arbitrary significand down onto the grid of subnormal exponents.
template <class F>
F MakeRegular(bool negative, std::int32_t exponent,
              typename F::Limbs limbs) {
  const std::int32_t lowest = F::kEmin - (F::kPrecision - 1);
  if (exponent < lowest) exponent = lowest;
  EXPECT_LE(exponent, F::kEmax);
  if (exponent < F::kEmin) {
    const int clear = F::kEmin - exponent;
    detail::ShiftRight(limbs, clear +
                                  static_cast<int>(F::kLimbs) * kLimbBits -
                                  F::kPrecision);
    detail::ShiftLeft(limbs, clear + static_cast<int>(F::kLimbs) * kLimbBits -
                                 F::kPrecision);
  }
  return F::FromRaw(negative, exponent, limbs);
}

template <class F>
std::int32_t RandomExponent(Random& random) {
  switch (random.Below(5)) {
    case 0:  // Close to one.
      return static_cast<std::int32_t>(random.Below(9)) - 4;
    case 1:  // Near overflow.
      return F::kEmax - static_cast<std::int32_t>(random.Below(
                            static_cast<limb_t>(F::kPrecision) + 4));
    case 2:  // Around the subnormal range.
      return F::kEmin + 4 -
             static_cast<std::int32_t>(
                 random.Below(static_cast<limb_t>(F::kPrecision) + 8));
    case 3:  // Half the range: products land near the extremes.
      return (random.Below(2) != 0 ? F::kEmax / 2 : F::kEmin / 2) + 4 -
             static_cast<std::int32_t>(
                 random.Below(static_cast<limb_t>(F::kPrecision) + 8));
    default: {
      // Covers every exponent from the smallest subnormal up to kEmax.
      const limb_t range = static_cast<limb_t>(F::kEmax) * 2 + F::kPrecision - 1;
      return F::kEmin - (F::kPrecision - 1) +
             static_cast<std::int32_t>(random.Below(range));
    }
  }
}

template <class F>
F RandomValue(Random& random) {
  const bool negative = random.Below(2) != 0;
  switch (random.Below(64)) {
    case 0:
      return F::Zero(negative);
    case 1:
      return F::Infinity(negative);
    case 2:
      return F::QuietNaN(negative);
    case 3:
      return F::SignalingNaN(negative);
    case 4:
      return F::MaxFinite(negative);
    case 5:
      return F::MinNormal(negative);
    case 6:
      return F::MinSubnormal(negative);
    default:
      return MakeRegular<F>(negative, RandomExponent<F>(random),
                            RandomSignificand<F>(random));
  }
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

template <class F>
void CheckRandomly(int iterations, limb_t seed) {
  Random random(seed);
  const int total = iterations * TestScale();
  coverage = {};
  for (int i = 0; i < total; ++i) {
    const F a = RandomValue<F>(random);
    const F b = RelatedValue<F>(random, a);
    for (const Op op : kOps) {
      for (const Rounding mode : kModes) {
        if (!Check(op, a, b, mode)) return;
      }
    }
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
  EXPECT_EQ(coverage.checks, 15L * total);
  EXPECT_GT(coverage.invalid, minimum);
  EXPECT_GT(coverage.overflow, minimum);
  EXPECT_GT(coverage.underflow, minimum);
  EXPECT_GT(coverage.inexact, minimum);
  EXPECT_GT(coverage.exact_regular, minimum);
  EXPECT_GT(coverage.exact_subnormal, minimum);
  EXPECT_GT(coverage.subnormal_results, minimum);
  EXPECT_GT(coverage.exact_cancellation, minimum);
}

constexpr int kIterations = 40000;

TEST(FloatOracle, RandomBinary128) { CheckRandomly<binary128>(kIterations, 1); }
TEST(FloatOracle, RandomBinary256) { CheckRandomly<binary256>(kIterations, 2); }
TEST(FloatOracle, RandomBinary512) { CheckRandomly<binary512>(kIterations, 3); }

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
}  // namespace hwymp
