// Copyright 2026 The hwyMP Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// The MPFR oracle shared by the differential tests. MPFR runs with its wide
// default exponent range, and the expected IEEE result is derived from it
// here: rounding directly to the precision available at the result's exponent
// avoids any double rounding, and overflow and subnormals are resolved
// explicitly. MPFR has no ties-away mode, so ties are detected by rounding
// towards zero with one extra bit.

#ifndef HWYMP_TESTS_ORACLE_H_
#define HWYMP_TESTS_ORACLE_H_

#include <gmp.h>
#include <gtest/gtest.h>
#include <mpfr.h>

#include <array>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>

#include "hwymp/float.h"

namespace hwymp::test {

constexpr Rounding kModes[] = {Rounding::kNearestEven, Rounding::kNearestAway,
                               Rounding::kTowardZero, Rounding::kUp,
                               Rounding::kDown};

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

// Rounds an exact value to the precision of `result` in any of the five
// modes. `apply(result, rnd)` must set `result` to the exact value rounded in
// direction rnd and return MPFR's ternary value.
template <class Apply>
int ApplyRounded(const Apply& apply, mpfr_ptr result, Rounding mode) {
  switch (mode) {
    case Rounding::kNearestEven:
      return apply(result, MPFR_RNDN);
    case Rounding::kTowardZero:
      return apply(result, MPFR_RNDZ);
    case Rounding::kUp:
      return apply(result, MPFR_RNDU);
    case Rounding::kDown:
      return apply(result, MPFR_RNDD);
    case Rounding::kNearestAway:
      break;
  }
  const long precision = mpfr_get_prec(result);
  Mpfr probe(precision + 1);
  const int ternary = apply(probe.get(), MPFR_RNDZ);
  const bool is_tie =
      ternary == 0 && mpfr_regular_p(probe.get()) &&
      static_cast<long>(mpfr_min_prec(probe.get())) == precision + 1;
  return apply(result, is_tie ? MPFR_RNDA : MPFR_RNDN);
}

// The IEEE result in F, with its flags, of an exact value given by `apply` as
// for ApplyRounded. A NaN result counts as an invalid operation.
template <class F, class Apply>
Expected<F> RoundExact(const Apply& apply, Rounding mode) {
  const long precision = F::kPrecision;

  // Rounding towards zero gives the true exponent of the exact result.
  Mpfr truncated(precision);
  const int truncated_ternary = apply(truncated.get(), MPFR_RNDZ);
  if (mpfr_nan_p(truncated.get())) {
    return {F::QuietNaN(), unsigned{kFlagInvalid}};
  }
  if (mpfr_inf_p(truncated.get())) {
    return {F::Infinity(mpfr_signbit(truncated.get()) != 0), 0u};
  }
  if (mpfr_zero_p(truncated.get()) && truncated_ternary != 0) {
    // Below MPFR's own exponent range, and so far below the smallest
    // subnormal of F: it rounds like any such tiny value.
    const bool negative = mpfr_signbit(truncated.get()) != 0;
    const bool round_up = (mode == Rounding::kUp && !negative) ||
                          (mode == Rounding::kDown && negative);
    return {round_up ? F::MinSubnormal(negative) : F::Zero(negative),
            unsigned{kFlagUnderflow} | kFlagInexact};
  }
  if (mpfr_zero_p(truncated.get())) {
    Mpfr exact(precision);
    ApplyRounded(apply, exact.get(), mode);
    return {FromMpfr<F>(exact.get()), 0u};
  }

  const bool negative = mpfr_signbit(truncated.get()) != 0;
  const long exponent = static_cast<long>(mpfr_get_exp(truncated.get())) - 1;

  // Rounded with an unbounded exponent, which decides overflow and tininess.
  Mpfr unbounded(precision);
  const int unbounded_ternary = ApplyRounded(apply, unbounded.get(), mode);
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
    const int ternary = ApplyRounded(apply, rounded.get(), mode);
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

int TestScale() {
  const char* text = std::getenv("HWYMP_TEST_SCALE");
  const int scale = text != nullptr ? std::atoi(text) : 1;
  return scale > 0 ? scale : 1;
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
F MakeRegular(bool negative, std::int32_t exponent, typename F::Limbs limbs) {
  const std::int32_t lowest = F::kEmin - (F::kPrecision - 1);
  if (exponent < lowest) exponent = lowest;
  EXPECT_LE(exponent, F::kEmax);
  if (exponent < F::kEmin) {
    const int clear = F::kEmin - exponent;
    detail::ShiftRight(
        limbs, clear + static_cast<int>(F::kLimbs) * kLimbBits - F::kPrecision);
    detail::ShiftLeft(
        limbs, clear + static_cast<int>(F::kLimbs) * kLimbBits - F::kPrecision);
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
      const limb_t range =
          static_cast<limb_t>(F::kEmax) * 2 + F::kPrecision - 1;
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

}  // namespace hwymp::test

#endif  // HWYMP_TESTS_ORACLE_H_
