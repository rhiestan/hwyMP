// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Boost.Multiprecision support: mpfix::boost_backend<P, E> is a backend for
// boost::multiprecision::number, so that Float can be used where code is
// written for Boost.Multiprecision types, as with cpp_bin_float or mpfr_float:
//
//   #include "mpfix/boost.h"
//
//   mpfix::boost_binary256 x = 2;
//   mpfix::boost_binary256 y = sqrt(x) + boost::multiprecision::erf(x);
//
// The aliases boost_binary128, boost_binary256 and boost_binary512 have
// expression templates off, as suits a fixed-size type; number<backend,
// et_on> works too. Arithmetic, sqrt, fma, the rounding functions, frexp,
// ldexp, fmod and remainder are exact or correctly rounded, as Float's are,
// and so are exp, exp2, log, log2, log10, the trigonometric and hyperbolic
// functions and their inverses, pow, atan2, hypot, cbrt, expm1, log1p, erf,
// erfc, tgamma and lgamma, which go to mpfix/math.h rather than to Boost's
// generic or Boost.Math implementations. Everything rounds to nearest.
//
// Boost.Multiprecision is not a dependency of mpfix: include this header only
// where Boost (1.79 or later, or the standalone Boost.Multiprecision) is
// available.

#ifndef MPFIX_BOOST_H_
#define MPFIX_BOOST_H_

#include <boost/multiprecision/number.hpp>
#include <cmath>
#include <cstddef>
#include <ios>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>

#include "mpfix/charconv.h"
#include "mpfix/float.h"
#include "mpfix/math.h"

namespace mpfix {

template <int P, int E>
class boost_backend {
 public:
  using value_type = Float<P, E>;
  // Boost converts the other built-in types to these first; all of them
  // convert exactly to the formats of 64 bits and more.
  using signed_types = std::tuple<long long>;
  using unsigned_types = std::tuple<unsigned long long>;
  using float_types = std::tuple<double, long double>;
  using exponent_type = int;

  boost_backend() noexcept : value_(value_type::Zero()) {}
  boost_backend(const value_type& value) noexcept  // NOLINT
      : value_(value) {}
  template <detail::Arithmetic T>
  boost_backend(T value) noexcept  // NOLINT
      : value_(static_cast<value_type>(value)) {}
  // From another format, explicit where it may round, as for Float.
  template <int P2, int E2>
    requires(P2 != P || E2 != E)
  explicit(!value_type::template kExactFrom<Float<P2, E2>>)
      boost_backend(const boost_backend<P2, E2>& other) noexcept
      : value_(static_cast<value_type>(other.value())) {}

  template <detail::Arithmetic T>
  boost_backend& operator=(T value) noexcept {
    value_ = static_cast<value_type>(value);
    return *this;
  }
  // Decimal or hexadecimal text, as from_string reads it, rounded to nearest.
  boost_backend& operator=(const char* text) {
    Context context;
    const value_type value = from_string<value_type>(text, context);
    if ((context.flags & kFlagInvalid) != 0) {
      throw std::runtime_error(std::string("mpfix: cannot parse \"") + text +
                               "\" as a number");
    }
    value_ = value;
    return *this;
  }

  void swap(boost_backend& other) noexcept { std::swap(value_, other.value_); }

  // As a stream with this precision and these flags would print the value.
  // With no digits, the shortest text that reads back.
  std::string str(std::streamsize digits, std::ios_base::fmtflags flags) const {
    const std::ios_base::fmtflags field = flags & std::ios_base::floatfield;
    if (digits <= 0 && field != std::ios_base::fixed) {
      std::chars_format fmt = std::chars_format::general;
      if (field == std::ios_base::scientific) {
        fmt = std::chars_format::scientific;
      }
      std::string text;
      if (field == (std::ios_base::fixed | std::ios_base::scientific)) {
        std::ostringstream out;
        out << std::hexfloat << value_;
        text = out.str();
      } else {
        text = detail::FormatShortest(value_, fmt,
                                      fmt == std::chars_format::general);
      }
      if ((flags & std::ios_base::showpos) != 0 && text[0] != '-') {
        text.insert(0, "+");
      }
      if ((flags & std::ios_base::uppercase) != 0) {
        for (char& c : text) {
          if (c >= 'a' && c <= 'z') c = static_cast<char>(c - 'a' + 'A');
        }
      }
      return text;
    }
    std::ostringstream out;
    out.flags(flags);
    out.precision(digits);
    out << value_;
    return out.str();
  }

  void negate() noexcept { value_ = -value_; }

  // Not called with NaNs: number compares those itself.
  int compare(const boost_backend& other) const noexcept {
    return value_ == other.value_ ? 0 : value_ < other.value_ ? -1 : 1;
  }
  // Compared exactly, in a format that holds both.
  template <detail::Arithmetic T>
  int compare(T other) const noexcept {
    using Common = Float<(P > 113 ? P : 113), (E > 15 ? E : 15)>;
    const Common a(value_);
    const Common b(other);
    return a == b ? 0 : a < b ? -1 : 1;
  }

  value_type& value() noexcept { return value_; }
  const value_type& value() const noexcept { return value_; }

 private:
  value_type value_;
};

using boost_binary128 =
    boost::multiprecision::number<boost_backend<113, 15>,
                                  boost::multiprecision::et_off>;
using boost_binary256 =
    boost::multiprecision::number<boost_backend<237, 19>,
                                  boost::multiprecision::et_off>;
using boost_binary512 =
    boost::multiprecision::number<boost_backend<489, 23>,
                                  boost::multiprecision::et_off>;

// --- The backend operations, found by argument-dependent lookup ---

template <int P, int E>
using BoostBackend = boost_backend<P, E>;

#define MPFIX_BOOST_BINARY(name, op)                                           \
  template <int P, int E>                                                      \
  void name(BoostBackend<P, E>& result, const BoostBackend<P, E>& a) noexcept { \
    result.value() = result.value() op a.value();                              \
  }                                                                            \
  template <int P, int E>                                                      \
  void name(BoostBackend<P, E>& result, const BoostBackend<P, E>& a,           \
            const BoostBackend<P, E>& b) noexcept {                            \
    result.value() = a.value() op b.value();                                   \
  }
MPFIX_BOOST_BINARY(eval_add, +)
MPFIX_BOOST_BINARY(eval_subtract, -)
MPFIX_BOOST_BINARY(eval_multiply, *)
MPFIX_BOOST_BINARY(eval_divide, /)
#undef MPFIX_BOOST_BINARY

// result = a b + c and result += a b, with one rounding.
template <int P, int E>
void eval_multiply_add(BoostBackend<P, E>& result, const BoostBackend<P, E>& a,
                       const BoostBackend<P, E>& b,
                       const BoostBackend<P, E>& c) noexcept {
  result.value() = fma(a.value(), b.value(), c.value());
}
template <int P, int E>
void eval_multiply_add(BoostBackend<P, E>& result, const BoostBackend<P, E>& a,
                       const BoostBackend<P, E>& b) noexcept {
  result.value() = fma(a.value(), b.value(), result.value());
}

template <class R, int P, int E>
  requires detail::Arithmetic<R>
void eval_convert_to(R* result, const BoostBackend<P, E>& value) noexcept {
  *result = static_cast<R>(value.value());
}

template <int P, int E>
bool eval_is_zero(const BoostBackend<P, E>& value) noexcept {
  return value.value().is_zero();
}
template <int P, int E>
int eval_get_sign(const BoostBackend<P, E>& value) noexcept {
  if (value.value().is_zero() || value.value().is_nan()) return 0;
  return value.value().signbit() ? -1 : 1;
}
template <int P, int E>
int eval_signbit(const BoostBackend<P, E>& value) noexcept {
  return value.value().signbit() ? 1 : 0;
}
template <int P, int E>
int eval_fpclassify(const BoostBackend<P, E>& value) noexcept {
  const Float<P, E>& x = value.value();
  if (x.is_nan()) return FP_NAN;
  if (x.is_inf()) return FP_INFINITE;
  if (x.is_zero()) return FP_ZERO;
  return x.is_subnormal() ? FP_SUBNORMAL : FP_NORMAL;
}

template <int P, int E, class I>
  requires std::is_integral_v<I>
void eval_frexp(BoostBackend<P, E>& result, const BoostBackend<P, E>& value,
                I* exponent) noexcept {
  int e = 0;
  result.value() = frexp(value.value(), &e);
  *exponent = static_cast<I>(e);
}

// Exponents beyond int's range over- or underflow as the largest ones do.
template <class I>
int ClampExponent(I n) noexcept {
  if constexpr (sizeof(I) > sizeof(int) || std::is_unsigned_v<I>) {
    constexpr I kHigh = static_cast<I>((std::numeric_limits<int>::max)());
    if (n > kHigh) return (std::numeric_limits<int>::max)();
    if constexpr (std::is_signed_v<I>) {
      constexpr I kLow = static_cast<I>((std::numeric_limits<int>::min)());
      if (n < kLow) return (std::numeric_limits<int>::min)();
    }
  }
  return static_cast<int>(n);
}

template <int P, int E, class I>
  requires std::is_integral_v<I>
void eval_ldexp(BoostBackend<P, E>& result, const BoostBackend<P, E>& value,
                I exponent) noexcept {
  result.value() = ldexp(value.value(), ClampExponent(exponent));
}
template <int P, int E, class I>
  requires std::is_integral_v<I>
void eval_scalbn(BoostBackend<P, E>& result, const BoostBackend<P, E>& value,
                 I exponent) noexcept {
  result.value() = scalbn(value.value(), ClampExponent(exponent));
}
template <int P, int E>
int eval_ilogb(const BoostBackend<P, E>& value) noexcept {
  return ilogb(value.value());
}

#define MPFIX_BOOST_UNARY(name, function)                                \
  template <int P, int E>                                                \
  void name(BoostBackend<P, E>& result, const BoostBackend<P, E>& value) { \
    result.value() = function(value.value());                            \
  }
MPFIX_BOOST_UNARY(eval_abs, abs)
MPFIX_BOOST_UNARY(eval_fabs, abs)
MPFIX_BOOST_UNARY(eval_sqrt, sqrt)
MPFIX_BOOST_UNARY(eval_floor, floor)
MPFIX_BOOST_UNARY(eval_ceil, ceil)
MPFIX_BOOST_UNARY(eval_trunc, trunc)
MPFIX_BOOST_UNARY(eval_round, round)
MPFIX_BOOST_UNARY(eval_logb, logb)
MPFIX_BOOST_UNARY(eval_exp, exp)
MPFIX_BOOST_UNARY(eval_exp2, exp2)
MPFIX_BOOST_UNARY(eval_log, log)
MPFIX_BOOST_UNARY(eval_log2, log2)
MPFIX_BOOST_UNARY(eval_log10, log10)
MPFIX_BOOST_UNARY(eval_sin, sin)
MPFIX_BOOST_UNARY(eval_cos, cos)
MPFIX_BOOST_UNARY(eval_tan, tan)
MPFIX_BOOST_UNARY(eval_asin, asin)
MPFIX_BOOST_UNARY(eval_acos, acos)
MPFIX_BOOST_UNARY(eval_atan, atan)
MPFIX_BOOST_UNARY(eval_sinh, sinh)
MPFIX_BOOST_UNARY(eval_cosh, cosh)
MPFIX_BOOST_UNARY(eval_tanh, tanh)
#undef MPFIX_BOOST_UNARY

// Round to nearest, ties to even, the only rounding mode here.
template <int P, int E>
void eval_rint(BoostBackend<P, E>& result, const BoostBackend<P, E>& value) {
  Context context{Rounding::kNearestEven};
  result.value() = round_to_integral(value.value(), context);
}
template <int P, int E>
void eval_nearbyint(BoostBackend<P, E>& result,
                    const BoostBackend<P, E>& value) {
  eval_rint(result, value);
}

#define MPFIX_BOOST_BINARY_FUNCTION(name, function)                        \
  template <int P, int E>                                                  \
  void name(BoostBackend<P, E>& result, const BoostBackend<P, E>& a,       \
            const BoostBackend<P, E>& b) {                                 \
    result.value() = function(a.value(), b.value());                       \
  }
MPFIX_BOOST_BINARY_FUNCTION(eval_pow, pow)
MPFIX_BOOST_BINARY_FUNCTION(eval_atan2, atan2)
MPFIX_BOOST_BINARY_FUNCTION(eval_hypot, hypot)
MPFIX_BOOST_BINARY_FUNCTION(eval_remainder, remainder)
#undef MPFIX_BOOST_BINARY_FUNCTION

// x - trunc(x / y) y, exactly, from the IEEE remainder: where that has the
// other sign than x, adding y with the sign of x gives the result, which is
// representable and so exact.
template <int P, int E>
void eval_fmod(BoostBackend<P, E>& result, const BoostBackend<P, E>& a,
               const BoostBackend<P, E>& b) {
  const Float<P, E>& x = a.value();
  const Float<P, E> r = remainder(x, b.value());
  if (r.is_nan() || r.is_zero()) {
    result.value() = r.is_zero() ? copysign(r, x) : r;
  } else if (r.signbit() == x.signbit()) {
    result.value() = r;
  } else {
    result.value() = r + copysign(b.value(), x);
  }
}

template <int P, int E>
std::size_t hash_value(const BoostBackend<P, E>& value) noexcept {
  // +0 and -0 compare equal, and so hash alike.
  const Float<P, E> x = value.value().is_zero() ? Float<P, E>::Zero()
                                                : value.value();
  std::size_t hash = 0;
  for (const limb_t limb : x.to_bits()) {
    hash ^= static_cast<std::size_t>(limb) + 0x9e3779b97f4a7c15u +
            (hash << 6) + (hash >> 2);
  }
  return hash;
}

}  // namespace mpfix

namespace boost {
namespace multiprecision {

template <int P, int E>
struct number_category<mpfix::boost_backend<P, E>>
    : public std::integral_constant<int, number_kind_floating_point> {};

template <int P, int E>
struct expression_template_default<mpfix::boost_backend<P, E>> {
  static constexpr expression_template_option value = et_off;
};

// The functions that Boost would otherwise compute with Boost.Math, from
// mpfix/math.h instead, correctly rounded.
#define MPFIX_BOOST_FUNCTION(name)                                       \
  template <int P, int E, expression_template_option ET>                 \
  number<mpfix::boost_backend<P, E>, ET> name BOOST_PREVENT_MACRO_SUBSTITUTION( \
      const number<mpfix::boost_backend<P, E>, ET>& x) {                 \
    return number<mpfix::boost_backend<P, E>, ET>(                       \
        mpfix::name(x.backend().value()));                               \
  }
MPFIX_BOOST_FUNCTION(asinh)
MPFIX_BOOST_FUNCTION(acosh)
MPFIX_BOOST_FUNCTION(atanh)
MPFIX_BOOST_FUNCTION(cbrt)
MPFIX_BOOST_FUNCTION(expm1)
MPFIX_BOOST_FUNCTION(log1p)
MPFIX_BOOST_FUNCTION(erf)
MPFIX_BOOST_FUNCTION(erfc)
MPFIX_BOOST_FUNCTION(tgamma)
MPFIX_BOOST_FUNCTION(lgamma)
#undef MPFIX_BOOST_FUNCTION

template <int P, int E, expression_template_option ET>
number<mpfix::boost_backend<P, E>, ET> nextafter BOOST_PREVENT_MACRO_SUBSTITUTION(
    const number<mpfix::boost_backend<P, E>, ET>& x,
    const number<mpfix::boost_backend<P, E>, ET>& y) {
  return number<mpfix::boost_backend<P, E>, ET>(
      nextafter(x.backend().value(), y.backend().value()));
}

}  // namespace multiprecision
}  // namespace boost

namespace std {

// The limits of Float<P, E>, as numbers.
template <int P, int E, boost::multiprecision::expression_template_option ET>
class numeric_limits<
    boost::multiprecision::number<mpfix::boost_backend<P, E>, ET>> {
  using L = numeric_limits<mpfix::Float<P, E>>;
  using N = boost::multiprecision::number<mpfix::boost_backend<P, E>, ET>;

 public:
  static constexpr bool is_specialized = true;
  static N(min)() noexcept { return N((L::min)()); }
  static N(max)() noexcept { return N((L::max)()); }
  static N lowest() noexcept { return N(L::lowest()); }
  static constexpr int digits = L::digits;
  static constexpr int digits10 = L::digits10;
  static constexpr int max_digits10 = L::max_digits10;
  static constexpr bool is_signed = true;
  static constexpr bool is_integer = false;
  static constexpr bool is_exact = false;
  static constexpr int radix = 2;
  static N epsilon() noexcept { return N(L::epsilon()); }
  static N round_error() noexcept { return N(L::round_error()); }
  static constexpr int min_exponent = L::min_exponent;
  static constexpr int min_exponent10 = L::min_exponent10;
  static constexpr int max_exponent = L::max_exponent;
  static constexpr int max_exponent10 = L::max_exponent10;
  static constexpr bool has_infinity = true;
  static constexpr bool has_quiet_NaN = true;
  static constexpr bool has_signaling_NaN = true;
  static constexpr float_denorm_style has_denorm = denorm_present;
  static constexpr bool has_denorm_loss = false;
  static N infinity() noexcept { return N(L::infinity()); }
  static N quiet_NaN() noexcept { return N(L::quiet_NaN()); }
  static N signaling_NaN() noexcept { return N(L::signaling_NaN()); }
  static N denorm_min() noexcept { return N(L::denorm_min()); }
  static constexpr bool is_iec559 = L::is_iec559;
  static constexpr bool is_bounded = true;
  static constexpr bool is_modulo = false;
  static constexpr bool traps = false;
  static constexpr bool tinyness_before = L::tinyness_before;
  static constexpr float_round_style round_style = round_to_nearest;
};

}  // namespace std

#endif  // MPFIX_BOOST_H_
