// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Complex numbers over Float: mpfix::Complex<F>, with the interface of
// std::complex and the special values of C99 Annex G for infinities and NaNs.
// std::complex<T> is unspecified for T other than the built-in floating-point
// types, and some standard libraries compute its functions in double.
//
// Addition, subtraction and multiplication are correctly rounded in each
// component, and so are abs (hypot), arg (atan2), norm, conj and proj.
// Division is within half an ulp plus 2^-(3P) ulp in each component. The
// functions sqrt, exp, log, log10, sin, cos, tan, sinh, cosh and tanh are
// evaluated with the correctly rounded real functions in P + 64 bits and
// formulas without cancellation, then rounded once: within half an ulp plus
// about 2^-60 ulp in each component, and correctly rounded where they reduce
// to real functions, as on the real axis. pow is exp(w log z), whose error is
// relative to |z^w| rather than to each component; integer powers up to 64
// multiply instead. Everything rounds to nearest; there is no Context.

#ifndef MPFIX_COMPLEX_H_
#define MPFIX_COMPLEX_H_

#include <algorithm>
#include <bit>
#include <ostream>

#include "mpfix/charconv.h"
#include "mpfix/float.h"
#include "mpfix/math.h"

namespace mpfix {

template <class F>
class Complex;

namespace detail {

template <class F>
constexpr int kComplexExponentBits = (std::min)(F::kExponentBits + 2, 29);

// Products of two components, exactly.
template <class F>
using ComplexExact = Float<2 * F::kPrecision, kComplexExponentBits<F>>;

// Sums of two such products, truncated, with a sticky bit.
template <class F>
using ComplexSum = Float<4 * F::kPrecision + 8, kComplexExponentBits<F>>;

// The format of the functions: 64 bits more, and exponents wide enough for
// F's range squared and for 2^(P + 72), which Kahan's tanh reaches, within the
// width the real functions support.
template <class F>
using ComplexWide =
    Float<F::kPrecision + 64,
          (std::min)(26, (std::max)(F::kExponentBits + 2,
                                    std::bit_width(static_cast<unsigned>(
                                        F::kPrecision + 72)) + 2))>;

// x y + s z, for s = +-1, correctly rounded to F. The products are exact, and
// their sum, truncated towards zero in 4P + 8 bits, carries a sticky bit.
template <class F>
F SumOfProducts(const F& x, const F& y, bool subtract, const F& w,
                const F& z) {
  using X = ComplexExact<F>;
  using S = ComplexSum<F>;
  const S first(X(x) * X(y));
  const S second(X(w) * X(z));
  Context context{Rounding::kTowardZero};
  const S sum = subtract ? sub(first, second, context)
                         : add(first, second, context);
  if (sum.is_zero() || !sum.is_finite()) return F(sum);
  return RoundTruncated<F>(sum.signbit(), abs(sum),
                           (context.flags & kFlagInexact) != 0, NearestEnv{});
}

}  // namespace detail

// A complex number with components of the Float type F.
template <class F>
class Complex {
 public:
  using value_type = F;

  Complex() noexcept : re_(F::Zero()), im_(F::Zero()) {}
  Complex(const F& re, const F& im = F::Zero()) noexcept  // NOLINT
      : re_(re), im_(im) {}

  F real() const noexcept { return re_; }
  F imag() const noexcept { return im_; }
  void real(const F& value) noexcept { re_ = value; }
  void imag(const F& value) noexcept { im_ = value; }

  Complex& operator+=(const Complex& z) noexcept {
    re_ = re_ + z.re_;
    im_ = im_ + z.im_;
    return *this;
  }
  Complex& operator-=(const Complex& z) noexcept {
    re_ = re_ - z.re_;
    im_ = im_ - z.im_;
    return *this;
  }
  Complex& operator*=(const Complex& z) { return *this = *this * z; }
  Complex& operator/=(const Complex& z) { return *this = *this / z; }
  Complex& operator+=(const F& x) noexcept {
    re_ = re_ + x;
    return *this;
  }
  Complex& operator-=(const F& x) noexcept {
    re_ = re_ - x;
    return *this;
  }
  Complex& operator*=(const F& x) noexcept {
    re_ = re_ * x;
    im_ = im_ * x;
    return *this;
  }
  Complex& operator/=(const F& x) noexcept {
    re_ = re_ / x;
    im_ = im_ / x;
    return *this;
  }

  friend Complex operator+(const Complex& z) noexcept { return z; }
  friend Complex operator-(const Complex& z) noexcept {
    return {-z.re_, -z.im_};
  }
  friend Complex operator+(Complex a, const Complex& b) noexcept {
    return a += b;
  }
  friend Complex operator-(Complex a, const Complex& b) noexcept {
    return a -= b;
  }
  friend Complex operator+(Complex a, const F& b) noexcept { return a += b; }
  friend Complex operator-(Complex a, const F& b) noexcept { return a -= b; }
  friend Complex operator*(Complex a, const F& b) noexcept { return a *= b; }
  friend Complex operator/(Complex a, const F& b) noexcept { return a /= b; }
  friend Complex operator+(const F& a, const Complex& b) noexcept {
    return {a + b.re_, b.im_};
  }
  friend Complex operator-(const F& a, const Complex& b) noexcept {
    return {a - b.re_, -b.im_};
  }
  friend Complex operator*(const F& a, const Complex& b) noexcept {
    return {a * b.re_, a * b.im_};
  }
  friend Complex operator/(const F& a, const Complex& b) {
    return Complex(a) / b;
  }

  // (a + bi)(c + di) = (ac - bd) + (ad + bc)i, each part rounded once, with
  // Annex G's recovery of infinities from NaN results.
  friend Complex operator*(const Complex& z, const Complex& w) {
    F a = z.re_, b = z.im_, c = w.re_, d = w.im_;
    if (a.is_finite() && b.is_finite() && c.is_finite() && d.is_finite()) {
      return {detail::SumOfProducts(a, c, true, b, d),
              detail::SumOfProducts(a, d, false, b, c)};
    }
    const F ac = a * c, bd = b * d, ad = a * d, bc = b * c;
    F x = ac - bd;
    F y = ad + bc;
    if (x.is_nan() && y.is_nan()) {
      const auto unit = [](const F& v) {
        return copysign(v.is_inf() ? F(1) : F::Zero(), v);
      };
      const auto zero_nan = [](F& v) {
        if (v.is_nan()) v = copysign(F::Zero(), v);
      };
      bool recalc = false;
      if (a.is_inf() || b.is_inf()) {
        a = unit(a);
        b = unit(b);
        zero_nan(c);
        zero_nan(d);
        recalc = true;
      }
      if (c.is_inf() || d.is_inf()) {
        c = unit(c);
        d = unit(d);
        zero_nan(a);
        zero_nan(b);
        recalc = true;
      }
      if (!recalc &&
          (ac.is_inf() || bd.is_inf() || ad.is_inf() || bc.is_inf())) {
        zero_nan(a);
        zero_nan(b);
        zero_nan(c);
        zero_nan(d);
        recalc = true;
      }
      if (recalc) {
        x = F::Infinity() * (a * c - b * d);
        y = F::Infinity() * (a * d + b * c);
      }
    }
    return {x, y};
  }

  // (a + bi)/(c + di) = ((ac + bd) + (bc - ad)i) / (c^2 + d^2), with c and d
  // scaled by 2^-k for k = logb(max(|c|, |d|)), as in Annex G, and the sums
  // truncated in 4P + 8 bits: each part within half an ulp plus 2^-(3P).
  friend Complex operator/(const Complex& z, const Complex& w) {
    F a = z.re_, b = z.im_, c = w.re_, d = w.im_;
    const F largest = (std::max)(abs(c), abs(d));
    if (a.is_finite() && b.is_finite() && largest.is_finite() &&
        !largest.is_zero()) {
      using S = detail::ComplexSum<F>;
      using X = detail::ComplexExact<F>;
      // Scaled in the exact format, whose exponent range keeps the smaller
      // part from underflowing.
      const int k = largest.raw_exponent();
      const X xa(a), xb(b);
      const X cs = scalbn(X(c), -k);
      const X ds = scalbn(X(d), -k);
      // Numerators and denominator, truncated: the products are exact.
      Context context{Rounding::kTowardZero};
      const auto sum = [&context](const X& p, const X& q, bool subtract,
                                  const X& r, const X& t) {
        const S first(p * q);
        const S second(r * t);
        return subtract ? sub(first, second, context)
                        : add(first, second, context);
      };
      const S denominator = sum(cs, cs, false, ds, ds);
      const S x = scalbn(sum(xa, cs, false, xb, ds) / denominator, -k);
      const S y = scalbn(sum(xb, cs, true, xa, ds) / denominator, -k);
      return {F(x), F(y)};
    }
    // Annex G.5.2.
    int ilogbw = 0;
    if (largest.is_finite() && !largest.is_zero()) {
      ilogbw = largest.raw_exponent();
      c = scalbn(c, -ilogbw);
      d = scalbn(d, -ilogbw);
    }
    const F denom = c * c + d * d;
    F x = scalbn((a * c + b * d) / denom, -ilogbw);
    F y = scalbn((b * c - a * d) / denom, -ilogbw);
    if (x.is_nan() && y.is_nan()) {
      const auto unit = [](const F& v) {
        return copysign(v.is_inf() ? F(1) : F::Zero(), v);
      };
      if (denom.is_zero() && (!a.is_nan() || !b.is_nan())) {
        x = copysign(F::Infinity(), c) * a;
        y = copysign(F::Infinity(), c) * b;
      } else if ((a.is_inf() || b.is_inf()) && c.is_finite() &&
                 d.is_finite()) {
        a = unit(a);
        b = unit(b);
        x = F::Infinity() * (a * c + b * d);
        y = F::Infinity() * (b * c - a * d);
      } else if (largest.is_inf() && a.is_finite() && b.is_finite()) {
        c = unit(c);
        d = unit(d);
        x = F::Zero() * (a * c + b * d);
        y = F::Zero() * (b * c - a * d);
      }
    }
    return {x, y};
  }

  friend bool operator==(const Complex& a, const Complex& b) noexcept {
    return a.re_ == b.re_ && a.im_ == b.im_;
  }
  friend bool operator==(const Complex& a, const F& b) noexcept {
    return a.re_ == b && a.im_ == F::Zero();
  }

  friend std::ostream& operator<<(std::ostream& out, const Complex& z) {
    return out << '(' << z.re_ << ',' << z.im_ << ')';
  }

 private:
  F re_;
  F im_;
};

template <class F>
F real(const Complex<F>& z) noexcept {
  return z.real();
}
template <class F>
F imag(const Complex<F>& z) noexcept {
  return z.imag();
}

// |z|, correctly rounded.
template <class F>
F abs(const Complex<F>& z) {
  return hypot(z.real(), z.imag());
}

// The argument in [-pi, pi], correctly rounded, with atan2's signed zeros.
template <class F>
F arg(const Complex<F>& z) {
  return atan2(z.imag(), z.real());
}

// |z|^2, correctly rounded.
template <class F>
F norm(const Complex<F>& z) {
  if (z.real().is_inf() || z.imag().is_inf()) return F::Infinity();
  if (!z.real().is_finite() || !z.imag().is_finite()) return z.real() + z.imag();
  return detail::SumOfProducts(z.real(), z.real(), false, z.imag(), z.imag());
}

template <class F>
Complex<F> conj(const Complex<F>& z) noexcept {
  return {z.real(), -z.imag()};
}

// The projection onto the Riemann sphere: infinities map to +infinity with
// the sign of the imaginary part kept on zero.
template <class F>
Complex<F> proj(const Complex<F>& z) noexcept {
  if (z.real().is_inf() || z.imag().is_inf()) {
    return {F::Infinity(), copysign(F::Zero(), z.imag())};
  }
  return z;
}

namespace detail {

template <class F>
using CW = ComplexWide<F>;

template <class F>
CW<F> Wide(const F& x) {
  return CW<F>(x);
}

// sinh(x + iy) for finite x and y, as sinh x cos y + i cosh x sin y.
template <class F>
Complex<F> SinhFinite(const F& x, const F& y) {
  if (y.is_zero()) return {sinh(x), y};
  if (x.is_zero()) return {x * copysign(F(1), cos(y)), sin(y)};
  const CW<F> wx = Wide(x), wy = Wide(y);
  return {F(sinh(wx) * cos(wy)), F(cosh(wx) * sin(wy))};
}

// cosh(x + iy) for finite x and y, as cosh x cos y + i sinh x sin y.
template <class F>
Complex<F> CoshFinite(const F& x, const F& y) {
  // The zero imaginary parts carry the sign of sinh x times that of sin y.
  if (y.is_zero()) return {cosh(x), y * copysign(F(1), x)};
  if (x.is_zero()) {
    const F s = sin(y);
    return {cos(y), x * copysign(F(1), s)};
  }
  const CW<F> wx = Wide(x), wy = Wide(y);
  return {F(cosh(wx) * cos(wy)), F(sinh(wx) * sin(wy))};
}

// tanh(x + iy) for finite x and y by Kahan's formula: with t = tan y,
// b = 1 + t^2, s = sinh x and r = sqrt(1 + s^2) = cosh x,
// tanh(x + iy) = (b r s + i t) / (1 + b s^2), which cancels nothing; beyond
// |x| = (P + 70) ln 2 / 2, the real part is +-1 to that precision and the
// imaginary part 4 t e^(-2|x|) / b, or sin(2y) 2 e^(-2|x|).
template <class F>
Complex<F> TanhFinite(const F& x, const F& y) {
  if (y.is_zero()) return {tanh(x), y};
  if (x.is_zero()) return {x, tan(y)};
  using W = CW<F>;
  const W wx = Wide(x), wy = Wide(y);
  const W one(1);
  if (abs(wx) > W((W::kPrecision + 6) * 0.3466)) {
    const W e = exp(scalbn(-abs(wx), 1));
    return {copysign(F(1), x), F(scalbn(sin(scalbn(wy, 1)) * e, 1))};
  }
  const W t = tan(wy);
  const W b = one + t * t;
  const W s = sinh(wx);
  const W r = sqrt(one + s * s);
  const W denominator = one + b * s * s;
  return {F(b * r * s / denominator), F(t / denominator)};
}

// The Annex G special values of sinh, cosh and tanh, for x + iy with a part
// infinite or NaN and x >= +0 or NaN; the callers apply the symmetries.
enum class HyperbolicKind { kSinh, kCosh, kTanh };

template <HyperbolicKind kKind, class F>
Complex<F> HyperbolicSpecial(const F& x, const F& y) {
  const F nan = F::QuietNaN();
  const F inf = F::Infinity();
  if (x.is_inf()) {
    if constexpr (kKind == HyperbolicKind::kTanh) {
      // 1 + i0 sin(2y), and 1 +- i0 for y infinite or NaN.
      if (!y.is_finite() || y.is_zero()) return {F(1), copysign(F::Zero(), y)};
      return {F(1), copysign(F::Zero(), sin(y) * cos(y))};
    } else {
      if (y.is_zero()) return {inf, y};
      if (y.is_finite()) return {inf * cos(y), inf * sin(y)};
      return {inf, nan};
    }
  }
  if (x.is_zero()) {
    // y is infinite or NaN here.
    if constexpr (kKind == HyperbolicKind::kCosh) return {nan, F::Zero()};
    return {x, nan};
  }
  if (x.is_nan() && y.is_zero()) return {nan, y};
  return {nan, nan};
}

}  // namespace detail

template <class F>
Complex<F> sinh(const Complex<F>& z) {
  const F x = z.real(), y = z.imag();
  if (x.is_finite() && y.is_finite()) return detail::SinhFinite(x, y);
  // sinh(-z) = -sinh(z).
  const bool flip = x.signbit() && !x.is_nan();
  const Complex<F> w = detail::HyperbolicSpecial<detail::HyperbolicKind::kSinh>(
      flip ? -x : x, flip ? -y : y);
  return flip ? -w : w;
}

template <class F>
Complex<F> cosh(const Complex<F>& z) {
  const F x = z.real(), y = z.imag();
  if (x.is_finite() && y.is_finite()) return detail::CoshFinite(x, y);
  // cosh(-z) = cosh(z).
  const bool flip = x.signbit() && !x.is_nan();
  return detail::HyperbolicSpecial<detail::HyperbolicKind::kCosh>(
      flip ? -x : x, flip ? -y : y);
}

template <class F>
Complex<F> tanh(const Complex<F>& z) {
  const F x = z.real(), y = z.imag();
  if (x.is_finite() && y.is_finite()) return detail::TanhFinite(x, y);
  // tanh(-z) = -tanh(z).
  const bool flip = x.signbit() && !x.is_nan();
  const Complex<F> w = detail::HyperbolicSpecial<detail::HyperbolicKind::kTanh>(
      flip ? -x : x, flip ? -y : y);
  return flip ? -w : w;
}

// sin z = -i sinh(iz), cos z = cosh(iz) and tan z = -i tanh(iz).
template <class F>
Complex<F> sin(const Complex<F>& z) {
  const Complex<F> w = sinh(Complex<F>(-z.imag(), z.real()));
  return {w.imag(), -w.real()};
}
template <class F>
Complex<F> cos(const Complex<F>& z) {
  return cosh(Complex<F>(-z.imag(), z.real()));
}
template <class F>
Complex<F> tan(const Complex<F>& z) {
  const Complex<F> w = tanh(Complex<F>(-z.imag(), z.real()));
  return {w.imag(), -w.real()};
}

// e^z = e^x (cos y + i sin y).
template <class F>
Complex<F> exp(const Complex<F>& z) {
  const F x = z.real(), y = z.imag();
  const F nan = F::QuietNaN();
  if (x.is_finite() && y.is_finite()) {
    if (y.is_zero()) return {exp(x), y};
    using W = detail::CW<F>;
    const W e = exp(W(x));
    const W wy(y);
    return {F(e * cos(wy)), F(e * sin(wy))};
  }
  if (x.is_inf()) {
    if (x.signbit()) {
      if (y.is_finite()) {
        return {F::Zero() * cos(y), F::Zero() * sin(y)};
      }
      return {F::Zero(), copysign(F::Zero(), y)};
    }
    if (y.is_zero()) return {x, y};
    if (y.is_finite()) return {x * cos(y), x * sin(y)};
    return {x, nan};
  }
  if (x.is_nan() && y.is_zero()) return {x, y};
  return {nan, nan};
}

// log z = log |z| + i arg z. Near |z| = 1, log |z| = log1p(x^2 + y^2 - 1) / 2
// with x^2 + y^2 - 1 computed without rounding first.
template <class F>
Complex<F> log(const Complex<F>& z) {
  const F x = z.real(), y = z.imag();
  const F angle = atan2(y, x);
  if (x.is_inf() || y.is_inf()) return {F::Infinity(), angle};
  if (x.is_nan() || y.is_nan()) return {F::QuietNaN(), F::QuietNaN()};
  if (x.is_zero() && y.is_zero()) return {-F::Infinity(), angle};
  using W = detail::CW<F>;
  const F ax = abs(x), ay = abs(y);
  const F largest = (std::max)(ax, ay);
  if (largest >= F(0.5) && largest <= F(2)) {
    // x^2 - 1 is exact in the wide sum format, and so the whole sum up to its
    // truncation.
    using X = detail::ComplexExact<F>;
    using S = detail::ComplexSum<F>;
    Context context{Rounding::kTowardZero};
    const S s = add(sub(S(X(x) * X(x)), S(1), context), S(X(y) * X(y)),
                    context);
    if (s.is_zero()) return {F::Zero(), angle};
    return {F(scalbn(log1p(W(s)), -1)), angle};
  }
  return {F(log(hypot(W(x), W(y)))), angle};
}

template <class F>
Complex<F> log10(const Complex<F>& z) {
  using W = detail::CW<F>;
  const F x = z.real(), y = z.imag();
  if (!x.is_finite() || !y.is_finite() || (x.is_zero() && y.is_zero())) {
    const Complex<F> w = log(z);
    return {w.real(), w.imag() / F(2.302585092994045684)};
  }
  const W inverse = W(1) / log(W(10));
  const F ax = abs(x), ay = abs(y);
  const F largest = (std::max)(ax, ay);
  W re;
  if (largest >= F(0.5) && largest <= F(2)) {
    using X = detail::ComplexExact<F>;
    using S = detail::ComplexSum<F>;
    Context context{Rounding::kTowardZero};
    const S s = add(sub(S(X(x) * X(x)), S(1), context), S(X(y) * X(y)),
                    context);
    re = s.is_zero() ? W::Zero() : scalbn(log1p(W(s)), -1);
  } else {
    re = log(hypot(W(x), W(y)));
  }
  return {F(re * inverse), F(atan2(W(y), W(x)) * inverse)};
}

// The principal square root, by Kahan's method: t = sqrt((|x| + |z|) / 2),
// and sqrt z = t + i y / (2t) for x >= 0, |y| / (2t) + i copysign(t, y)
// otherwise.
template <class F>
Complex<F> sqrt(const Complex<F>& z) {
  const F x = z.real(), y = z.imag();
  const F inf = F::Infinity();
  if (y.is_inf()) return {inf, y};
  if (x.is_inf()) {
    if (x.signbit()) {
      if (y.is_nan()) return {y, copysign(inf, y)};
      return {F::Zero(), copysign(inf, y)};
    }
    return {x, y.is_nan() ? y : copysign(F::Zero(), y)};
  }
  if (x.is_nan() || y.is_nan()) return {F::QuietNaN(), F::QuietNaN()};
  if (x.is_zero() && y.is_zero()) return {F::Zero(), y};
  if (y.is_zero()) {
    if (!x.signbit()) return {sqrt(x), y};
    return {F::Zero(), copysign(sqrt(-x), y)};
  }
  using W = detail::CW<F>;
  const W wx(x), wy(y);
  const W t = sqrt(scalbn(abs(wx) + hypot(wx, wy), -1));
  const W other = wy / scalbn(t, 1);
  if (!x.signbit()) return {F(t), F(other)};
  return {F(abs(other)), F(copysign(t, wy))};
}

namespace detail {

// z^n for an integer |n| <= 64 by squaring, in the wide format.
template <class F>
Complex<F> IntegerPower(const Complex<F>& z, int n) {
  using W = CW<F>;
  W re(1), im(0);
  W base_re(z.real()), base_im(z.imag());
  const bool invert = n < 0;
  unsigned m = static_cast<unsigned>(invert ? -n : n);
  while (m != 0) {
    if ((m & 1) != 0) {
      const W r = re * base_re - im * base_im;
      im = re * base_im + im * base_re;
      re = r;
    }
    m >>= 1;
    if (m != 0) {
      const W r = base_re * base_re - base_im * base_im;
      base_im = scalbn(base_re * base_im, 1);
      base_re = r;
    }
  }
  if (!invert) return {F(re), F(im)};
  const W denominator = re * re + im * im;
  return {F(re / denominator), F(-im / denominator)};
}

}  // namespace detail

// z^w = e^(w log z), with z^0 = 1 and 0^w = 0 for Re w > 0 and Im w = 0.
template <class F>
Complex<F> pow(const Complex<F>& z, const Complex<F>& w) {
  if (w.real().is_zero() && w.imag().is_zero()) return {F(1), F::Zero()};
  if (z.real().is_zero() && z.imag().is_zero() && w.imag().is_zero() &&
      w.real() > F::Zero()) {
    return {F::Zero(), F::Zero()};
  }
  const F n = w.real();
  if (w.imag().is_zero() && n.is_finite() && n == trunc(n) && abs(n) <= F(64) &&
      z.real().is_finite() && z.imag().is_finite()) {
    return detail::IntegerPower(z, static_cast<int>(n));
  }
  using W = detail::CW<F>;
  const W x(z.real()), y(z.imag());
  if (!z.real().is_finite() || !z.imag().is_finite() || !w.real().is_finite() ||
      !w.imag().is_finite()) {
    return exp(w * log(z));
  }
  // log z and the product in the wide format.
  const W lr = scalbn(log(x * x + y * y), -1);
  const W li = atan2(y, x);
  const W wr(w.real()), wi(w.imag());
  const W tr = wr * lr - wi * li;
  const W ti = wr * li + wi * lr;
  const W e = exp(tr);
  return {F(e * cos(ti)), F(e * sin(ti))};
}

template <class F>
Complex<F> pow(const Complex<F>& z, const F& w) {
  return pow(z, Complex<F>(w));
}
template <class F>
Complex<F> pow(const F& z, const Complex<F>& w) {
  return pow(Complex<F>(z), w);
}

// r (cos theta + i sin theta), each part rounded once from the wide format.
template <class F>
Complex<F> polar(const F& r, const F& theta = F::Zero()) {
  if (theta.is_zero()) return {r, r * theta};
  using W = detail::CW<F>;
  const W wr(r), wt(theta);
  return {F(wr * cos(wt)), F(wr * sin(wt))};
}

using complex128 = Complex<binary128>;
using complex256 = Complex<binary256>;
using complex512 = Complex<binary512>;

}  // namespace mpfix

#endif  // MPFIX_COMPLEX_H_
