// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Exact unsigned integers of any size, for the string conversions. They are
// only as fast as the conversions need: schoolbook multiplication with
// Karatsuba above a threshold, and long division after Knuth, The Art of
// Computer Programming, volume 2, section 4.3.1, algorithm D.

#ifndef MPFIX_BIGINT_H_
#define MPFIX_BIGINT_H_

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

#include "mpfix/limb_ops.h"

namespace mpfix::detail {

inline constexpr limb_t kTenToThe19 = 10000000000000000000u;

// Divides u1:u0 by a normalized d, given v = floor((2^128 - 1) / d) - 2^64,
// and requires u1 < d. N. Möller and T. Granlund, "Improved division by
// invariant integers", IEEE Transactions on Computers 60(2), 2011, algorithm 4.
MPFIX_ALWAYS_INLINE limb_t DivPreinverted(limb_t u1, limb_t u0, limb_t d,
                                          limb_t v, limb_t& remainder) {
  limb_t q1;
  limb_t q0 = mul_wide(v, u1, q1);
  carry_t carry = 0;
  q0 = adc(q0, u0, carry);
  q1 = adc(q1, u1, carry);
  ++q1;
  limb_t r = u0 - q1 * d;
  if (r > q0) {
    --q1;
    r += d;
  }
  if (r >= d) [[unlikely]] {
    ++q1;
    r -= d;
  }
  remainder = r;
  return q1;
}

// Divides the n-limb integer at `limbs`, least significant first, by 10^19 in
// place and returns the remainder. 10^19 has its top bit set, so it needs no
// normalization, and its reciprocal is a constant.
inline limb_t DivideLimbsBy10To19(limb_t* limbs, std::size_t n) {
  constexpr limb_t kInverse =
      0xd83c94fb6d2ac34au;  // (2^128 - 1) / 10^19 - 2^64
  limb_t remainder = 0;
  for (std::size_t i = n; i-- > 0;) {
    limbs[i] =
        DivPreinverted(remainder, limbs[i], kTenToThe19, kInverse, remainder);
  }
  return remainder;
}

// Divides the n-limb integer at `limbs`, least significant first, by a nonzero
// divisor in place and returns the remainder.
inline limb_t DivideLimbs(limb_t* limbs, std::size_t n, limb_t divisor) {
  // Shifting dividend and divisor alike leaves the quotient unchanged and
  // gives the division the normalized divisor it needs. The remainder of the
  // shifted division is the true one times 2^shift.
  const int shift = countl_zero(divisor);
  const limb_t d = divisor << shift;
  limb_t unused;
  const limb_t inverse = div_2by1(~d, ~limb_t{0}, d, unused);
  limb_t remainder = 0;
  if (shift != 0 && n != 0) remainder = limbs[n - 1] >> (kLimbBits - shift);
  for (std::size_t i = n; i-- > 0;) {
    const limb_t below =
        shift != 0 && i > 0 ? limbs[i - 1] >> (kLimbBits - shift) : 0;
    const limb_t digit = (limbs[i] << shift) | below;
    // The partial remainder stays below d, so each quotient limb fits.
    limbs[i] = DivPreinverted(remainder, digit, d, inverse, remainder);
  }
  return remainder >> shift;
}

// A vector of limbs that keeps up to kInline of them in place, so that the
// numbers of everyday conversions never touch the heap.
class LimbVector {
 public:
  static constexpr std::size_t kInline = 16;

  LimbVector() noexcept = default;
  LimbVector(const LimbVector& other) { CopyFrom(other); }
  LimbVector(LimbVector&& other) noexcept { MoveFrom(other); }
  LimbVector& operator=(const LimbVector& other) {
    if (this != &other) {
      size_ = 0;
      CopyFrom(other);
    }
    return *this;
  }
  LimbVector& operator=(LimbVector&& other) noexcept {
    if (this != &other) {
      heap_.reset();
      data_ = inline_;
      capacity_ = kInline;
      MoveFrom(other);
    }
    return *this;
  }

  std::size_t size() const noexcept { return size_; }
  bool empty() const noexcept { return size_ == 0; }
  limb_t* data() noexcept { return data_; }
  const limb_t* data() const noexcept { return data_; }
  limb_t* begin() noexcept { return data_; }
  limb_t* end() noexcept { return data_ + size_; }
  const limb_t* begin() const noexcept { return data_; }
  const limb_t* end() const noexcept { return data_ + size_; }
  limb_t& operator[](std::size_t i) noexcept { return data_[i]; }
  limb_t operator[](std::size_t i) const noexcept { return data_[i]; }
  limb_t& back() noexcept { return data_[size_ - 1]; }
  limb_t back() const noexcept { return data_[size_ - 1]; }

  void clear() noexcept { size_ = 0; }
  void pop_back() noexcept { --size_; }
  void push_back(limb_t value) {
    Reserve(size_ + 1);
    data_[size_++] = value;
  }
  // New limbs are set to `value`.
  void resize(std::size_t size, limb_t value = 0) {
    Reserve(size);
    for (std::size_t i = size_; i < size; ++i) data_[i] = value;
    size_ = size;
  }
  void assign(const limb_t* first, const limb_t* last) {
    const auto count = static_cast<std::size_t>(last - first);
    size_ = 0;
    Reserve(count);
    if (count != 0) std::memcpy(data_, first, count * sizeof(limb_t));
    size_ = count;
  }
  // Inserts `count` zero limbs at the bottom.
  void PrependZeros(std::size_t count) {
    Reserve(size_ + count);
    std::memmove(data_ + count, data_, size_ * sizeof(limb_t));
    std::fill(data_, data_ + count, limb_t{0});
    size_ += count;
  }
  // Removes the lowest `count` <= size() limbs.
  void DropLow(std::size_t count) noexcept {
    std::memmove(data_, data_ + count, (size_ - count) * sizeof(limb_t));
    size_ -= count;
  }

 private:
  void Reserve(std::size_t capacity) {
    if (capacity <= capacity_) return;
    const std::size_t grown = (std::max)(capacity, 2 * capacity_);
    std::unique_ptr<limb_t[]> storage(new limb_t[grown]);
    if (size_ != 0) std::memcpy(storage.get(), data_, size_ * sizeof(limb_t));
    heap_ = std::move(storage);
    data_ = heap_.get();
    capacity_ = grown;
  }

  void CopyFrom(const LimbVector& other) {
    Reserve(other.size_);
    if (other.size_ != 0) {
      std::memcpy(data_, other.data_, other.size_ * sizeof(limb_t));
    }
    size_ = other.size_;
  }

  // Requires this to be empty and in place.
  void MoveFrom(LimbVector& other) noexcept {
    if (other.heap_) {
      heap_ = std::move(other.heap_);
      data_ = heap_.get();
      capacity_ = other.capacity_;
    } else if (other.size_ != 0) {
      std::memcpy(inline_, other.inline_, other.size_ * sizeof(limb_t));
    }
    size_ = other.size_;
    other.data_ = other.inline_;
    other.capacity_ = kInline;
    other.size_ = 0;
  }

  limb_t inline_[kInline];
  limb_t* data_ = inline_;
  std::size_t size_ = 0;
  std::size_t capacity_ = kInline;
  std::unique_ptr<limb_t[]> heap_;
};

class BigUint {
 public:
  BigUint() = default;
  explicit BigUint(limb_t value) {
    if (value != 0) limbs_.push_back(value);
  }

  // From `count` limbs, least significant first.
  static BigUint FromLimbs(const limb_t* limbs, std::size_t count) {
    BigUint result;
    result.limbs_.assign(limbs, limbs + count);
    result.Trim();
    return result;
  }

  bool is_zero() const noexcept { return limbs_.empty(); }
  // Least significant first, with no zero limbs at the top.
  const LimbVector& limbs() const noexcept { return limbs_; }

  std::size_t bit_length() const noexcept {
    if (limbs_.empty()) return 0;
    return limbs_.size() * kLimbBits -
           static_cast<std::size_t>(countl_zero(limbs_.back()));
  }

  // The number of zero bits below the lowest set bit. Requires a nonzero value.
  std::size_t trailing_zeros() const noexcept {
    std::size_t i = 0;
    while (limbs_[i] == 0) ++i;
    return i * kLimbBits + static_cast<std::size_t>(countr_zero(limbs_[i]));
  }

  bool bit(std::size_t index) const noexcept {
    const std::size_t word = index / kLimbBits;
    return word < limbs_.size() && ((limbs_[word] >> (index % kLimbBits)) & 1);
  }

  // Whether any of the bits below `index` is set.
  bool any_bit_below(std::size_t index) const noexcept {
    const std::size_t word = index / kLimbBits;
    for (std::size_t i = 0; i < (std::min)(word, limbs_.size()); ++i) {
      if (limbs_[i] != 0) return true;
    }
    if (word >= limbs_.size()) return false;
    const int bits = static_cast<int>(index % kLimbBits);
    return bits != 0 && (limbs_[word] << (kLimbBits - bits)) != 0;
  }

  // *this = *this * factor + addend.
  void MulAdd(limb_t factor, limb_t addend) {
    limb_t carry = addend;
    for (limb_t& limb : limbs_) {
      limb_t hi;
      limb_t lo = mul_wide(limb, factor, hi);
      carry_t c = 0;
      lo = adc(lo, carry, c);
      limb = lo;
      carry = hi + c;
    }
    if (carry != 0) limbs_.push_back(carry);
    Trim();
  }

  // Divides by `divisor` in place and returns the remainder.
  limb_t DivBy10To19() {
    const limb_t remainder = DivideLimbsBy10To19(limbs_.data(), limbs_.size());
    Trim();
    return remainder;
  }

  limb_t DivSmall(limb_t divisor) {
    const limb_t remainder = DivideLimbs(limbs_.data(), limbs_.size(), divisor);
    Trim();
    return remainder;
  }

  void ShiftLeft(std::size_t count) {
    if (limbs_.empty() || count == 0) return;
    const std::size_t words = count / kLimbBits;
    const int bits = static_cast<int>(count % kLimbBits);
    if (bits != 0) {
      limb_t carry = 0;
      for (limb_t& limb : limbs_) {
        const limb_t next = limb >> (kLimbBits - bits);
        limb = (limb << bits) | carry;
        carry = next;
      }
      if (carry != 0) limbs_.push_back(carry);
    }
    limbs_.PrependZeros(words);
  }

  void ShiftRight(std::size_t count) {
    const std::size_t words = count / kLimbBits;
    if (words >= limbs_.size()) {
      limbs_.clear();
      return;
    }
    limbs_.DropLow(words);
    const int bits = static_cast<int>(count % kLimbBits);
    if (bits != 0) {
      for (std::size_t i = 0; i < limbs_.size(); ++i) {
        const limb_t above = i + 1 < limbs_.size() ? limbs_[i + 1] : 0;
        limbs_[i] = (limbs_[i] >> bits) | (above << (kLimbBits - bits));
      }
    }
    Trim();
  }

  void Add(const BigUint& other) {
    if (other.limbs_.size() > limbs_.size()) {
      limbs_.resize(other.limbs_.size(), 0);
    }
    carry_t carry = 0;
    for (std::size_t i = 0; i < limbs_.size(); ++i) {
      const limb_t addend = i < other.limbs_.size() ? other.limbs_[i] : 0;
      if (i >= other.limbs_.size() && carry == 0) break;
      limbs_[i] = adc(limbs_[i], addend, carry);
    }
    if (carry != 0) limbs_.push_back(1);
  }

  // Requires *this >= other.
  void Sub(const BigUint& other) {
    carry_t borrow = 0;
    for (std::size_t i = 0; i < limbs_.size(); ++i) {
      const limb_t subtrahend = i < other.limbs_.size() ? other.limbs_[i] : 0;
      if (i >= other.limbs_.size() && borrow == 0) break;
      limbs_[i] = sbb(limbs_[i], subtrahend, borrow);
    }
    Trim();
  }

  friend int Compare(const BigUint& a, const BigUint& b) noexcept {
    if (a.limbs_.size() != b.limbs_.size()) {
      return a.limbs_.size() < b.limbs_.size() ? -1 : 1;
    }
    for (std::size_t i = a.limbs_.size(); i-- > 0;) {
      if (a.limbs_[i] != b.limbs_[i]) return a.limbs_[i] < b.limbs_[i] ? -1 : 1;
    }
    return 0;
  }

  friend BigUint operator*(const BigUint& a, const BigUint& b) {
    BigUint result;
    if (a.is_zero() || b.is_zero()) return result;
    result.limbs_.resize(a.limbs_.size() + b.limbs_.size());
    Mul(a.limbs_.data(), a.limbs_.size(), b.limbs_.data(), b.limbs_.size(),
        result.limbs_.data());
    result.Trim();
    return result;
  }

  static BigUint Pow5(std::uint64_t exponent) {
    // 5^27 is the largest power of five in a limb.
    constexpr limb_t kPow5_27 = 7450580596923828125u;
    BigUint result(1);
    if (exponent >= 27) {
      // Square-and-multiply on 5^27, then the remaining small factor.
      const std::uint64_t chunks = exponent / 27;
      for (int bit = 63 - countl_zero(chunks); bit >= 0; --bit) {
        result = result * result;
        if ((chunks >> bit) & 1) result.MulAdd(kPow5_27, 0);
      }
    }
    limb_t factor = 1;
    for (std::uint64_t i = 0; i < exponent % 27; ++i) factor *= 5;
    result.MulAdd(factor, 0);
    return result;
  }

  // quotient = numerator / denominator and remainder = numerator % denominator.
  // Requires a nonzero denominator.
  static void DivRem(const BigUint& numerator, const BigUint& denominator,
                     BigUint& quotient, BigUint& remainder) {
    if (Compare(numerator, denominator) < 0) {
      quotient = BigUint();
      remainder = numerator;
      return;
    }
    if (denominator.limbs_.size() == 1) {
      quotient = numerator;
      remainder = BigUint(quotient.DivSmall(denominator.limbs_[0]));
      return;
    }

    // Normalize so that the divisor's top bit is set. The extra top limb of
    // u keeps the first quotient step like the others.
    const int shift = countl_zero(denominator.limbs_.back());
    BigUint v = denominator;
    v.ShiftLeft(static_cast<std::size_t>(shift));
    LimbVector u = numerator.limbs_;
    u.push_back(0);
    if (shift != 0) {
      for (std::size_t i = u.size() - 1; i > 0; --i) {
        u[i] = (u[i] << shift) | (u[i - 1] >> (kLimbBits - shift));
      }
      u[0] <<= shift;
    }

    const std::size_t n = v.limbs_.size();
    const std::size_t m = u.size() - n;
    const limb_t d1 = v.limbs_[n - 1];
    const limb_t d0 = v.limbs_[n - 2];
    const limb_t inverse = reciprocal_3by2(d1, d0);
    quotient.limbs_.clear();
    quotient.limbs_.resize(m, 0);

    for (std::size_t j = m; j-- > 0;) {
      // The top n limbs of the current window are below v, so the estimate
      // from the top three limbs is at most one too large.
      const limb_t u2 = u[j + n];
      const limb_t u1 = u[j + n - 1];
      const limb_t u0 = u[j + n - 2];
      limb_t q;
      if (u2 == d1 && u1 == d0) {
        q = ~limb_t{0};
      } else {
        limb_t r1;
        limb_t r0;
        q = div_3by2(u2, u1, u0, d1, d0, inverse, r1, r0);
      }

      // u[j .. j + n] -= q * v.
      limb_t carry = 0;
      carry_t borrow = 0;
      for (std::size_t i = 0; i < n; ++i) {
        limb_t hi;
        limb_t lo = mul_wide(q, v.limbs_[i], hi);
        carry_t c = 0;
        lo = adc(lo, carry, c);
        hi += c;
        u[j + i] = sbb(u[j + i], lo, borrow);
        carry = hi;
      }
      u[j + n] = sbb(u[j + n], carry, borrow);

      // Adding back corrects an estimate that was too large.
      while (borrow != 0) {
        --q;
        carry_t c = 0;
        for (std::size_t i = 0; i < n; ++i) {
          u[j + i] = adc(u[j + i], v.limbs_[i], c);
        }
        u[j + n] = adc(u[j + n], 0, c);
        // The carry out cancels the earlier borrow.
        borrow = c == 0 ? 1 : 0;
      }
      quotient.limbs_[j] = q;
    }
    quotient.Trim();

    u.resize(n);
    remainder.limbs_ = std::move(u);
    remainder.Trim();
    remainder.ShiftRight(static_cast<std::size_t>(shift));
  }

 private:
  static constexpr std::size_t kKaratsubaThreshold = 32;

  void Trim() noexcept {
    while (!limbs_.empty() && limbs_.back() == 0) limbs_.pop_back();
  }

  // r[0 .. an + bn) = a * b, for an, bn >= 1.
  static void MulBasecase(const limb_t* a, std::size_t an, const limb_t* b,
                          std::size_t bn, limb_t* r) {
    std::fill(r, r + an + bn, limb_t{0});
    for (std::size_t i = 0; i < an; ++i) {
      limb_t carry = 0;
      for (std::size_t j = 0; j < bn; ++j) {
        limb_t hi;
        limb_t lo = mul_wide(a[i], b[j], hi);
        carry_t c = 0;
        lo = adc(lo, r[i + j], c);
        hi += c;
        c = 0;
        lo = adc(lo, carry, c);
        hi += c;
        r[i + j] = lo;
        carry = hi;
      }
      r[i + bn] = carry;
    }
  }

  // Adds x[0 .. n) into r, propagating the carry up to r_end.
  static void AddInto(limb_t* r, const limb_t* r_end, const limb_t* x,
                      std::size_t n) {
    carry_t carry = 0;
    std::size_t i = 0;
    for (; i < n && r + i < r_end; ++i) r[i] = adc(r[i], x[i], carry);
    for (; carry != 0 && r + i < r_end; ++i) r[i] = adc(r[i], 0, carry);
  }

  // x[0 .. n) -= y[0 .. yn), for x >= y.
  static void SubFrom(limb_t* x, std::size_t n, const limb_t* y,
                      std::size_t yn) {
    carry_t borrow = 0;
    std::size_t i = 0;
    for (; i < yn; ++i) x[i] = sbb(x[i], y[i], borrow);
    for (; borrow != 0 && i < n; ++i) x[i] = sbb(x[i], 0, borrow);
  }

  // r[0 .. an + bn) = a * b.
  static void Mul(const limb_t* a, std::size_t an, const limb_t* b,
                  std::size_t bn, limb_t* r) {
    if (an < bn) {
      std::swap(a, b);
      std::swap(an, bn);
    }
    if (bn < kKaratsubaThreshold) {
      MulBasecase(a, an, b, bn, r);
      return;
    }
    if (an == bn) {
      Karatsuba(a, b, an, r);
      return;
    }
    // Unbalanced: multiply b by bn-limb slices of a.
    std::fill(r, r + an + bn, limb_t{0});
    std::vector<limb_t> partial(2 * bn);
    for (std::size_t offset = 0; offset < an; offset += bn) {
      const std::size_t length = (std::min)(bn, an - offset);
      Mul(a + offset, length, b, bn, partial.data());
      AddInto(r + offset, r + an + bn, partial.data(), length + bn);
    }
  }

  // r[0 .. 2n) = a * b for n-limb a and b.
  static void Karatsuba(const limb_t* a, const limb_t* b, std::size_t n,
                        limb_t* r) {
    const std::size_t low = n / 2;
    const std::size_t high = n - low;
    // z0 = a0 * b0 and z2 = a1 * b1 go straight into r.
    Mul(a, low, b, low, r);
    Mul(a + low, high, b + low, high, r + 2 * low);

    // (a0 + a1) * (b0 + b1) - z0 - z2 is the middle term.
    std::vector<limb_t> sum_a(high + 1, 0);
    std::vector<limb_t> sum_b(high + 1, 0);
    std::copy(a + low, a + n, sum_a.begin());
    std::copy(b + low, b + n, sum_b.begin());
    AddInto(sum_a.data(), sum_a.data() + high + 1, a, low);
    AddInto(sum_b.data(), sum_b.data() + high + 1, b, low);
    std::vector<limb_t> middle(2 * high + 2);
    Mul(sum_a.data(), high + 1, sum_b.data(), high + 1, middle.data());
    SubFrom(middle.data(), middle.size(), r, 2 * low);
    SubFrom(middle.data(), middle.size(), r + 2 * low, 2 * high);
    // The top limbs of the middle term are zero beyond the product's size.
    AddInto(r + low, r + 2 * n, middle.data(), middle.size());
  }

  LimbVector limbs_;
};

}  // namespace mpfix::detail

#endif  // MPFIX_BIGINT_H_
