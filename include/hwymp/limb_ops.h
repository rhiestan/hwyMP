// Copyright 2026 The hwyMP Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Lowest-level 64-bit limb primitives: add/subtract with carry, widening
// multiply, 128-by-64 division, shifts and bit scans.
//
// Every operation has a portable constexpr implementation in
// hwymp::detail::portable. The public functions forward to compiler intrinsics
// or inline assembly at run time. __int128 is deliberately not used here,
// because MSVC does not provide it.

#ifndef HWYMP_LIMB_OPS_H_
#define HWYMP_LIMB_OPS_H_

#include <bit>
#include <cstdint>
#include <type_traits>

#if defined(_M_X64) && !defined(_M_ARM64EC)
#define HWYMP_ARCH_X64 1
#elif defined(__x86_64__)
#define HWYMP_ARCH_X64 1
#else
#define HWYMP_ARCH_X64 0
#endif

#if defined(_M_ARM64) || defined(__aarch64__)
#define HWYMP_ARCH_ARM64 1
#else
#define HWYMP_ARCH_ARM64 0
#endif

// GNU-style inline assembly; also true for Clang targeting the MSVC ABI.
#if defined(__GNUC__) || defined(__clang__)
#define HWYMP_GNU_ASM 1
#else
#define HWYMP_GNU_ASM 0
#endif

#if defined(_MSC_VER)
#include <intrin.h>
#elif HWYMP_ARCH_X64
#include <immintrin.h>
#endif

#if defined(_MSC_VER) && !defined(__clang__)
#define HWYMP_ALWAYS_INLINE __forceinline
#else
#define HWYMP_ALWAYS_INLINE inline __attribute__((always_inline))
#endif

#if defined(__has_builtin)
#define HWYMP_HAS_BUILTIN(x) __has_builtin(x)
#else
#define HWYMP_HAS_BUILTIN(x) 0
#endif

namespace hwymp {

using limb_t = std::uint64_t;
// Carry or borrow flag, always 0 or 1. It matches the intrinsics' type so that
// chained adc calls compile to adc instructions without intermediate setc.
using carry_t = unsigned char;

inline constexpr int kLimbBits = 64;

namespace detail::portable {

constexpr limb_t adc(limb_t a, limb_t b, carry_t& carry) noexcept {
  const limb_t partial = a + b;
  const limb_t sum = partial + carry;
  carry = static_cast<carry_t>((partial < a) | (sum < partial));
  return sum;
}

constexpr limb_t sbb(limb_t a, limb_t b, carry_t& borrow) noexcept {
  const limb_t partial = a - b;
  const limb_t difference = partial - borrow;
  borrow = static_cast<carry_t>((a < b) | (partial < borrow));
  return difference;
}

constexpr limb_t mul_wide(limb_t a, limb_t b, limb_t& hi) noexcept {
  constexpr limb_t kMask = 0xFFFFFFFFu;
  const limb_t a_lo = a & kMask;
  const limb_t a_hi = a >> 32;
  const limb_t b_lo = b & kMask;
  const limb_t b_hi = b >> 32;

  const limb_t lo_lo = a_lo * b_lo;
  const limb_t hi_lo = a_hi * b_lo;
  const limb_t lo_hi = a_lo * b_hi;
  const limb_t hi_hi = a_hi * b_hi;

  const limb_t middle = (lo_lo >> 32) + (hi_lo & kMask) + lo_hi;
  hi = hi_hi + (hi_lo >> 32) + (middle >> 32);
  return (middle << 32) | (lo_lo & kMask);
}

constexpr int countl_zero(limb_t x) noexcept { return std::countl_zero(x); }

// Divides hi:lo by d using two 64-by-32 estimation steps with correction
// (Knuth, TAOCP Vol. 2, Algorithm D, in base 2^32). Requires hi < d.
constexpr limb_t div_2by1(limb_t hi, limb_t lo, limb_t d,
                          limb_t& remainder) noexcept {
  constexpr limb_t kBase = limb_t{1} << 32;
  constexpr limb_t kMask = kBase - 1;

  const int shift = countl_zero(d);
  const limb_t norm_d = d << shift;
  const limb_t norm_hi = shift == 0 ? hi : (hi << shift) | (lo >> (64 - shift));
  const limb_t norm_lo = lo << shift;

  const limb_t d1 = norm_d >> 32;
  const limb_t d0 = norm_d & kMask;
  const limb_t lo1 = norm_lo >> 32;
  const limb_t lo0 = norm_lo & kMask;

  limb_t q1 = norm_hi / d1;
  limb_t r = norm_hi - q1 * d1;
  while (q1 >= kBase || q1 * d0 > ((r << 32) | lo1)) {
    --q1;
    r += d1;
    if (r >= kBase) break;
  }

  const limb_t middle = ((norm_hi << 32) | lo1) - q1 * norm_d;
  limb_t q0 = middle / d1;
  r = middle - q0 * d1;
  while (q0 >= kBase || q0 * d0 > ((r << 32) | lo0)) {
    --q0;
    r += d1;
    if (r >= kBase) break;
  }

  remainder = (((middle << 32) | lo0) - q0 * norm_d) >> shift;
  return (q1 << 32) | q0;
}

}  // namespace detail::portable

namespace detail {

#if HWYMP_GNU_ASM && HWYMP_ARCH_X64 && !defined(_MSC_VER) && !defined(__BMI2__)
HWYMP_ALWAYS_INLINE limb_t mul_wide_asm(limb_t a, limb_t b,
                                        limb_t& hi) noexcept {
  limb_t lo;
  __asm__("mulq %3" : "=a"(lo), "=d"(hi) : "%0"(a), "r"(b) : "cc");
  return lo;
}
#elif HWYMP_GNU_ASM && HWYMP_ARCH_ARM64
HWYMP_ALWAYS_INLINE limb_t mul_wide_asm(limb_t a, limb_t b,
                                        limb_t& hi) noexcept {
  __asm__("umulh %0, %1, %2" : "=r"(hi) : "r"(a), "r"(b));
  return a * b;
}
#endif

#if HWYMP_GNU_ASM && HWYMP_ARCH_X64
HWYMP_ALWAYS_INLINE limb_t div_2by1_asm(limb_t hi, limb_t lo, limb_t d,
                                        limb_t& remainder) noexcept {
  limb_t quotient;
  __asm__("divq %4"
          : "=a"(quotient), "=d"(remainder)
          : "0"(lo), "1"(hi), "r"(d)
          : "cc");
  return quotient;
}
#endif

}  // namespace detail

// Returns a + b + carry; `carry` is both the incoming and outgoing carry.
HWYMP_ALWAYS_INLINE constexpr limb_t adc(limb_t a, limb_t b,
                                         carry_t& carry) noexcept {
  if (std::is_constant_evaluated()) return detail::portable::adc(a, b, carry);
#if HWYMP_ARCH_X64
  unsigned long long sum;
  carry = _addcarry_u64(carry, a, b, &sum);
  return sum;
#elif HWYMP_HAS_BUILTIN(__builtin_addcll)
  unsigned long long carry_out;
  const limb_t sum = __builtin_addcll(a, b, carry, &carry_out);
  carry = static_cast<carry_t>(carry_out);
  return sum;
#else
  return detail::portable::adc(a, b, carry);
#endif
}

// Returns a - b - borrow; `borrow` is both the incoming and outgoing borrow.
HWYMP_ALWAYS_INLINE constexpr limb_t sbb(limb_t a, limb_t b,
                                         carry_t& borrow) noexcept {
  if (std::is_constant_evaluated()) return detail::portable::sbb(a, b, borrow);
#if HWYMP_ARCH_X64
  unsigned long long difference;
  borrow = _subborrow_u64(borrow, a, b, &difference);
  return difference;
#elif HWYMP_HAS_BUILTIN(__builtin_subcll)
  unsigned long long borrow_out;
  const limb_t difference = __builtin_subcll(a, b, borrow, &borrow_out);
  borrow = static_cast<carry_t>(borrow_out);
  return difference;
#else
  return detail::portable::sbb(a, b, borrow);
#endif
}

// Returns the low 64 bits of a * b and stores the high 64 bits in `hi`.
HWYMP_ALWAYS_INLINE constexpr limb_t mul_wide(limb_t a, limb_t b,
                                              limb_t& hi) noexcept {
  if (std::is_constant_evaluated()) return detail::portable::mul_wide(a, b, hi);
#if defined(_MSC_VER) && HWYMP_ARCH_X64
  unsigned long long high;
  const limb_t low = _umul128(a, b, &high);
  hi = high;
  return low;
#elif HWYMP_ARCH_X64 && defined(__BMI2__)
  unsigned long long high;
  const limb_t low = _mulx_u64(a, b, &high);
  hi = high;
  return low;
#elif defined(_MSC_VER) && HWYMP_ARCH_ARM64
  hi = __umulh(a, b);
  return a * b;
#elif HWYMP_GNU_ASM && (HWYMP_ARCH_X64 || HWYMP_ARCH_ARM64)
  return detail::mul_wide_asm(a, b, hi);
#else
  return detail::portable::mul_wide(a, b, hi);
#endif
}

// Returns the number of leading zero bits; 64 when x == 0.
HWYMP_ALWAYS_INLINE constexpr int countl_zero(limb_t x) noexcept {
#if defined(_MSC_VER) && !defined(__clang__) && HWYMP_ARCH_X64
  // std::countl_zero adds a run-time CPU feature check on MSVC.
  if (!std::is_constant_evaluated()) {
    unsigned long index;
    return _BitScanReverse64(&index, x) ? 63 - static_cast<int>(index) : 64;
  }
#endif
  return std::countl_zero(x);
}

// Returns the number of trailing zero bits; 64 when x == 0.
HWYMP_ALWAYS_INLINE constexpr int countr_zero(limb_t x) noexcept {
#if defined(_MSC_VER) && !defined(__clang__) && HWYMP_ARCH_X64
  if (!std::is_constant_evaluated()) {
    unsigned long index;
    return _BitScanForward64(&index, x) ? static_cast<int>(index) : 64;
  }
#endif
  return std::countr_zero(x);
}

// Returns the high limb of (hi:lo) << count. Requires 0 < count < 64.
HWYMP_ALWAYS_INLINE constexpr limb_t shift_left_2(limb_t hi, limb_t lo,
                                                  int count) noexcept {
  return (hi << count) | (lo >> (kLimbBits - count));
}

// Returns the low limb of (hi:lo) >> count. Requires 0 < count < 64.
HWYMP_ALWAYS_INLINE constexpr limb_t shift_right_2(limb_t hi, limb_t lo,
                                                   int count) noexcept {
  return (lo >> count) | (hi << (kLimbBits - count));
}

// Returns hi:lo / d and stores the remainder. Requires hi < d, which also
// guarantees that d != 0 and that the quotient fits in one limb.
HWYMP_ALWAYS_INLINE constexpr limb_t div_2by1(limb_t hi, limb_t lo, limb_t d,
                                              limb_t& remainder) noexcept {
  if (std::is_constant_evaluated()) {
    return detail::portable::div_2by1(hi, lo, d, remainder);
  }
#if HWYMP_GNU_ASM && HWYMP_ARCH_X64
  return detail::div_2by1_asm(hi, lo, d, remainder);
#elif defined(_MSC_VER) && HWYMP_ARCH_X64
  unsigned long long rem;
  const limb_t quotient = _udiv128(hi, lo, d, &rem);
  remainder = rem;
  return quotient;
#else
  return detail::portable::div_2by1(hi, lo, d, remainder);
#endif
}

}  // namespace hwymp

#endif  // HWYMP_LIMB_OPS_H_
