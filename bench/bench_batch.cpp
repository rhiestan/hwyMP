// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Prototype for the batch decision: binary256 add and mul on Highway vectors,
// with each lane an independent number in structure-of-arrays form.
//
// A significand is nine limbs of radix 2^28 in 64-bit lanes, which leaves room
// to accumulate 32x32-bit products without intermediate carries. Results round
// to nearest-even. Regular in-range operands only: this measures the core of
// the algorithm, so the real thing would be somewhat slower. Every result is
// compared with the scalar Float.

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "hwy/aligned_allocator.h"
#include "hwy/highway.h"
#include "mpfix/float.h"

#if defined(_WIN32)
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

namespace {

namespace hn = hwy::HWY_NAMESPACE;

using F = mpfix::binary256;
using D = hn::ScalableTag<std::uint64_t>;
using Vec = hn::Vec<D>;

constexpr int kLimbs = 9;
constexpr int kRadixBits = 28;
constexpr std::uint64_t kMask = (std::uint64_t{1} << kRadixBits) - 1;
constexpr int kRadixTotalBits = kLimbs * kRadixBits;           // 252
constexpr int kSpareBits = kRadixTotalBits - F::kPrecision;    // 15
constexpr int kBatchSize = 1024;
constexpr int kExponentSpread = 8;

std::uint64_t Mix(std::uint64_t& state) {
  state += 0x9E3779B97F4A7C15u;
  std::uint64_t z = state;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9u;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBu;
  return z ^ (z >> 31);
}

std::vector<F> MakeValues(std::uint64_t seed) {
  std::uint64_t state = seed;
  std::vector<F> values;
  values.reserve(kBatchSize);
  for (int i = 0; i < kBatchSize; ++i) {
    F::Limbs limbs;
    for (auto& limb : limbs) limb = Mix(state);
    limbs[0] &= ~((std::uint64_t{1} << 19) - 1);
    limbs[3] |= std::uint64_t{1} << 63;
    const bool negative = (Mix(state) & 1) != 0;
    const int exponent =
        static_cast<int>(Mix(state) % (2 * kExponentSpread + 1)) -
        kExponentSpread;
    values.push_back(F::FromRaw(negative, exponent, limbs));
  }
  return values;
}

// Reads `count` bits at `position` from little-endian 64-bit limbs.
std::uint64_t ReadBits(const F::Limbs& limbs, int position, int count) {
  const int word = position / 64;
  const int offset = position % 64;
  std::uint64_t value = limbs[word] >> offset;
  if (offset + count > 64) value |= limbs[word + 1] << (64 - offset);
  return value & ((std::uint64_t{1} << count) - 1);
}

struct Soa {
  hwy::AlignedFreeUniquePtr<std::uint64_t[]> storage;
  std::uint64_t* limbs[kLimbs];
  std::uint64_t* exponent;  // Two's complement.
  std::uint64_t* negative;  // All ones when negative.

  Soa() {
    storage = hwy::AllocateAligned<std::uint64_t>(
        std::size_t{kLimbs + 2} * kBatchSize);
    for (int i = 0; i < kLimbs; ++i) {
      limbs[i] = storage.get() + std::ptrdiff_t{i} * kBatchSize;
    }
    exponent = storage.get() + std::ptrdiff_t{kLimbs} * kBatchSize;
    negative = storage.get() + std::ptrdiff_t{kLimbs + 1} * kBatchSize;
  }

  void Set(int index, const F& x) {
    for (int i = 0; i < kLimbs; ++i) {
      limbs[i][index] = ReadBits(x.raw_limbs(), 256 - kRadixTotalBits + i * kRadixBits,
                                 kRadixBits);
    }
    exponent[index] = static_cast<std::uint64_t>(std::int64_t{x.raw_exponent()});
    negative[index] = x.signbit() ? ~std::uint64_t{0} : 0;
  }

  F Get(int index) const {
    F::Limbs result{};
    for (int i = 0; i < kLimbs; ++i) {
      const int position = 256 - kRadixTotalBits + i * kRadixBits;
      const int word = position / 64;
      const int offset = position % 64;
      result[word] |= limbs[i][index] << offset;
      if (offset + kRadixBits > 64) {
        result[word + 1] |= limbs[i][index] >> (64 - offset);
      }
    }
    return F::FromRaw(negative[index] != 0,
                      static_cast<std::int32_t>(exponent[index]), result);
  }
};

HWY_INLINE Vec TestBit(Vec value, int bit) {
  const D d;
  return hn::And(hn::ShiftRightSame(value, bit), hn::Set(d, 1));
}

// Rounds to nearest-even. x[base] holds the low limb of the significand, and
// x[0..base) plus `sticky` hold the bits below it. x[base + 9] may hold an
// overflow into the next limb only after the increment here.
HWY_INLINE void RoundAndStore(Vec* HWY_RESTRICT x, int base, Vec sticky,
                              Vec exponent, Vec negative, Soa& out, int index) {
  const D d;
  const Vec low = x[base];
  const Vec round_bit = TestBit(low, kSpareBits - 1);
  const Vec odd = TestBit(low, kSpareBits);
  const Vec below =
      hn::And(low, hn::Set(d, (std::uint64_t{1} << (kSpareBits - 1)) - 1));
  for (int i = 0; i < base; ++i) sticky = hn::Or(sticky, x[i]);
  sticky = hn::Or(sticky, below);
  const auto sticky_set = hn::Ne(sticky, hn::Zero(d));
  const auto odd_set = hn::Ne(odd, hn::Zero(d));
  const Vec tie_or_more = hn::IfThenElseZero(
      hn::Or(sticky_set, odd_set), round_bit);

  Vec carry = hn::ShiftLeft<kSpareBits>(tie_or_more);
  carry = hn::Add(carry, hn::AndNot(hn::Set(d, (std::uint64_t{1} << kSpareBits) - 1), low));
  Vec limb = hn::And(carry, hn::Set(d, kMask));
  hn::Store(limb, d, out.limbs[0] + index);
  const Vec mask = hn::Set(d, kMask);
  for (int i = 1; i < kLimbs; ++i) {
    carry = hn::Add(x[base + i], hn::ShiftRight<kRadixBits>(carry));
    hn::Store(hn::And(carry, mask), d, out.limbs[i] + index);
  }
  const Vec overflow = hn::ShiftRight<kRadixBits>(carry);
  const auto overflowed = hn::Ne(overflow, hn::Zero(d));
  if (HWY_UNLIKELY(!hn::AllFalse(d, overflowed))) {
    const Vec top = hn::Load(d, out.limbs[kLimbs - 1] + index);
    hn::Store(hn::IfThenElse(overflowed, hn::Set(d, std::uint64_t{1} << 27), top),
              d, out.limbs[kLimbs - 1] + index);
  }
  hn::Store(hn::Add(exponent, overflow), d, out.exponent + index);
  hn::Store(negative, d, out.negative + index);
}

HWY_NOINLINE void MulBatch(const Soa& a, const Soa& b, Soa& out) {
  const D d;
  const hn::Repartition<std::uint32_t, D> d32;
  const size_t lanes = hn::Lanes(d);
  const Vec mask = hn::Set(d, kMask);

  for (int index = 0; index + static_cast<int>(lanes) <= kBatchSize;
       index += static_cast<int>(lanes)) {
    hn::Vec<decltype(d32)> va[kLimbs];
    hn::Vec<decltype(d32)> vb[kLimbs];
    for (int i = 0; i < kLimbs; ++i) {
      va[i] = hn::BitCast(d32, hn::Load(d, a.limbs[i] + index));
      vb[i] = hn::BitCast(d32, hn::Load(d, b.limbs[i] + index));
    }

    Vec column[2 * kLimbs - 1];
    for (Vec& sum : column) sum = hn::Zero(d);
    for (int i = 0; i < kLimbs; ++i) {
      for (int j = 0; j < kLimbs; ++j) {
        column[i + j] = hn::Add(column[i + j], hn::MulEven(va[i], vb[j]));
      }
    }

    Vec t[2 * kLimbs];
    Vec carry = hn::Zero(d);
    for (int i = 0; i < 2 * kLimbs - 1; ++i) {
      const Vec value = hn::Add(column[i], carry);
      t[i] = hn::And(value, mask);
      carry = hn::ShiftRight<kRadixBits>(value);
    }
    t[2 * kLimbs - 1] = carry;

    // The product is in [1, 4). Limbs 0..7 only ever contribute to sticky.
    const Vec high_bit = TestBit(t[2 * kLimbs - 1], kRadixBits - 1);
    const auto is_high = hn::Ne(high_bit, hn::Zero(d));
    Vec sticky = hn::Zero(d);
    for (int i = 0; i < kLimbs - 1; ++i) sticky = hn::Or(sticky, t[i]);
    for (int i = 2 * kLimbs - 1; i >= kLimbs; --i) {
      const Vec shifted = hn::And(
          hn::Or(hn::ShiftLeft<1>(t[i]), hn::ShiftRight<kRadixBits - 1>(t[i - 1])),
          mask);
      t[i] = hn::IfThenElse(is_high, t[i], shifted);
    }
    const Vec low_limb = hn::IfThenElse(
        is_high, t[kLimbs - 1],
        hn::And(hn::ShiftLeft<1>(t[kLimbs - 1]), mask));

    const Vec exponent =
        hn::Add(hn::Add(hn::Load(d, a.exponent + index),
                        hn::Load(d, b.exponent + index)),
                high_bit);
    const Vec negative = hn::Xor(hn::Load(d, a.negative + index),
                                 hn::Load(d, b.negative + index));
    t[kLimbs - 1] = low_limb;
    RoundAndStore(t + kLimbs - 1, 1, sticky, exponent, negative, out, index);
  }
}

// Shifts ext[] down by `amount` limbs where `active`, accumulating lost limbs.
template <int kAmount, int kSize>
HWY_INLINE void ShiftWords(Vec* HWY_RESTRICT ext, hn::Mask<D> active,
                           Vec& sticky) {
  const D d;
  Vec lost = hn::Zero(d);
  for (int i = 0; i < std::min(kAmount, kSize); ++i) {
    lost = hn::Or(lost, ext[i]);
  }
  sticky = hn::Or(sticky, hn::IfThenElseZero(active, lost));
  for (int i = 0; i < kSize; ++i) {
    const Vec moved = i + kAmount < kSize ? ext[i + kAmount] : hn::Zero(d);
    ext[i] = hn::IfThenElse(active, moved, ext[i]);
  }
}

// Shifts ext[] up by `amount` limbs where `active`.
template <int kAmount, int kSize>
HWY_INLINE void ShiftWordsUp(Vec* HWY_RESTRICT ext, hn::Mask<D> active) {
  const D d;
  for (int i = kSize - 1; i >= 0; --i) {
    const Vec moved = i >= kAmount ? ext[i - kAmount] : hn::Zero(d);
    ext[i] = hn::IfThenElse(active, moved, ext[i]);
  }
}

HWY_NOINLINE void AddBatch(const Soa& a, const Soa& b, Soa& out) {
  const D d;
  const hn::RebindToSigned<D> di;
  const int lanes = static_cast<int>(hn::Lanes(d));
  const Vec mask = hn::Set(d, kMask);
  // Two extra low limbs keep 56 bits below the significand.
  constexpr int kExtra = 2;
  constexpr int kSize = kLimbs + kExtra;

  for (int index = 0; index + lanes <= kBatchSize; index += lanes) {
    Vec big[kLimbs];
    Vec small[kLimbs];
    const Vec exp_a = hn::Load(d, a.exponent + index);
    const Vec exp_b = hn::Load(d, b.exponent + index);
    const Vec neg_a = hn::Load(d, a.negative + index);
    const Vec neg_b = hn::Load(d, b.negative + index);

    // Order the operands by magnitude.
    auto b_larger = hn::RebindMask(
        d, hn::Lt(hn::BitCast(di, exp_a), hn::BitCast(di, exp_b)));
    auto a_larger = hn::RebindMask(
        d, hn::Lt(hn::BitCast(di, exp_b), hn::BitCast(di, exp_a)));
    Vec limb_a[kLimbs];
    Vec limb_b[kLimbs];
    for (int i = kLimbs - 1; i >= 0; --i) {
      limb_a[i] = hn::Load(d, a.limbs[i] + index);
      limb_b[i] = hn::Load(d, b.limbs[i] + index);
      const auto undecided = hn::Not(hn::Or(a_larger, b_larger));
      b_larger = hn::Or(b_larger, hn::And(undecided, hn::Lt(limb_a[i], limb_b[i])));
      a_larger = hn::Or(a_larger, hn::And(undecided, hn::Lt(limb_b[i], limb_a[i])));
    }
    for (int i = 0; i < kLimbs; ++i) {
      big[i] = hn::IfThenElse(b_larger, limb_b[i], limb_a[i]);
      small[i] = hn::IfThenElse(b_larger, limb_a[i], limb_b[i]);
    }
    const Vec big_exponent = hn::IfThenElse(b_larger, exp_b, exp_a);
    const Vec small_exponent = hn::IfThenElse(b_larger, exp_a, exp_b);
    const Vec negative = hn::IfThenElse(b_larger, neg_b, neg_a);
    const auto subtract = hn::Ne(neg_a, neg_b);
    const Vec distance = hn::Sub(big_exponent, small_exponent);

    // Align the smaller operand.
    Vec ext[kSize];
    ext[0] = hn::Zero(d);
    ext[1] = hn::Zero(d);
    for (int i = 0; i < kLimbs; ++i) ext[i + kExtra] = small[i];
    Vec sticky = hn::Zero(d);
    const Vec limit = hn::Set(d, std::uint64_t{15} * kRadixBits);
    const Vec clamped = hn::Min(distance, limit);
    // Division by 28: the quotient only needs to be right below 16 * 28.
    const Vec words = hn::ShiftRight<32>(
        hn::MulEven(hn::BitCast(hn::Repartition<std::uint32_t, D>(), clamped),
                    hn::BitCast(hn::Repartition<std::uint32_t, D>(),
                                hn::Set(d, 153391690))));
    const Vec bits = hn::Sub(
        clamped, hn::Sub(hn::ShiftLeft<5>(words), hn::ShiftLeft<2>(words)));
    const auto has = [&](int bit) {
      return hn::Ne(hn::And(words, hn::Set(d, bit)), hn::Zero(d));
    };
    ShiftWords<8, kSize>(ext, has(8), sticky);
    ShiftWords<4, kSize>(ext, has(4), sticky);
    ShiftWords<2, kSize>(ext, has(2), sticky);
    ShiftWords<1, kSize>(ext, has(1), sticky);

    const Vec up = hn::Sub(hn::Set(d, kRadixBits), bits);
    sticky = hn::Or(sticky, hn::And(hn::Shl(ext[0], up), mask));
    for (int i = 0; i < kSize; ++i) {
      const Vec high = i + 1 < kSize ? ext[i + 1] : hn::Zero(d);
      ext[i] = hn::And(hn::Or(hn::Shr(ext[i], bits), hn::Shl(high, up)), mask);
    }
    // Everything lost is jammed into the lowest bit, far below the round bit.
    const auto any_lost = hn::Ne(sticky, hn::Zero(d));
    ext[0] = hn::IfThenElse(any_lost, hn::Or(ext[0], hn::Set(d, 1)), ext[0]);

    // Add or subtract limb by limb, then propagate signed carries.
    Vec result[kSize + 1];
    Vec carry = hn::Zero(d);
    for (int i = 0; i < kSize; ++i) {
      const Vec top = i >= kExtra ? big[i - kExtra] : hn::Zero(d);
      const Vec sum = hn::IfThenElse(subtract, hn::Sub(top, ext[i]),
                                     hn::Add(top, ext[i]));
      const Vec value = hn::Add(sum, carry);
      result[i] = hn::And(value, mask);
      carry = hn::BitCast(d, hn::ShiftRight<kRadixBits>(hn::BitCast(di, value)));
    }
    // After addition the carry is the overflow limb; after subtraction it is 0.
    Vec exponent = big_exponent;
    const auto overflowed = hn::Ne(carry, hn::Zero(d));
    if (HWY_LIKELY(!hn::AllFalse(d, overflowed))) {
      Vec previous = carry;
      for (int i = kSize - 1; i >= 0; --i) {
        const Vec shifted = hn::Or(hn::ShiftRight<1>(result[i]),
                                   hn::And(hn::ShiftLeft<kRadixBits - 1>(previous), mask));
        const Vec keep_low = i == 0 ? hn::And(result[0], hn::Set(d, 1)) : hn::Zero(d);
        previous = result[i];
        result[i] = hn::IfThenElse(overflowed, hn::Or(shifted, keep_low), result[i]);
      }
      exponent = hn::Add(exponent, hn::IfThenElseZero(overflowed, hn::Set(d, 1)));
    }

    // Normalize after cancellation: whole limbs first, then bits.
    Vec moved = hn::Zero(d);
    const auto zero_top = [&](int count) {
      Vec combined = hn::Zero(d);
      for (int i = 0; i < count; ++i) combined = hn::Or(combined, result[kSize - 1 - i]);
      return hn::Eq(combined, hn::Zero(d));
    };
    if (HWY_UNLIKELY(!hn::AllFalse(d, zero_top(1)))) {
      auto active = zero_top(8);
      moved = hn::IfThenElseZero(active, hn::Set(d, std::uint64_t{8} * kRadixBits));
      ShiftWordsUp<8, kSize>(result, active);
      active = zero_top(4);
      moved = hn::Add(moved, hn::IfThenElseZero(active, hn::Set(d, std::uint64_t{4} * kRadixBits)));
      ShiftWordsUp<4, kSize>(result, active);
      active = zero_top(2);
      moved = hn::Add(moved, hn::IfThenElseZero(active, hn::Set(d, std::uint64_t{2} * kRadixBits)));
      ShiftWordsUp<2, kSize>(result, active);
      active = zero_top(1);
      moved = hn::Add(moved, hn::IfThenElseZero(active, hn::Set(d, kRadixBits)));
      ShiftWordsUp<1, kSize>(result, active);
    }
    const Vec top = hn::Or(result[kSize - 1], hn::IfThenElseZero(hn::Eq(result[kSize - 1], hn::Zero(d)), hn::Set(d, 1)));
    const Vec leading = hn::Sub(hn::LeadingZeroCount(top), hn::Set(d, 64 - kRadixBits));
    const auto needs_shift = hn::Ne(leading, hn::Zero(d));
    if (HWY_LIKELY(!hn::AllFalse(d, needs_shift))) {
      const Vec complement = hn::Sub(hn::Set(d, kRadixBits), leading);
      for (int i = kSize - 1; i >= 0; --i) {
        const Vec lower = i > 0 ? result[i - 1] : hn::Zero(d);
        const Vec shifted = hn::And(
            hn::Or(hn::Shl(result[i], leading), hn::Shr(lower, complement)), mask);
        result[i] = hn::IfThenElse(needs_shift, shifted, result[i]);
      }
    }
    exponent = hn::Sub(exponent, hn::Add(moved, leading));

    RoundAndStore(result, kExtra, hn::Zero(d), exponent, negative, out, index);
  }
}

volatile std::uint64_t sink;

template <class Pass>
double BestNanoseconds(Pass pass, int operations) {
  using Clock = std::chrono::steady_clock;
  int repeats = 1;
  for (;;) {
    const auto start = Clock::now();
    for (int i = 0; i < repeats; ++i) pass();
    if (Clock::now() - start > std::chrono::milliseconds(20)) break;
    repeats *= 2;
  }
  double best = 1e300;
  for (int round = 0; round < 25; ++round) {
    const auto start = Clock::now();
    for (int i = 0; i < repeats; ++i) pass();
    const std::chrono::duration<double, std::nano> elapsed = Clock::now() - start;
    best = std::min(best, elapsed.count() / (static_cast<double>(repeats) * operations));
  }
  return best;
}

int Verify(const char* name, const Soa& out, const std::vector<F>& expected) {
  const int count = kBatchSize - kBatchSize % static_cast<int>(hn::Lanes(D()));
  int mismatches = 0;
  for (int i = 0; i < count; ++i) {
    if (expected[i].is_zero()) continue;  // The prototype skips exact zeros.
    if (out.Get(i).to_bits() != expected[i].to_bits()) {
      if (mismatches++ == 0) {
        const F::Bits got = out.Get(i).to_bits();
        const F::Bits want = expected[i].to_bits();
        std::printf("%s mismatch at %d: got %016llx %016llx want %016llx %016llx\n",
                    name, i, static_cast<unsigned long long>(got[3]),
                    static_cast<unsigned long long>(got[0]),
                    static_cast<unsigned long long>(want[3]),
                    static_cast<unsigned long long>(want[0]));
      }
    }
  }
  return mismatches;
}

}  // namespace

int main() {
#if defined(_WIN32)
  SetThreadAffinityMask(GetCurrentThread(), 1u << 2);
  SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
#endif
  const std::vector<F> a_values = MakeValues(1);
  const std::vector<F> b_values = MakeValues(2);
  Soa a;
  Soa b;
  Soa out;
  std::vector<F> scalar(kBatchSize);
  for (int i = 0; i < kBatchSize; ++i) {
    a.Set(i, a_values[i]);
    b.Set(i, b_values[i]);
  }

  std::printf("Highway target %s, %d lanes. Nanoseconds per number.\n",
              hwy::TargetName(HWY_TARGET), static_cast<int>(hn::Lanes(D())));

  const double scalar_mul = BestNanoseconds(
      [&] {
        for (int i = 0; i < kBatchSize; ++i) scalar[i] = a_values[i] * b_values[i];
        sink = scalar[kBatchSize / 2].raw_limbs()[0];
      },
      kBatchSize);
  const std::vector<F> expected_mul = scalar;
  const double batch_mul = BestNanoseconds(
      [&] {
        MulBatch(a, b, out);
        sink = out.limbs[0][kBatchSize / 2];
      },
      kBatchSize);
  const int mul_errors = Verify("mul", out, expected_mul);

  const double scalar_add = BestNanoseconds(
      [&] {
        for (int i = 0; i < kBatchSize; ++i) scalar[i] = a_values[i] + b_values[i];
        sink = scalar[kBatchSize / 2].raw_limbs()[0];
      },
      kBatchSize);
  const std::vector<F> expected_add = scalar;
  const double batch_add = BestNanoseconds(
      [&] {
        AddBatch(a, b, out);
        sink = out.limbs[0][kBatchSize / 2];
      },
      kBatchSize);
  const int add_errors = Verify("add", out, expected_add);

  // The cost of moving between Float storage and the SoA layout, per number,
  // which an array-of-structs container such as Eigen would pay for each of
  // the two inputs and the output.
  const double convert = BestNanoseconds(
      [&] {
        for (int i = 0; i < kBatchSize; ++i) out.Set(i, a_values[i]);
        sink = out.limbs[0][kBatchSize / 2];
      },
      kBatchSize);

  std::printf("%-6s %10s %10s %8s %8s\n", "op", "scalar", "batch", "speedup",
              "errors");
  std::printf("%-6s %10.2f %10.2f %8.2f %8d\n", "mul", scalar_mul, batch_mul,
              scalar_mul / batch_mul, mul_errors);
  std::printf("%-6s %10.2f %10.2f %8.2f %8d\n", "add", scalar_add, batch_add,
              scalar_add / batch_add, add_errors);
  std::printf("conversion to SoA: %.2f ns per number\n", convert);
  return 0;
}
