// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// libFuzzer target: add, sub, mul, div, sqrt, fma and remainder on operands
// of arbitrary bits, in every rounding mode and in four formats, against
// MPFR. Results and all five flags must match; operations on infinities and
// NaNs only must not crash, and give quiet NaNs.

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "mpfix/float.h"
#include "oracle.h"

namespace {

using namespace mpfix;        // NOLINT
using namespace mpfix::test;  // NOLINT

class Reader {
 public:
  Reader(const std::uint8_t* data, std::size_t size)
      : data_(data), size_(size) {}
  std::uint8_t Byte() { return position_ < size_ ? data_[position_++] : 0; }
  template <class F>
  F Value() {
    typename F::Bits bits{};
    for (auto& limb : bits) {
      for (int i = 0; i < 8; ++i) limb = limb << 8 | Byte();
    }
    return F::from_bits(bits);
  }

 private:
  const std::uint8_t* data_;
  std::size_t size_;
  std::size_t position_ = 0;
};

enum Op { kAdd, kSub, kMul, kDiv, kSqrt, kFma, kRemainder, kOps };

template <class F>
[[noreturn]] void Report(int op, Rounding mode, const F& a, const F& b,
                         const F& c, const F& actual, unsigned actual_flags,
                         const F& expected, unsigned expected_flags) {
  std::fprintf(stderr,
               "mismatch: P=%d op %d mode %d\n  a = %s\n  b = %s\n  c = %s\n"
               "  actual   %s flags %u\n  expected %s flags %u\n",
               F::kPrecision, op, static_cast<int>(mode), Hex(a).c_str(),
               Hex(b).c_str(), Hex(c).c_str(), Hex(actual).c_str(),
               actual_flags, Hex(expected).c_str(), expected_flags);
  std::abort();
}

template <class F>
void Check(Reader& reader) {
  const int op = reader.Byte() % kOps;
  const Rounding mode = kModes[reader.Byte() % 5];
  const F a = reader.Value<F>();
  const F b = reader.Value<F>();
  const F c = reader.Value<F>();

  Context context{mode};
  F actual;
  switch (op) {
    case kAdd: actual = add(a, b, context); break;
    case kSub: actual = sub(a, b, context); break;
    case kMul: actual = mul(a, b, context); break;
    case kDiv: actual = div(a, b, context); break;
    case kSqrt: actual = sqrt(a, context); break;
    case kFma: actual = fma(a, b, c, context); break;
    default: actual = remainder(a, b, context); break;
  }
  if (actual.is_signaling()) {
    Report(op, mode, a, b, c, actual, context.flags, F::QuietNaN(), 0);
  }
  // MPFR has no NaN payloads or signaling NaNs, and IEEE 754 leaves some of
  // their flags open: only finite operands are compared.
  const bool uses_b = op != kSqrt;
  const bool uses_c = op == kFma;
  if (!a.is_finite() || (uses_b && !b.is_finite()) ||
      (uses_c && !c.is_finite())) {
    return;
  }

  const View<F> va(a), vb(b), vc(c);
  Expected<F> expected;
  if (op == kDiv && b.is_zero() && !a.is_zero()) {
    expected = {F::Infinity(a.signbit() != b.signbit()),
                unsigned{kFlagDivByZero}};
  } else {
    expected = RoundExact<F>(
        [&](mpfr_ptr result, mpfr_rnd_t rnd) {
          switch (op) {
            case kAdd: return mpfr_add(result, va.get(), vb.get(), rnd);
            case kSub: return mpfr_sub(result, va.get(), vb.get(), rnd);
            case kMul: return mpfr_mul(result, va.get(), vb.get(), rnd);
            case kDiv: return mpfr_div(result, va.get(), vb.get(), rnd);
            case kSqrt: return CertifiedSqrt(result, va.get(), rnd);
            case kFma:
              return mpfr_fma(result, va.get(), vb.get(), vc.get(), rnd);
            default: return mpfr_remainder(result, va.get(), vb.get(), rnd);
          }
        },
        mode);
  }
  if (!SameResult(actual, expected.value) || context.flags != expected.flags) {
    Report(op, mode, a, b, c, actual, context.flags, expected.value,
           expected.flags);
  }
}

}  // namespace

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data,
                                      std::size_t size) {
  Reader reader(data, size);
  switch (reader.Byte() % 4) {
    case 0: Check<Float<11, 5>>(reader); break;
    case 1: Check<Float<53, 11>>(reader); break;
    case 2: Check<binary128>(reader); break;
    default: Check<binary256>(reader); break;
  }
  return 0;
}
