// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Single-threaded throughput of scalar operations on identical random inputs.
// All engines round to nearest. "MPFR" uses its native wide exponent range,
// which favours it; "MPFR IEEE" adds the emin/emax and subnormalize steps that
// are needed for IEEE-identical results.

#include <gmp.h>
#include <mpfr.h>

#include <algorithm>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "mpfix/float.h"

#ifdef MPFIX_BENCH_SOFTFLOAT
extern "C" {
#include "softfloat.h"
}
#endif

#if defined(_WIN32)
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

namespace {

constexpr int kCount = 1024;
constexpr int kExponentSpread = 8;

std::uint64_t Mix(std::uint64_t& state) {
  state += 0x9E3779B97F4A7C15u;
  std::uint64_t z = state;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9u;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBu;
  return z ^ (z >> 31);
}

struct Sample {
  bool negative;
  int exponent;
  std::vector<std::uint64_t> limbs;  // Normalized, least significant first.
};

std::vector<Sample> MakeSamples(int precision, std::uint64_t seed) {
  const int limb_count = (precision + 63) / 64;
  const int spare = limb_count * 64 - precision;
  std::uint64_t state = seed;
  std::vector<Sample> samples(kCount);
  for (Sample& sample : samples) {
    sample.negative = (Mix(state) & 1) != 0;
    sample.exponent = static_cast<int>(Mix(state) % (2 * kExponentSpread + 1)) -
                      kExponentSpread;
    sample.limbs.resize(limb_count);
    for (std::uint64_t& limb : sample.limbs) limb = Mix(state);
    if (spare != 0) sample.limbs[0] &= ~((std::uint64_t{1} << spare) - 1);
    sample.limbs.back() |= std::uint64_t{1} << 63;
  }
  return samples;
}

template <class F>
struct HwyEngine {
  using Value = F;
  static const char* Name() { return "mpfix"; }
  static void Init(Value&) {}
  static void Load(Value& x, const Sample& s) {
    typename F::Limbs limbs;
    std::copy(s.limbs.begin(), s.limbs.end(), limbs.begin());
    x = F::FromRaw(s.negative, s.exponent, limbs);
  }
  static void Add(Value& r, const Value& a, const Value& b) { r = a + b; }
  static void Mul(Value& r, const Value& a, const Value& b) { r = a * b; }
  static void Div(Value& r, const Value& a, const Value& b) { r = a / b; }
  static void Sqrt(Value& r, const Value& a) { r = sqrt(a); }
  static constexpr bool kHasFma = true;
  static void Fma(Value& r, const Value& a, const Value& b, const Value& c) {
    r = fma(a, b, c);
  }
  static std::uint64_t Fingerprint(const Value& x) { return x.raw_limbs()[0]; }
};

template <int P, int E, bool kIeee>
struct MpfrEngine {
  struct Value {
    mpfr_t x;
    Value() { mpfr_init2(x, P); }
    ~Value() { mpfr_clear(x); }
    Value(const Value&) = delete;
    Value& operator=(const Value&) = delete;
  };
  static const char* Name() { return kIeee ? "MPFR IEEE" : "MPFR"; }
  static void Init(Value&) {}
  static void Load(Value& v, const Sample& s) {
    std::vector<std::uint64_t> limbs = s.limbs;
    mpfr_t view;
    mpfr_custom_init_set(view, (s.negative ? -1 : 1) * MPFR_REGULAR_KIND,
                         s.exponent + 1, P, limbs.data());
    mpfr_set(v.x, view, MPFR_RNDN);
  }
  static void Finish(Value& r, int ternary) {
    if constexpr (kIeee) {
      ternary = mpfr_check_range(r.x, ternary, MPFR_RNDN);
      mpfr_subnormalize(r.x, ternary, MPFR_RNDN);
    }
  }
  static void Add(Value& r, const Value& a, const Value& b) {
    Finish(r, mpfr_add(r.x, a.x, b.x, MPFR_RNDN));
  }
  static void Mul(Value& r, const Value& a, const Value& b) {
    Finish(r, mpfr_mul(r.x, a.x, b.x, MPFR_RNDN));
  }
  static void Div(Value& r, const Value& a, const Value& b) {
    Finish(r, mpfr_div(r.x, a.x, b.x, MPFR_RNDN));
  }
  static void Sqrt(Value& r, const Value& a) {
    Finish(r, mpfr_sqrt(r.x, a.x, MPFR_RNDN));
  }
  static constexpr bool kHasFma = true;
  static void Fma(Value& r, const Value& a, const Value& b, const Value& c) {
    Finish(r, mpfr_fma(r.x, a.x, b.x, c.x, MPFR_RNDN));
  }
  static std::uint64_t Fingerprint(const Value& v) {
    return static_cast<std::uint64_t>(mpfr_get_exp(v.x));
  }
};

template <int P>
struct BoostEngine {
  using Value = boost::multiprecision::number<
      boost::multiprecision::cpp_bin_float<
          P, boost::multiprecision::backends::digit_base_2>,
      boost::multiprecision::et_off>;
  static const char* Name() { return "Boost"; }
  static void Init(Value&) {}
  static void Load(Value& x, const Sample& s) {
    x = 0;
    for (std::size_t i = s.limbs.size(); i-- > 0;) {
      x = ldexp(x, 64) + Value(s.limbs[i]);
    }
    x = ldexp(x, s.exponent - static_cast<int>(s.limbs.size()) * 64 + 1);
    if (s.negative) x = -x;
  }
  static void Add(Value& r, const Value& a, const Value& b) { r = a + b; }
  static void Mul(Value& r, const Value& a, const Value& b) { r = a * b; }
  static void Div(Value& r, const Value& a, const Value& b) { r = a / b; }
  static void Sqrt(Value& r, const Value& a) { r = sqrt(a); }
  // Boost has no correctly rounded fused multiply-add for this type.
  static constexpr bool kHasFma = false;
  static void Fma(Value&, const Value&, const Value&, const Value&) {}
  static std::uint64_t Fingerprint(const Value& x) {
    return static_cast<std::uint64_t>(x.backend().exponent());
  }
};

#ifdef MPFIX_BENCH_SOFTFLOAT
// Berkeley SoftFloat's binary128, rounding to nearest-even with all flags, as
// mpfix does.
struct SoftFloatEngine {
  using Value = float128_t;
  static const char* Name() { return "SoftFloat"; }
  static void Init(Value&) {}
  static void Load(Value& x, const Sample& s) {
    mpfix::binary128::Limbs limbs;
    std::copy(s.limbs.begin(), s.limbs.end(), limbs.begin());
    const auto bits =
        mpfix::binary128::FromRaw(s.negative, s.exponent, limbs).to_bits();
    x.v[0] = bits[0];
    x.v[1] = bits[1];
  }
  static void Add(Value& r, const Value& a, const Value& b) {
    r = f128_add(a, b);
  }
  static void Mul(Value& r, const Value& a, const Value& b) {
    r = f128_mul(a, b);
  }
  static void Div(Value& r, const Value& a, const Value& b) {
    r = f128_div(a, b);
  }
  static void Sqrt(Value& r, const Value& a) { r = f128_sqrt(a); }
  static constexpr bool kHasFma = true;
  static void Fma(Value& r, const Value& a, const Value& b, const Value& c) {
    r = f128_mulAdd(a, b, c);
  }
  static std::uint64_t Fingerprint(const Value& x) { return x.v[0]; }
};

// SoftFloat and mpfix agree bit for bit on the benchmark's inputs, which also
// shows that the values reach SoftFloat intact.
bool SoftFloatAgrees() {
  using F = mpfix::binary128;
  const std::vector<Sample> sa = MakeSamples(113, 1);
  const std::vector<Sample> sb = MakeSamples(113, 2);
  const std::vector<Sample> sc = MakeSamples(113, 3);
  const auto same = [](const F& x, const float128_t& y) {
    const auto bits = x.to_bits();
    return bits[0] == y.v[0] && bits[1] == y.v[1];
  };
  for (int i = 0; i < kCount; ++i) {
    F a, b, c;
    HwyEngine<F>::Load(a, sa[i]);
    HwyEngine<F>::Load(b, sb[i]);
    HwyEngine<F>::Load(c, sc[i]);
    float128_t x, y, z;
    SoftFloatEngine::Load(x, sa[i]);
    SoftFloatEngine::Load(y, sb[i]);
    SoftFloatEngine::Load(z, sc[i]);
    if (!same(a + b, f128_add(x, y)) || !same(a * b, f128_mul(x, y)) ||
        !same(a / b, f128_div(x, y)) || !same(sqrt(abs(a)), f128_sqrt(
                                            {{x.v[0], x.v[1] & ~(1ull << 63)}})) ||
        !same(fma(a, b, c), f128_mulAdd(x, y, z))) {
      return false;
    }
  }
  return true;
}
#endif

volatile std::uint64_t sink;

// Repeats `pass`, which performs `operations` operations, and returns the best
// time per operation in nanoseconds, which is the least disturbed measurement.
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
    const std::chrono::duration<double, std::nano> elapsed =
        Clock::now() - start;
    best = std::min(best, elapsed.count() / (double(repeats) * operations));
  }
  return best;
}

template <class Engine>
struct Times {
  double add;
  double mul;
  double div;
  double sqrt;
  double fma;
  double dot;
};

template <class Engine>
Times<Engine> Measure(int precision) {
  using Value = typename Engine::Value;
  const std::vector<Sample> sample_a = MakeSamples(precision, 1);
  const std::vector<Sample> sample_b = MakeSamples(precision, 2);
  const std::vector<Sample> sample_c = MakeSamples(precision, 3);
  auto a = std::make_unique<Value[]>(kCount);
  auto b = std::make_unique<Value[]>(kCount);
  auto c = std::make_unique<Value[]>(kCount);
  auto positive = std::make_unique<Value[]>(kCount);
  auto r = std::make_unique<Value[]>(kCount);
  for (int i = 0; i < kCount; ++i) {
    Engine::Load(a[i], sample_a[i]);
    Engine::Load(b[i], sample_b[i]);
    Engine::Load(c[i], sample_c[i]);
    Sample unsigned_sample = sample_a[i];
    unsigned_sample.negative = false;
    Engine::Load(positive[i], unsigned_sample);
    Engine::Load(r[i], sample_a[i]);
  }
  Value sum;
  Value product;
  Value next;
  Engine::Load(product, sample_a[0]);
  Engine::Load(next, sample_a[0]);

  Times<Engine> times;
  times.add = BestNanoseconds(
      [&] {
        for (int i = 0; i < kCount; ++i) Engine::Add(r[i], a[i], b[i]);
        sink = Engine::Fingerprint(r[kCount / 2]);
      },
      kCount);
  times.mul = BestNanoseconds(
      [&] {
        for (int i = 0; i < kCount; ++i) Engine::Mul(r[i], a[i], b[i]);
        sink = Engine::Fingerprint(r[kCount / 2]);
      },
      kCount);
  times.div = BestNanoseconds(
      [&] {
        for (int i = 0; i < kCount; ++i) Engine::Div(r[i], a[i], b[i]);
        sink = Engine::Fingerprint(r[kCount / 2]);
      },
      kCount);
  times.sqrt = BestNanoseconds(
      [&] {
        for (int i = 0; i < kCount; ++i) Engine::Sqrt(r[i], positive[i]);
        sink = Engine::Fingerprint(r[kCount / 2]);
      },
      kCount);
  times.fma = -1;
  if constexpr (Engine::kHasFma) {
    times.fma = BestNanoseconds(
        [&] {
          for (int i = 0; i < kCount; ++i) Engine::Fma(r[i], a[i], b[i], c[i]);
          sink = Engine::Fingerprint(r[kCount / 2]);
        },
        kCount);
  }
  // A dot product over pairs, restarted so that the sum stays in range. Counted
  // per multiply-add step.
  times.dot = BestNanoseconds(
      [&] {
        Engine::Load(sum, sample_a[0]);
        for (int i = 0; i < kCount; ++i) {
          Engine::Mul(product, a[i], b[i]);
          Engine::Add(next, sum, product);
          Engine::Add(sum, next, product);
        }
        sink = Engine::Fingerprint(sum);
      },
      kCount * 3);
  return times;
}

template <class Engine>
void Report(const char* format, int precision) {
  const Times<Engine> t = Measure<Engine>(precision);
  std::printf("%-12s %-11s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n", format,
              Engine::Name(), t.add, t.mul, t.div, t.sqrt, t.fma, t.dot);
  std::fflush(stdout);
}

template <int P, int E>
void ReportFormat(const char* format) {
  mpfr_set_emin(mpfr_get_emin_min());
  mpfr_set_emax(mpfr_get_emax_max());
  Report<HwyEngine<mpfix::Float<P, E>>>(format, P);
  Report<MpfrEngine<P, E, false>>(format, P);
  // MPFR's exponent is one larger than IEEE's, and its emin refers to the
  // smallest subnormal.
  const int emax = (1 << (E - 1)) - 1;
  mpfr_set_emin((1 - emax) - (P - 1) + 1);
  mpfr_set_emax(emax + 1);
  Report<MpfrEngine<P, E, true>>(format, P);
  mpfr_set_emin(mpfr_get_emin_min());
  mpfr_set_emax(mpfr_get_emax_max());
  Report<BoostEngine<P>>(format, P);
#ifdef MPFIX_BENCH_SOFTFLOAT
  if constexpr (P == 113 && E == 15) {
    if (!SoftFloatAgrees()) {
      std::printf("SoftFloat disagrees with mpfix; not measured\n");
    } else {
      Report<SoftFloatEngine>(format, P);
    }
  }
#endif
}

}  // namespace

int main() {
#if defined(_WIN32)
  SetThreadAffinityMask(GetCurrentThread(), 1u << 2);
  SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
#endif
  std::printf("GMP %s, MPFR %s. Nanoseconds per operation; lower is better.\n",
              gmp_version, mpfr_get_version());
  std::printf("%-12s %-11s %9s %9s %9s %9s %9s %9s\n", "format", "library",
              "add", "mul", "div", "sqrt", "fma", "dot step");
  ReportFormat<113, 15>("binary128");
  ReportFormat<128, 15>("full 128");
  ReportFormat<237, 19>("binary256");
  ReportFormat<256, 19>("full 256");
  ReportFormat<489, 23>("binary512");
  ReportFormat<512, 23>("full 512");
  return 0;
}
