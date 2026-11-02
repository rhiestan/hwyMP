// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Single-threaded throughput of the exponentials and logarithms, on identical
// arguments for every library, all rounding to nearest. Exponentials take
// arguments in (-20, 20); logarithms take positive values with binary
// exponents in [-20, 20); expm1 and log1p take arguments below 1/4 in
// magnitude. MPFR runs at the format's precision; Boost.Multiprecision's
// cpp_bin_float has exp and log only.

#include <gmp.h>
#include <mpfr.h>

#include <algorithm>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <vector>

#include "mpfix/charconv.h"
#include "mpfix/math.h"

namespace {

constexpr int kCount = 256;

std::uint64_t Mix(std::uint64_t& state) {
  state += 0x9E3779B97F4A7C15u;
  std::uint64_t z = state;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9u;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBu;
  return z ^ (z >> 31);
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
  for (int round = 0; round < 9; ++round) {
    const auto start = Clock::now();
    for (int i = 0; i < repeats; ++i) pass();
    const std::chrono::duration<double, std::nano> elapsed =
        Clock::now() - start;
    best = std::min(best, elapsed.count() / (double(repeats) * operations));
  }
  return best;
}

enum class Range { kExp, kLog, kSmall };

template <class F>
std::vector<F> Arguments(Range range) {
  std::uint64_t state = static_cast<std::uint64_t>(range) + 1;
  std::vector<F> values(kCount);
  for (F& value : values) {
    typename F::Limbs limbs;
    for (auto& limb : limbs) limb = Mix(state);
    limbs[F::kLimbs - 1] |= std::uint64_t{1} << 63;
    const bool negative = (Mix(state) & 1) != 0;
    switch (range) {
      case Range::kExp:  // Below 2^4.32 = 20.
        value = scalbn(F::FromRaw(negative, 0, limbs),
                       static_cast<int>(Mix(state) % 5));
        value = value * F(0.6);
        break;
      case Range::kLog:
        value = F::FromRaw(false, static_cast<int>(Mix(state) % 40) - 20,
                           limbs);
        break;
      case Range::kSmall:
        value = F::FromRaw(negative, -3 - static_cast<int>(Mix(state) % 20),
                           limbs);
        break;
    }
  }
  return values;
}

// The first calls build the constants and tables, once per format.
template <int P, int E>
void ReportFirstCalls(const char* format) {
  using F = mpfix::Float<P, E>;
  using Clock = std::chrono::steady_clock;
  const auto start = Clock::now();
  sink = mpfix::exp(F(0.5)).raw_limbs()[0];
  const auto middle = Clock::now();
  sink = mpfix::log(F(3)).raw_limbs()[0];
  const auto end = Clock::now();
  const std::chrono::duration<double, std::micro> exp_us = middle - start;
  const std::chrono::duration<double, std::micro> log_us = end - middle;
  std::printf("%-10s first exp %.0f us, first log %.0f us\n", format,
              exp_us.count(), log_us.count());
}

template <int P, int E>
void ReportFormat(const char* format) {
  using F = mpfix::Float<P, E>;
  using Boost = boost::multiprecision::number<
      boost::multiprecision::cpp_bin_float<
          P, boost::multiprecision::backends::digit_base_2>,
      boost::multiprecision::et_off>;

  struct Case {
    const char* name;
    Range range;
    F (*mpfix_fn)(const F&);
    int (*mpfr_fn)(mpfr_ptr, mpfr_srcptr, mpfr_rnd_t);
    Boost (*boost_fn)(const Boost&);
  };
  const Case cases[] = {
      {"exp", Range::kExp, [](const F& x) { return mpfix::exp(x); }, mpfr_exp,
       [](const Boost& x) { return Boost(exp(x)); }},
      {"exp2", Range::kExp, [](const F& x) { return mpfix::exp2(x); },
       mpfr_exp2, nullptr},
      {"expm1", Range::kSmall, [](const F& x) { return mpfix::expm1(x); },
       mpfr_expm1, nullptr},
      {"log", Range::kLog, [](const F& x) { return mpfix::log(x); }, mpfr_log,
       [](const Boost& x) { return Boost(log(x)); }},
      {"log2", Range::kLog, [](const F& x) { return mpfix::log2(x); },
       mpfr_log2, nullptr},
      {"log10", Range::kLog, [](const F& x) { return mpfix::log10(x); },
       mpfr_log10, nullptr},
      {"log1p", Range::kSmall, [](const F& x) { return mpfix::log1p(x); },
       mpfr_log1p, nullptr},
  };

  for (const Case& c : cases) {
    const std::vector<F> values = Arguments<F>(c.range);
    std::vector<__mpfr_struct> mpfr_values(kCount);
    std::vector<Boost> boost_values(kCount);
    mpfr_t out;
    mpfr_init2(out, P);
    for (int i = 0; i < kCount; ++i) {
      mpfr_init2(&mpfr_values[i], P);
      const std::string text = mpfix::to_string(values[i]);
      mpfr_set_str(&mpfr_values[i], text.c_str(), 10, MPFR_RNDN);
      boost_values[i] = Boost(text);
    }

    const double mpfix_ns = BestNanoseconds(
        [&] {
          std::uint64_t acc = 0;
          for (const F& x : values) acc += c.mpfix_fn(x).raw_limbs()[0];
          sink = acc;
        },
        kCount);
    const double mpfr_ns = BestNanoseconds(
        [&] {
          std::uint64_t acc = 0;
          for (auto& x : mpfr_values) {
            c.mpfr_fn(out, &x, MPFR_RNDN);
            acc += static_cast<std::uint64_t>(mpfr_get_exp(out));
          }
          sink = acc;
        },
        kCount);
    double boost_ns = 0;
    if (c.boost_fn != nullptr) {
      boost_ns = BestNanoseconds(
          [&] {
            std::uint64_t acc = 0;
            for (const Boost& x : boost_values) {
              acc += static_cast<std::uint64_t>(c.boost_fn(x).backend().exponent());
            }
            sink = acc;
          },
          kCount);
    }
    if (boost_ns > 0) {
      std::printf("%-10s %-6s %10.1f %10.1f %10.1f %8.2fx\n", format, c.name,
                  mpfix_ns, mpfr_ns, boost_ns, mpfr_ns / mpfix_ns);
    } else {
      std::printf("%-10s %-6s %10.1f %10.1f %10s %8.2fx\n", format, c.name,
                  mpfix_ns, mpfr_ns, "-", mpfr_ns / mpfix_ns);
    }
    for (auto& x : mpfr_values) mpfr_clear(&x);
    mpfr_clear(out);
  }
}

}  // namespace

int main() {
  ReportFirstCalls<113, 15>("binary128");
  ReportFirstCalls<237, 19>("binary256");
  ReportFirstCalls<489, 23>("binary512");
  std::printf("%-10s %-6s %10s %10s %10s %9s\n", "format", "fn", "mpfix ns",
              "MPFR ns", "Boost ns", "MPFR/mpfix");
  ReportFormat<113, 15>("binary128");
  ReportFormat<237, 19>("binary256");
  ReportFormat<489, 23>("binary512");
  return 0;
}
