// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Single-threaded throughput of decimal string conversions, on identical
// values and strings for every library, all rounding to nearest. "print" gives
// the round-trip number of significant digits in scientific notation, and
// "parse" reads the shortest strings back. Values come near 1, or spread over
// decimal exponents from -300 to 300.

#include <gmp.h>
#include <mpfr.h>

#include <algorithm>
#include <boost/multiprecision/cpp_bin_float.hpp>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "mpfix/charconv.h"

#if defined(_WIN32)
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

namespace {

constexpr int kCount = 512;

std::uint64_t Mix(std::uint64_t& state) {
  state += 0x9E3779B97F4A7C15u;
  std::uint64_t z = state;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9u;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBu;
  return z ^ (z >> 31);
}

volatile std::size_t sink;

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
  for (int round = 0; round < 15; ++round) {
    const auto start = Clock::now();
    for (int i = 0; i < repeats; ++i) pass();
    const std::chrono::duration<double, std::nano> elapsed =
        Clock::now() - start;
    best = std::min(best, elapsed.count() / (double(repeats) * operations));
  }
  return best;
}

template <int P, int E>
void ReportFormat(const char* format, bool wide) {
  using F = mpfix::Float<P, E>;
  using Boost = boost::multiprecision::number<
      boost::multiprecision::cpp_bin_float<
          P, boost::multiprecision::backends::digit_base_2>,
      boost::multiprecision::et_off>;
  // Enough significant digits to read back.
  const int digits = static_cast<int>(P * 0.30103) + 2;

  std::uint64_t state = 1;
  std::vector<F> values(kCount);
  std::vector<std::string> strings(kCount);
  for (int i = 0; i < kCount; ++i) {
    typename F::Limbs limbs;
    for (auto& limb : limbs) limb = Mix(state);
    constexpr int kSpare = static_cast<int>(F::kLimbs) * 64 - P;
    if constexpr (kSpare != 0) limbs[0] &= ~((std::uint64_t{1} << kSpare) - 1);
    limbs.back() |= std::uint64_t{1} << 63;
    const int spread = wide ? 1000 : 8;
    const int exponent =
        static_cast<int>(Mix(state) % (2 * spread + 1)) - spread;
    values[i] = F::FromRaw((Mix(state) & 1) != 0, exponent, limbs);
    strings[i] = mpfix::to_string(values[i]);
  }
  struct MpfrValue {
    mpfr_t x;
  };
  std::vector<MpfrValue> mpfr_values(kCount);
  std::vector<Boost> boost_values(kCount);
  for (int i = 0; i < kCount; ++i) {
    mpfr_init2(mpfr_values[i].x, P);
    mpfr_strtofr(mpfr_values[i].x, strings[i].c_str(), nullptr, 10, MPFR_RNDN);
    boost_values[i] = Boost(strings[i]);
  }

  char buffer[1024];
  const double print_ours = BestNanoseconds(
      [&] {
        std::size_t total = 0;
        for (const F& x : values) {
          const auto result =
              mpfix::to_chars(buffer, buffer + sizeof(buffer), x,
                              std::chars_format::scientific, digits - 1);
          total += static_cast<std::size_t>(result.ptr - buffer);
        }
        sink = total;
      },
      kCount);
  const double shortest_ours = BestNanoseconds(
      [&] {
        std::size_t total = 0;
        for (const F& x : values) {
          const auto result =
              mpfix::to_chars(buffer, buffer + sizeof(buffer), x);
          total += static_cast<std::size_t>(result.ptr - buffer);
        }
        sink = total;
      },
      kCount);
  const double parse_ours = BestNanoseconds(
      [&] {
        std::size_t total = 0;
        for (const std::string& s : strings) {
          F x;
          mpfix::from_chars(s.data(), s.data() + s.size(), x);
          total += x.raw_limbs()[0];
        }
        sink = total;
      },
      kCount);

  const double print_mpfr = BestNanoseconds(
      [&] {
        std::size_t total = 0;
        for (auto& value : mpfr_values) {
          mpfr_exp_t exponent;
          mpfr_get_str(buffer, &exponent, 10, static_cast<std::size_t>(digits),
                       value.x, MPFR_RNDN);
          total += static_cast<std::size_t>(exponent);
        }
        sink = total;
      },
      kCount);
  mpfr_t parsed;
  mpfr_init2(parsed, P);
  const double parse_mpfr = BestNanoseconds(
      [&] {
        std::size_t total = 0;
        for (const std::string& s : strings) {
          mpfr_strtofr(parsed, s.c_str(), nullptr, 10, MPFR_RNDN);
          total += static_cast<std::size_t>(mpfr_get_exp(parsed));
        }
        sink = total;
      },
      kCount);
  mpfr_clear(parsed);

  const double print_boost = BestNanoseconds(
      [&] {
        std::size_t total = 0;
        for (const Boost& x : boost_values) {
          total += x.str(digits, std::ios_base::scientific).size();
        }
        sink = total;
      },
      kCount);
  const double parse_boost = BestNanoseconds(
      [&] {
        std::size_t total = 0;
        for (const std::string& s : strings) {
          const Boost x(s);
          total += static_cast<std::size_t>(x.backend().exponent());
        }
        sink = total;
      },
      kCount);

  for (auto& value : mpfr_values) mpfr_clear(value.x);
  const char* range = wide ? "wide" : "near 1";
  std::printf("%-10s %-7s %-6s %9.0f %9.0f %9.0f\n", format, range, "mpfix",
              print_ours, shortest_ours, parse_ours);
  std::printf("%-10s %-7s %-6s %9.0f %9s %9.0f\n", format, range, "MPFR",
              print_mpfr, "", parse_mpfr);
  std::printf("%-10s %-7s %-6s %9.0f %9s %9.0f\n", format, range, "Boost",
              print_boost, "", parse_boost);
  std::fflush(stdout);
}

}  // namespace

int main() {
#if defined(_WIN32)
  SetThreadAffinityMask(GetCurrentThread(), 1u << 2);
  SetPriorityClass(GetCurrentProcess(), HIGH_PRIORITY_CLASS);
#endif
  std::printf("GMP %s, MPFR %s. Nanoseconds per conversion; lower is better.\n",
              gmp_version, mpfr_get_version());
  std::printf("%-10s %-7s %-6s %9s %9s %9s\n", "format", "values", "library",
              "print", "shortest", "parse");
  for (const bool wide : {false, true}) {
    ReportFormat<113, 15>("binary128", wide);
    ReportFormat<237, 19>("binary256", wide);
    ReportFormat<489, 23>("binary512", wide);
  }
  return 0;
}
