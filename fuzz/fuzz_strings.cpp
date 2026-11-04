// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// libFuzzer target for the conversions to and from text, in double's format
// and in binary128. Arbitrary text must parse without crashing, and plain
// decimal numbers must round as MPFR rounds them, in every mode, with the
// same flags. The bytes, read as a value, must print in shortest form and in
// hexadecimal so that the text reads back exactly, and in scientific notation
// with a precision exactly as MPFR prints it.

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <string_view>

#include "mpfix/charconv.h"
#include "oracle.h"

namespace {

using namespace mpfix;        // NOLINT
using namespace mpfix::test;  // NOLINT

[[noreturn]] void Fail(const char* what, int precision, std::string_view text,
                       const std::string& detail) {
  std::fprintf(stderr, "%s, P=%d: \"%.*s\" %s\n", what, precision,
               static_cast<int>(text.size()), text.data(), detail.c_str());
  std::abort();
}

// A decimal number that MPFR and from_chars read alike: an optional minus
// sign, digits with at most one point, and an exponent of up to six digits.
bool IsPlainDecimal(std::string_view text) {
  std::size_t i = 0;
  if (i < text.size() && text[i] == '-') ++i;
  int digits = 0;
  bool point = false;
  for (; i < text.size(); ++i) {
    if (text[i] >= '0' && text[i] <= '9') {
      ++digits;
    } else if (text[i] == '.' && !point) {
      point = true;
    } else {
      break;
    }
  }
  if (digits == 0) return false;
  if (i == text.size()) return true;
  if (text[i] != 'e' && text[i] != 'E') return false;
  ++i;
  if (i < text.size() && (text[i] == '-' || text[i] == '+')) ++i;
  const std::size_t exponent_digits = text.size() - i;
  if (exponent_digits == 0 || exponent_digits > 6) return false;
  for (; i < text.size(); ++i) {
    if (text[i] < '0' || text[i] > '9') return false;
  }
  return true;
}

template <class F>
void CheckParse(std::string_view text, Rounding mode) {
  F value = F::Zero();
  Context context{mode};
  const auto [end, ec] = from_chars(text.data(), text.data() + text.size(),
                                    value, std::chars_format::general,
                                    context);
  if (ec != std::errc{}) return;
  const std::string parsed(
      text.substr(0, static_cast<std::size_t>(end - text.data())));
  if (!IsPlainDecimal(parsed)) return;
  const Expected<F> expected = RoundExact<F>(
      [&](mpfr_ptr result, mpfr_rnd_t rnd) {
        return mpfr_strtofr(result, parsed.c_str(), nullptr, 10, rnd);
      },
      mode);
  if (!SameResult(value, expected.value) || context.flags != expected.flags) {
    Fail("parse", F::kPrecision, parsed,
         Hex(value) + " flags " + std::to_string(context.flags) + ", MPFR " +
             Hex(expected.value) + " flags " + std::to_string(expected.flags));
  }
}

template <class F>
void CheckPrint(const F& x, int precision) {
  if (!x.is_finite()) return;
  char buffer[2048];
  const auto reads_back = [&](std::chars_format fmt, const char* what) {
    const auto [end, ec] =
        fmt == std::chars_format::general
            ? to_chars(buffer, buffer + sizeof buffer, x)
            : to_chars(buffer, buffer + sizeof buffer, x, fmt);
    if (ec != std::errc{}) Fail(what, F::kPrecision, "", "too long");
    F back;
    const auto parsed = from_chars(buffer, end, back, fmt);
    if (parsed.ec != std::errc{} || parsed.ptr != end ||
        back.to_bits() != x.to_bits()) {
      Fail(what, F::kPrecision, std::string_view(buffer, end - buffer),
           "does not read back as " + Hex(x));
    }
  };
  reads_back(std::chars_format::general, "shortest");
  reads_back(std::chars_format::hex, "hex");

  const auto [end, ec] = to_chars(buffer, buffer + sizeof buffer, x,
                                  std::chars_format::scientific, precision);
  if (ec != std::errc{}) Fail("scientific", F::kPrecision, "", "too long");
  const View<F> view(x);
  char expected[2048];
  mpfr_snprintf(expected, sizeof expected, "%.*Re", precision, view.get());
  if (std::string_view(buffer, end - buffer) != expected) {
    Fail("scientific", F::kPrecision, std::string_view(buffer, end - buffer),
         std::string("MPFR prints ") + expected);
  }
}

template <class F>
F ValueFrom(const std::uint8_t* data, std::size_t size) {
  typename F::Bits bits{};
  std::size_t position = 0;
  for (auto& limb : bits) {
    for (int i = 0; i < 8; ++i) {
      limb = limb << 8 | (position < size ? data[position++] : 0);
    }
  }
  return F::from_bits(bits);
}

}  // namespace

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data,
                                      std::size_t size) {
  if (size == 0) return 0;
  const Rounding mode = kModes[data[0] % 5];
  const std::string_view text(reinterpret_cast<const char*>(data + 1),
                              size - 1 < 400 ? size - 1 : 400);
  CheckParse<Float<53, 11>>(text, mode);
  CheckParse<binary128>(text, mode);
  // Up to 60 digits after the point.
  const int precision = data[0] % 61;
  CheckPrint(ValueFrom<Float<53, 11>>(data + 1, size - 1), precision);
  CheckPrint(ValueFrom<binary128>(data + 1, size - 1), precision);
  return 0;
}
