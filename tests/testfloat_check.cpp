// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Checks mpfix against Berkeley TestFloat: runs testfloat_gen for one
// function in every rounding mode it takes, with and without exact where that
// applies, and compares each result and its five flags with SoftFloat's.
//
//   testfloat_check <path to testfloat_gen> <function> [level]
//
// f16, f32, f64 and f128 are Float<11, 5>, Float<24, 8>, Float<53, 11> and
// binary128. Where IEEE 754 leaves a choice, either answer passes: the payload
// and sign of NaN results, the integer result of an invalid conversion, and
// invalid for fma(0, inf, qNaN). Tininess is detected after rounding, as
// mpfix does.

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <string_view>
#include <type_traits>
#include <vector>

#include "mpfix/float.h"

#if defined(_WIN32)
#define popen _popen
#define pclose _pclose
#endif

namespace {

using namespace mpfix;  // NOLINT

using F16 = Float<11, 5>;
using F32 = Float<24, 8>;
using F64 = Float<53, 11>;
using F128 = binary128;

std::uint64_t ParseHex(std::string_view text) {
  std::uint64_t value = 0;
  for (const char c : text) {
    const int digit = c <= '9' ? c - '0' : (c | 0x20) - 'a' + 10;
    value = value << 4 | static_cast<std::uint64_t>(digit);
  }
  return value;
}

template <class F>
F ReadFloat(std::string_view text) {
  typename F::Bits bits{};
  // 64 bits per limb, the least significant limb last in the text.
  std::size_t end = text.size();
  for (std::size_t i = 0; i < bits.size() && end > 0; ++i) {
    const std::size_t begin = end > 16 ? end - 16 : 0;
    bits[i] = ParseHex(text.substr(begin, end - begin));
    end = begin;
  }
  return F::from_bits(bits);
}

template <class F>
std::string Hex(const F& x) {
  std::string text;
  const auto bits = x.to_bits();
  constexpr int kDigits = (F::kPrecision + F::kExponentBits) / 4;
  for (int i = kDigits - 1; i >= 0; --i) {
    const auto limb = bits[static_cast<std::size_t>(i / 16)];
    text += "0123456789ABCDEF"[(limb >> (4 * (i % 16))) & 15];
  }
  return text;
}

template <class I>
I ReadInt(std::string_view text) {
  return static_cast<I>(ParseHex(text));
}

// TestFloat's flag bits: invalid, infinite (divide by zero), overflow,
// underflow and inexact.
unsigned TestFloatFlags(unsigned flags) {
  unsigned result = 0;
  if ((flags & kFlagInvalid) != 0) result |= 0x10;
  if ((flags & kFlagDivByZero) != 0) result |= 0x08;
  if ((flags & kFlagOverflow) != 0) result |= 0x04;
  if ((flags & kFlagUnderflow) != 0) result |= 0x02;
  if ((flags & kFlagInexact) != 0) result |= 0x01;
  return result;
}

template <class F>
bool SameFloat(const F& actual, const F& expected) {
  // Any quiet NaN stands for a NaN result.
  if (expected.is_nan()) return actual.is_nan() && !actual.is_signaling();
  return actual.to_bits() == expected.to_bits();
}

enum class Kind { kFloat, kInt, kBool };

struct Function {
  std::string name;
  bool rounds = true;     // Takes a rounding mode.
  bool has_exact = false;  // Takes -exact and -notexact.
};

struct Checker {
  std::string function;
  Rounding mode = Rounding::kNearestEven;
  bool exact = false;
  const char* mode_name = "";
  long cases = 0;
  long failures = 0;

  void Fail(const std::vector<std::string_view>& fields,
            const std::string& actual, unsigned actual_flags) {
    ++failures;
    if (failures > 20) return;
    std::string line;
    for (const auto field : fields) {
      line += field;
      line += ' ';
    }
    std::printf("FAIL %s %s%s: %s-> %s %02X\n", function.c_str(), mode_name,
                exact ? " -exact" : "", line.c_str(), actual.c_str(),
                actual_flags);
  }

  // Compares a float result and flags, the last two fields.
  template <class F>
  // The flags are read here, after the result is computed: an argument
  // `context.flags` next to the call that sets them could be read first.
  void Expect(const std::vector<std::string_view>& fields, const F& actual,
              const Context& context, bool invalid_optional = false) {
    const unsigned flags = context.flags;
    const F expected = ReadFloat<F>(fields[fields.size() - 2]);
    unsigned expected_flags =
        static_cast<unsigned>(ParseHex(fields[fields.size() - 1]));
    unsigned actual_flags = TestFloatFlags(flags);
    if (invalid_optional) {
      expected_flags &= ~0x10u;
      actual_flags &= ~0x10u;
    }
    if (!SameFloat(actual, expected) || actual_flags != expected_flags) {
      Fail(fields, Hex(actual), actual_flags);
    }
  }

  template <class I>
  void ExpectInt(const std::vector<std::string_view>& fields, I actual,
                 const Context& context) {
    const unsigned flags = context.flags;
    const I expected = ReadInt<I>(fields[1]);
    const unsigned expected_flags =
        static_cast<unsigned>(ParseHex(fields[2]));
    const unsigned actual_flags = TestFloatFlags(flags);
    // The value of an invalid conversion is not specified.
    const bool value_ok = (expected_flags & 0x10) != 0 || actual == expected;
    if (!value_ok || actual_flags != expected_flags) {
      Fail(fields, std::to_string(actual), actual_flags);
    }
  }

  void ExpectBool(const std::vector<std::string_view>& fields, bool actual,
                  const Context& context) {
    const unsigned flags = context.flags;
    const bool expected = fields[2] == "1";
    const unsigned expected_flags =
        static_cast<unsigned>(ParseHex(fields[3]));
    const unsigned actual_flags = TestFloatFlags(flags);
    if (actual != expected || actual_flags != expected_flags) {
      Fail(fields, actual ? "1" : "0", actual_flags);
    }
  }
};

template <class F>
bool FloatOperation(Checker& check, std::string_view op,
                    const std::vector<std::string_view>& f) {
  Context context{check.mode};
  const auto a = [&] { return ReadFloat<F>(f[0]); };
  const auto b = [&] { return ReadFloat<F>(f[1]); };
  if (op == "add") {
    check.Expect(f, add(a(), b(), context), context);
  } else if (op == "sub") {
    check.Expect(f, sub(a(), b(), context), context);
  } else if (op == "mul") {
    check.Expect(f, mul(a(), b(), context), context);
  } else if (op == "div") {
    check.Expect(f, div(a(), b(), context), context);
  } else if (op == "rem") {
    check.Expect(f, remainder(a(), b(), context), context);
  } else if (op == "sqrt") {
    check.Expect(f, sqrt(a(), context), context);
  } else if (op == "mulAdd") {
    const F x = a(), y = b(), z = ReadFloat<F>(f[2]);
    // IEEE 754 leaves invalid for 0 * inf + qNaN to the implementation.
    const bool zero_times_inf = (x.is_zero() && y.is_inf()) ||
                                (x.is_inf() && y.is_zero());
    const bool optional = zero_times_inf && z.is_nan() && !z.is_signaling();
    check.Expect(f, fma(x, y, z, context), context, optional);
  } else if (op == "roundToInt") {
    const F x = a();
    const F z = check.exact ? round_to_integral_exact(x, context)
                            : round_to_integral(x, context);
    check.Expect(f, z, context);
  } else if (op == "eq" || op == "le" || op == "lt" || op == "eq_signaling" ||
             op == "le_quiet" || op == "lt_quiet") {
    const bool signaling =
        op == "le" || op == "lt" || op == "eq_signaling";
    const auto order = signaling ? compare_signaling(a(), b(), context)
                                 : compare_quiet(a(), b(), context);
    bool result;
    if (op.substr(0, 2) == "eq") {
      result = order == 0;
    } else if (op.substr(0, 2) == "le") {
      result = order <= 0;
    } else {
      result = order < 0;
    }
    check.ExpectBool(f, result, context);
  } else {
    return false;
  }
  return true;
}

// Calls `body` with a value of the float type named `name`.
template <class Body>
bool WithFloatType(std::string_view name, Body&& body) {
  if (name == "f16") return body(F16{}), true;
  if (name == "f32") return body(F32{}), true;
  if (name == "f64") return body(F64{}), true;
  if (name == "f128") return body(F128{}), true;
  return false;
}

template <class Body>
bool WithIntType(std::string_view name, Body&& body) {
  if (name == "i32") return body(std::int32_t{}), true;
  if (name == "i64") return body(std::int64_t{}), true;
  if (name == "ui32") return body(std::uint32_t{}), true;
  if (name == "ui64") return body(std::uint64_t{}), true;
  return false;
}

// Checks one line; false for a function it does not know.
bool CheckLine(Checker& check, const std::vector<std::string_view>& f) {
  const std::string_view name = check.function;
  const auto underscore = name.find('_');
  const std::string_view from = name.substr(0, underscore);
  const std::string_view rest = name.substr(underscore + 1);
  if (rest.substr(0, 3) == "to_") {
    const std::string_view to = rest.substr(3);
    bool known = false;
    // Float to float.
    WithFloatType(from, [&](auto source) {
      using From = decltype(source);
      WithFloatType(to, [&](auto target) {
        using To = decltype(target);
        Context context{check.mode};
        check.Expect(f, convert<To>(ReadFloat<From>(f[0]), context),
                     context);
        known = true;
      });
      // Float to integer.
      WithIntType(to, [&](auto target) {
        using I = decltype(target);
        Context context;
        const From x = ReadFloat<From>(f[0]);
        const I z = check.exact ? to_integer_exact<I>(x, check.mode, context)
                                : to_integer<I>(x, check.mode, context);
        check.ExpectInt(f, z, context);
        known = true;
      });
    });
    // Integer to float.
    WithIntType(from, [&](auto source) {
      using I = decltype(source);
      WithFloatType(to, [&](auto target) {
        using To = decltype(target);
        Context context{check.mode};
        check.Expect(f, convert<To>(ReadInt<I>(f[0]), context),
                     context);
        known = true;
      });
    });
    return known;
  }
  bool known = false;
  WithFloatType(from, [&](auto value) {
    known = FloatOperation<decltype(value)>(check, rest, f);
  });
  return known;
}

Function Describe(const std::string& name) {
  Function function{name};
  const auto has = [&](std::string_view part) {
    return name.find(part) != std::string::npos;
  };
  const bool to_int = has("_to_i") || has("_to_ui");
  if (has("_eq") || has("_le") || has("_lt") || has("_rem")) {
    function.rounds = false;
  }
  if (to_int || has("_roundToInt")) function.has_exact = true;
  return function;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr, "usage: %s <testfloat_gen> <function> [level]\n",
                 argv[0]);
    return 2;
  }
  const std::string generator = argv[1];
  const Function function = Describe(argv[2]);
  const std::string level = argc > 3 ? argv[3] : "1";

  struct Mode {
    const char* option;
    Rounding rounding;
  };
  const Mode modes[] = {{"-rnear_even", Rounding::kNearestEven},
                        {"-rnear_maxMag", Rounding::kNearestAway},
                        {"-rminMag", Rounding::kTowardZero},
                        {"-rmin", Rounding::kDown},
                        {"-rmax", Rounding::kUp}};
  long cases = 0;
  long failures = 0;
  for (const Mode& mode : modes) {
    for (const bool exact : {false, true}) {
      if (exact && !function.has_exact) continue;
      std::string command = "\"" + generator + "\" -level " + level +
                            " -tininessafter " + mode.option +
                            (exact ? " -exact " : " -notexact ") +
                            function.name;
      std::FILE* pipe = popen(command.c_str(), "r");
      if (pipe == nullptr) {
        std::fprintf(stderr, "cannot run %s\n", command.c_str());
        return 2;
      }
      Checker check{function.name, mode.rounding, exact, mode.option};
      char buffer[512];
      while (std::fgets(buffer, sizeof buffer, pipe) != nullptr) {
        std::vector<std::string_view> fields;
        std::string_view line(buffer);
        while (!line.empty()) {
          const auto start = line.find_first_not_of(" \r\n");
          if (start == std::string_view::npos) break;
          line.remove_prefix(start);
          const auto end = line.find_first_of(" \r\n");
          fields.push_back(line.substr(0, end));
          if (end == std::string_view::npos) break;
          line.remove_prefix(end);
        }
        if (fields.size() < 3) continue;
        ++check.cases;
        if (!CheckLine(check, fields)) {
          std::fprintf(stderr, "unknown function %s\n", function.name.c_str());
          pclose(pipe);
          return 2;
        }
      }
      if (pclose(pipe) != 0 || check.cases == 0) {
        std::fprintf(stderr, "%s failed or gave no cases\n", command.c_str());
        return 2;
      }
      cases += check.cases;
      failures += check.failures;
      if (!function.rounds) break;
    }
    if (!function.rounds) break;
  }
  std::printf("%s: %ld cases, %ld failures\n", function.name.c_str(), cases,
              failures);
  return failures == 0 ? 0 : 1;
}
