#!/usr/bin/env python3
# Copyright 2026 The hwyMP Authors
# SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause
"""Checks that the oracle tests catch deliberately injected bugs.

Each mutation changes one place in a temporary copy of the headers and rebuilds
an oracle test against it: float_oracle_test for the arithmetic, and
convert_oracle_test for the conversions. A mutant that survives means the tests
have a blind spot. Windows only for now, as it links the vcpkg libraries
directly.
Usage: mutation_check.py <build directory> [compiler] [jobs]
"""

import pathlib
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

ROOT = pathlib.Path(__file__).resolve().parent.parent

# (name, text to find, replacement). Each text must occur exactly once.
#
# Not listed, because no input is known to reach it: dropping the second borrow
# in NextRootLimb's correction. That needs an estimate two or more too high with
# a remainder of exactly zero above limb 2j. The code stays for correctness.
# Nor the overflow check in RoundSignificand: Finish overflows correctly by
# itself unless the exponent exceeds 32 bits, which MPFR on Windows cannot
# parse. charconv_test's ParseSyntax covers that case.
MUTATIONS = [
    ("ties round down", "return round_bit & (sticky | is_odd);",
     "return round_bit & sticky;"),
    ("ties always up", "return round_bit & (sticky | is_odd);",
     "return round_bit;"),
    ("ties away ignores round bit", "      return round_bit;\n    case Rounding::kTowardZero:",
     "      return round_bit && sticky;\n    case Rounding::kTowardZero:"),
    ("upward ignores sticky", "return !negative & (round_bit | sticky);",
     "return !negative & round_bit;"),
    ("downward uses wrong sign", "return negative & (round_bit | sticky);",
     "return !negative & (round_bit | sticky);"),
    ("rest ignored when rounding",
     "sticky = ((low & (half - 1)) | value.rest) != 0;",
     "sticky = (low & (half - 1)) != 0;"),
    ("carry lost after rounding", "if (carry != 0) value.limbs[kLimbs - 1] = kTopBit;\n    return carry != 0;",
     "if (carry != 0) value.limbs[kLimbs - 1] = kTopBit;\n    return false;"),
    ("overflow after rounding missed",
     "if (exponent > kEmax) [[unlikely]] return Overflow(negative, env);",
     "if (exponent > kEmax + 1) [[unlikely]] return Overflow(negative, env);"),
    ("overflow to infinity in all modes",
     "return to_infinity ? Infinity(negative) : MaxFinite(negative);",
     "return Infinity(negative);"),
    ("tininess before rounding",
     "tiny = !RoundToPrecision(negative, copy,\n                               detail::SilentEnv{env.rounding()});",
     "tiny = true;"),
    ("underflow raised when exact",
     "if (round_bit || sticky) {\n      env.raise(kFlagInexact);\n      if (tiny) env.raise(kFlagUnderflow);\n    }",
     "env.raise(kFlagInexact);\n    if (tiny) env.raise(kFlagUnderflow);"),
    ("subnormal carry keeps exponent",
     "if (round_up && IncrementAtBit(limbs, count)) ++exponent;",
     "if (round_up) IncrementAtBit(limbs, count);"),
    ("subnormal boundary off by one", "if (discard >= kTotalBits) {",
     "if (discard > kTotalBits) {"),
    ("alignment drops lost words",
     "for (int i = 0; i + 1 < words; ++i) {\n      lost |= lesser[static_cast<std::size_t>(i)];\n    }",
     ""),
    ("alignment drops lost bits",
     "lost |= (wide[0] << 1) << complement;",
     ""),
    ("far alignment threshold too small",
     "if (distance >= (kWords + 1) * kLimbBits) [[unlikely]] {",
     "if (distance >= kWords * kLimbBits) [[unlikely]] {"),
    ("addition carry drops shifted bit",
     "sum.rest = ((sum.limbs[0] << 1) << complement) | (sum.rest >> shift) |\n                 (sum.rest & 1);",
     "sum.rest = (sum.rest >> shift) | (sum.rest & 1);"),
    ("addition carry drops sticky",
     "sum.rest = ((sum.limbs[0] << 1) << complement) | (sum.rest >> shift) |\n                 (sum.rest & 1);",
     "sum.rest = ((sum.limbs[0] << 1) << complement) | (sum.rest >> shift);"),
    ("subtraction ignores rest borrow",
     "difference.rest = sbb(0, difference.rest, borrow);",
     "borrow = 0;"),
    ("cancellation sign not flipped",
     "lesser = temporary;\n        negative = !negative;",
     "lesser = temporary;"),
    ("single-bit cancellation exponent wrong",
     "      exponent -= shift;\n    }\n    return Finish(negative, exponent, difference, env);",
     "    }\n    return Finish(negative, exponent, difference, env);"),
    ("single-bit cancellation drops rest",
     "(value.limbs[0] << shift) | ((value.rest >> 1) >> complement);",
     "(value.limbs[0] << shift);"),
    ("deep cancellation taken for single-bit",
     "difference.limbs[kLimbs - 1] >> (kLimbBits - 2) == 0",
     "false"),
    ("cancellation exponent wrong",
     "exponent -= NormalizeDifference(difference);",
     "exponent -= NormalizeDifference(difference) - 1;"),
    ("exact zero sign in downward mode",
     "if (order == 0) return Zero(env.rounding() == Rounding::kDown);",
     "if (order == 0) return Zero(false);"),
    ("column carry lost",
     "column_high += carry;",
     ""),
    ("product low limbs not sticky",
     "for (std::size_t i = 0; i + 1 < kLimbs; ++i) lost |= product[i];",
     ""),
    ("row carry lost", "product[i + kLimbs] = carry;", ""),
    ("word-aligned shift drops low bits",
     "result.rest = (lesser[0] << 1) << complement;",
     "result.rest = 0;"),
    ("multiplication exponent off by one",
     "a.exponent_ + b.exponent_ + 1 - shift;",
     "a.exponent_ + b.exponent_ - shift;"),
    ("multiplying zero by infinity is valid",
     "if (a.is_zero() || b.is_zero()) return InvalidOperation(env);\n      return Infinity(negative);",
     "return Infinity(negative);"),
    ("infinities of opposite sign add",
     "if (b.is_inf() && a.negative_ != b_sign) return InvalidOperation(env);",
     ""),
    ("quiet NaN raises invalid", "if (a.is_signaling() || b.is_signaling()) env.raise(kFlagInvalid);\n    Float result",
     "env.raise(kFlagInvalid);\n    Float result"),
    ("NaN stays signaling", "SetQuietBit(result.limbs_);\n    return result;",
     "return result;"),
    ("zeros of opposite sign in downward mode",
     "return Zero(env.rounding() == Rounding::kDown);\n    }\n    if (a.is_zero()) {",
     "return Zero(false);\n    }\n    if (a.is_zero()) {"),
    ("subtrahend sign ignored for zero",
     "result.negative_ = b_sign;\n      return result;",
     "return result;"),
    # Division.
    ("division integer part inverted exponent",
     "a.exponent_ - b.exponent_ - 1 + integer_part;",
     "a.exponent_ - b.exponent_ + integer_part;"),
    ("division ignores remainder",
     "value.rest |= detail::IsZero(r) ? 0 : 1;",
     ""),
    ("division stops one limb short",
     "kSpareBits >= 1 ? kLimbs : kLimbs + 1;",
     "kLimbs;"),
    ("division skips the add-back",
     "while (over < 0) [[unlikely]] {\n        --q;",
     "while (false) {\n        --q;"),
    ("division add-back keeps quotient",
     "while (over < 0) [[unlikely]] {\n        --q;",
     "while (over < 0) [[unlikely]] {"),
    ("division ignores equal top limbs",
     "if (u2 == d1 && u1 == d0) [[unlikely]] {",
     "if (false) {"),
    ("division loses carry for equal top limbs",
     "over = carry;\n    } else {",
     "} else {"),
    ("division ignores low divisor limbs",
     "r[i] = sbb(r[i], low, borrow);",
     ""),
    ("division by zero raises nothing",
     "env.raise(kFlagDivByZero);\n      return Infinity(negative);",
     "return Infinity(negative);"),
    ("zero over zero is valid",
     "if (a.is_zero()) return InvalidOperation(env);\n      env.raise",
     "env.raise"),
    # Square root.
    ("root of negative is valid",
     "if (x.negative_ != 0) return InvalidOperation(env);",
     ""),
    ("root loses sign of negative zero",
     "if (x.is_zero()) return x;\n      if (x.negative_ != 0)",
     "if (x.is_zero()) return Zero(false);\n      if (x.negative_ != 0)"),
    ("root ignores odd exponent",
     "const int shift = (x.exponent_ & 1) ^ 1;",
     "const int shift = 1;"),
    ("root exponent truncates instead of floors",
     "Finish(false, x.exponent_ >> 1, value, env);",
     "Finish(false, x.exponent_ / 2, value, env);"),
    ("root rounds without half test",
     "value.rest = (above_half ? kTopBit : 0) | (detail::IsZero(rem) ? 0 : 1);",
     "value.rest = detail::IsZero(rem) ? 0 : 1;"),
    ("root ignores remainder",
     "value.rest = (above_half ? kTopBit : 0) | (detail::IsZero(rem) ? 0 : 1);",
     "value.rest = above_half ? kTopBit : 0;"),
    ("root half test ignores high limbs",
     "above_half = above_half || rem[i] != 0;",
     ""),
    ("root limb never corrected",
     "while (negative > 0) [[unlikely]] {",
     "while (false) {"),
    ("root estimate not clamped",
     "if (u2 < d1 || (u2 == d1 && u1 < d0)) {",
     "if (u2 <= d1) {"),
    ("root keeps a stale reciprocal",
     "if (j + 3 == kLimbs) {",
     "if (false) {"),
    ("integer root not corrected downward",
     "while (negative) {\n    // (s - 1)^2",
     "while (false) {\n    // (s - 1)^2"),
    # Fused multiply-add.
    ("fma sticky lost in alignment",
     "lesser[0] |= detail::IsZero(lost) ? 0 : 1;",
     ""),
    ("fma far operand vanishes",
     "lesser.fill(0);\n      lesser[0] = 1;",
     "lesser.fill(0);"),
    ("fma sum carry drops sticky",
     "total[0] |= sticky;",
     ""),
    ("fma sum carry keeps exponent",
     "total[kWide - 1] |= kTopBit;\n        ++exponent;",
     "total[kWide - 1] |= kTopBit;"),
    ("fma rest ignores low limbs",
     "for (std::size_t i = 0; i < kLimbs; ++i) lost |= total[i];",
     ""),
    ("fma keeps sign after swap",
     "return Finish(big_negative, exponent, value, env);",
     "return Finish(product_negative, exponent, value, env);"),
    ("fma exact zero ignores mode",
     "if (detail::IsZero(total)) return Zero(env.rounding() == Rounding::kDown);",
     "if (detail::IsZero(total)) return Zero(false);"),
    ("fma compares exponents only",
     "(big_exponent == small_exponent &&\n         detail::CompareLimbs(big, lesser) < 0)",
     "false"),
    ("fma zero times infinity is valid",
     "if (a.is_zero() || b.is_zero()) return InvalidOperation(env);\n      if (c.is_inf()",
     "if (c.is_inf()"),
    ("fma ignores opposite infinity",
     "if (c.is_inf() && (c.negative_ != 0) != product_negative) {\n        return InvalidOperation(env);\n      }",
     ""),
    ("fma ignores signaling third operand",
     "if (a.is_signaling() || b.is_signaling() || c.is_signaling()) {",
     "if (a.is_signaling() || b.is_signaling()) {"),
    ("fma zero sum ignores mode",
     "if ((c.negative_ != 0) == product_negative) return c;\n      return Zero(env.rounding() == Rounding::kDown);",
     "if ((c.negative_ != 0) == product_negative) return c;\n      return Zero(false);"),
    # Remainder.
    ("remainder ties round up",
     "if (order > 0 || (order == 0 && odd)) {",
     "if (order >= 0) {"),
    ("remainder ties round down",
     "if (order > 0 || (order == 0 && odd)) {",
     "if (order > 0) {"),
    ("remainder keeps sign when rounding up",
     "result_negative = !negative;",
     ""),
    ("remainder of half distance ignored",
     "if (detail::CompareLimbs(x.limbs_, divisor) <= 0) return x;",
     "return x;"),
    ("remainder half distance tie rounds up",
     "if (detail::CompareLimbs(x.limbs_, divisor) <= 0) return x;",
     "if (detail::CompareLimbs(x.limbs_, divisor) < 0) return x;"),
    ("remainder drops overflowed bit",
     "odd = overflow || detail::CompareLimbs(r, divisor) >= 0;",
     "odd = detail::CompareLimbs(r, divisor) >= 0;"),
    ("remainder parity from wrong step",
     "odd = (q & 1) != 0;",
     ""),
    ("remainder of zero loses sign",
     "if (detail::IsZero(r)) return Zero(negative);",
     "if (detail::IsZero(r)) return Zero(false);"),
    ("remainder ignores high top bit",
     "if ((r[kLimbs - 1] & kTopBit) == 0) {\n      Limbs doubled = r;",
     "{\n      Limbs doubled = r;"),
    # Rounding to integers, and exponents.
    ("integral rounds below one to zero",
     "const bool round_bit = x.exponent_ == -1;",
     "const bool round_bit = false;"),
    ("integral half is sticky",
     "const bool sticky = x.exponent_ < -1 || x.limbs_ != TopBitOnly();",
     "const bool sticky = true;"),
    ("integral below one loses sign",
     "return Zero(negative);\n    }\n\n    constexpr int kTotalBits",
     "return Zero(false);\n    }\n\n    constexpr int kTotalBits"),
    ("integral below one is exact",
     "const bool sticky = x.exponent_ < -1 || x.limbs_ != TopBitOnly();\n      env.raise(kFlagInexact);",
     "const bool sticky = x.exponent_ < -1 || x.limbs_ != TopBitOnly();"),
    ("integral inexact not raised",
     "if (!round_bit && !sticky) return x;\n    env.raise(kFlagInexact);",
     "if (!round_bit && !sticky) return x;"),
    ("integral silent form raises inexact",
     "context.flags |= flags & ~unsigned{kFlagInexact};",
     "context.flags |= flags;"),
    ("integral carry keeps exponent",
     "++exponent;\n      // Only formats",
     "// Only formats"),
    ("integral keeps low bits",
     "detail::ShiftRight(limbs, count);\n    detail::ShiftLeft(limbs, count);\n    std::int32_t exponent = x.exponent_;",
     "std::int32_t exponent = x.exponent_;"),
    ("scalbn ignores shift", "std::int64_t exponent = std::int64_t{x.exponent_} + n;",
     "std::int64_t exponent = std::int64_t{x.exponent_};"),
    ("scalbn quiets nothing",
     "if (x.is_nan()) return PropagateNaN(x, x, env);\n    if (!IsRegular(x)) return x;\n    // Clamping",
     "if (!IsRegular(x)) return x;\n    // Clamping"),
]


# Mutations of the conversions, caught by convert_oracle_test: (name, header,
# text to find, replacement).
CONVERSION_MUTATIONS = [
    ("rounding drops the next limb", "float.h",
     "value.rest = limbs[count - kLimbs - 1];", "value.rest = 0;"),
    ("rounding ignores lower limbs", "float.h",
     "sticky = sticky || limbs[i] != 0;", "(void)limbs[i];"),
    ("rounding ignores sticky", "float.h",
     "if (sticky) value.rest |= 1;", ""),
    ("signaling NaN converts silently", "float.h",
     "if (x.is_signaling()) env.raise(kFlagInvalid);\n      // Both",
     "// Both"),
    ("converted NaN not quieted", "float.h",
     "SetQuietBit(limbs);\n      return FromParts(negative, kExpNaN, limbs);",
     "return FromParts(negative, kExpNaN, limbs);"),
    ("negative integers lose their sign", "float.h",
     "if (negative) magnitude = limb_t{0} - magnitude;", ""),
    ("most negative integer out of range", "float.h",
     "static_cast<limb_t>((Limits::max)()) + (negative ? 1 : 0);",
     "static_cast<limb_t>((Limits::max)());"),
    ("convertToInteger raises inexact", "float.h",
     "if (exact) env.raise(local.flags & kFlagInexact);",
     "env.raise(local.flags & kFlagInexact);"),
    ("negative to unsigned allowed", "float.h",
     "if ((negative && !Limits::is_signed) || magnitude > limit) {",
     "if (magnitude > limit) {"),
    ("integer of NaN raises nothing", "float.h",
     "if (is_nan()) {\n      env.raise(kFlagInvalid);\n      return 0;",
     "if (is_nan()) {\n      return 0;"),
    ("fast parse rounds its upper bound down", "charconv.h",
     "const T high = mul(d_high, p_high, up);",
     "const T high = mul(d_high, p_high, down);"),
    ("fast parse always inexact", "charconv.h",
     "if (!EqualsDecimal(result, digits, t)) env.raise(kFlagInexact);",
     "env.raise(kFlagInexact);"),
    ("fast parse tiny without underflow", "charconv.h",
     "env.raise(kFlagInexact | kFlagUnderflow);",
     "env.raise(kFlagInexact);"),
    ("fast parse misses overflow", "charconv.h",
     "if (magnitude > Unbounded(F::MaxFinite())) {",
     "if (magnitude > Unbounded(F::MaxFinite()) && false) {"),
    ("trailing fraction zeros not counted", "charconv.h",
     "--end[1];\n    ++result.exponent;", "--end[1];"),
    ("eight digits misread", "charconv.h",
     "chars = chars * 10 + (chars >> 8);", "chars = chars * 10 + (chars >> 16);"),
    ("full chunk not flushed", "charconv.h",
     "// A full chunk takes no more digits.\n        if (chunk_digits == 19) {",
     "if (false) {"),
    ("parsed exponent sign lost", "charconv.h",
     "if (exponent_negative) exponent = -exponent;", ""),
    ("exact output ties round up", "charconv.h",
     "round_bit = half >= 0;", "round_bit = half > 0;"),
    ("scaled exactness off by one", "charconv.h",
     "if (t < 0) return e - t >= 0;", "if (t < 0) return e - t > 0;"),
    ("fast shortest ignores the closer lower end", "charconv.h",
     "const bool lower_closer = power_of_two && x.raw_exponent() > F::kEmin;",
     "const bool lower_closer = false;"),
    ("fast shortest ties round up", "charconv.h",
     "order = (floor_prefix.back() - '0') % 2 == 0 ? -1 : 1;", "order = 1;"),
    ("hexadecimal rounding parity", "charconv.h",
     ": fraction.bit(0);", ": !fraction.bit(0);"),
    ("Karatsuba forgets the low product", "bigint.h",
     "SubFrom(middle.data(), middle.size(), r, 2 * low);", ""),
    ("powers of five lose a factor", "bigint.h",
     "for (std::uint64_t i = 0; i < exponent % 27; ++i) factor *= 5;",
     "for (std::uint64_t i = 1; i < exponent % 27; ++i) factor *= 5;"),
]

HEADERS = ["float.h", "limb_ops.h", "bigint.h", "charconv.h"]


def run_mutant(mutation, sources, build_flags, compiler, work):
    name, header, old, new = mutation
    if sources[header].count(old) != 1:
        return name, "BAD PATTERN"
    directory = work / str(abs(hash(name)))
    (directory / "hwymp").mkdir(parents=True)
    for file in HEADERS:
        text = sources[file]
        if file == header:
            text = text.replace(old, new)
        (directory / "hwymp" / file).write_text(text, encoding="utf8")
    executable = directory / "mutant.exe"
    command = [compiler, "-std=c++20", "-O1", f"-I{directory}", *build_flags,
               "-o", str(executable)]
    compiled = subprocess.run(command, capture_output=True, text=True)
    if compiled.returncode != 0:
        return name, "COMPILE ERROR\n" + compiled.stderr[-600:]
    try:
        # Generous, since a mutant that only slows down must not count as
        # caught.
        result = subprocess.run([str(executable)], capture_output=True,
                                text=True, timeout=900)
    except subprocess.TimeoutExpired:
        return name, "caught (timeout)"
    return name, "caught" if result.returncode != 0 else "SURVIVED"


def main():
    build = pathlib.Path(sys.argv[1]).resolve()
    compiler = sys.argv[2] if len(sys.argv) > 2 else "clang++"
    jobs = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    installed = build.parent / "vcpkg_installed" / "x64-windows-static-md-v143"
    deps = build / "_deps" / "googletest-src" / "googletest"
    if not deps.exists():
        deps = ROOT.parent / "does-not-exist"
    libraries = build / "lib"
    common = [
        f"-I{deps / 'include'}", f"-I{installed / 'include'}",
        "-D_CRT_SECURE_NO_WARNINGS", "-D_DLL", "-Xlinker", "/ignore:4099",
        str(libraries / "gtest_main.lib"), str(libraries / "gtest.lib"),
        str(installed / "lib" / "mpfr.lib"), str(installed / "lib" / "gmp.lib"),
        "-lmsvcrt", "-lucrt", "-lvcruntime", "-fuse-ld=lld",
    ]
    arithmetic = [str(ROOT / "tests" / "float_oracle_test.cpp"), *common]
    conversions = [str(ROOT / "tests" / "convert_oracle_test.cpp"), *common]
    sources = {file: (ROOT / "include" / "hwymp" / file).read_text(
        encoding="utf8") for file in HEADERS}
    runs = [((name, "float.h", old, new), arithmetic)
            for name, old, new in MUTATIONS]
    runs += [(mutation, conversions) for mutation in CONVERSION_MUTATIONS]
    with tempfile.TemporaryDirectory() as temporary:
        work = pathlib.Path(temporary)
        # Unchanged builds must pass, or "caught" would mean nothing.
        anchor = MUTATIONS[0][1]
        for index, flags in enumerate((arithmetic, conversions)):
            control = run_mutant(("control", "float.h", anchor, anchor),
                                 sources, flags, compiler,
                                 work / f"control{index}")
            if control[1] != "SURVIVED":
                print(f"Control build did not pass: {control[1]}")
                return 1
        with ThreadPoolExecutor(jobs) as pool:
            results = list(pool.map(
                lambda run: run_mutant(run[0], sources, run[1], compiler, work),
                runs))
    failures = 0
    for name, outcome in results:
        print(f"{outcome:20} {name}")
        failures += not outcome.startswith("caught")
    print(f"{len(results) - failures}/{len(results)} mutants caught")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
