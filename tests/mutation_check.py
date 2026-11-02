#!/usr/bin/env python3
# Copyright 2026 The hwyMP Authors
# SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause
"""Checks that the oracle tests catch deliberately injected bugs.

Each mutation changes one line of a temporary copy of float.h and rebuilds the
oracle test against it. A mutant that survives means the tests have a blind
spot. Windows only for now, as it links the vcpkg libraries directly.
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
     "limbs[kLimbs - 1] = kTopBit;\n        ++exponent;",
     "limbs[kLimbs - 1] = kTopBit;"),
    ("subnormal boundary off by one", "if (discard >= kTotalBits) {",
     "if (discard > kTotalBits) {"),
    ("alignment drops lost words",
     "for (int i = 0; i + 1 < words; ++i) {\n      lost |= small[static_cast<std::size_t>(i)];\n    }",
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
     "small = temporary;\n        negative = !negative;",
     "small = temporary;"),
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
     "result.rest = (small[0] << 1) << complement;",
     "result.rest = 0;"),
    ("multiplication exponent off by one",
     "a.exponent_ + b.exponent_ + 1 - shift;",
     "a.exponent_ + b.exponent_ - shift;"),
    ("multiplying zero by infinity is valid",
     "if (a.is_zero() || b.is_zero()) return InvalidOperation(env);",
     ""),
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
]


def run_mutant(mutation, source, build_flags, compiler, work):
    name, old, new = mutation
    if source.count(old) != 1:
        return name, "BAD PATTERN"
    directory = work / str(abs(hash(name)))
    (directory / "hwymp").mkdir(parents=True)
    (directory / "hwymp" / "float.h").write_text(source.replace(old, new),
                                                 encoding="utf8")
    shutil.copy(ROOT / "include" / "hwymp" / "limb_ops.h", directory / "hwymp")
    executable = directory / "mutant.exe"
    command = [compiler, "-std=c++20", "-O1", f"-I{directory}", *build_flags,
               "-o", str(executable)]
    compiled = subprocess.run(command, capture_output=True, text=True)
    if compiled.returncode != 0:
        return name, "COMPILE ERROR\n" + compiled.stderr[-600:]
    try:
        result = subprocess.run([str(executable)], capture_output=True,
                                text=True, timeout=120)
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
    flags = [
        str(ROOT / "tests" / "float_oracle_test.cpp"),
        f"-I{deps / 'include'}", f"-I{installed / 'include'}",
        "-D_CRT_SECURE_NO_WARNINGS", "-D_DLL", "-Xlinker", "/ignore:4099",
        str(libraries / "gtest_main.lib"), str(libraries / "gtest.lib"),
        str(installed / "lib" / "mpfr.lib"), str(installed / "lib" / "gmp.lib"),
        "-lmsvcrt", "-lucrt", "-lvcruntime", "-fuse-ld=lld",
    ]
    source = (ROOT / "include" / "hwymp" / "float.h").read_text(encoding="utf8")
    with tempfile.TemporaryDirectory() as temporary:
        work = pathlib.Path(temporary)
        # An unchanged build must pass, or "caught" would mean nothing.
        anchor = MUTATIONS[0][1]
        control = run_mutant(("control", anchor, anchor), source, flags,
                             compiler, work)
        if control[1] != "SURVIVED":
            print(f"Control build did not pass: {control[1]}")
            return 1
        with ThreadPoolExecutor(jobs) as pool:
            results = list(pool.map(
                lambda m: run_mutant(m, source, flags, compiler, work),
                MUTATIONS))
    failures = 0
    for name, outcome in results:
        print(f"{outcome:20} {name}")
        failures += not outcome.startswith("caught")
    print(f"{len(results) - failures}/{len(results)} mutants caught")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
