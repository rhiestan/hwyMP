// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Confirms that GMP and MPFR link and work. They are test and benchmark
// dependencies only: the mpfix library itself must never use them.

#include <gmp.h>
#include <gtest/gtest.h>
#include <mpfr.h>

#include <cstdio>

namespace {

TEST(Oracle, ReportsVersions) {
  std::printf("GMP %s (limb bits %d), MPFR %s\n", gmp_version, mp_bits_per_limb,
              mpfr_get_version());
  EXPECT_EQ(mp_bits_per_limb, 64);
}

TEST(Oracle, Sqrt2RoundsCorrectlyAtBinary256Precision) {
  mpfr_t x, square;
  mpfr_init2(x, 237);
  mpfr_init2(square, 474);

  const int ternary = mpfr_sqrt_ui(x, 2, MPFR_RNDN);
  EXPECT_NE(ternary, 0);

  // The exact square of a 237-bit value fits in 474 bits.
  EXPECT_EQ(mpfr_sqr(square, x, MPFR_RNDN), 0);
  mpfr_sub_ui(square, square, 2, MPFR_RNDN);
  mpfr_abs(square, square, MPFR_RNDN);
  EXPECT_LT(mpfr_get_exp(square), -230);

  mpfr_clear(square);
  mpfr_clear(x);
  mpfr_free_cache();
}

}  // namespace
