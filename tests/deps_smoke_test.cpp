// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Confirms that Eigen is wired into the build on every platform.

#include <gtest/gtest.h>

#include <Eigen/Dense>

namespace {

TEST(Dependencies, EigenSolve) {
  Eigen::Matrix2d matrix;
  matrix << 2, 0, 0, 4;
  const Eigen::Vector2d solution =
      matrix.partialPivLu().solve(Eigen::Vector2d(2, 8));
  EXPECT_DOUBLE_EQ(solution(0), 1.0);
  EXPECT_DOUBLE_EQ(solution(1), 2.0);
}

}  // namespace
