// Copyright 2026 The hwyMP Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Confirms that Highway and Eigen are wired into the build on every platform.

#include <gtest/gtest.h>

#include <Eigen/Dense>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "hwy/highway.h"

namespace {

TEST(Dependencies, HighwayVectorAddition) {
  namespace hn = hwy::HWY_NAMESPACE;
  const hn::ScalableTag<std::uint64_t> d;
  const std::size_t lanes = hn::Lanes(d);
  ASSERT_GE(lanes, 1u);

  std::vector<std::uint64_t> values(lanes, 20);
  const auto sum = hn::Add(hn::LoadU(d, values.data()), hn::Set(d, 22));
  hn::StoreU(sum, d, values.data());
  for (const std::uint64_t value : values) EXPECT_EQ(value, 42u);
}

TEST(Dependencies, EigenSolve) {
  Eigen::Matrix2d matrix;
  matrix << 2, 0, 0, 4;
  const Eigen::Vector2d solution =
      matrix.partialPivLu().solve(Eigen::Vector2d(2, 8));
  EXPECT_DOUBLE_EQ(solution(0), 1.0);
  EXPECT_DOUBLE_EQ(solution(1), 2.0);
}

}  // namespace
