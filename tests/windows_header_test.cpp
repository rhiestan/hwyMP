// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// The headers must compile after <windows.h> without NOMINMAX or
// WIN32_LEAN_AND_MEAN, which define macros such as min, max, near, far and
// small.

#if defined(_WIN32)
#include <windows.h>
#endif

#include <gtest/gtest.h>

#include "mpfix/charconv.h"
#include "mpfix/float.h"

namespace mpfix {
namespace {

TEST(WindowsHeader, Compiles) {
  const binary128 x = from_string<binary128>("0.1");
  EXPECT_EQ(to_string(x), "0.1");
  EXPECT_EQ(static_cast<int>(binary128(2.5)), 2);
}

}  // namespace
}  // namespace mpfix
