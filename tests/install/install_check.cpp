// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// Uses each public header of an installed mpfix.

#include <cstdio>

#include "mpfix/charconv.h"
#include "mpfix/complex.h"
#include "mpfix/math.h"

int main() {
  const mpfix::binary128 two(2);
  const mpfix::complex128 z(two, two);
  const std::string text = mpfix::to_string(sqrt(two));
  std::printf("%s %s\n", text.c_str(),
              mpfix::to_string(abs(z) * mpfix::exp(two)).c_str());
  return text == "1.414213562373095048801688724209698" ? 0 : 1;
}
