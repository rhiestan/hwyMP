// Copyright 2026 The mpfix Authors
// SPDX-License-Identifier: Apache-2.0 OR BSD-3-Clause

// A main for the fuzz targets without libFuzzer, as on MSVC: with file
// arguments it runs each file once, to replay inputs that libFuzzer found;
// otherwise it runs random inputs, half of them text of the characters that
// numbers are made of. MPFIX_TEST_SCALE multiplies their number.
//
//   fuzz_arithmetic_random [iterations] | fuzz_arithmetic_random file...

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

extern "C" int LLVMFuzzerTestOneInput(const std::uint8_t* data,
                                      std::size_t size);

int main(int argc, char** argv) {
  if (argc > 1 && std::string(argv[1]).find_first_not_of("0123456789") !=
                      std::string::npos) {
    for (int i = 1; i < argc; ++i) {
      std::ifstream file(argv[i], std::ios::binary);
      const std::vector<char> bytes((std::istreambuf_iterator<char>(file)),
                                    std::istreambuf_iterator<char>());
      LLVMFuzzerTestOneInput(
          reinterpret_cast<const std::uint8_t*>(bytes.data()), bytes.size());
    }
    return 0;
  }
  long iterations = argc > 1 ? std::atol(argv[1]) : 10000;
  if (const char* scale = std::getenv("MPFIX_TEST_SCALE")) {
    iterations *= std::atol(scale);
  }
  std::mt19937_64 random(12345);
  static constexpr char kNumberChars[] = "0123456789000000.eE-+";
  std::vector<std::uint8_t> data;
  for (long i = 0; i < iterations; ++i) {
    const std::size_t size = 1 + random() % 120;
    data.resize(size);
    const bool text = (i & 1) != 0;
    for (std::size_t j = 0; j < size; ++j) {
      data[j] = text && j > 0
                    ? static_cast<std::uint8_t>(
                          kNumberChars[random() % (sizeof kNumberChars - 1)])
                    : static_cast<std::uint8_t>(random());
    }
    LLVMFuzzerTestOneInput(data.data(), data.size());
  }
  std::printf("%ld inputs\n", iterations);
  return 0;
}
