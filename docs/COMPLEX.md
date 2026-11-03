# Complex numbers

mpfix has its own complex type, `mpfix::Complex<F>` in `mpfix/complex.h`, with
the interface of `std::complex` and the special values of C99 Annex G. This
page says why, and what `std::complex<Float>` gives instead with each standard
library.

## mpfix::Complex

- Addition, subtraction and multiplication are correctly rounded in each
  component, and so are `abs` (hypot), `arg` (atan2) and `norm`.
- Division is within half an ulp plus 2^-(3P) ulp in each component.
- `sqrt`, `exp`, `log`, `log10`, `sin`, `cos`, `tan`, `sinh`, `cosh`, `tanh`,
  `asin`, `acos`, `atan`, `asinh`, `acosh` and `atanh` are evaluated in P + 64
  bits with formulas that do not cancel, then rounded once: within half an ulp
  plus about 2^-60 ulp in each component, and correctly rounded where they
  reduce to real functions, as on the axes. Against MPC they match in all but
  1 component in 1500 to 7500, depending on the format, which is off by one
  ulp.
- `pow` is exp(w log z), accurate relative to |z^w|; integer powers up to 64
  multiply instead.
- Everything rounds to nearest; there is no `Context`.

## std::complex\<Float\>

`std::complex<T>` is unspecified for T other than float, double and long
double, so each library does what its generic code happens to do.
`tests/std_complex_test.cpp` records it for `std::complex<binary128>`:

| Library | Compiles | Precision |
| --- | --- | --- |
| libstdc++ (GCC 16, and Clang on Linux) | everything | full, except `acos`: about 64 bits |
| libc++ (Xcode 16.4) | `+`, `-`, `*`, `/`, `conj`, construction | full for those |
| MSVC STL | deprecated (STL4037) | `log` does not compile; `exp`, `cos`, `sinh`, `atan2` and others compute in double |

- **libstdc++** finds mpfix's `sqrt`, `exp`, `atan2` and the rest by
  argument-dependent lookup. Its `acos` is π/2 − `asin` z with π/2 a
  `long double` constant, which caps its real part at about 64 bits.
- **libc++** calls `std::abs`, `std::hypot`, `std::exp`, `std::isinf` and the
  other `<cmath>` functions qualified, and those accept only built-in types:
  `abs`, `arg`, `norm`, `polar`, `sqrt`, `exp` and `pow` do not compile, and
  the other functions build on them.
- **MSVC's STL** deprecates `std::complex` of other types (warning STL4037),
  needs constexpr arithmetic that Float does not have for `log`, and its
  generic code converts to double for `exp`, `cos`, `sinh`, `atan2`, `isinf`
  and `signbit`, losing precision silently.

Use `mpfix::Complex` for portable, accurate results.

## With Eigen

`mpfix/eigen.h` makes `mpfix::Complex<F>` an Eigen complex scalar whose real
type is F, so real and complex scalars and matrices mix in expressions.
Matrices of `Complex<F>` work with the LU, QR, Cholesky and SVD
decompositions, `SelfAdjointEigenSolver` on Hermitian matrices, `ComplexSchur`
and `ComplexEigenSolver`; Eigen assigns integer literals to the scalar, so
`Complex` converts implicitly from built-in numbers.

Eigen hard-codes `std::complex<RealScalar>` as the complex type of real
matrices: `EigenSolver`, `RealSchur`'s and `MatrixBase::eigenvalues()`, the
matrix functions and `GeneralizedEigenSolver` return `std::complex<F>`. That
works with libstdc++, where `EigenSolver` on binary128 is tested, but not with
libc++ or MSVC's STL. For the eigenvalues of a real matrix there, cast it to
`Complex<F>` and use `ComplexEigenSolver`:

```cpp
Eigen::Matrix<mpfix::complex256, Eigen::Dynamic, Eigen::Dynamic> c =
    a.cast<mpfix::complex256>();
Eigen::ComplexEigenSolver<decltype(c)> solver(c);
```
