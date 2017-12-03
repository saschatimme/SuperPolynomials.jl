# SuperPolynomials

Exploration of evaluating multivariate Polynomials using the `@generated` macro.

```julia
import SuperPolynomials
const SP = SuperPolynomials
import FixedPolynomials
const FP = FixedPolynomials
using BenchmarkTools

# We have 24 variables
N = 24
# with 20 terms
M = 20
# generate random exponents with degree between 0 and 5 (per variable)
exponents = round.(Int, 5 * rand(N, M))
# some random coefficients
coefficients = rand(M)
g = SP.SuperPolynomial(coefficients, exponents)
w = rand(N)

f = FP.Polynomial(exponents, coefficients)
cfg = FP.GradientConfig(f)

@btime SP.evaluate($g, $w) # 130.837 ns (0 allocations: 0 bytes)
@btime FP.evaluate($f, $w, $cfg) # 1.580 μs (0 allocations: 0 bytes)
```
