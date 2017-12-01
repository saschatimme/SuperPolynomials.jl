using SuperPolynomials
const SP = SuperPolynomials
using BenchmarkTools
using Base.Test
import FixedPolynomials
const FP = FixedPolynomials
import DynamicPolynomials: @polyvar
import Homotopies

@polyvar x y z

F = Homotopies.randomsystem(Float64, 5, 4, density=0.2)
f = F[1]
# f = FP.Polynomial{Float64}((5x^2*y^3+2y^2+3x^5+x-3))
g = SP.SuperPolynomial(f.coefficients, f.exponents)
w = rand(FP.nvariables(f))
cfg = FP.GradientConfig(f)

SP.horner6(g, w)
SP.evaluate(g, w)
FP.evaluate(f, w)
@btime SP.evaluate($g, $w)
@btime SP.horner8($g, $w)
@btime FP.evaluate(f, w)

function horner_impl(f::SP.SuperPolynomial{T, NVars, NTerms, Val{Exponents}}, x) where {T, NVars, NTerms, Exponents}
    matrix = SP.convert_to_matrix(NVars, NTerms, Exponents)
    E = [matrix[:,j] for j=1:NTerms]
    P = sortperm(E, lt=SP.revlexless, rev=true)
    I = SP.lower_set(E)

    d_I = NVars
    M = length(I)

    _xs = []
    for k=1:d_I
        xk = Symbol("x", k)
        push!(_xs, :($xk = x[$k]))
    end
    push!(_xs, :(@inbounds r0 = f.coefficients[$(P[1])]))
    r0_zero = false
    last_k = 0
    l = 2
    r = [Symbol("r", k) for k=1:d_I]
    rs_zero = trues(d_I)

    for k=1:d_I
        push!(_xs, :($(r[k]) = zero($T)))
    end

    for n=2:M
        # k = max{1≤j≤d_I : α(n)_j ≠ α(n-1)_j}
        k = 1
        for j = d_I:-1:1
            if I[n][j] != I[n-1][j]
                k = j
                break
            end
        end
        non_zero_r = r[filter(i -> !rs_zero[i], 1:k)]
        if r0_zero && length(non_zero_r) > 1
            summand = Expr(:call, :+, non_zero_r...)
        elseif r0_zero
            summand = :($(non_zero_r[1]))
        elseif isempty(non_zero_r)
            summand = :r0
        else
            summand = Expr(:call, :+, :r0, non_zero_r...)
        end
        xk = Symbol("x", k)
        push!(_xs, :($(r[k]) = $xk*($(summand))))
        rs_zero[k] = false
        if l ≤ NTerms && I[n] == E[P[l]]
            push!(_xs, :(@inbounds r0 = f.coefficients[$(P[l])]))
            r0_zero = false
            l += 1
        else
            if n == M
                push!(_xs, :(r0 = zero($T)))
            end
            r0_zero = true
        end
        for i = 1:k-1
            if n == M
                push!(_xs, :($(r[i]) = zero($T)))
            end
            rs_zero[i] = true
        end
    end

    Expr(:block,
        _xs...,
        Expr(:call, :+, :r0, r...))
end

@time lower_set(E)
