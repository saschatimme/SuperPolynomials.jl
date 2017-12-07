export PolynomialSystem, genevaluate!, genjacobian!

"""
    PolynomialSystem([T, ], F::Vector{<:MP.AbstractPolynomialLike})

Construct a system of SuperPolynomials.

    PolynomialSystem([T, ], F::Vector{<:FP.Polynomial})
"""
struct PolynomialSystem{T}
    polynomials::Vector{Polynomial{T, N, NTerms, Exponents} where Exponents where NTerms where N}
end

function PolynomialSystem(ps::Vector{P}) where {T, N, P<:Polynomial{T, N}}
    @assert !isempty(ps) "The system cannot be empty"
    @assert length(ps) < 129 "Currently only systems of up to 128 polynomials are supported."
    PolynomialSystem{T}(ps)
end

function PolynomialSystem(ps::Vector{<:MP.AbstractPolynomialLike})
    variables = sort!(union(Iterators.flatten(MP.variables.(ps))), rev=true)
    PolynomialSystem([Polynomial(p, variables) for p in ps])
end

function PolynomialSystem(ps::Vector{<:FP.Polynomial})
    PolynomialSystem([Polynomial(p) for p in ps])
end

function PolynomialSystem(::Type{T}, ps::Vector{<:MP.AbstractPolynomialLike}) where T
    variables = sort!(union(Iterators.flatten(MP.variables.(ps))), rev=true)
    PolynomialSystem([Polynomial(T, p, variables) for p in ps])
end

for N = 1:128
    @eval begin
        function genevaluate!($([Symbol("p", i) for i=1:N]...))
            (u, x) -> begin
                $(Expr(:block, [:(u[$i] = evaluate($(Symbol("p", i)), x)) for i=1:N]...))
                u
            end
        end
    end

    @eval begin
        function genjacobian!($([Symbol("p", i) for i=1:N]...))
            (U, x) -> begin
                $(Expr(:block, [:(gradient!(U, $(Symbol("p", i)), x, $i)) for i=1:N]...))
                U
            end
        end
    end
end

"""
    genevaluate!(F::PolynomialSystem)

Generate an evaluation function `(u, x) ->  u .= F(x)` for the polynomial system `F`.
"""
function genevaluate!(ps::PolynomialSystem)
    genevaluate!(ps.polynomials...)
end

"""
    genevaluate!(F::PolynomialSystem)

Generate an evaluation function `(u, x) ->  u .= J_F(x)` for the Jacobian of the
polynomial system `F`.
"""
function genjacobian!(ps::PolynomialSystem)
    genjacobian!(ps.polynomials...)
end
