export SuperPolynomialSystem, genevaluate!, genjacobian!

struct SuperPolynomialSystem{T}
    polynomials::Vector{SuperPolynomial{T, N, NTerms, Exponents} where Exponents where NTerms where N}
end

function SuperPolynomialSystem(ps::Vector{P}) where {T, N, P<:SuperPolynomial{T, N}}
    @assert !isempty(ps) "The system cannot be empty"
    @assert length(ps) < 129 "Currently only systems of up to 128 polynomials are supported."
    SuperPolynomialSystem{T}(ps)
end

function SuperPolynomialSystem(ps::Vector{<:MP.AbstractPolynomialLike})
    variables = sort!(union(Iterators.flatten(MP.variables.(ps))), rev=true)
    SuperPolynomialSystem([SuperPolynomial(p, variables) for p in ps])
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
    genevaluate!(F::SuperPolynomialSystem)

Generate an evaluation function `(u, x) ->  u .= F(x)` for the polynomial system `F`.
"""
function genevaluate!(ps::SuperPolynomialSystem)
    genevaluate!(ps.polynomials...)
end

"""
    genevaluate!(F::SuperPolynomialSystem)

Generate an evaluation function `(u, x) ->  u .= J_F(x)` for the Jacobian of the
polynomial system `F`.
"""
function genjacobian!(ps::SuperPolynomialSystem)
    genjacobian!(ps.polynomials...)
end
