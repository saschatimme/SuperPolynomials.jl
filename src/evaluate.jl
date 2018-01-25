export evaluate

function evaluate_impl(::Type{S}, f::Type{Polynomial{T, NVars, NTerms, Val{Exponents}}}) where {S, T, NVars, NTerms, Exponents}
    M = convert_to_matrix(NVars, NTerms, Exponents)

    G = ComputeGraph()
    pow_instr = all_powers_of_variables!(G, M)

    j_prods = products_to_evaluate(M)
    prods = last.(j_prods)
    res = group_products!(G, prods)

    graph_instr = construct_instructions(G, j_prods) do j_prod, instr
        (:(@inbounds out = muladd(c[$(j_prod[1])], $instr, out)))
    end

    Expr(:block,
        :(out = zero($(promote_type(S, T)))),
        :(c = f.coefficients),
        pow_instr...,
        graph_instr...,
        :(out))
end

"""
    products_to_evaluate(exponents)

Compute a list of all products to evaluate in the form
    (j, factors)

where `i` is the partial derivative, `j` is the corresponding term, `k` is the factor from the derivative
and `factors` is a vector of `Int`s with the occuring monomials.
"""
function products_to_evaluate(exponents)
    out = Vector{Tuple{Int, Vector{Symbol}}}()
    m, n = size(exponents)
    for j=1:n
        factors = Symbol[]
        for l=1:m
            k = exponents[l,j]
            if k > 0
                push!(factors, x_((l, k)))
            end
        end
        push!(out, (j, factors))
    end
    out
end


@generated function evaluate(f::Polynomial{T, NVars, NTerms, Val{Exponents}}, x::AbstractVector{S}) where {S, T, NVars, NTerms, Exponents}
    evaluate_impl(S, f)
end

(f::Polynomial)(x) = evaluate(f, x)
