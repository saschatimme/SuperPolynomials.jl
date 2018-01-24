export evaluate

function evaluate_impl(f::Type{Polynomial{T, NVars, NTerms, Val{Exponents}}}) where {T, NVars, NTerms, Exponents}
    M = convert_to_matrix(NVars, NTerms, Exponents)

    grouped_powers = group_powers(M)

    as = [:(out = zero($T))]
    for j = 1:NTerms
        factors = []
        push!(factors, :(f.coefficients[$j]))
        nmultiplications = sum(M[:, j] .> 0)
        mult_counter = 0
        if nmultiplications == 0
            push!(as, :(@inbounds out += $(factors[1])))
        else
            for i=1:NVars
                k = M[i, j]
                if k > 0
                    mult_counter += 1
                    xik = x_((i, k))
                    if mult_counter == nmultiplications
                        a = batch_arithmetic_ops(:*, factors)
                        push!(as, :(@inbounds out = muladd($a, $xik, out)))
                    else
                        push!(factors, :($xik))
                    end
                end
            end
        end
    end

    Expr(:block,
        grouped_powers...,
        as...,
        :(out))
end

@generated function evaluate(f::Polynomial{T, NVars, NTerms, Val{Exponents}}, x) where {T, NVars, NTerms, Exponents}
    evaluate_impl(f)
end

@inline (f::Polynomial)(x) = evaluate(f, x)
