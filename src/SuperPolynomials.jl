__precompile__(true)

module SuperPolynomials
    export SuperPolynomial

    struct SuperPolynomial{T, NVars, NTerms, Exponents<:Val}
        coefficients::Vector{T}
    end

    function SuperPolynomial(coefficients::AbstractVector{T}, exponents::Matrix{<:Integer}) where T
        NVars = size(exponents, 1)
        @assert length(coefficients) == size(exponents, 2)

        exps, coeffs = normalize_exponents_coeffs(exponents, coefficients)
        NTerms = length(coeffs)
        Exponents = Val{ntuple(i -> exps[i], length(exps))}
        SuperPolynomial{T, NVars, NTerms, Exponents}(coeffs)
    end

    include("utilities.jl")
    include("evaluate.jl")
    include("horner.jl")
    include("gradient.jl")
end
