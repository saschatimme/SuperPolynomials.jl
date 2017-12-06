__precompile__(true)

module SuperPolynomials
    export SuperPolynomial

    import MultivariatePolynomials
    const MP = MultivariatePolynomials

    import FixedPolynomials
    const FP = FixedPolynomials

    struct SuperPolynomial{T, NVars, NTerms, Exponents<:Val}
        coefficients::Vector{T}
    end

    function SuperPolynomial(coefficients::Vector{T}, exponents::Vector{<:Integer}, NVars) where T
        @assert length(coefficients)*NVars == length(exponents)
        NTerms = length(coefficients)
        Exponents = Val{ntuple(i -> exponents[i], length(exponents))}
        SuperPolynomial{T, NVars, NTerms, Exponents}(coefficients)
    end

    function SuperPolynomial(coefficients::Vector{T}, exponents::Matrix{<:Integer}) where T
        NVars = size(exponents, 1)
        exps, coeffs = normalize_exponents_coeffs(exponents, coefficients)
        SuperPolynomial(coeffs, exps, NVars)
    end

    function SuperPolynomial(p::MP.AbstractPolynomialLike, variables = MP.variables(p))
        exps, coeffs = mp_exponents_coefficients(p, variables)
        SuperPolynomial(coeffs, exps, length(variables))
    end

    SuperPolynomial(p::FP.Polynomial) = SuperPolynomial(p.coefficients, p.exponents)


    include("utilities.jl")
    include("evaluate.jl")
    include("horner.jl")
    include("gradient.jl")


    export exponents

    """
        exponents(f)

    Get the exponents of `f` as a matrix where each column represents a monomial of the
    polynomial.
    """
    function exponents(::SuperPolynomial{T, NVars, NTerms, Val{Exponents}}) where {T, NVars, NTerms, Exponents}
        convert_to_matrix(NVars, NTerms, Exponents)
    end

    """
        coefficients(f)

    Get the coefficients of `f`.
    """
    function coefficients(f::SuperPolynomial{T, NVars, NTerms, Val{Exponents}}) where {T, NVars, NTerms, Exponents}
        f.coefficients
    end

    """
        scale_coefficients!(f, λ)

    Get the coefficients of `f`.
    """
    function scale_coefficients!(f::SuperPolynomial{T, NVars, NTerms, Val{Exponents}}, λ) where {T, NVars, NTerms, Exponents}
        scale!(f.coefficients, λ)
        f
    end
end
