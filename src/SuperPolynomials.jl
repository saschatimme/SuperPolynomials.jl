__precompile__(true)

module SuperPolynomials
    import MultivariatePolynomials
    const MP = MultivariatePolynomials

    import FixedPolynomials
    const FP = FixedPolynomials

    using Nullables


    export Polynomial
    """
        Polynomial([T, ] f::MP.AbstractPolynomial, [variables])

    Construct a Polynomial from f.

        Polynomial(f::FixedPolynomials.Polynomial)
    """
    struct Polynomial{T, NVars, NTerms, Exponents<:Val}
        coefficients::Vector{T}
    end

    function Polynomial(coefficients::Vector{T}, exponents::Vector{<:Integer}, NVars) where T
        @assert length(coefficients)*NVars == length(exponents)
        NTerms = length(coefficients)
        Exponents = Val{ntuple(i -> exponents[i], length(exponents))}
        Polynomial{T, NVars, NTerms, Exponents}(coefficients)
    end

    function Polynomial(coefficients::Vector{T}, exponents::Matrix{<:Integer}) where T
        NVars = size(exponents, 1)
        exps, coeffs = normalize_exponents_coeffs(exponents, coefficients)
        Polynomial(coeffs, exps, NVars)
    end

    function Polynomial(p::MP.AbstractPolynomialLike, variables = MP.variables(p))
        exps, coeffs = mp_exponents_coefficients(p, variables)
        Polynomial(coeffs, exps, length(variables))
    end

    function Polynomial(::Type{T}, p::MP.AbstractPolynomialLike, variables = MP.variables(p)) where T
        exps, coeffs = mp_exponents_coefficients(p, variables)
        Polynomial(convert.(T, coeffs), exps, length(variables))
    end

    Polynomial(p::FP.Polynomial) = Polynomial(p.coefficients, p.exponents)

    include("promotion_conversion.jl")
    include("utilities.jl")
    include("evaluate.jl")
    include("horner.jl")
    include("gradient_helpers.jl")
    include("gradient.jl")
    include("system.jl")

    export exponents, coefficients, scale_coefficients!

    """
        exponents(f)

    Get the exponents of `f` as a matrix where each column represents a monomial of the
    polynomial.
    """
    function exponents(::Polynomial{T, NVars, NTerms, Val{Exponents}}) where {T, NVars, NTerms, Exponents}
        convert_to_matrix(NVars, NTerms, Exponents)
    end

    """
        coefficients(f)

    Get the coefficients of `f`.
    """
    coefficients(f::Polynomial) = f.coefficients

    """
        scale_coefficients!(f, λ)

    Scale the coefficients of `f` with the factor `λ`.
    """
    function scale_coefficients!(f::Polynomial, λ)
        scale!(f.coefficients, λ)
        f
    end
    function scale_coefficients!(F::PolynomialSystem, λ)
        for f in F.polynomials
            scale!(f.coefficients, λ)
        end
        F
    end
end
