export gradient!


"""
    gradient!(u, f, x)

Compute the gradient of `f` at `x`, i.e. `∇f(x)`, and store the result in `u`.

    gradient!(U::Matrix, f, x, i)

Compute the gradient of `f` at `x`, i.e. `∇f(x)`, and store the result in the
`i`-th row of `U`.
"""
@generated function gradient!(u::AbstractMatrix, f::Polynomial{T, NVars, NTerms, Val{Exponents}}, x, row) where {T, NVars, NTerms, Exponents}
    gradient_impl(f, true)
end

@generated function gradient!(u::AbstractVector, f::Polynomial{T, NVars, NTerms, Val{Exponents}}, x) where {T, NVars, NTerms, Exponents}
    gradient_impl(f, false)
end


function gradient_impl(f::Type{Polynomial{T, NVars, NTerms, Val{Exponents}}}, matrix_assignement=false) where {T, NVars, NTerms, Exponents}
    M = convert_to_matrix(NVars, NTerms, Exponents)
    M_reduced = max.(M .- 1, 0)

    expressions = group_powers(M_reduced)

    for i=1:NVars
        push!(expressions, :($(u_(i)) = zero($T)))
    end
    coefficients_assignment = :(coefficients = f.coefficients)

    term_bases = gradient_term_bases(:coefficients, M_reduced)

    partial_derivs = partial_derivatives(M, T)

    assignments = Expr[]
    for i=1:size(M, 1)
        if matrix_assignement
            push!(assignments, :(u[row, $i] = $(u_(i))))
        else
            push!(assignments, :(u[$i] = $(u_(i))))
        end
    end
    Expr(:block,
        expressions...,
        coefficients_assignment,
        term_bases...,
        partial_derivs...,
        assignments...,
        :(u)
        )
end
