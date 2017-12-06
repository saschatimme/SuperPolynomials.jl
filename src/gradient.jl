export gradient!

function grouped_derivative_factors!(as, M, NVars, j)
    is_first = true
    last_xi = :default
    for i=1:NVars-1
        if M[i, j] > 0
            _xi = Symbol("_x", i)
            xi = Symbol("x", i)
            if is_first
                if i == 1
                    push!(as, :($_xi = $xi))
                else
                    factors = Symbol[]
                    for l=1:i
                         if M[l, j] > 0
                             push!(factors, Symbol("x", l))
                         end
                    end
                    e = batch_arithmetic_ops(:*, factors)
                    push!(as, :($_xi = $e))
                end
                is_first = false
            else
                push!(as, :($_xi = $last_xi * $xi))
            end
            last_xi = _xi
        end
    end

    is_first = true
    last_yi = :default
    for i=NVars:-1:2
        if M[i, j] > 0
            _yi = Symbol("_y", i)
            xi = Symbol("x", i)
            if is_first
                if i == NVars
                    push!(as, :($_yi = $xi))
                else
                    factors = Symbol[]
                    for l=NVars:-1:i
                         if M[l, j] > 0
                             push!(factors, Symbol("x", l))
                         end
                    end
                    e = batch_arithmetic_ops(:*, factors)
                    push!(as, :($_yi = $e))
                end
                is_first = false
            else
                push!(as, :($_yi = $last_yi * $xi))
            end
            last_yi = _yi
        end
    end
    return as
end

function gradient_impl(f::Type{SuperPolynomial{T, NVars, NTerms, Val{Exponents}}}, matrix_assignement) where {T, NVars, NTerms, Exponents}
    M = convert_to_matrix(NVars, NTerms, Exponents)
    M_reduced = max.(M .- 1, 0)

    grouped_powers = group_powers(M_reduced, NVars)

    as = []
    for i=1:NVars
        ui = Symbol("u", i)
        push!(as, :($ui = zero($T)))
    end
    # push!(as, :(u .= zero($T)))

    for j=1:NTerms
        ops = []
        push!(ops, :(f.coefficients[$j]))
        # we compute the common factor to all partial derivatives
        for i = 1:NVars
            k = M_reduced[i, j]
            if k > 0
                push!(ops, :($(Symbol("x", i, "_", k))))
            end
        end
        push!(as, :(c = $(batch_arithmetic_ops(:*, ops))))

        grouped_derivative_factors!(as, M, NVars, j)


        nonzero_term_exps = Vector{NTuple{2, Int}}()
        for i=1:NVars
            k = M[i, j]
            if k > 0
                push!(nonzero_term_exps, (i, k))
            end
        end
        for l = 1:length(nonzero_term_exps)
            (i, k) = nonzero_term_exps[l]
            ops = Any[:c]
            if k > 1
                push!(ops, :($k))
            end
            if l > 1
                _xi = Symbol("_x", nonzero_term_exps[l-1][1])
                push!(ops, _xi)
            end
            if l < length(nonzero_term_exps)
                _yi = Symbol("_y", nonzero_term_exps[l+1][1])
                push!(ops, _yi)
            end
            ui = Symbol("u", i)
            if length(ops) == 2
                push!(as, :($ui = muladd($(ops[1]), $(ops[2]), $ui)))
            elseif length(ops) == 1
                push!(as, :($ui += $(ops[1])))
            else
                a = batch_arithmetic_ops(:*, ops[1:end-1])
                b = ops[end]
                push!(as, :($ui = muladd($a, $b, $ui)))
            end
        end
    end

    vecassignments = []
    rowassignments = []
    for i=1:NVars
        ui = Symbol("u", i)
        push!(vecassignments, :(u[$i] = $ui))
        push!(rowassignments, :(u[row, $i] = $ui))
    end


    if matrix_assignement
        assignmentblock = Expr(:block, rowassignments...)
    else
        assignmentblock = Expr(:block, vecassignments...)
    end

    return Expr(:block, grouped_powers..., as..., assignmentblock, :u)
end


"""
    gradient!(u, f, x)

Compute the gradient of `f` at `x`, i.e. `∇f(x)`, and store the result in `u`.

    gradient!(U::Matrix, f, x, i)

Compute the gradient of `f` at `x`, i.e. `∇f(x)`, and store the result in the
`i`-th row of ``.
"""
@generated function gradient!(u::AbstractMatrix, f::SuperPolynomial{T, NVars, NTerms, Val{Exponents}}, x, row) where {T, NVars, NTerms, Exponents}
    gradient_impl(f, true)
end

@generated function gradient!(u::AbstractVector, f::SuperPolynomial{T, NVars, NTerms, Val{Exponents}}, x) where {T, NVars, NTerms, Exponents}
    gradient_impl(f, false)
end
