function gradient_term_bases(coeff_name, reduced_exponents)
    term_bases = []
    ops = []
    # we compute the common factor to all partial derivatives
    for j = 1:size(reduced_exponents, 2)
        push!(ops, :($(coeff_name)[$j]))
        for i = 1:size(reduced_exponents, 1)
            k = reduced_exponents[i, j]
            if k > 0
                push!(ops, :($(x_((i, k)))))
            end
        end
        cj = c_(j)
        push!(term_bases, :(@inbounds $cj = $(batch_arithmetic_ops(:*, ops))))
        empty!(ops)
    end
    term_bases
end



function partial_derivatives(exponents, ::Type{<:Number})
    list = partial_derivatives_helper(exponents)

    exprs = Expr[]
    partial_derivatives_subfactors!(exprs, list)

    sort!(list, lt=((a, b) -> first(a) < first(b)))
    partials = Symbol[]
    last_i = 1
    for (i, j, k, factors) in list
        if i != last_i
            if !isempty(partials)
                push!(exprs, :($(u_(last_i)) = $(batch_arithmetic_ops(:+, partials))))
                empty!(partials)
            end
        end
        partial = u_(i, j)
        if k == 1 && isempty(factors)
            expr = :($(partial) = $(c_(j)))
        elseif k == 1
            expr = :($(partial) = $(c_(j)) * $(x_(factors)))
        elseif isempty(factors)
            expr = :($(partial) = $k * $(c_(j)))
        else
            expr = :($(partial) = $k * $(c_(j)) * $(x_(factors)))
        end
        last_i = i
        push!(partials, partial)
        push!(exprs, expr)
    end

    if !isempty(partials)
        push!(exprs, :($(u_(last_i)) = $(batch_arithmetic_ops(:+, partials))))
        empty!(partials)
    end
    exprs
end

function partial_derivatives(exponents, ::Type{<:Union{Float32, Float64}})
    list = partial_derivatives_helper(exponents)

    exprs = Expr[]
    partial_derivatives_subfactors!(exprs, list)

    # partials = [Vector{Symbol}() for i=1:size(exponents, 1)]

    for (i, j, k, factors) in sort(list, lt=((a, b) -> a[2] < b[2]))
        ui = u_(i)
        if k == 1 && isempty(factors)
            expr = :($ui += $(c_(j)))
        elseif k == 1
            expr = :($ui = muladd($(c_(j)), $(x_(factors)), $ui))
        elseif isempty(factors)
            expr = :($ui = muladd($k, $(c_(j)), $ui))
        else
            expr = :($ui = muladd($k * $(c_(j)), $(x_(factors)), $ui))
        end
        push!(exprs, expr)
    end
    exprs
end
function partial_derivatives_subfactors!(instructions, list)
    computed_values = Vector{IntSet}()
    values = sort!(unique(last.(list)), lt=((a, b) -> length(a) < length(b)))
    for (k, v) in enumerate(values)
        subsets, to_compute_subsets = find_partition(computed_values, values, k)
        append!(computed_values, to_compute_subsets)
        push!(computed_values, v)
        if length(subsets) == 1 && isempty(to_compute_subsets)
            continue
        end
        operands = x_.(subsets)
        for sub in to_compute_subsets
            prod = batch_arithmetic_ops(:*, x_.(collect(sub)))
            operand = x_(sub)
            push!(instructions, :($(operand) = $prod))
            push!(operands, operand)
        end
        if length(operands) > 1
            push!(instructions, :($(x_(v)) = $(batch_arithmetic_ops(:*, operands))))
        end
    end

    instructions
end


"""
    partial_derivatives_helper(exponents)

Compute a list of all occuring partial derivatives in the form
    (i, j, k, factors)

where `i` is the partial derivative, `j` is the corresponding term, `k` is the factor from the derivative
and `factors` is a vector of `Int`s with the occuring monomials.
"""
function partial_derivatives_helper(exponents)
    out = Vector{Tuple{Int, Int, Int, IntSet}}()
    m, n = size(exponents)
    for i=1:m, j=1:n
        k = exponents[i, j]
        if k > 0
            factors = IntSet()
            for l=1:m
                if l != i && exponents[l,j] > 0
                    push!(factors, l)
                end
            end
            # factors = [l for l=1:m if l != i && exponents[l,j] > 0]
            push!(out, (i, j, k, factors))
        end
    end
    out
end
