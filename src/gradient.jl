export gradient!

function new_gradient_impl(::Type{T}, f::Type{Polynomial{S, NVars, NTerms, Val{Exponents}}}) where {T, S, NVars, NTerms, Exponents}
    exponents = convert_to_matrix(NVars, NTerms, Exponents)
    G = ComputeGraph()
    reduced_exponents = max.(exponents .- 1, 0)
    pow_instr = all_powers_of_variables!(G, reduced_exponents)

    ijk_prods, prods = gradient_products_to_evaluate!(G, exponents, reduced_exponents)
    res = group_products!(G, prods)

    init = Expr[]
    for i=1:size(exponents, 1)
        push!(init, :($(u_(i)) = zero($T)))
    end
    push!(init, :(c = f.coefficients))
    for j=1:size(exponents, 2)
        push!(init, :(@inbounds $(Symbol("c", j)) = c[$j]))
    end

    sort!(ijk_prods, lt=((a, b) -> (a[1], a[2]) < (b[1], b[2])))

    summations = Vector{Symbol}()
    computed_vertices = Set{Symbol}()
    last_i = 1
    graph_instr = Expr[]
    for i=1:NVars
        sub_ijk_prods = filter(x -> x[1] == i, ijk_prods)
        construct_instructions!(graph_instr, computed_vertices, G, sub_ijk_prods) do ijk_prod, instr
            i, j, k, prod = ijk_prod
            if isempty(prod)
                :(0)
            else
                :($(x_(prod)) = $instr)
            end
        end
        for ijk_prod in sub_ijk_prods
            i, j, k, prod = ijk_prod
            push!(graph_instr, :(@inbounds $(u_(i, j)) = $(k) * $(x_(prod))))
            push!(summations, u_(i, j))
        end
        if !isempty(summations)
            push!(graph_instr, :($(u_(i)) = $(batch_arithmetic_ops(:+, summations))))
            empty!(summations)
        end
    end

    assignments = [:(u[$i] = $(u_(i))) for i=1:size(exponents, 1)]

    Expr(:block,
        init...,
        pow_instr...,
        graph_instr...,
        assignments...,
        :(u))
end



function gradient_products_to_evaluate!(G, exponents, reduced_exponents)
    out = Vector{Tuple{Int, Int, Int, Vector{Symbol}}}()
    out_factors = Vector{Vector{Symbol}}()
    ops = Symbol[]
    # terms = Symbol[]
    m, n = size(reduced_exponents)
    for j = 1:n
        # push!(ops, Symbol("c", j))
        # construct the bases for the partial derivatives
        for i = 1:m
            k = reduced_exponents[i, j]
            if k > 0
                push!(ops, x_((i, k)))
            end
        end
        base_term = x_(ops)
        for op in ops
            addedge!(G, op, base_term)
        end
        cj = Symbol("c", j)
        push!(ops, cj)
        base_term_old = base_term
        base_term = x_(ops)
        addedge!(G, cj, base_term)
        addedge!(G, base_term_old, base_term)

        # now we have to handle all partial derivatives
        for i=1:m
            k = exponents[i, j]
            if k > 0
                factors = Symbol[]
                partial_derivatives = Symbol[]
                for l=1:m
                    kk = exponents[l,j]
                    if kk > 0
                        if l != i
                            push!(partial_derivatives, x_((l, kk)))
                            push!(factors, x_(l))
                        elseif kk > 1
                            push!(partial_derivatives, x_((l, kk - 1)))
                        end
                    end
                end
                push!(partial_derivatives, Symbol("c", j))
                push!(out, (i, j, k, partial_derivatives))

                # now we still have to construct the compute graph
                # the factors will be added separetly to the graph (we want to make optimizations)
                if !(isempty(partial_derivatives) || isempty(ops) || isempty(factors))
                    x_part = x_(partial_derivatives)
                    if !hasvertex(G, x_part)
                         x_facts = x_(factors)

                         addedge!(G, base_term, x_part)
                         # we have to introduce a fake node
                        if base_term == x_facts
                            mod_x_facts = Symbol("MOD$(i)_$(j)_$(x_facts)")
                            addedge!(G, base_term, mod_x_facts)
                            addedge!(G, mod_x_facts, x_part)
                        else
                            addedge!(G, x_facts, x_part)
                        end
                    end
                end
                push!(out_factors, factors)
            end
        end


        empty!(ops)
    end
    # terms
    out, out_factors
end

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

"""
    gradient!(u, f, x)

Compute the gradient of `f` at `x`, i.e. `∇f(x)`, and store the result in `u`.

    gradient!(U::Matrix, f, x, i)

Compute the gradient of `f` at `x`, i.e. `∇f(x)`, and store the result in the
`i`-th row of `U`.
"""
@generated function gradient!(u::AbstractMatrix{T}, f::Polynomial{S, NVars, NTerms, Val{Exponents}}, x, row) where {T, S, NVars, NTerms, Exponents}
    gradient_impl(T, f, true)
end

@generated function gradient!(u::AbstractVector{T}, f::Polynomial{S, NVars, NTerms, Val{Exponents}}, x) where {T, S, NVars, NTerms, Exponents}
    gradient_impl(T, f, false)
end

@generated function new_gradient!(u::AbstractVector{T}, f::Polynomial{S, NVars, NTerms, Val{Exponents}}, x) where {T, S, NVars, NTerms, Exponents}
    new_gradient_impl(T, f)
end
function gradient_impl(::Type{T}, f::Type{Polynomial{S, NVars, NTerms, Val{Exponents}}}, matrix_assignement=false) where {T, S, NVars, NTerms, Exponents}
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
