"""
    all_powers_of_variables!(G::ComputeGraph, M::Matrix)

Create the most efficient way to compute all occuring powers. Returns the instructions
and adds the occuring powers as nodes to `G`.
"""
function all_powers_of_variables!(G, M::Matrix)
    xs = Expr[]
    for i = 1:size(M, 1)
        exps = occuring_exponents(M, i)
        push!(xs, :($(x_(i)) = x[$i]))
        addvertex!(G, x_(i))

        if isempty(exps)
            continue
        end

        last_k = exps[1]
        last_xik = x_((i, last_k))
        p = static_pow(x_(i), last_k)
        push!(xs, :($last_xik = $p))
        addvertex!(G, last_xik)
        for j=2:length(exps)
            k = exps[j]
            xik = x_((i,k))
            if _pow_already_computed(exps, k - last_k)
                push!(xs, :($xik = $last_xik * $(x_((i, k - last_k)))))
            else
                p = static_pow(x_(i), k - last_k)
                push!(xs, :($xik = $last_xik * $p))
            end
            addvertex!(G, xik)
            last_k = k
            last_xik = xik
        end
    end
    return xs
end


"""
    group_products!(G, products)

Group `products` together such the necesarry amount of multiplications is low and
store the information the in the `ComputeGraph` `G`.
"""
function group_products!(G, products)
    computed_values = Vector{eltype(products)}()
    values = sort(products, lt=((a, b) -> length(a) < length(b)))
    for (k, v) in enumerate(values)
        if length(v) == 1
            continue
        end
        subsets, to_compute_subsets = find_partition(computed_values, values, k)
        append!(computed_values, to_compute_subsets)
        push!(computed_values, v)

        if length(subsets) == 1 && isempty(to_compute_subsets)
            continue
        end

        xv = x_(v)
        for w in to_compute_subsets
            xw = x_(w)
            for w_i in w
                addedge!(G, w_i, xw)
            end
            if xv != xw
                addedge!(G, xw, xv)
            end
        end
        for w in subsets
            xw1 = x_(w)
            if xv != xw1
                addedge!(G, xw1, xv)
            end
        end
    end
    values
end


function find_partition(computed_values, values, n)
    a = values[n]

    subsets, rest = find_partition_from_computed(computed_values, a)
    to_compute_subsets = Vector{eltype(subsets)}()
    if !isnull(rest)
        find_partition_from_future!(subsets, to_compute_subsets, values, get(rest), n)
    end

    subsets, to_compute_subsets
end

function find_partition_from_computed!(subsets, computed_values, a, n)
    if length(a) == 0
        return Nullable(a, false)
    end
    if length(a) == 1
        push!(subsets, a)
        return Nullable(a, false)
    end

    min_val = Nullable(a)
    min_subsets = Vector{eltype(subsets)}()
    for i=n:-1:1
        b = computed_values[i]
        if b ⊆ a
            c = setdiff(a, b)
            c_subsets = [b]
            val = find_partition_from_computed!(c_subsets, computed_values, c, i - 1)
            if isnull(val) || length(get(val)) < length(get(min_val))
                min_val = val
                min_subsets = c_subsets
            end
        end
        if isnull(min_val)
            break
        end
    end
    append!(subsets, min_subsets)
    # push!(subsets, (a, :to_compute))
    return min_val
end

function find_partition_from_computed(computed_values, a)
    subsets = Vector{eltype(computed_values)}()
    to_compute = find_partition_from_computed!(
        subsets, computed_values, a, length(computed_values))
    (subsets, to_compute)
end

function find_partition_from_future!(subsets, future_subsets, values, a, n)
    if length(a) == 0
        return nothing
    end
    if length(a) == 1
        push!(subsets, a)
        return nothing
    end

    for i=n:length(values)
        b = values[i]
        a_b = a ∩ b
        if 1 < length(a_b) < length(a)
            push!(future_subsets, a_b)
            c = setdiff(a, a_b)
            find_partition_from_future!(subsets, future_subsets, values, c, i + 1)
            return nothing
        end
    end
    push!(future_subsets, a)
    return nothing
end


"""
    construct_instructions(G::ComputeGraph, vertex_names)

Construct the necessary instructions to compute all vertices
"""
function construct_instructions(f, G::ComputeGraph, vertex_names)
    computed_vertices = Set{Symbol}()
    instructions = Expr[]
    for (i, vertex_name) in enumerate(vertex_names)
        last_instruction = construct_sub_instructions!(instructions, computed_vertices, G, vertex_name)
        push!(instructions, f(i, vertex_name, last_instruction))
    end
    instructions
end

function construct_sub_instructions!(instructions, computed_vertices, G::ComputeGraph, vertex_name::Symbol)
    nbs = in_neighbors(G, vertex_name)
    if isempty(nbs) || vertex_name ∈ computed_vertices
        return :($vertex_name)
    end

    operands = Symbol[]
    for nb in nbs
        if nb.name ∉ computed_vertices
            construct_sub_instructions!(instructions, computed_vertices, G, nb.name)
            push!(computed_vertices, nb.name)
        end
        push!(operands, nb.name)
    end

    prod = batch_arithmetic_ops(:*, operands)
    # we are in the recursion start
    if !isempty(out_neighbors(G, vertex_name))
        push!(instructions, :($vertex_name = $prod))
    end

    return prod
end
