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

    function normalize_exponents_coeffs(exponents::Matrix, coefficients::AbstractVector{T}) where T
        E = [exponents[:,j] for j=1:size(exponents, 2)]
        P = sortperm(E, lt=revlexless, rev=true)
        f_exponents = Vector{Int}()
        f_coefficients = Vector{T}()
        k = 1
        while k ≤ length(P)
            v = E[P[k]]
            l = k
            while l < length(P) && E[P[l+1]] == v
                l += 1
            end
            append!(f_exponents, v)
            if l == k
                push!(f_coefficients, coefficients[P[k]])
            else
                c = sum(i -> coefficients[P[i]], k:l)
                push!(f_coefficients, c)

            end
            k = l + 1
        end
        f_exponents, f_coefficients
    end

    function revlexless(a, b)
        n = length(a)
        for j=n:-1:1
            if a[j] > b[j]
                return false
            elseif a[j] < b[j]
                return true
            end
        end
        return false
    end


    function exponents(::SuperPolynomial{T, NVars, NTerms, Val{Exponents}}) where {T, NVars, NTerms, Exponents}
        convert_to_matrix(NVars, NTerms, Exponents)
    end

    function convert_to_matrix(nvars, nterms, exponents)
        [exponents[nvars*(j - 1) + i] for i=1:nvars, j=1:nterms]
    end

    @inline pow(x::AbstractFloat, k::Integer) = Base.FastMath.pow_fast(x, k)
    # simplified from Base.power_by_squaring
    @inline function pow(x::Number, p::Integer)
        if p == 1
            return copy(x)
        elseif p == 0
            return one(x)
        elseif p == 2
            return x*x
        end
        t = trailing_zeros(p) + 1
        p >>= t
        while (t -= 1) > 0
            x *= x
        end
        y = x
        while p > 0
            t = trailing_zeros(p) + 1
            p >>= t
            while (t -= 1) >= 0
                x *= x
            end
            y *= x
        end
        return y
    end


    function group_powers(M::Matrix, NVars)
        _xs = []
        for i=1:NVars
            i_exps = sort!(unique(M[i,:]))
            if first(i_exps) == 0
                shift!(i_exps)
            end
            if !isempty(i_exps)
                last_k = first(i_exps)
                xi = Symbol("_x", i)
                push!(_xs, :($xi = x[$i]))
                last_xik = Symbol("_x", i, "_", last_k)
                if last_k > 1
                    push!(_xs, :($last_xik = pow($xi, $last_k)))
                else
                    push!(_xs, :($last_xik = $xi))
                end
                for k in Iterators.drop(i_exps, 1)
                    xik = Symbol("_x", i, "_",k)
                    if (k - last_k) > 1
                        push!(_xs, :($xik = $last_xik * pow($xi, $(k - last_k))))
                    else
                        push!(_xs, :($xik = $last_xik * $xi))
                    end
                    last_k = k
                    last_xik = xik
                end
            end
        end
        return _xs
    end

    function gradient_impl(u, f::Type{SuperPolynomial{T, NVars, NTerms, Val{Exponents}}}, x) where {T, NVars, NTerms, Exponents}
        M = convert_to_matrix(NVars, NTerms, Exponents)

        _xs = []
        for i=1:NVars
            i_unique_exps = unique(M[i,:])
            i_exps = sort!(unique(max.(i_unique_exps .- 1, 0)))
            if first(i_exps) == 0
                shift!(i_exps)
            end
            if !isempty(i_exps)
                last_k = first(i_exps)
                xi = Symbol("x", i)
                push!(_xs, :($xi = x[$i]))
                last_xik = Symbol("x", i, "_", last_k)
                if last_k > 1
                    push!(_xs, :($last_xik = pow($xi, $last_k)))
                else
                    push!(_xs, :($last_xik = $xi))
                end
                for k in Iterators.drop(i_exps, 1)
                    xik = Symbol("x", i, "_",k)
                    if (k - last_k) > 1
                        push!(_xs, :($xik = $last_xik * pow($xi, $(k - last_k))))
                    else
                        push!(_xs, :($xik = $last_xik * $xi))
                    end
                    last_k = k
                    last_xik = xik
                end
            end
        end

        as = []
        push!(as, :(u .= zero($T)))
        # end
        for j=1:NTerms
            c = :(f.coefficients[$j])

            E = M[:, j]
            # we compute the common factor to all partial derivatives
            for (i, k) in enumerate(max.(E .- 1, 0))
                if k > 0
                    xik = Symbol("x", i, "_", k)
                    c = :($c * $xik)
                end
            end
            push!(as, :(c = $c))
            for i=1:NVars
                if M[i, j] > 0
                    f_i = :(c * $(M[i, j]))

                    for l=1:NVars
                        if l == i || M[l, j] == 0
                            continue
                        end
                        xl = Symbol("x", l)
                        f_i = :($f_i * $xl)
                    end
                    push!(as, :(u[$i] += $f_i))
                end
            end
        end

        return Expr(:block, _xs..., as..., :u)
    end


    function eval_and_gradient_impl(u, f::Type{SuperPolynomial{T, NVars, NTerms, Val{Exponents}}}, x) where {T, NVars, NTerms, Exponents}
        M = convert_to_matrix(NVars, NTerms, Exponents)

        _xs = []
        for i=1:NVars
            i_unique_exps = unique(M[i,:])
            i_exps = sort!(unique(max.(i_unique_exps .- 1, 0)))
            if first(i_exps) == 0
                shift!(i_exps)
            end
            if !isempty(i_exps)
                last_k = first(i_exps)
                xi = Symbol("x", i)
                push!(_xs, :($xi = x[$i]))
                last_xik = Symbol("x", i, "_", last_k)
                if last_k > 1
                    push!(_xs, :($last_xik = pow($xi, $last_k)))
                else
                    push!(_xs, :($last_xik = $xi))
                end
                for k in Iterators.drop(i_exps, 1)
                    xik = Symbol("x", i, "_",k)
                    if (k - last_k) > 1
                        push!(_xs, :($xik = $last_xik * pow($xi, $(k - last_k))))
                    else
                        push!(_xs, :($xik = $last_xik * $xi))
                    end
                    last_k = k
                    last_xik = xik
                end
            end
        end

        as = []
        push!(as, :(out = zero($T)))
        push!(as, :(u .= zero($T)))
        # end
        for j=1:NTerms
            c = :(f.coefficients[$j])

            E = M[:, j]
            # we compute the common factor to all partial derivatives
            for (i, k) in enumerate(max.(E .- 1, 0))
                if k > 0
                    xik = Symbol("x", i, "_", k)
                    c = :($c * $xik)
                end
            end
            push!(as, :(c = $c))
            push!(as, :(out_term = c))
            for i=1:NVars

                if M[i, j] > 0
                    xi = Symbol("x", i)
                    push!(as, :(out_term *= $xi))
                    f_i = :(c * $(M[i, j]))

                    for l=1:NVars
                        if l == i || M[l, j] == 0
                            continue
                        end
                        xl = Symbol("x", l)
                        f_i = :($f_i * $xl)
                    end
                    push!(as, :(u[$i] += $f_i))
                end
            end
            push!(as, :(out += out_term))
        end

        return Expr(:block, _xs..., as..., :out)
    end


    @generated function gradient5!(u::AbstractVector, f::SuperPolynomial{T, NVars, NTerms, Val{Exponents}}, x) where {T, NVars, NTerms, Exponents}
        gradient_impl(u, f, x)
    end

    @generated function eval_and_gradient!(u::AbstractVector, f::SuperPolynomial{T, NVars, NTerms, Val{Exponents}}, x) where {T, NVars, NTerms, Exponents}
        eval_and_gradient_impl(u, f, x)
    end

    function gradient_alloc_impl(u, f::Type{SuperPolynomial{T, NVars, NTerms, Val{Exponents}}}, x) where {T, NVars, NTerms, Exponents}
        M = convert_to_matrix(NVars, NTerms, Exponents)

        _xs = []
        for i=1:NVars
            i_unique_exps = unique(M[i,:])
            i_exps = sort!(unique(max.(i_unique_exps .- 1, 0)))
            if first(i_exps) == 0
                shift!(i_exps)
            end
            xi = Symbol("x", i)
            push!(_xs, :($xi = x[$i]))
            if !isempty(i_exps)
                last_k = first(i_exps)
                last_xik = Symbol("x", i, "_", last_k)
                if last_k > 1
                    push!(_xs, :($last_xik = pow($xi, $last_k)))
                else
                    push!(_xs, :($last_xik = $xi))
                end
                for k in Iterators.drop(i_exps, 1)
                    xik = Symbol("x", i, "_",k)
                    if (k - last_k) > 1
                        push!(_xs, :($xik = $last_xik * pow($xi, $(k - last_k))))
                    else
                        push!(_xs, :($xik = $last_xik * $xi))
                    end
                    last_k = k
                    last_xik = xik
                end
            end
        end

        as = []
        push!(as, :(u .= zero($T)))
        # end
        for j=1:NTerms
            ops = []
            push!(ops, :(f.coefficients[$j]))


            E = M[:, j]
            # we compute the common factor to all partial derivatives
            for (i, k) in enumerate(max.(E .- 1, 0))
                if k > 0
                    push!(ops, :($(Symbol("x", i, "_", k))))
                end
            end
            if length(ops) > 1
                c = Expr(:call, :*, ops...)
            else
                c = :($(first(ops)))
            end
            push!(as, :(c = $c))


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
                            if length(factors) > 1
                                e = Expr(:call, :*, factors...)
                                push!(as, :($_xi = $e))
                            else
                                push!(as, :($_xi = $(factors[1])))
                            end

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
                            if length(factors) > 1
                                e = Expr(:call, :*, factors...)
                                push!(as, :($_yi = $e))
                            else
                                push!(as, :($_yi = $(factors[1])))
                            end

                        end
                        is_first = false
                    else
                        push!(as, :($_yi = $last_yi * $xi))
                    end
                    last_yi = _yi
                end
            end

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
                if length(ops) == 2
                    push!(as, :(u[$i] = muladd($(ops[1]), $(ops[2]), u[$i])))
                elseif length(ops) == 1
                    push!(as, :(u[$i] += $(ops[1])))
                else
                    a = Expr(:call, :*, ops[1:end-1]...)
                    b = ops[end]
                    push!(as, :(u[$i] = muladd($a, $b, u[$i]))) # if M[i, j] > 0
                end
            end
        end

        return Expr(:block, _xs..., as..., :u)
    end


    @generated function gradient_alloc!(u::AbstractVector, f::SuperPolynomial{T, NVars, NTerms, Val{Exponents}}, x) where {T, NVars, NTerms, Exponents}
        gradient_alloc_impl(u, f, x)
    end



    function evaluate_impl(f::Type{SuperPolynomial{T, NVars, NTerms, Val{Exponents}}}, x) where {T, NVars, NTerms, Exponents}
        M = convert_to_matrix(NVars, NTerms, Exponents)

        _xs = group_powers(M, NVars)

        as = []
        for j = 1:NTerms
            # aj = Symbol("a", j)
            a = :(f.coefficients[$j])
            nmultiplications = sum(M[:, j] .> 0)
            mult_counter = 0
            if nmultiplications == 0
                push!(as, :(sum += $a))
            else
                for i=1:NVars
                    k = M[i, j]
                    if k > 0
                        mult_counter += 1
                        xik = Symbol("_x", i, "_", k)
                        if mult_counter == nmultiplications
                            push!(as, :(@inbounds sum = muladd($a, $xik, sum)))
                        else
                            a = :($a * $xik)
                        end
                    end
                end
            end
        end

        Expr(:block,
            _xs...,
            :(sum = zero($T)),
            as...,
            :(sum))
    end

    @generated function evaluate(f::SuperPolynomial{T, NVars, NTerms, Val{Exponents}}, x) where {T, NVars, NTerms, Exponents}
        evaluate_impl(f, x)
    end


    function revlexless(a, b)
        n = length(a)
        for j=n:-1:1
            if a[j] > b[j]
                return false
            elseif a[j] < b[j]
                return true
            end
        end
        return false
    end


    function decrease_lex(a)
        b = copy(a)
        for i=1:length(a)
            if a[i] > 0
                b[i] -= 1
                return b
            end
        end
        b
    end


    function lower_set(E)
        sorted = sort(E, lt=revlexless, rev=true)
        zero_el = zero(last(sorted))
        if zero_el != last(sorted)
            push!(sorted, zero_el)
        end
        set = [shift!(sorted)]
        while !isempty(sorted)
            next = shift!(sorted)
            next_lower = decrease_lex(last(set))
            while next_lower != next && !revlexless(next_lower, next)
                push!(set, next_lower)
                next_lower = decrease_lex(next_lower)
            end
            push!(set, next)
        end
        set
    end

    function next_k(I, n, d_I)
        k = 1
        for j = d_I:-1:1
            if I[n][j] != I[n-1][j]
                return j
            end
        end
        k
    end

    # multivariat horner scheme
    @generated function horner(f::SuperPolynomial{T, NVars, NTerms, Val{Exponents}}, x) where {T, NVars, NTerms, Exponents}
        matrix = convert_to_matrix(NVars, NTerms, Exponents)
        E = [matrix[:,j] for j=1:NTerms]
        P = sortperm(E, lt=revlexless, rev=true)
        I = lower_set(E)

        d_I = NVars
        M = length(I)

        # we save all powers computed to precompute them later
        xpowers = [Vector{Int}() for _=1:d_I]

        _xs = []
        push!(_xs, :(r0 = f.coefficients[$(P[1])]))
        r0_zero = false
        last_k = 0
        l = 2
        r = [Symbol("r", k) for k=1:d_I]
        rs_zero = trues(d_I)

        for k=1:d_I
            push!(_xs, :($(r[k]) = zero($T)))
        end

        n = 2
        while n ≤ M
            # k = max{1≤j≤d_I : α(n)_j ≠ α(n-1)_j}
            k = next_k(I, n, d_I)
            n_i = n + 1

            # We want to find out if we can write instead of
            #
            # r1 = x1 * (r0 + r1)
            # r0 = zero(Float64)
            # r1 = x1 * (r0 + r1)
            #
            # r1 = x1_2 * (r0 + r1)
            if l ≤ NTerms && I[n] != E[P[l]]
                while n_i ≤ M
                    n_k = next_k(I, n_i, d_I)
                    # @show n_k
                    if n_k != k
                        break
                    end
                    if l ≤ NTerms && (I[n_i] == E[P[l]])
                        n_i += 1
                        break
                    end
                    n_i += 1
                end
            end

            non_zero_r = r[filter(i -> !rs_zero[i], 1:k)]
            if r0_zero && length(non_zero_r) > 1
                summand = Expr(:call, :+, non_zero_r...)
            elseif r0_zero && isempty(non_zero_r)
                summand = :(zero($T))
            elseif r0_zero
                summand = :($(non_zero_r[1]))
            elseif isempty(non_zero_r)
                summand = :r0
            else
                summand = Expr(:call, :+, :r0, non_zero_r...)
            end
            if n_i - n > 1
                push!(xpowers[k], n_i - n)
                xkj = Symbol("x", k, "_", (n_i - n))
                push!(_xs, :($(r[k]) = $(xkj)*($(summand))))
                n = n_i - 1
            else
                xk = Symbol("x", k)
                push!(_xs, :($(r[k]) = $xk*($(summand))))
            end
            rs_zero[k] = false
            if l ≤ NTerms && I[n] == E[P[l]]
                push!(_xs, :(r0 = f.coefficients[$(P[l])]))
                r0_zero = false
                l += 1
            else
                if n == M
                    push!(_xs, :(r0 = zero($T)))
                end
                r0_zero = true
            end
            for i = 1:k-1
                if n == M || !rs_zero[i]
                    push!(_xs, :($(r[i]) = zero($T)))
                end
                rs_zero[i] = true
            end
            n += 1
        end

        # We now still have to compute the powers
        _ps = []
        for (i, x_i_powers) in enumerate(xpowers)
            i_exps = sort!(unique(x_i_powers))
            if !isempty(i_exps)
                last_k = first(i_exps)
                xi = Symbol("x", i)
                push!(_ps, :($xi = x[$i]))
                last_xik = Symbol("x", i, "_", last_k)
                if last_k > 1
                    push!(_ps, :($last_xik = pow($xi, $last_k)))
                else
                    push!(_ps, :($last_xik = $xi))
                end
                for k in Iterators.drop(i_exps, 1)
                    xik = Symbol("x", i, "_",k)
                    if (k - last_k) > 1
                        push!(_ps, :($xik = $last_xik * pow($xi, $(k - last_k))))
                    else
                        push!(_ps, :($xik = $last_xik * $xi))
                    end
                    last_k = k
                    last_xik = xik
                end
            else
                xi = Symbol("x", i)
                push!(_ps, :($xi = x[$i]))
            end
        end


        Expr(:block,
            _ps...,
            _xs...,
            Expr(:call, :+, :r0, r...))
    end


end # module
