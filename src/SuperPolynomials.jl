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

    function evaluate_impl(f::SuperPolynomial{T, NVars, NTerms, Val{Exponents}}, x) where {T, NVars, NTerms, Exponents}
        M = convert_to_matrix(NVars, NTerms, Exponents)

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
        M = convert_to_matrix(NVars, NTerms, Exponents)

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
                    xik = Symbol("_x", i, "_", k)
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
                            push!(as, :(sum = muladd($a, $xik, sum)))
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

end # module
