# function eval_and_gradient_impl(u, f::Type{SuperPolynomial{T, NVars, NTerms, Val{Exponents}}}, x) where {T, NVars, NTerms, Exponents}
#     M = convert_to_matrix(NVars, NTerms, Exponents)
#
#     _xs = []
#     for i=1:NVars
#         i_unique_exps = unique(M[i,:])
#         i_exps = sort!(unique(max.(i_unique_exps .- 1, 0)))
#         if first(i_exps) == 0
#             shift!(i_exps)
#         end
#         if !isempty(i_exps)
#             last_k = first(i_exps)
#             xi = Symbol("x", i)
#             push!(_xs, :($xi = x[$i]))
#             last_xik = Symbol("x", i, "_", last_k)
#             if last_k > 1
#                 push!(_xs, :($last_xik = pow($xi, $last_k)))
#             else
#                 push!(_xs, :($last_xik = $xi))
#             end
#             for k in Iterators.drop(i_exps, 1)
#                 xik = Symbol("x", i, "_",k)
#                 if (k - last_k) > 1
#                     push!(_xs, :($xik = $last_xik * pow($xi, $(k - last_k))))
#                 else
#                     push!(_xs, :($xik = $last_xik * $xi))
#                 end
#                 last_k = k
#                 last_xik = xik
#             end
#         end
#     end
#
#     as = []
#     push!(as, :(out = zero($T)))
#     push!(as, :(u .= zero($T)))
#     # end
#     for j=1:NTerms
#         c = :(f.coefficients[$j])
#
#         E = M[:, j]
#         # we compute the common factor to all partial derivatives
#         for (i, k) in enumerate(max.(E .- 1, 0))
#             if k > 0
#                 xik = Symbol("x", i, "_", k)
#                 c = :($c * $xik)
#             end
#         end
#         push!(as, :(c = $c))
#         push!(as, :(out_term = c))
#         for i=1:NVars
#
#             if M[i, j] > 0
#                 xi = Symbol("x", i)
#                 push!(as, :(out_term *= $xi))
#                 f_i = :(c * $(M[i, j]))
#
#                 for l=1:NVars
#                     if l == i || M[l, j] == 0
#                         continue
#                     end
#                     xl = Symbol("x", l)
#                     f_i = :($f_i * $xl)
#                 end
#                 push!(as, :(u[$i] += $f_i))
#             end
#         end
#         push!(as, :(out += out_term))
#     end
#
#     return Expr(:block, _xs..., as..., :out)
# end
