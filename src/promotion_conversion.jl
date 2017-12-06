function promote_rule(::Type{S}, ::Type{SuperPolynomial{T, NV, NT, E}}) where {S<:Number, T, NV, NT, E}
    SuperPolynomial{promote_type(S, T), NV, NT, E}
end
function promote_rule(::Type{SuperPolynomial{S, NV, NT, E}}, ::Type{SuperPolynomial{T, NV, NT, E}}) where {S<:Number, T, NV, NT, E}
    SuperPolynomial{promote_type(S, T), NV, NT, E}
end

function convert(::Type{SuperPolynomial{S, NV, NT, E}}, f::SuperPolynomial{T, NV, NT, E}) where {S, T, NV, NT, E}
    SuperPolynomial{S, NV, NT, E}(convert.(S, f.coefficients))
end
