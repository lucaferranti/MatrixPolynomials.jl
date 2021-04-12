function companion(p::Pol{T}) where {T<:Number}
    p /= lc(p)
    C = zeros(T, deg(p), deg(p))
    C[2:end, 1:end-1] .= Matrix(I, deg(p)-1, deg(p)-1)
    C[:, end] = -p.coeffs[1:end-1]
    return C
end

function roots(p::Pol{T}) where {T<:Number}
    C = companion(p)
    return eigvals(C)
end
