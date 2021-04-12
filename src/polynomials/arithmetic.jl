Base.:+(p::Pol) = p
Base.:-(p::Pol) = Pol([-c for c in p.coeffs])


function Base.:+(p1::Pol, p2::Pol)
    maxdeg = max(deg(p1), deg(p2))
    coeffs = [p1[i]+p2[i] for i in 1:maxdeg+1]
    return Pol(coeffs)
end

Base.:+(p1::Union{Pol, Number}, p2::Union{Pol, Number}) = +(promote(p1, p2)...)
Base.:-(p1::Union{Pol, Number}, p2::Union{Pol, Number}) = p1 + (-p2)

function Base.:*(p1::Pol{T}, p2::Pol{T}) where {T<:Number}
    coeffs = zeros(T, deg(p1)+deg(p2)+1)
    @inbounds for (i, c1) in enumerate(p1.coeffs)
        @inbounds for (j, c2) in enumerate(p2.coeffs)
            coeffs[i+j-1] += c1*c2
        end
    end
    return Pol(coeffs)
end

Base.:*(p1::Union{Pol, Number}, p2::Union{Pol, Number}) = *(promote(p1, p2)...)

Base.:/(p::Pol, n::Number)  = Pol([c/n for c in p.coeffs])

function Base.:^(p::Pol, n::Integer)
    iszero(n) && return Pol([1])
    if iseven(n)
        return (p*p)^(n÷2)
    else
        return p*(p*p)^((n-1)÷2)
    end
end

function Base.divrem(a::Pol, b::Pol)
    r = a
    q = zero(a)
    d = deg(b)
    c = lc(b)
    while !iszero(r) && deg(r) >= d
        s = lc(r)/c
        coeffs = zeros(typeof(s), deg(r)-d+1)
        coeffs[end] = s
        s = Pol(coeffs)
        q, r = q + s, r - s*b
    end
    return q, r
end

function Base.div(a::Pol, b::Pol)
    q, _ = divrem(a, b)
    return q
end

function Base.rem(a::Pol, b::Pol)
    _, r = divrem(a, b)
    return r
end

function Base.gcd(a::Pol, b::Pol)
    r0 = a
    r1 = b
    while !iszero(r1)
        r1, r0 = rem(r0, r1), r1
    end
    return r0/lc(r0)
end

# comparisons

function Base.:(==)(p1::Pol, p2::Pol)
    deg(p1) == deg(p2) || return false
    @inbounds for i in 1:length(p1.coeffs)
        p1.coeffs[i] == p2.coeffs[i] || return false
    end
    return true
end

for op in (:>, :<, :>=, :<=)
    @eval Base.$op(p1::Pol, p2::Pol) = $op(deg(p1), deg(p2))
end

Base.iszero(p::Pol) = iszero(p.coeffs)
