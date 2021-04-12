struct Pol{T<:Number}
    coeffs::Vector{T}
    function Pol{T}(coeffs) where {T <: Number}
        length(coeffs) == 1 || removeTailZeroes!(coeffs)
        new(coeffs)
    end
end

Pol(coeffs::Vector{T}) where {T<:Number} = Pol{T}(coeffs)
Pol() = Pol([0])
Pol(n::Number) = Pol([n])

rev(p::Pol) = Pol(p.coeffs[end:-1:1])

# Get properties
coeffs(p::Pol) = p.coeffs
coeff(p::Pol, d::Integer) = p[d+1]

lc(p::Pol) = p.coeffs[end]

function lm(p::Pol{T}) where {T<:Number}
    coeffs = zeros(T, length(p.coeffs))
    coeffs[end] = p.coeffs[end]
    return Pol(coeffs)
end

deg(p::Pol) = length(p.coeffs)-1

function variable(x::String, t=Int)
    global _varname = x
    return one(t)*Pol([0, 1])
end

Base.length(p::Pol) = length(p.coeffs)

function subs(p::Pol{T}, x::S) where {T,S<:Number}
    res = lc(p)
    @inbounds for i in length(p)-1:-1:1
        res = x*res+p[i]
    end
    return res
end

# PROMOTION AND CONVERSION
Base.promote_rule(::Type{Pol{T}}, ::Type{S}) where {T, S <: Number} = Pol{promote_type(T, S)}
Base.promote_rule(::Type{Pol{T}}, ::Type{Pol{S}}) where {T, S <: Number} = Pol{promote_type(T, S)}
Base.convert(::Type{Pol{T}}, x::S) where {T,S<:Number} = Pol{promote_type(T, S)}([x])
Base.convert(::Type{Pol{T}}, p::Pol{S}) where {T, S <: Number} = Pol{promote_type(T, S)}(p.coeffs)

Base.zero(p::Pol{T}) where {T<:Number} = Pol(zero(T))
Base.zero(::Type{Pol{T}}) where {T<:Number} = Pol(zero(T))
Base.zero(::Type{Pol}) = Pol(zero(Int))
Base.one(::Type{Pol{T}}) where {T <: Number} = Pol(one(T))

# indexing
function Base.getindex(p::Pol{T}, i::Integer) where {T<:Number}
    i < 1 && throw(BoundsError(p, i))
    i > length(p.coeffs) && return zero(T)
    return p.coeffs[i]
end

# DISPLAYING
function Base.show(io::IO, p::Pol)
    print(io, string(p))
end

function Base.string(p::Pol)
    if iszero(p)
        return "0"
    else
        str = ""
        @inbounds for (i, c) in Iterators.reverse(enumerate(p.coeffs))
            iszero(c) && continue
            str *= c>0 ? "+" : "-"
            str *= (abs(c) != 1 || i==1) ? string(abs(c)) : ""
            str *= i==1 ? "" : _varname;
            str *= i<=2 ? "" : toSuperscript(i-1)
        end
        str = str[1]=='+' ? str[2:end] : str
        return str
    end
end

function toSuperscript(n::Integer)
    exps = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹']
    dig = reverse(digits(n))
    return join([exps[d+1] for d in dig])
end
