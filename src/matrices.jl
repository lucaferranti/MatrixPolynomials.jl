

deg(A::AbstractArray{Pol{T}, N}) where {T<:Number, N} = maximum(deg.(A))
deg(A::AbstractArray{T, N}) where {T<:Number, N} = 0


lc(P::AbstractArray{Pol{T}, N}) where {T<:Number, N} = coeff(P, deg(P))
lc(P::AbstractArray{T, N}) where {T<:Number, N} = P

coeff(A::AbstractArray{Pol{T}, N}, k::Int) where {T<:Number, N} = coeff.(A, k)

function rev(P::AbstractArray{Pol{T}, N}) where {T<:Number, N}
    res = zero(P)
    for i in 0:deg(P)
        res += Pol([0, 1])^(deg(P)-i)*coeff(P, i)
    end
    return res
end

subs(P::AbstractArray{Pol{T}, N}, val::S) where {T,S<:Number, N} = subs.(P, val)

Base.:+(A::AbstractArray{Pol{T}, N}, b::Number) where {T<:Number, N} = A .+ b
Base.:+(b::Number, A::AbstractArray{Pol{T}, N}) where {T<:Number, N} = b .+ A

Base.:-(A::AbstractArray{Pol{T}, N}, b::Number) where {T<:Number, N} = A .- b
Base.:-(b::Number, A::AbstractArray{Pol{T}, N}) where {T<:Number, N} = b .- A

for op in (:+, :-, :*)
    @eval Base.$op(A::AbstractArray, b::Pol) = [$op(a, b) for a in A]
    @eval Base.$op(b::Pol, A::AbstractArray) = [$op(b, a) for a in A]
end

function toPencil(P::AbstractArray{Pol{T}, N}) where {T<:Number, N}
    k = deg(P)

    m = size(P, 1)
    n = size(P, 2)

    # patological case of polynomial of degree zero
    if iszero(k)
        B = zeros(T, m, n)
        A = hcat(-coeff(P, 0)) # ensure array is 2-dimensional
        return A, B
    end

    A = zeros(T, (k-1)*n+m, k*n)
    B = zeros(T, (k-1)*n+m, k*n)
    B[1:m, 1:n] = lc(P)
    B[m+1:end, n+1:end] = Matrix(I, (k-1)*n, (k-1)*n)
    A[m+1:end, 1:(k-1)*n] = Matrix(I, (k-1)*n, (k-1)*n)
    @inbounds for i in 1:k
        A[1:m, (i-1)*n+1:i*n] = -coeff(P, k-i)
    end
    return A, B
end

function toPencil(P::AbstractArray{T, N}) where {T<:Number, N}
    m = size(P, 1)
    n = size(P, 2)
    B = zeros(T, m, n)
    A = hcat(-P)
    return A, B
end

function polrank(M::AbstractArray{Pol{T}, N}, tol=0.0) where {T<:Number, N}
    rnk  = 0
    imax = deg(M)*min(size(M)...)
    for i=1:imax
        rnk = max(rnk, rank(subs(M, i), atol=tol))
    end
    return rnk
end
