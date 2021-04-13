"""
computes the eigenvalues of the polynomial matrix P
"""
function eig(P::AbstractArray{Pol{T}, N}, tol=0) where {T<:Number, N}

    e, d, Ar, Br = colNullities(P, tol)
    infMult = isempty(d) ? 0 : sum(collect(1:length(d)).*d)

    e1, d1, Ar, Br = rowNullities(Ar, Br, tol)
    e2, _, _, _ = rowNullities(P, tol)
    (e1 == e2 && isempty(d1)) || return @error "could not compute the eigenvalues. Try using a higher tolerance"
    eigs = eigvals(Ar, Br)
    append!(eigs, [Inf for _ in 1:infMult])

    return eigs
end

"""
computes the column indices of the polynomial matrix P
"""
function colIndices(P, tol=0)
    e, _, _, _ = colNullities(P, tol)
    _, _, Ar, Br = rowNullities(P, tol)
    e1, _, _, _ = colNullities(Ar, Br, tol)
    e1 == e || return @error "could not compute the column indices. Try using a higher tolerance"
    ϵ = null2indices(e)
    k = max(deg(P), 1)
    return [i-k+1 for i in ϵ]
end

"""
computes the row indices of the polynomial matrix P
"""
function rowIndices(P, tol=0)
    e, _, _, _ = rowNullities(P, tol)
    _, _, Ar, Br = colNullities(P, tol)
    e1, _, _, _ = rowNullities(Ar, Br, tol)
    e == e1 || return @error "could not compute the row indices. Try using a higher tolerance"

    return null2indices(e)
end

"""
computes how many times ∞ is an eigenvalue of the polynomial matrix P. That is, how many
times 0 is an eigenvalue of the reversal of P
"""
function infMultiplicities(P, tol=0)
    _, d, _, _ = colNullities(P, tol)
    _, d1, _, _ = rowNullities(P, tol)
    d == d1 || return @error "could not compute the infinite eigenvalues. Try using a higher tolerance"
    return inf2indices(d)
end

"""
computes the multiplicity of the eigenvalue λ of P.
"""
function eigMultiplicities(λ::Number, P, tol=0)
    e = eigNullities(λ, P, tol)
    return inf2indices(e)
end

"""
computes the multiplicity of each eigenvalue in arr.
"""
function eigMultiplicities(arr, P, tol=0)
    d = minimum(abs.(diff(arr)))
    d < 1e-6*maximum(abs.(arr)) && @warn "Some eigenvalues are close to each others (d=$d). The multiplicities might be incorrect"
    nulls = Dict{T, Array{Int, 1}}()
    @inbounds for λ in arr
        e = eigMultiplicities(λ, P, tol)
        push!(nulls, λ => e)
    end
    return nulls
end
