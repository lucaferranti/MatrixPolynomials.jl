function blkdiag(blocks::Array{Array{Pol{T}, 2}, 1}) where {T<:Number}
    m = sum(size.(blocks, 1))
    n = sum(size.(blocks, 2))
    J = zeros(Pol{T}, m, n)
    offset1 = 0
    offset2 = 0
    for blk in blocks
        J[offset1+1:offset1+size(blk, 1), offset2+1:offset2+size(blk, 2)] = blk
        offset1 += size(blk, 1)
        offset2 += size(blk, 2)
    end
    return J
end

"""
creates the matrix block of size n×(n+1) corresponding to column index n
"""
function colBlock(n::Int)
    var = Pol([0, 1])
    J = zeros(Pol{Int}, n, n+1)
    for ii in 1:n
        J[ii, ii] = var
        J[ii, ii+1] = -1
    end
    return J
end

colBlock(arr::AbstractArray{Int, 1}) = blkdiag(colBlock.(arr))

"""
creates the matrix block of size (n+1)×n corresponding to row index n
"""
function rowBlock(n::Int)
    var = Pol([0, 1])
    J = zeros(Pol{Int}, n+1, n)
    @inbounds for ii in 1:n
        J[ii, ii] = -1
        J[ii+1, ii] = var
    end
    return J
end

rowBlock(arr::AbstractArray{Int, 1}) = blkdiag(rowBlock.(arr))

"""
creates the block matrix corresponding to infinite eigenvalue of multiplicity n
"""
function infBlock(n::Int)
    var = Pol([0, 1])
    J = zeros(Pol{Int}, n, n)
    @inbounds for ii in 1:n-1
        J[ii, ii] = -1
        J[ii+1, ii] = var
    end
    J[n, n] = -1
    return J
end

infBlock(arr::AbstractArray{Int, 1}) = blkdiag(infBlock.(arr))

"""
creates the block matrix of size n×n corresponding to eigenvalue λ
"""
function eigBlock(λ::T, n::Int) where {T<:Number}
    var = Pol([0, 1])
    J = zeros(Pol{T}, n, n)
    @inbounds for ii in 1:n-1
        J[ii, ii] = var - λ
        J[ii+1, ii] = -1
    end
    J[n, n] = var - λ
    return J
end

eigBlock(λ::Number, arr::AbstractArray{Int, 1}) = blkdiag(eigBlock.(λ, arr))
