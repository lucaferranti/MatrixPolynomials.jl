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
