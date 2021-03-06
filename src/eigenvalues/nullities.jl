function colNullities(P::AbstractArray, tol=0)
    A, B = toPencil(P)
    return colNullities(A, B, tol)
end

"""
    e, d, Ar, Br = colNullities(A, B, tol=0)
computes the column nullities of the pencil Bx-A.
e = column nullities
d =  nullities of infinite eigenvalue
x*Br - Ar = reduced pencil containing only finite eigenvalues and row indices
"""
function colNullities(A::AbstractArray{T, N}, B::AbstractArray{T, N}, tol=0) where {T<:Number, N}

    r = rank(B, atol=tol)
    s = size(B, 2) - r
    blocks  = [s]
    ranks = Int[]

    while s != 0
        F = svd(B, full=true);
        B = B*F.V
        A = A*F.V

        F = svd(A[:, r+1:end], full=true)
        Ua = F.U
        ra = rank(A[:, r+1:end], atol=tol)
        push!(ranks, ra)
        P = Matrix(I, size(Ua))
        P = P[[ra+1:end; 1:ra], :]
        A = P*Ua'*A
        B = P*Ua'*B

        A = A[1:end-ra, 1:r]
        B = B[1:end-ra, 1:r]
        r = rank(B, atol=tol)
        s = size(B, 2) - r
        push!(blocks, s)
    end
    e = blocks[1:end-1] - ranks
    d = ranks - blocks[2:end]

    removeTailZeroes!(e)
    removeTailZeroes!(d)
    iszero(d) && (d = Int[])
    iszero(e) && (e = Int[])
    return e, d, A, B
end

function rowNullities(P::AbstractArray, tol=0)
    A, B = toPencil(P)
    return rowNullities(A, B, tol)
end

"""
    e, d, Ar, Br = colNullities(A, B, tol=0)
computes the row nullities of the pencil Bx-A.
e = row nullities
d =  nullities of infinite eigenvalue
x*Br - Ar = reduced pencil containing only finite eigenvalues and column indices
"""
function rowNullities(A::AbstractArray{T, N}, B::AbstractArray{T, N}, tol=0) where {T<:Number, N}

    r = rank(B, atol=tol)
    s = size(B, 1) - r
    blocks = [s]
    ranks = Int[]
    while s != 0
        F = svd(B, full=true)
        Ub = F.U
        P = Matrix(I, size(Ub))
        P = P[[r+1:end; 1:r], :]
        B = P*Ub'*B
        A = P*Ub'*A
        F = svd(A[1:s, :], full=true)
        ra = rank(A[1:s, :], atol=tol)
        push!(ranks, ra)
        A = A*F.V
        B = B*F.V
        A = A[s+1:end, ra+1:end]
        B = B[s+1:end, ra+1:end]
        r = rank(B, atol=tol)
        s = size(B, 1) - r
        push!(blocks, s)
    end

    e = blocks[1:end-1] - ranks
    d = ranks - blocks[2:end]

    removeTailZeroes!(e)
    removeTailZeroes!(d)
    iszero(d) && (d = Int[])
    iszero(e) && (e = Int[])
    return e, d, A, B
end

function eigNullities(??::Number, A::AbstractArray{T, N}, B::AbstractArray{T, N}, tol=0) where {T<:Number, N}
    A = A - ??*B
    r = rank(A, atol=tol)
    s = size(A, 2) - r

    blocks  = [s]

    while s != 0
        F = svd(A, full=true);
        B = B*F.V
        A = A*F.V

        F = svd(B[:, r+1:end], full=true)
        Ua = F.U

        P = Matrix(I, size(Ua))
        P = P[[s+1:end; 1:s], :]
        A = P*Ua'*A
        B = P*Ua'*B

        A = A[1:r, 1:r]
        B = B[1:r, 1:r]
        r = rank(A, atol=tol)
        s = size(A, 2) - r
        push!(blocks, s)
    end
    e = blocks[1:end-1] - blocks[2:end]

    return e
end

function eigNullities(??::Number, P::AbstractArray, tol=0)
    _, _, Ar, Br = colNullities(P, tol)
    _, _, Ar, Br = rowNullities(Ar, Br, tol)
    return eigNullities(??, Ar, Br, tol)
end
