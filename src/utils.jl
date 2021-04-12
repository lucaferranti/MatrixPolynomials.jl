function removeTailZeroes!(arr)
    while length(arr)>1 && iszero(arr[end])
        pop!(arr)
    end
    nothing
end

function null2indices(arr)
    res = Int[]
    isempty(arr) && (return res)

    @inbounds for (i, x) in enumerate(arr)
        append!(res, [i-1 for _ in 1:x])
    end
    return res
end

function inf2indices(arr)
    res = Int[]
    isempty(arr) && return res

    @inbounds for (i, x) in enumerate(arr)
        append!(res, [i for _ in 1:x])
    end
    return res
end

function uniquecomplex(arr::Array{T, 1}, tol=0) where {T<:Number}

    iszero(tol) && (tol = maximum(eps.(abs.(arr))))
    res = T[]
    while true
        for (i, c) in enumerate(arr)
            dist = abs.(arr .- c)
            idx = findall(x -> x > tol, dist)
            push!(res, c)
            arr = arr[idx]
            break
        end
        isempty(arr) && break
    end
    return res
end
