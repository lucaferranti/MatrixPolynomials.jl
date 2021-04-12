function wilkinson(n::Integer)
    var = Pol([0, 1])
    res = var-1
    @inbounds for i in 2:n
        res *= var-i
    end
    return res
end
