module MatrixPolynomials

using LinearAlgebra

export
    Pol, variable, rev,
    lc, lm, subs,
    zero, one, izero, polrank,
    coeffs, coeff, deg,
    companion, roots,
    toPencil, eig, colNullities, rowNullities, colIndices, rowIndices, infMultiplicities, eigNullities, eigMultiplicities,
    wilkinson,
    colBlock, blkdiag, rowBlock, infBlock, eigBlock


include("utils.jl")
include("polynomials/polynomials.jl")
include("polynomials/arithmetic.jl")
include("polynomials/specialPolynomials.jl")
include("polynomials/roots.jl")
include("matrices.jl")
include("eigenvalues/eigenvalues.jl")
include("eigenvalues/nullities.jl")
include("kronecker/blocks.jl")
end # module
