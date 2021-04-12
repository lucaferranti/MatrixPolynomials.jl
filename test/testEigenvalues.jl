using MatrixPolynomials, LinearAlgebra
using Test

x = variable("x")

@testset "column nullities" begin
    J1 = [x -1]
    e, d, Ar, Br = colNullities(J1)
    @test e == [0, 1]
    @test colIndices(J1) == [1]
    @test isempty(d)
    @test isempty(Ar)
    @test isempty(Br)

    J2 = [x -1 0;
          0  x -1]

    e, d, Ar, Br = colNullities(J2)
    @test e == [0, 0, 1]
    @test colIndices(J2) == [2]
    @test isempty(d)
    @test isempty(Ar)
    @test isempty(Br)

    J3 = [zeros(2, 2) J2]
    e, d, Ar, Br = colNullities(J3)
    @test e == [2, 0, 1]
    @test colIndices(J3) == [0, 0, 2]
    @test isempty(d)
    @test isempty(Ar)
    @test isempty(Br)

    J = [x]
    e, d, Ar, Br = colNullities(J)
    @test isempty(e)
    @test isempty(d)
    @test Ar == hcat([0])
    @test Br == hcat([1])

    J = [-1]
    e, d, Ar, Br = colNullities(J)
    @test isempty(e)
    @test d == [1]
    @test isempty(Ar)
    @test isempty(Br)

    J = [x -1 0 0 0;
        0 0 -1 0 0;
        0 0 0 x-1 0]

    e, d, Ar, Br = colNullities(J)
    @test e == [1, 1]
    @test colIndices(J) == [0, 1]
    @test d == [1]

    J = [x -1 0 0 0 0 0;
         0 0 x -1 0 0 0
         0 0 0 0 -1 0 0;
         0 0 0 0 0 x-1 0]

    e, d, Ar, Br = colNullities(J)
    @test e == [1, 2]
    @test colIndices(J) == [0, 1, 1]
    @test d == [1]

    J = [x -1 0 0 0 0  0 0  0 0 0;
         0 0 x -1 0 0  0 0  0 0 0;
         0 0 0 0 -1 0  0 0  0 0 0;
         0 0 0 0 0 x-1 0 0  0 0 0;
         0 0 0 0 0 0   0 0 -1 0 0;
         0 0 0 0 0 0   0 0  x 0 0;
         0 0 0 0 0 0   0 0  0 -1 0
         0 0 0 0 0 0   0 0  0  x -1]

    e, d, Ar, Br = colNullities(J)
    @test e == [2, 2]
    @test colIndices(J) == [0, 0, 1, 1]
    @test d == [1, 1]
end

@testset "row nullities" begin
    J1 = [-1;x]
    e, d, Ar, Br = rowNullities(J1)
    @test e == [0, 1]
    @test rowIndices(J1) == [1]
    @test isempty(d)
    @test isempty(Ar)
    @test isempty(Br)

    J2 = [-1 0;
        x -1;
        0 x]

    e, d, Ar, Br = rowNullities(J2)
    @test e == [0, 0, 1]
    @test rowIndices(J2) == [2]
    @test isempty(d)
    @test isempty(Ar)
    @test isempty(Br)

    J3 = [zeros(2, 2); J2]
    e, d, Ar, Br = rowNullities(J3)
    @test e == [2, 0, 1]
    @test rowIndices(J3) == [0, 0, 2]
    @test isempty(d)
    @test isempty(Ar)
    @test isempty(Br)

    J = [x]
    e, d, Ar, Br = rowNullities(J)
    @test isempty(e)
    @test isempty(d)
    @test Ar == hcat([0])
    @test Br == hcat([1])

    J = [-1]
    e, d, Ar, Br = rowNullities(J)
    @test isempty(e)
    @test d == [1]
    @test isempty(Ar)
    @test isempty(Br)

    J = [-1 0 0;
        x 0 0;
        0 -1 0;
        0 0 x-1;
        0 0 0]

    e, d, Ar, Br = rowNullities(J)
    @test e == [1, 1]
    @test rowIndices(J) == [0, 1]
    @test d == [1]

    J = [-1 0 0 0;
        x 0 0 0;
        0 -1 0 0;
        0 x 0 0;
        0 0 -1 0;
        0 0 0 0]

    e, d, Ar, Br = rowNullities(J)
    @test e == [1, 2]
    @test rowIndices(J) == [0, 1, 1]
    @test d == [1]

    J = [-1 0 0 0 0;
          x 0 0 0 0;
          0 -1 0 0 0;
          0 x 0 0 0;
          0 0 -1 0 0;
          0 0 0 -1 0;
          0 0 0 x -1;
          0 0 0 0 0]

    e, d, Ar, Br = rowNullities(J)
    @test e == [1, 2]
    @test rowIndices(J) == [0, 1, 1]
    @test d == [1, 1]
end

@testset "eigenvalues computations" begin
    J = [x-1 0;
         0   x-2]
    @test eig(J) ≈ [1, 2]
    @test isempty(infMultiplicities(J))

    A = [x^2 1 0; 0 0 0;0 0 x]
    @test eig(A) ≈ [0, Inf]
    @test infMultiplicities(A) == [1]

    A = [1 0 0;0 x 0;0 0 x^2-x]
    @test eig(A) ≈ [0, 0, 1, Inf, Inf, Inf]
    @test infMultiplicities(A) == [1, 2]

    J = [x -1 0 0 0 0  0 0  0 0 0;
    0 0 x -1 0 0  0 0  0 0 0;
    0 0 0 0 -1 0  0 0  0 0 0;
    0 0 0 0 0 x-1 0 0  0 0 0;
    0 0 0 0 0 0   0 0 -1 0 0;
    0 0 0 0 0 0   0 0  x 0 0;
    0 0 0 0 0 0   0 0  0 -1 0
    0 0 0 0 0 0   0 0  0  x -1]
    @test eig(J) ≈ [1, Inf, Inf, Inf]
    @test infMultiplicities(J) == [1, 2]

    A = [x-1.1 0 0
         0    -1 0
         0     x -1]
    @test eig(A) ≈ [1.1, Inf, Inf]
    @test infMultiplicities(A) == [2]

    A = zeros(Pol{Int}, (14, 16))
    A[1,3:4] = [x -1]
    A[2:3, 5:7] = [x -1 0;0 x -1]
    A[5:8, 8:10] = [-1 0 0;x -1 0;0 x -1;0 0 x]
    A[9, 11] = -1
    A[10:11, 12:13] = [-1 0;x -1]
    A[12, 14] = x-2
    A[13:14, 15:16] = [x-3 0;-1 x-3]

    @test eig(A) ≈ [2, 3, 3, Inf, Inf, Inf]
    @test colIndices(A) == [0, 0, 1, 2]
    @test rowIndices(A) == [0, 3]
    @test infMultiplicities(A) == [1, 2]

    F = svd(rand(size(A)...), full=true)
    A1 = F.U*A*F.Vt
    @test eig(A1, 1e-6) ≈ [2, 3, 3, Inf, Inf, Inf]
end
