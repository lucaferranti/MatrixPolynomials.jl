using MatrixPolynomials, LinearAlgebra
using Test

x = variable("x")

@testset "matrix operations" begin
    A = [x x-1;1 0]
    @test Matrix(I, size(A))*A == A
    @test A+A == 2*A
    @test A+A == [2x 2x-2;2 0]
    @test A*A == A^2
    @test A^2 == [x^2+x-1 x^2-x;x x-1]
    @test iszero(A-A)

    @test x*A == [x^2 x^2-x;x 0]
    @test A*x == x*A
    @test A+2 == [x+2 x+1;3 2]
    @test 2+A == A+2
    @test A-1 == [x-1 x-2;0 -1]
    @test 1-A == [1-x 2-x;0 1]

    B = [1 1;1 1]
    @test x+B == [x+1 x+1; x+1 x+1]
    @test B+x == x+B
    @test x^2*B == [x^2 x^2;x^2 x^2]
    @test B*x^2 == x^2*B

    @test deg(A) == 1
    @test deg(A^2) == 2
    @test lc(A) == [1 1;0 0]
    @test coeff(A, deg(A)) == lc(A)
    @test coeff(A, 0) == [0 -1;1 0]
    @test rev(A) == [1 -x+1;x 0]

    A1, B1 = toPencil(A)
    @test A1 == -coeff(A, 0)
    @test B1 == coeff(A, 1)

    A = [x^2 1 0;0 0 0;0 0 x]
    A1, B1 = toPencil(A)
    @test A1 == [0 0 0 0 -1 0; 0 0 0 0 0 0; 0 0 -1 0 0 0; 1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0]
    @test B1 == [1 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]
end
