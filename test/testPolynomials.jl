using MatrixPolynomials
using Test

x = variable("x")

@testset "Polynomials definition" begin

    p1 = x^2+x+1
    @test p1 == Pol([1, 1, 1])

    p2 = x^2+2x+3
    @test p2 == Pol([3, 2, 1])

    @test rev(p1) == p1
    @test rev(p2) == Pol([1, 2, 3])

    @test deg(p1) == 2
    @test deg(Pol([3])) == 0

    p2 = x+x+x^2+x^2
    @test p2 == 2x^2+2x
    @test lc(p2) == 2
    @test lm(p2) == 2x^2

    @test (x+1)^2 == x^2+2x+1
end

@testset "Polynomial opereations" begin
    p1 = x^2+x+1

    @test p1 + 0 == p1
    @test 0 + p1 == p1
    @test iszero(p1-p1)
    @test length(p1-p1) == 1
    @test deg(p1-p1) == 0

    p2 = Pol(1.0)
    @test p1+p2 == 1.0x^2+1.0x+2.0

    @test (x+1)^2 == x^2+2x+1
    @test (x+1)*(x+2) == x^2+3x+2
    @test p1*1 == p1
    @test 2*p1 == 2x^2+2x+2

    p1 = x^2+2x+1
    p2 = x+1
    q, r = divrem(p1, p2)
    @test iszero(r)
    @test q == x+1

end
