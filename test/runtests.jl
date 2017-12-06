using SuperPolynomials
using Base.Test
import FixedPolynomials
const FP = FixedPolynomials
import DynamicPolynomials: @polyvar
# write your own tests here

@testset "Evaluate, horner, gradient!" begin
    for T in [Float64, Complex128]
        for N=2:8
            M = 4*N

            exponents = round.(Int, 4 * rand(N, M))
            # some random coefficients
            coefficients = rand(T, M)
            g = SuperPolynomial(coefficients, exponents)
            w = rand(T, N)
            f = FP.Polynomial(exponents, coefficients)
            cfg = FP.GradientConfig(f)

            @test FP.evaluate(f, w, cfg) ≈ evaluate(g, w)
            @test FP.evaluate(f, w, cfg) ≈ evalhorner(g, w)

            u = zeros(T, N)
            v = zeros(T, N)
            V = zeros(T, 2, N)
            FP.gradient!(u, f, w, cfg)
            gradient!(v, g, w)
            gradient!(V, g, w, 2)
            @test all(u .≈ v)
            @test all(u .≈ V[2,:])
            @test all(0 .≈ V[1,:])
        end
    end
end

@testset "Constructors" begin
    @polyvar x y z

    f = SuperPolynomial(3.0 * x * y * z^3 - 2x^2*z)
    @test exponents(f) == [1 2; 1 0; 3 1]

    g = SuperPolynomial(3.0 * x * z^3 - 2x^2*z, [x, y, z])
    @test exponents(g) == [1 2; 0 0; 3 1]

    h = SuperPolynomial(FP.Polynomial{Complex128}(3.0 * x * y * z^3 - 2x^2*z))
    @test exponents(h) == [1 2; 1 0; 3 1]
end
