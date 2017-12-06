using SuperPolynomials
using Base.Test
import FixedPolynomials
const FP = FixedPolynomials
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
