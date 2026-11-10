using Test
using Capsule: Capsule
using LinearAlgebra: issuccess
using RadiiPolynomial: inf

# each script is a proof written as one reads it, with its own names: it runs in a module of its own
script(file) = @eval Main module $(gensym(file)) include($(joinpath(@__DIR__, file))) end
proved(cert, bound) = issuccess(cert) && all(Capsule.radius(cert) .≤ bound)

@testset "Capsule" begin

    @testset "finite-dimensional problem, operators with an error" begin
        m = script("play.jl")
        @test proved(m.cert, 1e-16)
    end

    @testset "periodic orbit of Lorenz: the explicit proof and `validate` agree" begin
        m = script("lorenz_po_validate.jl")
        @test m.proved
        @test proved(m.cert, [3e-11, 4e-10, 3e-10, 2e-15])
        @test Capsule.radius(m.cert) ≈ inf.(m.r) rtol = 1e-8
    end

    @testset "initial value problems and a steady state" begin
        m = script("logistic_ivp.jl")
        @test proved(m.cert, 3e-6)       # logistic, Taylor
        @test proved(m.cert_c, 2e-13)    # logistic, Chebyshev
        @test proved(m.cert_lv, 6e-6)    # Lotka-Volterra
        @test proved(m.cert_ch, 2e-14)   # Cahn-Hilliard
    end

    @testset "lifts of finite-dimensional problems and of an initial value problem" begin
        m = script("lift_validate.jl")
        @test proved(m.cert₁, 1e-11)     # x² - 2 + p, Chebyshev in t
        @test proved(m.cert₂, 1e-9)      # exp(x) - 2 - p
        @test proved(m.cert₃, 1e-12)     # x² - 2 + p, Taylor in t
        @test proved(m.cert₄, 1e-11)     # logistic, λ ∈ [0.9, 1.1]
    end

    @testset "the prototype of the lift" begin
        m = script("lift_prototype.jl")
        @test m.cert.success && m.cert_t.success && m.cert₂.success
    end

    @testset "lift of a steady state: a branch of Cahn-Hilliard" begin
        m = script("cahn_hilliard_continuation.jl")
        @test proved(m.cert, 1e-5)
    end

    @testset "quasilinear problem: a steady state of SKT and a branch through it" begin
        m = script("skt_steady_state.jl")
        @test proved(m.cert, 1e-12)
        @test proved(m.cert_t, 1e-9)
    end

    @testset "lift of a problem with a scalar unknown: a branch of periodic orbits of Lorenz" begin
        m = script("lorenz_po_continuation.jl")
        @test proved(m.cert_t, [3e-10, 8e-10, 7e-10, 2e-12])
    end

    @testset "pseudo-arclength: an arc of SKT, at a node and along the arclength" begin
        m = script("skt_pseudo_arclength.jl")
        @test proved(m.cert, 1e-13)
        @test proved(m.cert_t, 1e-8)
    end

end
