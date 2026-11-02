module CahnHilliard

# Two-parameter continuation of Cahn–Hilliard steady states (paper, Section 4.1):
# F(v, c, β) = Δv + β(v - v³ - c) on cosine series (Neumann on (0, 1)), over the
# U-shaped region Λ parametrised by Coons' map θ of [-1, 1]² (two choices of
# boundary curves: "no pinch" and "pinch"). The bounds are those of statement.md
# (Section 4.1, Step 3, and Appendix A).
# From cahn_hilliard.jl of OlivierHnt/ChebyshevContinuation.jl, Step 2 (Newton
# continuation) removed: its output v_grid is an input.

using RadiiPolynomial, LinearAlgebra

export ch_certificate, curves_no_pinch, curves_with_pinch, N₁, N₂, K, θ, F, DF, X_CH, X_Λ





# Step 1: Defining the problem

# Chebyshev orders in s = (s₁, s₂) and cosine order of ū

const N₁, N₂ = 20, 20

const K = 20

# ℓ¹ norms (paper, eq. (ell1_multi) and Section 4) with ω ≡ 1 and χ ≡ 1: the sums run
# over n ∈ Z² and k ∈ Z, so every nonzero Chebyshev or cosine index has weight 2
# (one `IdentityWeight()` per factor of the tensor space)

const X_Λ  = Ell1(IdentityWeight(), IdentityWeight())                   # ℓ^{1,⊗2}(R): Chebyshev ⊗ Chebyshev
const X_CH = Ell1(IdentityWeight(), IdentityWeight(), IdentityWeight()) # ℓ^{1,⊗2}(𝒳): (Chebyshev ⊗ Chebyshev ⊗ cosine)_sym

# Coons' parameterization (paper, eq. (coons) of Section 3.2)
#
# The boundary curves γ₁, …, γ₄ : [-1, 1] → ∂Λ and the map θ are written in interval
# arithmetic with the exact constants of Section 4.1 (c_max = 1/10, g(c) = 79/2 + 150c²).
# They accept a number s (point evaluation) or a Chebyshev sequence s (polynomial
# arithmetic): in the latter case θ(curves, s₁, s₂) is an interval Chebyshev polynomial
# enclosing the exact θ, which is a polynomial of degree ≤ 2 in each variable
# (statement.md, Example 3).

struct BoundaryCurves{T₁,T₂,T₃,T₄}
    γ₁ :: T₁
    γ₂ :: T₂
    γ₃ :: T₃
    γ₄ :: T₄
end

α(s) = (s + exact(1)) / exact(2) # [-1,1] → [0,1]

const c_max = interval(1) / interval(10)

g(c) = interval(79) / interval(2) + exact(150) * c^2 # lower boundary β = g(c) of Λ

function θ(curves::BoundaryCurves, s₁, s₂)
    P₁, P₂ = curves.γ₁(interval(-1)), curves.γ₁(interval(1))
    P₃, P₄ = curves.γ₃(interval(-1)), curves.γ₃(interval(1))
    return θ_edge(s₁, s₂, curves.γ₁, curves.γ₃) .+
           θ_edge(-s₂, s₁, curves.γ₄, curves.γ₂) .-
           θ_corners(s₁, s₂, P₁, P₂, P₃, P₄)
end

θ_edge(s₁, s₂, γᵢ, γⱼ) = γᵢ(s₁) .+ α(s₂) .* (γⱼ(-s₁) .- γᵢ(s₁))

θ_corners(s₁, s₂, P₁, P₂, P₃, P₄) =
    P₁ .+ α(s₁) .* (P₂ .- P₁) .+ α(s₂) .* (P₄ .- P₁ .+ α(s₁) .* (P₃ .- P₄ .- P₂ .+ P₁))

# the corner identities (paper, eq. (bdy_rel)) that make θ(∂[-1,1]²) = ∂Λ:
# γ₁(-1) = γ₄(1), γ₁(1) = γ₂(-1), γ₂(1) = γ₃(-1), γ₃(1) = γ₄(-1)
corner_mismatches(curves::BoundaryCurves) = reduce(vcat, [
    curves.γ₁(interval(-1)) .- curves.γ₄(interval(1)),
    curves.γ₁(interval(1))  .- curves.γ₂(interval(-1)),
    curves.γ₂(interval(1))  .- curves.γ₃(interval(-1)),
    curves.γ₃(interval(1))  .- curves.γ₄(interval(-1))])

#- boundary curves with no pinch (paper, eq. (parabolacurves1))

γ_no_pinch₁(s) = [(s - exact(3)) / exact(4) * c_max, g((s - exact(3)) / exact(4) * c_max)]
γ_no_pinch₂(s) = [s / exact(2) * c_max,              g(s / exact(2) * c_max)]
γ_no_pinch₃(s) = [(s + exact(3)) / exact(4) * c_max, g((s + exact(3)) / exact(4) * c_max)]
γ_no_pinch₄(s) = [-s * c_max, interval(41)]

const curves_no_pinch = BoundaryCurves(γ_no_pinch₁, γ_no_pinch₂, γ_no_pinch₃, γ_no_pinch₄)

#- boundary curves with pinch (paper, eq. (parabolacurves2))

γ_with_pinch₁(s) = [-c_max, interval(41)]
γ_with_pinch₂(s) = [s * c_max, g(s * c_max)]
γ_with_pinch₃(s) = [c_max, interval(41)]
γ_with_pinch₄(s) = [-s * c_max, interval(41)]

const curves_with_pinch = BoundaryCurves(γ_with_pinch₁, γ_with_pinch₂, γ_with_pinch₃, γ_with_pinch₄)

#

F(v, c, β) = Laplacian() * v + β * (v - v^3 - c)

DF(v, c, β) = Laplacian() + Multiplication(β * (exact(1) - exact(3) * v^2))

# utilities

struct PseudoInverseLaplacian <: AbstractDiagonalOperator end
RadiiPolynomial.getcoefficient(::PseudoInverseLaplacian, (codom, i)::Tuple{SymmetricSpace{<:Fourier},Integer}, (dom, j)::Tuple{SymmetricSpace{<:Fourier},Integer}) =
    (i == j) & !(i == j == 0) ? inv(-(frequency(dom) * exact(i))^2) : zero(frequency(dom))

# the Chebyshev nodes -cos(π j/N), j = 0, …, N, of the paper (eq. (multi_chebyshev_nodes))
nodes(N) = [-cospi(interval(j) / interval(N)) for j ∈ 0:N]

# number of Chebyshev nodes on which a polynomial of degree n is sampled: at least
# n + 1, rounded up to 2^k + 1 so that the FFT of RadiiPolynomial is radix-2
fast_nodes(n) = nextpow(2, n) + 1

# samples of a Chebyshev polynomial on the m₁ × m₂ grid of RadiiPolynomial (nodes
# ordered from s = 1 down to s = -1 in each variable), and the interpolant of such
# samples (exact when the degree is at most m - 1 in each variable)
sample(a, m₁, m₂) = to_grid(a, (m₁, m₂))
interpolate(x_grid, m₁, m₂) = real(to_coef(x_grid, Chebyshev(m₁ - 1) ⊗ Chebyshev(m₂ - 1)))





# Step 2 (the part that does not search): Coons' map, finite part of the approximate
# inverse at the nodes, Chebyshev interpolants.
# `v_grid[j₁+1, j₂+1]` is the approximate zero at θ(curves, -cos(π j₁/N₁), -cos(π j₂/N₂)).

function ch_certificate(curves, v_grid::Matrix{Sequence})
    size(v_grid) == (N₁+1, N₂+1) || throw(ArgumentError("expected $(N₁+1)×$(N₂+1) nodes"))

    # Coons' map at the nodes (interval enclosures of the exact values)
    θ_grid = [θ(curves, s₁, s₂) for s₁ ∈ nodes(N₁), s₂ ∈ nodes(N₂)]

    # finite part of the approximate inverse at the nodes, in binary64
    A_finite_grid = Matrix{LinearOperator}(undef, N₁+1, N₂+1)

    Π = Projection(space(v_grid[1,1]))

    for j₂ ∈ 1:N₂+1, j₁ ∈ 1:N₁+1
        cⱼ, βⱼ = mid.(θ_grid[j₁,j₂])

        A_finite_grid[j₁,j₂] = inv(Π * DF(v_grid[j₁,j₂], cⱼ, βⱼ) * Π)
    end

    # retrieve the Chebyshev interpolants (binary64 FFT, then thin intervals)

    v_cheb = interval(real(to_coef(reverse(v_grid; dims=(1,2)), Chebyshev(N₁) ⊗ Chebyshev(N₂))))

    A_finite_cheb = interval(real(to_coef(reverse(A_finite_grid; dims=(1,2)), Chebyshev(N₁) ⊗ Chebyshev(N₂))))

    # θ̄ = (c̄, β̄): Coons' map as an interval Chebyshev polynomial, by polynomial
    # arithmetic on the coordinate functions s₁, s₂ (no interpolation)

    s₁ = zeros(Interval{Float64}, Chebyshev(1) ⊗ Chebyshev(0)); s₁[(1,0)] = interval(0.5) # s₁ = T₁(s₁)
    s₂ = zeros(Interval{Float64}, Chebyshev(0) ⊗ Chebyshev(1)); s₂[(0,1)] = interval(0.5) # s₂ = T₁(s₂)

    c_cheb, β_cheb = θ(curves, s₁, s₂)

    N_θ = max(maximum(order(space(c_cheb))), maximum(order(space(β_cheb)))) # Chebyshev order of θ̄ (2)

    # consistency of θ̄ with the pointwise evaluation of θ (both enclose the exact θ,
    # so they must intersect at every node) and the corner identities of Coons' map
    θ_poly_grid = [reverse(sample(a, N₁ + 1, N₂ + 1); dims=(1,2)) for a ∈ (c_cheb, β_cheb)]
    θ_consistent = all(!isempty_interval(intersect_interval(θ_grid[j₁,j₂][i], real(θ_poly_grid[i][j₁,j₂])))
                       for j₁ ∈ 1:N₁+1, j₂ ∈ 1:N₂+1, i ∈ 1:2)
    corners = corner_mismatches(curves)
    θ_ok = θ_consistent && all(in_interval(0, d) for d ∈ corners)
    θ_mismatch = maximum(mag, corners)





    # Step 3: Estimating the bounds (paper, Section 4.1, Step 3; derivation in statement.md)

    Δ⁻¹ = PseudoInverseLaplacian()

    # the degree counts below hold in each variable separately (N = N₁ or N₂)

    #- Y bound: Y = ‖𝒜 𝓕̃(ū)‖, 𝓕̃(ū) = Δū + β̄(ū - ū³ - c̄)
    #  A has degree N, β̄ū³ degree 3N + N_θ: the product has degree 4N + N_θ

    m₁_Y, m₂_Y = fast_nodes(4N₁ + N_θ), fast_nodes(4N₂ + N_θ)

    c_grid_Y = sample(c_cheb, m₁_Y, m₂_Y);
    β_grid_Y = sample(β_cheb, m₁_Y, m₂_Y);

    v_grid_Y = sample(v_cheb, m₁_Y, m₂_Y);

    A_grid_Y = sample(A_finite_cheb, m₁_Y, m₂_Y) .+ Δ⁻¹ * (interval(I) - interval(Π)) # A(s) = A_{≤K}(s)Π_{≤K} + Δ₀⁻¹Π_{>K}

    Y_grid = A_grid_Y .* F.(v_grid_Y, c_grid_Y, β_grid_Y)
    Y = norm(interpolate(Y_grid, m₁_Y, m₂_Y), X_CH)

    #- Z₁ bound: Z₁ = max(‖Π_{≤q} - 𝒜 D𝓕̃(ū) Π_{≤q}‖, ‖β̄(1 - 3ū²)‖ / (π(K+1))²), q = 3K
    #  (the second term bounds the columns > q, on which
    #  𝒜 D𝓕̃(ū) Π_{>q} = Δ₀⁻¹Π_{>K} (Δ + β̄(1 - 3ū²)) Π_{>q}; see statement.md)
    #  A D𝓕̃(ū) has degree N + 2N + N_θ; β̄(1 - 3ū²) has degree 2N + N_θ

    m₁_Z, m₂_Z = fast_nodes(3N₁ + N_θ), fast_nodes(3N₂ + N_θ)

    c_grid_Z = sample(c_cheb, m₁_Z, m₂_Z);
    β_grid_Z = sample(β_cheb, m₁_Z, m₂_Z);

    v_grid_Z = sample(v_cheb, m₁_Z, m₂_Z);

    A_grid_Z = sample(A_finite_cheb, m₁_Z, m₂_Z) .+ Δ⁻¹ * (interval(I) - interval(Π))

    q = 3K
    Π_q = interval(Projection(evensym(Fourier(q, π))))

    Z₁_grid = Π_q .- A_grid_Z .* (DF.(v_grid_Z, c_grid_Z, β_grid_Z) .* Π_q)
    Z₁_finite = opnorm(interpolate(Z₁_grid, m₁_Z, m₂_Z), X_CH)

    tail_grid = β_grid_Z .* (exact(1) .- exact(3) .* v_grid_Z .^ 2)
    Z₁_tail = norm(interpolate(tail_grid, m₁_Z, m₂_Z), X_CH) / (interval(π) * interval(K + 1))^2

    Z₁ = max(Z₁_finite, Z₁_tail)

    #- Z₂ bound: Z₂ = 3‖β̄‖ ‖𝒜‖ (2‖ū‖ + r⋆), ‖𝒜‖ = max(‖𝒜_{≤K}‖, 1/(π(K+1))²)

    r_star = 10sup(Y)

    normA = max(opnorm(A_finite_cheb, X_CH), interval(1) / (interval(π) * interval(K + 1))^2)

    Z₂ = exact(3) * norm(β_cheb, X_Λ) * normA * (exact(2) * norm(v_cheb, X_CH) + exact(r_star))

    #- radii polynomial (Theorem 1.1: Y + Z₁ r + Z₂ r²/2 < r and Z₁ + Z₂ r < 1, the
    #  convention of `interval_of_existence`)

    ie, proved = interval_of_existence(Y, Z₁, Z₂, r_star)

    return (; Y, Z₁, Z₂, r_star, ie, proved, Z₁_finite, Z₁_tail, normA,
              θ_ok, θ_consistent, θ_mismatch, N_θ, m_Y = (m₁_Y, m₂_Y), m_Z = (m₁_Z, m₂_Z),
              c_cheb, β_cheb, v_cheb)
end

end
