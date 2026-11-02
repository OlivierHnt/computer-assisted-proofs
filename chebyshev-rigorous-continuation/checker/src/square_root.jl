module SquareRoot

# Square root example (paper, Section 3.1.2): F(x, λ) = x² - λ, pseudo-arclength
# continuation over λ ∈ [0, 1], lifted to (ℓ¹)² (weights ω ≡ 1).
# From square_root.jl of OlivierHnt/ChebyshevContinuation.jl, Step 2 (Newton
# continuation) removed: its outputs u_grid (points) and u̇_grid (tangents) are inputs.

using RadiiPolynomial, LinearAlgebra

export square_root_certificate





# Step 1: Defining the problem

F(x, λ) = x^2 - λ

DλF(x, λ) = exact(-1)

DxF(x, λ) = exact(2) * x

F_pa(u, u̇, w) = [F(u[1], u[2])
                 sum(u̇ .* (u - w))]

DF_pa(u, u̇) = [DxF(u[1], u[2]) DλF(u[1], u[2])
               u̇[1]            u̇[2]]

const N = 2^5

const arclength = sqrt(5) + 0.5 * log(2 + sqrt(5))





# Step 2 (the part that does not search): approximate inverse at the nodes and
# Chebyshev interpolants.

function square_root_certificate(u_grid::Vector{Vector{Float64}}, u̇_grid::Vector{Vector{Float64}})
    length(u_grid) == length(u̇_grid) == N+1 || throw(ArgumentError("expected $(N+1) nodes"))

    A_grid = Vector{Matrix{Float64}}(undef, N+1)

    A_grid[1] = inv(DF_pa(u_grid[1], u̇_grid[1]))

    for j ∈ 2:N+1
        A_grid[j] = inv(DF_pa(u_grid[j], u̇_grid[j]))
    end

    # retrieve the Chebyshev interpolants

    u_cheb = [interval(real(to_coef(getindex.(u_grid, i), Chebyshev(N)))) for i ∈ 1:2]

    u̇_cheb = [interval(real(to_coef(getindex.(u̇_grid, i), Chebyshev(N)))) for i ∈ 1:2]

    A_cheb = [interval(real(to_coef(getindex.(A_grid, i, j), Chebyshev(N)))) for i ∈ 1:2, j ∈ 1:2]





    # Step 3: Estimating the bounds

    #- Y bound

    Y = norm(norm.(A_cheb * F_pa(u_cheb, u̇_cheb, u_cheb), 1), 1)

    #- Z₁ bound

    Z₁ = opnorm(norm.(exact(I(2)) - A_cheb * DF_pa(u_cheb, u̇_cheb), 1), 1)

    #- Z₂ bound: Z₂ = 2(‖A₁₁‖ + ‖A₂₁‖), since D²F_pa only involves the first column of A

    r_star = Inf

    Z₂ = exact(2) * (norm(A_cheb[1, 1], 1) + norm(A_cheb[2, 1], 1))

    #-

    ie, proved = interval_of_existence(Y, Z₁, Z₂, r_star)

    return (; Y, Z₁, Z₂, r_star, ie, proved)
end

end
