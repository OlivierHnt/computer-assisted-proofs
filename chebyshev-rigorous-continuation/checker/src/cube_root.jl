module CubeRoot

# Cube root example (paper, Section 2.2.3): F(x, λ) = x³ - λ - 2 on λ ∈ [-1, 1],
# lifted to ℓ¹ (weights ω ≡ 1) in the Chebyshev variable.
# From cube_root.jl of OlivierHnt/ChebyshevContinuation.jl, Step 2 (Newton
# continuation) removed: its output x_grid is an input.

using RadiiPolynomial

export cube_root_certificate





# Step 1: Defining the problem

F(x, λ) = x^3 - λ - exact(2)

DF(x) = exact(3) * x^2

const N = 2^3

λ_grid() = [-cospi(j/N) for j ∈ 0:N]





# Step 2 (the part that does not search): approximate inverse at the nodes and
# Chebyshev interpolants. `x_grid[j]` is the approximate zero at λ_grid()[j].

function cube_root_certificate(x_grid::Vector{Float64})
    length(x_grid) == N+1 || throw(ArgumentError("expected $(N+1) nodes"))

    A_grid = Vector{Float64}(undef, N+1)

    A_grid[1] = inv(DF(x_grid[1]))

    for j ∈ 2:N+1
        A_grid[j] = inv(DF(x_grid[j]))
    end

    # retrieve the Chebyshev interpolants

    λ_cheb = interval(Sequence(Chebyshev(1), [0.0, 0.5]))

    x_cheb = interval(real(to_coef(reverse(x_grid), Chebyshev(N))))

    A_cheb = interval(real(to_coef(reverse(A_grid), Chebyshev(N))))





    # Step 3: Estimating the bounds

    #- Y bound

    Y = norm(A_cheb * F(x_cheb, λ_cheb), 1)

    #- Z₁ bound

    Z₁ = norm(exact(1) - A_cheb * DF(x_cheb), 1)

    #- Z₂ bound

    r_star = 1e-6

    Z₂ = exact(6) * norm(A_cheb, 1) * (norm(x_cheb, 1) + exact(r_star))

    #-

    ie, proved = interval_of_existence(Y, Z₁, Z₂, r_star)

    return (; Y, Z₁, Z₂, r_star, ie, proved)
end

end
