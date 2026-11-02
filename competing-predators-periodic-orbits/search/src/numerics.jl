# Numerical search for the approximate family: ODE integration of one periodic
# orbit, Newton refinement, and continuation in κ over the Chebyshev nodes
# (construct_tube). From src/numerics.jl of OlivierHnt/BHWPeriodicCoexistence.jl,
# commit 6effd7d.

# Fourier coefficients of f sampled on nextpow(2, 2K+1) equispaced points
_sample(f, s::Fourier) = (m = nextpow(2, dimension(s)); to_coef([complex(f(2π * j / m)) for j ∈ 0:m-1], s))

function vectorfield!(du, u, params, t)
    x₁, x₂, s = u
    δ₁, δ₂, α₁, α₂, λ₁, λ₂, τ = params

    du[1] = τ * δ₁ * x₁ * (s - λ₁) / (α₁ + s)
    du[2] = τ * δ₂ * x₂ * (s - λ₂) / (α₂ + s)
    du[3] = τ * (1 - s - x₁ / (α₁ + s) - x₂ / (α₂ + s)) * s

    return du
end

function vectorfield_lin!(du, u, params, t)
    J = params[1]
    J_at_t = real(evaluate.(J, t))
    mul!(view(du, 1:3), J_at_t, view(u, 1:3))
    mul!(view(du, 4:6), J_at_t, view(u, 4:6))
    mul!(view(du, 7:9), J_at_t, view(u, 7:9))
    return du
end

#

function refine_po!(x_po, params; F_po_cache = zeros(eltype(x_po), space(x_po)), DF_po_cache = zeros(eltype(x_po), space(x_po), space(x_po)))
    approx_tang = differentiate(component(x_po, 2))
    newton!(x_po, F_po_cache, DF_po_cache; verbose = false, tol = 1e-15) do F, DF, x
        F_po!(F, x, params, approx_tang)
        DF_po!(DF, x, params, approx_tang)
        F, DF
    end; conjugacy_symmetry!(x_po)

    approx_tang = differentiate(component(x_po, 2))
    _, success = newton!(x_po, F_po_cache, DF_po_cache; verbose = false, tol = 1e-15) do F, DF, x
        F_po!(F, x, params, approx_tang)
        DF_po!(DF, x, params, approx_tang)
        F, DF
    end; conjugacy_symmetry!(x_po)

    success || println("Periodic orbit failure")

    return x_po, success
end

function refine_lin!(x_lin, J; F_lin_cache = zeros(eltype(x_lin), space(x_lin)), DF_lin_cache = zeros(eltype(x_lin), space(x_lin), space(x_lin)))
    _, success = newton!(x_lin, F_lin_cache, DF_lin_cache; verbose = false) do F, DF, x
        F_lin!(F, x, J)
        DF_lin!(DF, x, J)
        F, DF
    end; conjugacy_symmetry!(x_lin)

    success || println("Floquet normal form failure")

    return x_lin, success
end

#

function get_po(u0, τ, params, K) # get periodic solution using numerical integration
    prob = ODEProblem(vectorfield!, u0, (0.0, 2π), [params ; τ])
    sol = solve(prob, Tsit5(); abstol = 1e-12, reltol = 1e-12)

    x_po = zeros(ComplexF64, ScalarSpace()^3 × Fourier(K, 1.0)^3)
    x_po[1] = τ
    x_po[2] = sol(0; idxs = 1)
    x_po[3] = sol(0; idxs = 2)
    for i ∈ 1:3
        coefficients(component(component(x_po, 2), i)) .= coefficients(_sample(t -> sol(t; idxs = i), Fourier(K, 1.0)))
    end
    coefficients(component(component(x_po, 2), 1)) ./= x_po[2]
    coefficients(component(component(x_po, 2), 2)) ./= x_po[3]
    conjugacy_symmetry!(x_po)

    refine_po!(x_po, params)

    return x_po
end

function get_lin(x_po, params, K_lin) # get linearized solution using numerical integration
    J = τ∂ᵤf(x_po, params)
    prob = ODEProblem(vectorfield_lin!, float(I[1:3,1:3]), (0.0, 2π), (J,))
    sol = solve(prob, Tsit5(); abstol = 1e-12, reltol = 1e-12)

    monodromy = reshape(sol(2π), 3, 3)
    C = log( monodromy ) ./ (2π)

    x_lin = zeros(ComplexF64, ScalarSpace()^9 × Fourier(K_lin, 1.0)^9)
    coefficients(component(x_lin, 1)) .= C[:]
    for i ∈ 1:3, j ∈ 1:3
        project!(component(component(x_lin, 2), i+3*(j-1)), _sample(t -> (reshape(sol(t), 3, 3) * exp(-C * t))[i,j], space(space(x_lin)[2])))
    end

    refine_lin!(x_lin, J)

    return x_lin
end

function construct_tube(u_init, τ_init, κ_init, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂); increment, K, K_lin, N, round)
    # validate the first periodic orbit
    params_init = [(m₁-d₁)/γ, (m₂-d₂)/γ, a₁ / κ_init, a₂ / κ_init, a₁ / κ_init * d₁/(m₁-d₁), a₂ / κ_init * d₂/(m₂-d₂)]
    x_po = get_po(u_init, τ_init, params_init, K)

    # reach boundary Q_2, {X₁ = 0}
    κ₂ = reach_boundary!(x_po, κ_init, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂), increment, 2)

    # reach boundary Q_1, {X₂ = 0}
    κ₁ = reach_boundary!(x_po, κ₂, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂), -increment, 3)

    if round == true
        κ₂ = ceil(κ₂)
        κ₁ = floor(κ₁)
    end

    # Floquet normal form
    params = [(m₁-d₁)/γ, (m₂-d₂)/γ, a₁ / κ₁, a₂ / κ₁, a₁ / κ₁ * d₁/(m₁-d₁), a₂ / κ₁ * d₂/(m₂-d₂)]
    x_lin = get_lin(x_po, params, K_lin)

    # whole family
    x_po_grid, x_lin_grid = continuation!(x_po, x_lin, κ₁, κ₂, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂), N)

    #

    return x_po_grid, x_lin_grid, κ₁, κ₂
end

function reach_boundary!(x_po, κ, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂), increment, i)
    println("Trying to reach boundary Q_$(i-1)...")
    κ_opt = κ
    success = true
    while real(x_po[i]) ≥ -0.01 # 0
        params = ((m₁-d₁)/γ, (m₂-d₂)/γ, a₁ / κ_opt, a₂ / κ_opt, a₁ / κ_opt * d₁/(m₁-d₁), a₂ / κ_opt * d₂/(m₂-d₂))
        _, success = refine_po!(x_po, params)
        if !success
            println("Failure to reach boundary Q_$(i-1)")
            break
        end
        κ_opt += increment
    end
    success && println("Successfully reached boundary Q_$(i-1)")
    return κ_opt
end

function continuation!(x_po, x_lin, κ₁, κ₂, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂), N)
    N_fft = nextpow(2, 2N + 1) # size of the discrete Chebyshev transform (64 for N = 30); the grid has N_fft ÷ 2 + 1 nodes
    npts = N_fft ÷ 2 + 1
    κ⁻¹_grid = [(κ₁ + κ₂ - cospi(2j/N_fft)*(κ₂ - κ₁))/(2κ₁*κ₂) for j ∈ npts-1:-1:0]

    x_po_grid = Vector{typeof(x_po)}(undef, npts)
    F_po_cache = zeros(eltype(x_po), space(x_po))
    DF_po_cache = zeros(eltype(x_po), space(x_po), space(x_po))

    x_lin_grid = Vector{typeof(x_lin)}(undef, npts)
    F_lin_cache = zeros(eltype(x_lin), space(x_lin))
    DF_lin_cache = zeros(eltype(x_lin), space(x_lin), space(x_lin))

    for j ∈ 1:npts
        params = ((m₁-d₁)/γ, (m₂-d₂)/γ, a₁ * κ⁻¹_grid[j], a₂ * κ⁻¹_grid[j], a₁ * κ⁻¹_grid[j] * d₁/(m₁-d₁), a₂ * κ⁻¹_grid[j] * d₂/(m₂-d₂))

        # periodic orbit
        refine_po!(x_po, params; F_po_cache = F_po_cache, DF_po_cache = DF_po_cache)
        x_po_grid[j] = copy(x_po)

        # Floquet normal form
        J = τ∂ᵤf(x_po, params)
        refine_lin!(x_lin, J; F_lin_cache = F_lin_cache, DF_lin_cache = DF_lin_cache)
        x_lin_grid[j] = copy(x_lin)
    end

    return x_po_grid, x_lin_grid
end
