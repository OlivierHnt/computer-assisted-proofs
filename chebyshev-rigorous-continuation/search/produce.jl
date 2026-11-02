# SEARCH — not part of the proof.
#
# Runs the binary64 Newton continuations of the four examples (Step 2 of the
# scripts of OlivierHnt/ChebyshevContinuation.jl) and writes their outputs, the
# approximate solutions on the Chebyshev nodes (and, for the pseudo-arclength
# examples, the approximate tangents), to data/. The approximate inverses are
# not written: the checker recomputes them. Any other method of producing a
# candidate would do; the checker does not care how the data were obtained.
#
# Usage (from the package root; this is the `produce` command of entry.json):
#   julia --project=search search/produce.jl [outdir] [examples...]
# (examples: any subset of cube_root, square_root, cahn_hilliard, skt; default all)
# then record the hashes:
#   julia --project=checker tools/hash_data.jl
#
# The SKT continuation starts from the initial guesses in
# search/skt_initial_guesses.bin (skt_data.jld2 of that repository: 14 vectors
# of 163 binary64, in the order of SKT.ARC_NAMES). The Cahn–Hilliard initial
# guess contains rand(); the generator is seeded so that the run is reproducible.

import Pkg
Pkg.instantiate()

using RadiiPolynomial, LinearAlgebra, Random
include(joinpath(@__DIR__, "..", "checker", "src", "ChebyshevContinuationCertificate.jl"))
using .ChebyshevContinuationCertificate

BLAS.set_num_threads(min(2, Sys.CPU_THREADS))

outdir = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "..", "data")
mkpath(outdir)
examples = length(ARGS) ≥ 2 ? ARGS[2:end] : ["cube_root", "square_root", "cahn_hilliard", "skt"]

# ------------------------------------------------------------------------
# cube root

if "cube_root" ∈ examples
let F = CubeRoot.F, DF = CubeRoot.DF, N = CubeRoot.N

λ_grid = [-cospi(j/N) for j ∈ 0:N]

x_grid = Vector{Float64}(undef, N+1)

# initialize

x_bar = cbrt(λ_grid[1] + 2)

x_grid[1] = x_bar

# run continuation scheme

for j ∈ 2:N+1
    λ = λ_grid[j]

    predictor = x_grid[j-1]

    x_barⱼ, convergedⱼ = newton(predictor) do x
        return F(x, λ), DF(x)
    end
    convergedⱼ || error("Newton failed")

    x_grid[j] = x_barⱼ
end

write_cube_root(joinpath(outdir, "cube_root.bin"), x_grid)
println("cube root: written")
end
end

# ------------------------------------------------------------------------
# square root

if "square_root" ∈ examples
let F = SquareRoot.F, DxF = SquareRoot.DxF, DλF = SquareRoot.DλF, F_pa = SquareRoot.F_pa, DF_pa = SquareRoot.DF_pa,
    N = SquareRoot.N, arclength = SquareRoot.arclength

arclength_grid = [0.5 * (1 - cospi(j/N)) * arclength for j ∈ 0:N]

u_grid = Vector{Vector{Float64}}(undef, N+1)
u̇_grid = Vector{Vector{Float64}}(undef, N+1)

# initialize

λ_bar = 1.0

x_init = 1.0
x_bar, converged = newton(x_init) do x
    return F(x, λ_bar), DxF(x, λ_bar)
end
converged || error("Newton failed")

u_bar = [x_bar, λ_bar]

u_grid[1] = u_bar

u̇_grid[1] = vec(nullspace([DxF(u_grid[1][1], u_grid[1][2]) DλF(u_grid[1][1], u_grid[1][2])]))
if u̇_grid[1][end] > 0 # decrease the parameter
    u̇_grid[1] .*= -1
end

# run continuation scheme

for j ∈ 2:N+1
    δ = arclength_grid[j] - arclength_grid[j-1]

    predictor = u_grid[j-1] + δ * u̇_grid[j-1]

    u_barⱼ, convergedⱼ = newton(predictor) do u
        return F_pa(u, u̇_grid[j-1], predictor), DF_pa(u, u̇_grid[j-1])
    end
    convergedⱼ || error("Newton failed")

    u_grid[j] = u_barⱼ

    u̇_grid[j] = vec(nullspace([DxF(u_grid[j][1], u_grid[j][2]) DλF(u_grid[j][1], u_grid[j][2])]))
    if sum(u̇_grid[j-1] .* u̇_grid[j]) < 0 # keep the same direction
        u̇_grid[j] .*= -1
    end
end

write_square_root(joinpath(outdir, "square_root.bin"), u_grid, u̇_grid)
println("square root: written")
end
end

# ------------------------------------------------------------------------
# Cahn–Hilliard (both parameterizations)

# Coons' parameterization in binary64 (the checker has its own interval version
# with the exact constants c_max = 1/10, g(c) = 79/2 + 150c²)

α(s) = (s + 1) / 2 # [-1,1] → [0,1]

θ_binary64(curves, s₁, s₂) =
    θ_edge(s₁, s₂, curves.γ₁, curves.γ₃) +
    θ_edge(-s₂, s₁, curves.γ₄, curves.γ₂) -
    θ_corners(s₁, s₂; P₁ = curves.γ₁(-1), P₂ = curves.γ₁(1), P₃ = curves.γ₃(-1), P₄ = curves.γ₃(1))

θ_edge(s₁, s₂, γᵢ, γⱼ) = γᵢ(s₁) + α(s₂) * (γⱼ(-s₁) - γᵢ(s₁))

θ_corners(s₁, s₂; P₁, P₂, P₃, P₄) =
    P₁ + α(s₁) * (P₂ - P₁) + α(s₂) * (P₄ - P₁ + α(s₁) * (P₃ - P₄ - P₂ + P₁))

γ_no_pinch₁(s) = [ 0.1 * (s - 3) / 4, 39.5 + 150 * (0.1 * (s - 3) / 4)^2]
γ_no_pinch₂(s) = [ 0.1 * s       / 2, 39.5 + 150 * (0.1 * s       / 2)^2]
γ_no_pinch₃(s) = [ 0.1 * (s + 3) / 4, 39.5 + 150 * (0.1 * (s + 3) / 4)^2]
γ_no_pinch₄(s) = [-0.1 * s, 41]

γ_with_pinch₁(s) = [-0.1, 41]
γ_with_pinch₂(s) = [ 0.1 * s, 39.5 + 150 * (0.1 * s)^2]
γ_with_pinch₃(s) = [ 0.1, 41]
γ_with_pinch₄(s) = [-0.1 * s, 41]

curves_no_pinch_binary64   = CahnHilliard.BoundaryCurves(γ_no_pinch₁, γ_no_pinch₂, γ_no_pinch₃, γ_no_pinch₄)
curves_with_pinch_binary64 = CahnHilliard.BoundaryCurves(γ_with_pinch₁, γ_with_pinch₂, γ_with_pinch₃, γ_with_pinch₄)

if "cahn_hilliard" ∈ examples
Random.seed!(20260922)

let θ = θ_binary64, F = CahnHilliard.F, DF = CahnHilliard.DF

for (case, curves) ∈ (("no_pinch", curves_no_pinch_binary64), ("pinch", curves_with_pinch_binary64))

N₁, N₂ = CahnHilliard.N₁, CahnHilliard.N₂

θ_grid = [θ(curves, -cospi(j₁/N₁), -cospi(j₂/N₂)) for j₁ ∈ 0:N₁, j₂ ∈ 0:N₂]

v_grid = Matrix{Sequence}(undef, N₁+1, N₂+1)

# initialize

c, β = θ_grid[1,1]

K = CahnHilliard.K
v_init = Sequence(evensym(Fourier(K, π)), [-0.1022666446473428 ; 0 ; 0.044213651292169816 ; [rand()/5.6^k for k = 3:K]])
v_bar, converged = newton(v_init) do v
    return F(v, c, β), DF(v, c, β)
end
converged || error("Newton failed")

v_grid[1,1] = v_bar

# run continuation scheme

for j₂ ∈ 1:N₂+1, j₁ ∈ 1:N₁+1
    if !(j₁ == j₂ == 1)
        cⱼ, βⱼ = θ_grid[j₁,j₂]

        predictor = j₁ == 1 ? v_grid[1,j₂-1] : v_grid[j₁-1,j₂]

        v_barⱼ, convergedⱼ = newton(predictor) do v
            return F(v, cⱼ, βⱼ), DF(v, cⱼ, βⱼ)
        end
        convergedⱼ || error("Newton failed")

        v_grid[j₁,j₂] = v_barⱼ
    end
end

write_cahn_hilliard(joinpath(outdir, "cahn_hilliard_$case.bin"), v_grid)
println("Cahn–Hilliard $case: written")
end
end
end

# ------------------------------------------------------------------------
# SKT

if "skt" ∈ examples
let F = SKT.F, DvF = SKT.DvF, F_pa = SKT.F_pa, DF_pa = SKT.DF_pa, tangent = SKT.tangent, K = SKT.K

guesses = reshape(read_floats(joinpath(@__DIR__, "skt_initial_guesses.bin"), 163 * length(SKT.ARC_NAMES)), 163, :)

# The last arc of branches 2 and 3 starts from the last node and tangent of the
# preceding arc, so the junction data coincide bit for bit (checker condition 4.J).
chained_to = Dict("u_bar2d" => "u_bar2c", "u_bar3d" => "u_bar3c")

for (i, name) ∈ enumerate(SKT.ARC_NAMES)

arclength = SKT.arclength_of(name)

N = SKT.arc_N(name)

arclength_grid = [0.5 * (1 - cospi(j/N)) * arclength for j ∈ 0:N]

u_grid = Vector{Sequence}(undef, N+1)
u̇_grid = Vector{Sequence}(undef, N+1)

# initialize

if haskey(chained_to, name)
    prev_grid, prev_tangent = read_skt(joinpath(outdir, "skt_$(chained_to[name]).bin"), chained_to[name])
    u_grid[1] = prev_grid[end]
    u̇_grid[1] = prev_tangent[end]
else
    data_name = guesses[:,i] # = skt_data.jld2[name]

    λ_bar = data_name[end]

    v_init = Sequence(evensym(Fourier(K, π)) × evensym(Fourier(K, π)), data_name[1:end-1])
    v_bar, converged = newton(v_init) do v
        return F(unpack(v), λ_bar), DvF(unpack(v), λ_bar)
    end
    converged || error("Newton failed")

    u_grid[1] = Sequence(evensym(Fourier(K, π)) × evensym(Fourier(K, π)) × ScalarSpace(),
        [coefficients(v_bar) ; λ_bar])

    direction = SKT.direction_of(name)
    u̇_grid[1] = tangent(u_grid[1])
    if u̇_grid[1][end] * direction < 0
        u̇_grid[1] *= -1
    end
end

# run continuation scheme

for j ∈ 2:N+1
    δ = arclength_grid[j] - arclength_grid[j-1]

    predictor = u_grid[j-1] + δ * u̇_grid[j-1]

    u_barⱼ, convergedⱼ = newton(predictor) do u
        return F_pa(unpack(u), unpack(u̇_grid[j-1]), unpack(predictor)),
               DF_pa(unpack(u), unpack(u̇_grid[j-1]))
    end
    convergedⱼ || error("Newton failed")

    u_grid[j] = u_barⱼ

    u̇_grid[j] = tangent(u_grid[j])
    if sum(coefficients(u̇_grid[j-1]) .* coefficients(u̇_grid[j])) < 0 # keep the same direction
        u̇_grid[j] *= -1
    end
end

write_skt(joinpath(outdir, "skt_$name.bin"), u_grid, u̇_grid)
println("SKT $name: written")
end
end
end

println("\nInputs written to $(abspath(outdir)).")
println("Now record the hashes:  julia --project=checker tools/hash_data.jl")
