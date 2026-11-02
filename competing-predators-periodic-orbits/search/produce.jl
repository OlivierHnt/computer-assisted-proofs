# SEARCH — not part of the proof.
#
# Finds the numerical approximation of the family of periodic orbits and of
# its Floquet normal form on the Chebyshev nodes (ODE integration, Newton
# refinement, continuation in κ: `construct_tube`), then writes them to data/.
#
# Everything here may be replaced by any other method of producing an
# approximation; the checker does not care how the data were obtained.
#
# Usage (from the package root; this is the `produce` command of entry.json):
#   julia --project=search search/produce.jl
# then record the hashes with
#   julia --project=checker tools/hash_data.jl

import Pkg
Pkg.instantiate()

module CompetingPredatorsSearch
    using RadiiPolynomial
    using LinearAlgebra
    using OrdinaryDiffEq
include(joinpath(@__DIR__, "..", "checker", "src", "zero_finding_problems.jl"))
include(joinpath(@__DIR__, "src", "numerics.jl"))
    export construct_tube
end

using .CompetingPredatorsSearch, RadiiPolynomial, LinearAlgebra

BLAS.set_num_threads(min(2, Sys.CPU_THREADS))

# parameters of Theorem 1.1, starting point of the continuation, truncation orders

γ = m₁ = m₂ = y₁ = y₂ = 1.0
d₁, d₂, a₁, a₂ = 0.8, 0.5, 10.0, 41.0

κ_init = 95.0
u_init = [0.32090923676599215, 0.0428131951743777, 0.3426766250219656]
τ_init = 3.5

K = 20
K_lin = 20
N = 30

x_po_grid, x_lin_grid, κ₁, κ₂ = construct_tube(u_init, τ_init, κ_init, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂); increment = 0.1, K = K, K_lin = K_lin, N = N, round = true)

# the checker hard-codes the range [κ₁, κ₂] = [92, 129] of Theorem 1.1
(κ₁, κ₂) == (92.0, 129.0) || @warn "the continuation gave [κ₁, κ₂] = [$κ₁, $κ₂], not [92, 129]: the checker will evaluate these data on [92, 129]"

# raw little-endian binary64, node by node (η = -1 first), each coefficient as (re, im)

datadir = joinpath(@__DIR__, "..", "data")
mkpath(datadir)
for (f, grid) ∈ (("x_po_grid.bin", x_po_grid), ("x_lin_grid.bin", x_lin_grid))
    open(joinpath(datadir, f), "w") do io
        for x ∈ grid
            write(io, htol.(reinterpret(Float64, coefficients(x))))
        end
    end
    println("wrote data/", f, " (", length(grid), " nodes × ", length(grid[1]), " complex coefficients)")
end
