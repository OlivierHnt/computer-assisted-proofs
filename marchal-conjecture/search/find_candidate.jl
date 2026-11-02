# Producer, not part of the proof. Finds the candidate x̄ on the Chebyshev
# nodes by Newton continuation in floating point and writes it to data/ as
# intervals. Any other way of producing the file would do.
#
# Usage (from the package root; this is the `produce` command of entry.json):
#   julia --project=search search/find_candidate.jl

import Pkg
Pkg.instantiate()

using RadiiPolynomial, LinearAlgebra
include(joinpath(@__DIR__, "..", "checker", "src", "MarchalCertificate.jl"))
using .MarchalCertificate

BLAS.set_num_threads(Sys.CPU_THREADS)

K = 70 # truncation order for the Fourier series
N = 20 # truncation order for the Chebyshev series
npts = N + 1

X = marchal_space(K, 1.0) # symmetric space of x = (a, β, α, u, v, w)

# initial data: Lagrange triangle at Ω = 1
# (stored coefficient = mode k ≥ 0 of the cosine or sine series)

x = zeros(ComplexF64, X)
u, v, w = component(x, 4), component(x, 5), component(x, 6)

component(u, 1)[2] =       3^(-1/6) / 2
component(u, 2)[2] = -im * 3^(-1/6) / 2
component(u, 3)[1] = 1 / 2

component(v, 1)[2] =  im * 3^(-1/6)
component(v, 2)[2] =       3^(-1/6)
component(v, 3)[1] =  im / 2

w[0] = 3^(-1/3)

# continuation in Ω over the Chebyshev–Lobatto nodes of [0, 1], from Ω = 1 down to Ω = 0

Ω_grid = [0.5 + 0.5cospi(j/N) for j ∈ 0:N]
x_grid = Vector{typeof(x)}(undef, npts)

F_cache  = zeros(ComplexF64, X)
DF_cache = zeros(ComplexF64, X, X)

for j ∈ 1:npts
    x_grid[j] = copy(j == 1 ? x : x_grid[j-1])

    _, success = newton!(x_grid[j], F_cache, DF_cache; tol = 1e-15, verbose = true) do F, DF, x
        F_marchal!(F, x, Ω_grid[j])
        DF_marchal!(DF, x, Ω_grid[j])
        return F, DF
    end
    success || error("Newton failed at node $j (Ω = $(Ω_grid[j]))")
    println("Success: Newton's method converged at the Chebyshev node Ωⱼ = $(Ω_grid[j])\n")

    # clean up the numerics: real-valued functions, β = α = 0
    impose_reality!(x_grid[j])
    x_grid[j][2] = 0
    x_grid[j][3] = 0
end

# ---- node values as thin intervals, exact data at the end nodes -----------

Xi = interval(X)
x_nodes = map(xⱼ -> Sequence(Xi, interval.(coefficients(xⱼ))), x_grid) # thin intervals

# at Ω = 1, the periodic orbit is the Lagrange triangle

coefficients(x_nodes[1]) .= interval(0)
u, v, w = component(x_nodes[1], 4), component(x_nodes[1], 5), component(x_nodes[1], 6)

component(u, 1)[2] =                 inv(sqrt(cbrt(interval(3)))) / interval(2)
component(u, 2)[2] = interval(-im) * inv(sqrt(cbrt(interval(3)))) / interval(2)
component(u, 3)[1] = interval(1) / interval(2)

component(v, 1)[2] = interval( im) * inv(sqrt(cbrt(interval(3))))
component(v, 2)[2] =                 inv(sqrt(cbrt(interval(3))))
component(v, 3)[1] = interval( im) / interval(2)

w[0] = inv(cbrt(interval(3)))

# at Ω = 0, the orbit lies in the yz-plane: u_x = v_x = 0 exactly
# (the figure-eight shape, Lemma 2.2 (iv), is proved by the checker)

coefficients(component(component(x_nodes[npts], 4), 1)) .= interval(0)
coefficients(component(component(x_nodes[npts], 5), 1)) .= interval(0)

# ---- write the frozen inputs ----------------------------------------------

datadir = joinpath(@__DIR__, "..", "data")
write_witness(datadir, x_nodes)

println("\nInputs written to $(abspath(datadir)).")
println("Now record the hashes:  julia --project=checker tools/hash_data.jl")
