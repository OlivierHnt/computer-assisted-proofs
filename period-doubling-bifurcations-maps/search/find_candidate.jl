# SEARCH — not part of the proof.
#
# Finds the numerical approximation X̄ = (x̄, v̄, λ̄, w̄, μ̄, h̄) of a zero of the
# period-doubling map F_pd by floating-point Newton iterations and writes it to
# data/X_bar.bin. Any other way of producing X̄ would do.
# Verbatim from OlivierHnt/PeriodDoublingMaps.jl/kot_schaffer.jl, lines 32–45.
#
# Usage (from the package root; this is the `produce` command of entry.json):
#   julia --project=search search/find_candidate.jl
# then re-record the hash with
#   julia --project=checker tools/hash_data.jl

import Pkg
Pkg.instantiate()

using LinearAlgebra
BLAS.set_num_threads(min(2, Sys.CPU_THREADS))

include(joinpath(@__DIR__, "..", "checker", "src", "kot_schaffer_map.jl"))

# nontrivial fixed point at λ ≈ 2.7, critical eigenvectors, and h from (2.4)

λ_init = 2.7
x_init, _ = newton(x -> (f(x, λ_init, β, κ) - x, Df(x, λ_init, β, κ) - I), [0.9 ; zeros(n-1)])
ν, V = eigen(Df(x_init, λ_init, β, κ))
k = argmin(abs.(ν .+ 1))
v_init = V[:,k] / norm(V[:,k])
w_init = inv(V)[k,:] * norm(V[:,k])
h_init = -(Df(x_init, λ_init, β, κ) - I) \ (B(v_init, λ_init, β, κ) * v_init)

# Newton on F_pd

X_bar, converged = newton(X -> (F(X, β, κ), DF(X, β, κ)), [x_init ; v_init ; λ_init ; w_init ; -1.0 ; h_init]; verbose = true)
converged || error("Newton's method did not converge")

out = joinpath(@__DIR__, "..", "data", "X_bar.bin")
mkpath(dirname(out))
write_Xbar(out, X_bar)
println("wrote ", normpath(out), " (", filesize(out), " bytes)")
