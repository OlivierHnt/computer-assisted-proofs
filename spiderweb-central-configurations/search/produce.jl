# PRODUCER — not part of the proof.
#
# Replays the continuation loop of the original `spiderweb.proof_𝑛ℓ` (2019
# code, src/spiderweb.jl) for λ = -1, m₀ = 0, mᵢ = 1/l, n ∈ 5:100,
# l ∈ 18:200, and writes the Newton output x̄(n, l) of every case to
# data/radii.bin. The only change to the loop is that each call
# `radii(x, m0, m, λ, n, l)` (the proof, done by checker/) is replaced by
# recording a copy of x. Everything else, including the order of the loops,
# the initial guesses and the in-place updates of x (which also act on the
# saved vector x2, since `x2 = x` aliases it), is kept as in the original
# (the identifiers 𝑛, ℓ of the 2019 code are written n, l).
#
# Usage (from the package root):
#   julia --project=search search/produce.jl
#
# Byte layout of data/radii.bin: raw little-endian IEEE binary64, no header;
# for l = 18, …, 200 (outer), for n = 5, …, 100 (inner): x̄(n, l)[1:n].

import Pkg
Pkg.instantiate()

using LinearAlgebra
BLAS.set_num_threads(1) # fixed, so that the LAPACK solves in Newton's method are reproducible

include(joinpath(@__DIR__, "src", "spiderweb.jl"))
using .spiderweb: newton

const N_RANGE = collect(5:1:100)
const L_RANGE = collect(18:1:200)

function produce_nl(λ::Float64, n_range::Vector, l_range::Vector)
	sols = Dict{Tuple{Int,Int},Vector{Float64}}()
	record(x, n, l) = (sols[(n, l)] = copy(x))
	# --- body of spiderweb.proof_𝑛ℓ, radii(...) ↦ record(...) ---
	m0 = 0. # set central mass to zero
	x = float.(1:n_range[1])
	x2 = x
	for l ∈ l_range
	    for n ∈ n_range
	        m = ones(Float64, n)./l # set m_i = 1/l for i > 0
	        if n == n_range[1] && l == l_range[1] # initial guess for Newton's method
	            x = newton(x, m0, m, λ, n, l)
	            record(x, n, l)
	            x2 = x # saving solution with the smallest number of circles
			elseif n > n_range[1]
				push!(x, x[n-1]^2/x[n-2])
				x .*= range(.98, stop = .95, length = n) # correction for better convergence
				x = newton(x, m0, m, λ, n, l)
				record(x, n, l)
	        elseif n == n_range[1] && l > l_range[1]
	            x = x2[1:n]
	            x .*= 0.8
	            x = newton(x, m0, m, λ, n, l)
	            record(x, n, l)
	            x2 = x # saving solution with the smallest number of circles
	        end
	    end
	    println("l = ", l, " done (n = ", n_range[1], ":", n_range[end], ")"); flush(stdout)
	end
	return sols
end

# optional arguments for a partial run: n_max l_max output_path
n_max = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 100
l_max = length(ARGS) ≥ 2 ? parse(Int, ARGS[2]) : 200
out   = length(ARGS) ≥ 3 ? ARGS[3] : joinpath(@__DIR__, "..", "data", "radii.bin")
n_range = filter(≤(n_max), N_RANGE)
l_range = filter(≤(l_max), L_RANGE)

t = @elapsed sols = produce_nl(-1., n_range, l_range)

open(out, "w") do io
    for l ∈ l_range, n ∈ n_range
        x = sols[(n, l)]
        length(x) == n || error("wrong length for (n, l) = ($n, $l)")
        write(io, htol.(x))
    end
end
println("wrote ", out, " (", filesize(out), " bytes) in ", round(t; digits = 1), " s")
