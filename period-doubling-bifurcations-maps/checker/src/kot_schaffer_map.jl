# The 10-mode Galerkin projection of the Kot–Schaffer model and its
# period-doubling map F_pd, verbatim from kot_schaffer.jl of
# OlivierHnt/PeriodDoublingMaps.jl (lines 1–30). Shared by checker/verify.jl
# and search/find_candidate.jl. Every function is generic in the number type:
# called with Float64 data (search; approximate inverse) or with
# Interval{Float64} data (checker).
#
# Unknown X = (x, v, λ, w, μ, h) ∈ ℝ^(4n+2):
#   x = X[1:n], v = X[n+1:2n], λ = X[2n+1], w = X[2n+2:3n+1], μ = X[3n+2], h = X[3n+3:4n+2]

using RadiiPolynomial, LinearAlgebra

n = 10
S = evensym(Fourier(n-1, 1.0))
M(κ, x) = coefficients(project(Multiplication(κ * Sequence(space(κ), x)), space(κ), space(κ)))

f(x, λ, β, κ) = λ * (β .* (x - M(κ, x) * x))
Df(x, λ, β, κ) = λ * Diagonal(β) * (exact(I) - exact(2) * M(κ, x))
B(u, λ, β, κ) = -exact(2) * λ * Diagonal(β) * M(κ, u)

parts(X) = X[1:n], X[n+1:2n], X[2n+1], X[2n+2:3n+1], X[3n+2], X[3n+3:4n+2]

function F(X, β, κ)
    x, v, λ, w, μ, h = parts(X)
    A = Df(x, λ, β, κ)
    return [f(x, λ, β, κ) - x ; A * v + v ; v' * v - exact(1) ; A' * w - μ * w ; w' * v - exact(1) ; A * h - h + B(v, λ, β, κ) * v]
end

function DF(X, β, κ)
    x, v, λ, w, μ, h = parts(X)
    A, Bv, Bh = Df(x, λ, β, κ), B(v, λ, β, κ), B(h, λ, β, κ)
    N = reduce(hcat, [B(e, λ, β, κ)' * w for e = eachcol(Matrix(one(λ) * exact(I), n, n))])
    O, o, z = zeros(eltype(X), n, n), zeros(eltype(X), n), zero(λ)
    return [A - exact(I)   O               f(x, λ, β, κ) / λ        O                 o     O
            Bv             A + exact(I)    A * v / λ                O                 o     O
            o'             exact(2) * v'   z                        o'                z     o'
            N              O               A' * w / λ               A' - μ * exact(I) -w    O
            o'             w'              z                        v'                z     o'
            Bh             exact(2) * Bv   (A * h + Bv * v) / λ     O                 o     A - exact(I)]
end

# Model data (paper, Section 6): β_k = 2^-k, κ₀ = 0.8, κ₁ = -0.2, κ_k = 0 (k ≥ 2).
β = [2.0^-k for k = 0:n-1]
κ = Sequence(S, [0.8 ; -0.2 ; zeros(n-2)])

# ------------------------------------------------------------------
# Input I/O. X_bar.bin: the 4n+2 = 42 entries of X̄ = (x, v, λ, w, μ, h) in
# that order, raw little-endian IEEE binary64, no header.
# ------------------------------------------------------------------

write_Xbar(path, X::Vector{Float64}) = open(io -> write(io, htol.(X)), path, "w")

function read_Xbar(path)
    filesize(path) == (4n+2) * 8 || error("$path has the wrong size")
    return ltoh.(reinterpret(Float64, read(path)))
end
