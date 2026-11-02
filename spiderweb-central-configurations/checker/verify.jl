# Checker for one case (n, l). Reads x̄(n, l) from data/radii.bin and checks
# in interval arithmetic the two conditions of statement.md: the ball lies in
# the domain, and the radii polynomial of Theorem 6.1 (Y₀, Z₀ of Section 6.1,
# Z₂ of Lemma 6.3). Nothing is searched for or solved here.
#
# Usage (from the package root):
#   julia --project=checker checker/verify.jl <n> <l>     # n ∈ 5:100, l ∈ 18:200
# Prints each condition prefixed pass or fail; exit code 0 iff both hold.
# The theorem needs every case (entry.json: "run" with "cases").

import Pkg
Pkg.instantiate()

using RadiiPolynomial, LinearAlgebra, SHA, JSON
include(joinpath(@__DIR__, "src", "spiderweb.jl"))
using .Spiderweb
include(joinpath(@__DIR__, "src", "radii_check.jl"))
include(joinpath(@__DIR__, "src", "common.jl"))

length(ARGS) == 2 || error("usage: julia --project=checker checker/verify.jl <n> <l>")
n, l = parse(Int, ARGS[1]), parse(Int, ARGS[2])
(n ∈ N_RANGE && l ∈ L_RANGE) || error("(n, l) = ($n, $l) is outside n ∈ $N_RANGE, l ∈ $L_RANGE")

x̄ = load_inputs(joinpath(@__DIR__, "..", "data"))
println()

ok = report_case(n, l, x̄(n, l))
println()
println(ok ? "CASE (n,l) = ($n,$l) PASSES" : "CASE (n,l) = ($n,$l) FAILS",
        " — the theorem needs every n ∈ $N_RANGE, l ∈ $L_RANGE (entry.json: run × cases)")
exit(ok ? 0 : 1)
