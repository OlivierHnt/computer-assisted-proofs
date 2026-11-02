# CHECKER, one case — the proof of the case (n, l).
#
# Reads the frozen Newton output x̄(n, l) (data/radii.bin) and evaluates in
# interval arithmetic the conditions of statement.md for that case: the ball
# lies in the domain, and the radii polynomial of Theorem 6.1 (bounds Y₀, Z₀
# of Section 6.1, Z₂ of Lemma 6.3) certified by RadiiPolynomial. Nothing is
# searched for or solved here: only evaluated. See src/radii_check.jl
# (`check_case`) and src/common.jl.
#
# Usage (from the package root):
#   julia --project=checker checker/verify.jl <n> <l>     # n ∈ 5:100, l ∈ 18:200
# Prints the conditions of the case, prefixed pass or fail. Exit code 0 iff
# both hold. The theorem is the conjunction over all cases: entry.json
# ("run" with "cases") runs this script once per (n, l) ∈ 5:100 × 18:200.

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
