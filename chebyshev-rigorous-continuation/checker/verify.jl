# CHECKER — the proof.
#
# Reads the frozen inputs (the floating-point approximate solutions produced by
# the continuation, on the Chebyshev nodes), recomputes the approximate inverses
# in floating point exactly as the original scripts do, encloses everything in
# interval arithmetic, and evaluates the conditions of statement.md: the bounds
# Y, Z₁, Z₂ of Theorem 1.1 (applied to the lifted map on ℓ¹ in the Chebyshev
# variable, Theorem 2.11) and the radii-polynomial inequalities, for the four
# examples of the paper. Nothing is searched for or solved here: only evaluated
# (the approximate inverses are a matrix inversion of the frozen data; the
# conditions hold whatever their exact value).
#
# Usage (from the package root):
#   julia --project=checker checker/verify.jl <example>...   one or several of
#                        cube_root, square_root, cahn_hilliard, skt (the `cases`
#                        of entry.json); none means all four
# Prints each condition of the examples run, in the order of statement.md,
# prefixed pass or fail. Exit code 0 iff every one of them holds. The four
# results are independent (one example each); all four are the proof.

import Pkg
Pkg.instantiate()

using RadiiPolynomial, LinearAlgebra, SHA, JSON
include(joinpath(@__DIR__, "src", "ChebyshevContinuationCertificate.jl"))
using .ChebyshevContinuationCertificate

BLAS.set_num_threads(min(2, Sys.CPU_THREADS))
# the interval matrix products use IntervalArithmetic's own OpenBLAS, which defaults to
# every CPU thread on Apple silicon; cap it like BLAS
RadiiPolynomial.IntervalArithmetic.configure(nthreads = min(2, Sys.CPU_THREADS))

verdicts = Bool[]
function condition(label, ok, enclosure)
    push!(verdicts, ok)
    println(ok ? "pass " : "fail ", label, ": ", enclosure)
    flush(stdout)
    return ok
end

# intervals printed at full binary64 precision (shortest round-trip decimal of each endpoint)
fmt(x) = "[" * repr(inf(x)) * ", " * repr(sup(x)) * "]"

# the four conditions of one application of Theorem 1.1
function report(prefix, b)
    condition("$prefix.1 $(b.label): Y is a finite bound", isbounded(b.Y) && inf(b.Y) ≥ 0, "Y ∈ $(fmt(b.Y))")
    condition("$prefix.2 $(b.label): Z₁ is a finite bound", isbounded(b.Z₁) && inf(b.Z₁) ≥ 0, "Z₁ ∈ $(fmt(b.Z₁))")
    condition("$prefix.3 $(b.label): Z₂ is a finite bound (r⋆ = $(b.r_star))", isbounded(b.Z₂) && inf(b.Z₂) ≥ 0, "Z₂ ∈ $(fmt(b.Z₂))")
    condition("$prefix.4 $(b.label): radii polynomial, interval of existence non-empty in [0, r⋆]", b.proved && !isempty_interval(b.ie),
        "interval of existence $(fmt(b.ie)); error bound r = $(inf(b.ie))")
end

# ---- load and authenticate the frozen inputs --------------------------------

root    = joinpath(@__DIR__, "..")
datadir = joinpath(root, "data")
hashes  = JSON.parsefile(joinpath(datadir, "hashes.json"))

for f ∈ sort(collect(keys(hashes)))
    p = joinpath(datadir, f)
    isfile(p) || error("$f is missing: run `julia --project=checker tools/fetch_data.jl` first")
    bytes2hex(open(sha256, p)) == hashes[f]["sha256"] || error("sha256 mismatch for $f: not the recorded input")
    println("pass sha256 ", f)
end
println()

EXAMPLES = ["cube_root", "square_root", "cahn_hilliard", "skt"] # the `cases` of entry.json
examples = isempty(ARGS) ? EXAMPLES : ARGS
all(x -> x ∈ EXAMPLES, examples) || error("unknown example in $ARGS; use ", join(EXAMPLES, ", "))
t_start = time()

# ------------------------------------------------------------------------
# 1. Cubic root, Section 2.2.3

if "cube_root" ∈ examples; let
    x_grid = read_cube_root(joinpath(datadir, "cube_root.bin"))
    b = CubeRoot.cube_root_certificate(x_grid)
    report("1", (; b..., label = "cube root"))
end; end

# ------------------------------------------------------------------------
# 2. Square root, pseudo-arclength, Section 3.1.2

if "square_root" ∈ examples; let
    u_grid, u̇_grid = read_square_root(joinpath(datadir, "square_root.bin"))
    b = SquareRoot.square_root_certificate(u_grid, u̇_grid)
    report("2", (; b..., label = "square root"))
end; end

# ------------------------------------------------------------------------
# 3. Cahn–Hilliard, two-parameter continuation, Section 4.1

if "cahn_hilliard" ∈ examples
    for (i, case, curves) ∈ ((1, "no_pinch", CahnHilliard.curves_no_pinch), (2, "pinch", CahnHilliard.curves_with_pinch))
        label = "Cahn–Hilliard, $(replace(case, "_" => " "))"
        v_grid = read_cahn_hilliard(joinpath(datadir, "cahn_hilliard_$case.bin"))
        b = CahnHilliard.ch_certificate(curves, v_grid)
        # θ̄ = (c̄, β̄) encloses the exact Coons map θ of Λ: the corner identities of Coons'
        # construction hold and θ̄ agrees with the pointwise evaluation of θ at the nodes
        condition("3.$i.0 $label: θ̄ encloses Coons' map θ (corner identities, node values)", b.θ_ok,
            "corner mismatch ≤ $(repr(b.θ_mismatch)); Chebyshev order of θ̄ = $(b.N_θ); grids $(b.m_Y) (Y), $(b.m_Z) (Z₁)")
        report("3.$i", (; b..., label))
    end
end

# ------------------------------------------------------------------------
# 4. SKT system, pseudo-arclength, fourteen arcs, Section 4.2

if "skt" ∈ examples
    arcs = Dict{String,Any}()
    for (i, name) ∈ enumerate(SKT.ARC_NAMES)
        u_grid, u̇_grid = read_skt(joinpath(datadir, "skt_$name.bin"), name)
        b = SKT.skt_certificate(name, u_grid, u̇_grid)
        report("4.$i", (; b..., label = "SKT $name"))
        arcs[name] = (; b..., u_grid, u̇_grid)
    end

    # Connection of consecutive arcs of the same branch (statement.md, Result 4): the
    # last node of the arc a (s = -1) and the first node of the arc b (s = 1) carry the
    # same point and the same tangent, so F_pa(·, -1) of a and F_pa(·, 1) of b are the
    # same map G on 𝒳_pa; u⋆_b(1) is a zero of G within r_b of that point and u⋆_a(-1)
    # is the unique zero of G within r_a⁺ (Theorem 2.11): they coincide if r_b ≤ r_a⁺
    # (or, symmetrically, r_a ≤ r_b⁺).
    branches = [["u_bar1a", "u_bar1b", "u_bar1c"], ["u_bar2a", "u_bar2b", "u_bar2c", "u_bar2d"],
                ["u_bar3a", "u_bar3b", "u_bar3c", "u_bar3d"], ["u_bar4a", "u_bar4b", "u_bar4c"]]
    junctions = [(a, b) for branch ∈ branches for (a, b) ∈ zip(branch[1:end-1], branch[2:end])]
    for (j, (a, b)) ∈ enumerate(junctions)
        A, B = arcs[a], arcs[b]
        same_point   = coefficients(A.u_grid[end]) == coefficients(B.u_grid[1])
        same_tangent = coefficients(A.u̇_grid[end]) == coefficients(B.u̇_grid[1])
        r_a, r_a⁺ = inf(A.ie), sup(A.ie)
        r_b, r_b⁺ = inf(B.ie), sup(B.ie)
        connected = A.proved && B.proved && same_point && same_tangent && (r_b ≤ r_a⁺ || r_a ≤ r_b⁺)
        condition("4.J$j SKT $a → $b: same node data at the junction and the balls nest", connected,
            "same point: $same_point, same tangent: $same_tangent; r_$(a[6:end]) = $r_a, r_$(a[6:end])⁺ = $r_a⁺, r_$(b[6:end]) = $r_b, r_$(b[6:end])⁺ = $r_b⁺")
    end
end

# ---- verdict --------------------------------------------------------------

all_pass = all(verdicts)
println()
println(all_pass ? "ALL $(length(verdicts)) CONDITIONS PASS" : "$(count(!, verdicts)) OF $(length(verdicts)) CONDITIONS FAILED")
println("examples run: ", join(examples, ", "), " (the results need every example of `cases` in entry.json)")
println("checker time: ", round(time() - t_start; digits = 1), " s")
exit(all_pass ? 0 : 1) # 0 iff every condition of the examples run holds
