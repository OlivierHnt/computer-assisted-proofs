# CHECKER — the proof.
#
# Reads the frozen inputs (data/*.bin: for each of the four solutions, the
# approximate solution x̄, the equation parameters, the equilibrium c and the
# eigen-data Λ, Γ as intervals, and the constants of the proof) and evaluates
# the conditions of statement.md: the hypotheses on the inputs (0a–0d), the
# three inequalities of Proposition 2.2 (1a–1c, plus 1d with a ψ̂ including |c|
# and |Γ|), the bounds Y, Z₁, Z₂ and the radii polynomial of Theorem 3.3 (2–5),
# and the C⁰ error of Corollary 3.5 against the value stated in Theorems
# 1.1–1.3 (6). Nothing is searched for or solved: the only auxiliary objects
# computed in floating point are the approximate inverses A = inv(mid(DF)) (as
# in the original code) and the 4×4 A_eig of Appendix B.
#
# Usage (from the package root):
#   julia --project=checker checker/verify.jl <case>...    one or several of
#                        kg0, kg1, sh, fhn (the `cases` of entry.json); none
#                        means all four
# Case names: kg0 (Klein–Gordon m = 0), kg1 (Klein–Gordon m = 1),
# sh (Swift–Hohenberg), fhn (FitzHugh–Nagumo). kg1 and fhn need a large server
# (see README.md).
# Prints each condition of the cases run, in the order of statement.md,
# prefixed pass or fail. Exit code 0 iff every one of them holds. The cases are
# independent; the theorems need all four.

import Pkg
Pkg.instantiate()

using RadiiPolynomial, SHA, JSON
include(joinpath(@__DIR__, "src", "EllipticPDELocalizedRadialSol.jl"))
using .EllipticPDELocalizedRadialSol

verdicts = Bool[]
function condition(label, ok, enclosure)
    push!(verdicts, ok)
    println(ok ? "pass " : "fail ", label, ": ", enclosure)
    flush(stdout)
    return ok
end

# ---- BLAS and directed rounding ----------------------------------------------
# IntervalArithmetic 1.0 encloses interval matrix products (its default `:fast`
# mode, Rump's midpoint–radius algorithm) with BLAS products computed under
# directed rounding, in a dedicated OpenBLAS built to propagate the rounding
# mode to its threads (OpenBLASConsistentFPCSR_jll). Condition 0 tests that the
# configured number of threads gives the same enclosures as one thread, and
# that these enclose exact products. CHECKER_BLAS_THREADS=n sets the number of
# threads of both BLAS libraries (default: IntervalArithmetic's default).

using LinearAlgebra: BLAS
if haskey(ENV, "CHECKER_BLAS_THREADS")
    IntervalArithmetic.configure(nthreads = parse(Int, ENV["CHECKER_BLAS_THREADS"]))
    BLAS.set_num_threads(parse(Int, ENV["CHECKER_BLAS_THREADS"]))
end
blas_threads = IntervalArithmetic.configuration_options.nthreads

function blas_rounding_test(nt)
    n = 1024
    A = [((i*7919 + j*104729) % 1009) / 3 - 168.0 for i ∈ 1:n, j ∈ 1:n]
    B = [((i*104723 + j*7907) % 1013) / 7 - 72.0 for i ∈ 1:n, j ∈ 1:n]
    prod_(t) = (IntervalArithmetic.configure(nthreads = t); interval.(A) * interval.(B))
    C₁ = prod_(1)
    Cₙ = prod_(nt)
    same = all(isequal_interval.(C₁, Cₙ))
    # exact products of a few entries (BigFloat, 256 bits, exact for these data)
    encl = setprecision(BigFloat, 256) do
        all(in_interval(sum(big(A[i,k]) * big(B[k,j]) for k ∈ 1:n), Cₙ[i,j]) for (i, j) ∈ ((1, 1), (17, 513), (512, 3), (1024, 1024), (700, 300)))
    end
    return same && encl,
        "interval product of two 1024×1024 matrices: $nt-thread enclosure bitwise equal to 1-thread: $same; encloses the exact entries sampled: $encl"
end
ok₀, msg₀ = blas_rounding_test(blas_threads)
condition("0 BLAS honours directed rounding with $blas_threads thread(s)", ok₀, msg₀)

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

# ---- the four cases ----------------------------------------------------------
# number, short name, title, theorem, C⁰ distance stated in the theorem (as a
# rational p/q, compared rigorously)

META = Dict(
    "klein_gordon_m_0" => (k = 1, short = "kg0", title = "Klein–Gordon m=0", thm = "Theorem 1.1", stated = (29, 10^8)),   # 2.9e-7
    "klein_gordon_m_1" => (k = 2, short = "kg1", title = "Klein–Gordon m=1", thm = "Theorem 1.1", stated = (55, 10^7)),   # 5.5e-6
    "swift_hohenberg"  => (k = 3, short = "sh",  title = "Swift–Hohenberg",  thm = "Theorem 1.2", stated = (24, 10^6)),   # 2.4e-5
    "fitzhugh_nagumo"  => (k = 4, short = "fhn", title = "FitzHugh–Nagumo",  thm = "Theorem 1.3", stated = (98, 10^8)),   # 9.8e-7
)

selected = isempty(ARGS) ? [META[c.name].short for c ∈ CASES] : ARGS
all(s -> s ∈ [m.short for m ∈ values(META)], selected) || error("unknown case in $ARGS; use kg0, kg1, sh, fhn")

_full(x) = EllipticPDELocalizedRadialSol._full(x)

for case ∈ CASES
    m = META[case.name]
    m.short ∈ selected || continue
    println()
    data = read_case(joinpath(datadir, case.name * ".bin"), case)
    eq = data.equation
    println("== case $(m.k): $(m.title) on ℝ^$(case.d) ($(m.thm)); n_T = $(case.n_T_num) + $(case.n_T_pad), n_C = $(case.n_C_num) + $(case.n_C_pad)")
    tag(l) = replace(l, r"^(\S+) " => SubstitutionString("\\1 $(m.title): "))
    cond(l, ok, e) = condition(tag(l), ok, e)
    t_case = time()

    # 0. hypotheses on the inputs (added; see statement.md)
    local ok, msg
    ok, msg = check_parameters(data);         cond("$(m.k).0a constants satisfy the hypotheses of Theorem 3.3", ok, msg)
    ok, msg = check_symmetry(data);           cond("$(m.k).0b x̄ = S(x̄) (Lemma 3.4)", ok, msg)
    ok, msg = check_equilibrium(eq, data);    cond("$(m.k).0c parameters, equilibrium c and eigen-data (Λ, Γ) enclose the exact ones", ok, msg)

    # 1–5. the original proof (Proposition 2.2 and Theorem 3.3)
    after_LP = (LP, aux) -> begin
        ok, msg = check_LP_strict(eq, data.c, data.Λ, data.Γ, aux.λ̂, aux.δ, aux.μ, interval(data.ℒ_x), interval(data.ℒ_y))
        cond("$(m.k).1d Proposition 2.2 with ψ̂ including |c| and |Γ|", ok, msg)
    end
    res = prove(data, (condition = cond, after_LP = after_LP), "$(m.k)")

    # 6. Corollary 3.5: C⁰ distance, compared with the value stated in the theorem
    stated = interval(m.stated[1]) / interval(m.stated[2])
    cond("$(m.k).6 C⁰ error (Corollary 3.5) ≤ $(inf(stated)) ($(m.thm))", res.error_sup ≤ inf(stated),
        "max(ρ̄, |Γ|(ρ̄ + (1+ℒ_y)μ)) ≤ $(res.error_sup) with ρ̄ = $(isempty_interval(res.ρ) ? "∅" : inf(res.ρ)) (original code prints $(res.error))")

    println("info case $(m.k) timings (s): ", join(["$k $(round(v; sigdigits = 3))" for (k, v) ∈ res.timings], ", "),
        "; total $(round(time() - t_case; sigdigits = 3)); max RSS so far $(round(Sys.maxrss() / 2^30; sigdigits = 3)) GiB")
end

# ---- verdict --------------------------------------------------------------

all_pass = all(verdicts)
println()
println(all_pass ? "ALL $(length(verdicts)) CONDITIONS PASS" : "$(count(!, verdicts)) OF $(length(verdicts)) CONDITIONS FAILED")
println("cases run: ", join(selected, ", "), " (the theorems need every case of `cases` in entry.json)")
exit(all_pass ? 0 : 1) # 0 iff every condition of the cases run holds
