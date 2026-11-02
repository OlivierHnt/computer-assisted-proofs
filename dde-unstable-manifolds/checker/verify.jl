# CHECKER — the proof.
#
# For each of the three examples: reads the frozen approximate eigenvalue λ̄,
# proves an enclosure of the unstable eigenvalue (radii polynomial theorem for
# the scalar characteristic function), computes the interval Taylor
# coefficients of the unstable manifold up to order N from the homological
# equations, and evaluates the bounds C_N, Y, Z₁, Z₂ of Lemma 3.3 and the radii
# polynomial of Theorem 3.1 (paper: Hénot, Lessard, Mireles James, JDDE 2022,
# revised version). Everything is evaluated in binary64 interval arithmetic
# (RadiiPolynomial / IntervalArithmetic); the only input is λ̄.
#
# Usage (from the package root):
#   julia --project=checker checker/verify.jl
# Prints each condition, in the order of statement.md, prefixed pass or fail.
# Exit code 0 iff every condition holds.

import Pkg
Pkg.instantiate()

using RadiiPolynomial, SHA, JSON

include(joinpath(@__DIR__, "src", "io.jl"))
include(joinpath(@__DIR__, "src", "vanderpol.jl"))
include(joinpath(@__DIR__, "src", "seir.jl"))
include(joinpath(@__DIR__, "src", "mackeyglass.jl"))

verdicts = Bool[]
function condition(label, ok, enclosure)
    push!(verdicts, ok)
    println(ok ? "pass " : "fail ", label, ": ", enclosure)
    return ok
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

# ---- the conditions, shared by the three cases ------------------------------
#
# Every enclosure must be bounded and carry IntervalArithmetic's guarantee flag
# (no non-interval number was mixed into its computation).

finite(x) = isbounded(x) && inf(x) ≥ 0 && isguaranteed(x)

# Theorem 3.1 with p(r) = Z₂ r²/2 + (Z₁ − 1) r + Y. RadiiPolynomial's
# interval_of_existence(Y, Z₁, Z₂, R) certifies p(r) ≤ 0 and Z₁ + Z₂ r < 1 for
# every r in the interval it returns, a subset of [0, R]; its radii polynomial
# is Y + (Z₁ − 1) r + Z₂ r²/2, the paper's convention, so Z₂ is passed as is.
# The Z₂ bound holds on the ball of radius R. The radius retained is
# r₀ = inf of the interval; the hypotheses of Theorem 3.1 are re-evaluated at
# r₀ as a direct check. For r ≥ 0, p(r) and Z₁ + Z₂ r increase with Y, Z₁, Z₂,
# so their suprema over the enclosures are their values at the upper endpoints;
# these and r₀ are binary64 numbers, and the check evaluates them exactly in
# rational arithmetic. (An interval evaluation would not do: r₀ is the root of p
# up to a few ulps, so p(r₀) ≈ −10⁻³² is below the rounding error of the
# evaluation when Y ∼ 10⁻¹⁶.) Used twice: for the eigenvalue (R = r_λ, giving ε)
# and for the tail of the Taylor series (R = r*, giving r₀).
function radii_polynomial(Y, Z₁, Z₂, R)
    I_r, success = interval_of_existence(Y, Z₁, Z₂, R)
    success && !isnai(I_r) && isbounded(I_r) || return false, I_r, NaN
    r₀ = inf(I_r)
    y, z₁, z₂, r = Rational{BigInt}.((sup(Y), sup(Z₁), sup(Z₂), r₀))
    p = z₂ * r^2 / 2 + (z₁ - 1) * r + y
    ok = isguaranteed(I_r) && p ≤ 0 && z₁ + z₂ * r < 1 && r₀ ≤ R
    return ok, I_r, r₀
end

# conditions 0a–0c: the eigenvalue enclosure λ ∈ {|λ − λ̄| ≤ ε} (see
# eigenvalue_bounds in src/common.jl); returns the enclosure as an interval box
function eigenvalue_conditions(id, M, λ̄)
    Y, Z₁, Z₂ = M.eigenvalue_bounds(λ̄)
    condition("$id.0a eigenvalue: Y = |A f(λ̄)| is a finite bound", finite(Y), "Y ∈ $(enc(Y))")
    condition("$id.0b eigenvalue: Z₁ = |1 − A f′(λ̄)| and Z₂ = |A| sup|f″| on |λ − λ̄| ≤ r_λ = $(M.r_λ) are finite bounds",
        finite(Z₁) && finite(Z₂), "Z₁ ∈ $(enc(Z₁)), Z₂ ∈ $(enc(Z₂))")
    ok, I_r, ε = radii_polynomial(Y, Z₁, Z₂, M.r_λ)
    λ = M.ball(λ̄, ε)
    ok &= inf(real(λ)) > 0
    condition("$id.0c eigenvalue (radii polynomial): f has exactly one zero λ in |λ − λ̄| ≤ ε, Re λ > 0", ok,
        ok ? "r ∈ $(enc(I_r)), ε = $ε (paper: $(M.ε_paper)), λ ∈ $(enc(λ))" : "interval_of_existence failed")
    return λ
end

function bounds_conditions(id, C, Y, Z₁, Z₂, r_star)
    condition("$id.1 C_N (Prop. A.1): denominator positive, ‖Ψ(⟨α,Λ⟩)⁻¹‖ ≤ C for |α| > N", inf(C) > 0 && isbounded(C) && isguaranteed(C), "C ∈ $(enc(C))")
    condition("$id.2 Y (Lemma 3.3, eq. (30)) is a finite bound", finite(Y), "Y ∈ $(enc(Y))")
    condition("$id.3 Z₁ (Lemma 3.3, eq. (31)) is a finite bound", finite(Z₁), "Z₁ ∈ $(enc(Z₁))")
    condition("$id.4 Z₂ (Lemma 3.3, eq. (32), r* = $r_star) is a finite bound", finite(Z₂), "Z₂ ∈ $(enc(Z₂))")
    ok, I_r, r₀ = radii_polynomial(Y, Z₁, Z₂, r_star)
    condition("$id.5 radii polynomial (Thm 3.1): p(r) ≤ 0 and Z₁ + Z₂ r < 1 on an interval ⊂ [0, r*]", ok,
        ok ? "r ∈ $(enc(I_r)), tail bound r₀ = $r₀" : "interval_of_existence failed")
    return r₀
end

# the largest coefficient of order N−2 ≤ |α| ≤ N (information only)
top_size(A, N) = maximum(sup(abs(x[α])) for x ∈ A for α ∈ indices(space(x)) if N - 2 ≤ sum(α) ≤ N)

# ---- I. delayed van der Pol (Section 4.1) -----------------------------------

let M = VanDerPol
    λ̄ = read_complex(joinpath(datadir, "vdp_lambda.bin"))
    println("I. delayed van der Pol, N = $(M.N), scaling $(round(mid(M.s); digits = 6)), r* = $(M.r_star), λ̄ = $(mid(real(λ̄))) + im*$(mid(imag(λ̄)))")
    λ = eigenvalue_conditions("I", M, λ̄)
    condition("I.0d eigenvector: Ψ(λ) ξ ∋ 0 for ξ = s (1, λ)", M.check_eigenvector(λ), "s = $(enc(M.s))")
    A, nonresonant, min_abs_f = M.solve_homological_equations(λ)
    condition("I.0e non-resonance: f(⟨α,Λ⟩) ≠ 0 for 2 ≤ |α| ≤ $(M.N) (Π_N(A) solved from (25) order by order)", nonresonant,
        "min |f(⟨α,Λ⟩)| ≥ $min_abs_f, max |A_α| over N−2 ≤ |α| ≤ N: $(top_size(A, M.N))")
    bounds_conditions("I", M.bounds(A, λ)..., M.r_star)
    println()
end

# ---- II. SEIR (Section 4.2) -------------------------------------------------

let M = SEIR
    λ̄ = read_real(joinpath(datadir, "seir_lambda.bin"))
    println("II. SEIR, N = $(M.N), scaling $(round(mid(M.s); digits = 6)), r* = $(M.r_star), λ̄ = $(mid(λ̄))")
    λ = eigenvalue_conditions("II", M, λ̄)
    condition("II.0d eigenvector: Ψ(λ) ξ ∋ 0 for ξ = s (ν, −(λ + μ))", M.check_eigenvector(λ), "s = $(enc(M.s))")
    A, nonresonant, min_abs_f = M.solve_homological_equations(λ)
    condition("II.0e non-resonance: f(kλ) ≠ 0 for 2 ≤ k ≤ $(M.N) (Π_N(A) solved from (10) order by order)", nonresonant,
        "min |f(kλ)| ≥ $min_abs_f, max |A_k| over N−2 ≤ k ≤ N: $(top_size(A, M.N))")
    bounds_conditions("II", M.bounds(A, λ)..., M.r_star)
    println()
end

# ---- III. Mackey–Glass, ρ = 3 (Section 4.3) ---------------------------------

let M = MackeyGlass
    λ̄ = read_real(joinpath(datadir, "mg_lambda.bin"))
    println("III. Mackey–Glass ρ = $(M.ρ), N = $(M.N), scaling $(round(mid(M.s); digits = 6)), r* = $(M.r_star), λ̄ = $(mid(λ̄))")
    λ = eigenvalue_conditions("III", M, λ̄)
    condition("III.0d eigenvector: Ψ(λ) ξ ∋ 0 for ξ = s (1, 0)", M.check_eigenvector(λ), "s = $(enc(M.s))")
    A, nonresonant, min_abs_f = M.solve_homological_equations(λ)
    condition("III.0e non-resonance: f(kλ) ≠ 0 for 2 ≤ k ≤ $(M.N) (Π_N(A) solved from (10) order by order)", nonresonant,
        "min |f(kλ)| ≥ $min_abs_f, max |A_k| over N−2 ≤ k ≤ N: $(top_size(A, M.N))")
    bounds_conditions("III", M.bounds(A, λ)..., M.r_star)
    println()
end

# ---- verdict --------------------------------------------------------------

all_pass = all(verdicts)
println(all_pass ? "ALL $(length(verdicts)) CONDITIONS PASS" : "$(count(!, verdicts)) OF $(length(verdicts)) CONDITIONS FAILED")
exit(all_pass ? 0 : 1)
