# CHECKER — the proof.
#
# Reads the frozen numerical approximation X̄ of a zero of the period-doubling
# map F_pd, and evaluates in interval arithmetic the conditions of
# statement.md: the Newton–Kantorovich bounds Y, Z (Theorem 5.1 of the paper),
# the Gershgorin conditions (Theorem 5.2), the normal-form coefficient c and
# the crossing speed ν̇_c (Theorem 5.3, Theorem 6.1). Nothing is searched for
# or solved here: only evaluated. The arithmetic is that of kot_schaffer.jl
# (OlivierHnt/PeriodDoublingMaps.jl, lines 47–74), unchanged; only the input
# and the printing differ.
#
# Usage (from the package root):
#   julia --project=checker checker/verify.jl
# Prints each condition, in the order of statement.md, prefixed pass or fail.
# Exit code 0 iff every condition holds.

import Pkg
Pkg.instantiate()

using LinearAlgebra
BLAS.set_num_threads(min(2, Sys.CPU_THREADS))

using RadiiPolynomial, SHA, JSON
include(joinpath(@__DIR__, "src", "kot_schaffer_map.jl")) # n, S, f, Df, B, F, DF, β, κ

verdicts = Bool[]
function condition(label, ok, enclosure)
    push!(verdicts, ok)
    println(ok ? "pass " : "fail ", label, ": ", enclosure)
    return ok
end

# intervals printed as [inf, sup] with Julia's shortest round-trip Float64
# representation, i.e. bit-exactly
fmt(x::Interval) = "[$(inf(x)), $(sup(x))]"

# exact comparison of a Float64 with a decimal number p / 10^q
exact_ge(x::Float64, p, q) = Rational{BigInt}(x) ≥ big(p) // big(10)^q
exact_le(x::Float64, p, q) = Rational{BigInt}(x) ≤ big(p) // big(10)^q

# ---- load and authenticate the frozen input ---------------------------------

root    = joinpath(@__DIR__, "..")
datadir = joinpath(root, "data")
hashes  = JSON.parsefile(joinpath(datadir, "hashes.json"))

for file ∈ sort(collect(keys(hashes)))
    p = joinpath(datadir, file)
    isfile(p) || error("$file is missing: run `julia --project=checker tools/fetch_data.jl` first")
    bytes2hex(open(sha256, p)) == hashes[file]["sha256"] || error("sha256 mismatch for $file: not the recorded input")
    println("pass sha256 ", file)
end

X_bar = read_Xbar(joinpath(datadir, "X_bar.bin"))
println("loaded X̄ ∈ ℝ^$(length(X_bar)) (n = $n), λ̄ = $(X_bar[2n+1]), μ̄ = $(X_bar[3n+2])")
println()

# ---- auxiliary object: the approximate inverse Q ----------------------------
# Any matrix will do for Theorem 5.1; it is recomputed here in floating point,
# as in the original script, and then used as a thin interval matrix.

Q = inv(DF(X_bar, β, κ))

# ---- enclosures of the data --------------------------------------------------
# β_k = 2^-k is exact in binary64; κ₀ = 0.8 and κ₁ = -0.2 are enclosed.

β_i = interval.(β)
κ_i = Sequence(interval(S), [I"0.8" ; I"-0.2" ; zeros(Interval{Float64}, n-2)])
X_i = interval.(X_bar)
Q_i = interval.(Q)

# ---- Theorem 5.1 (Newton–Kantorovich), step (1) of Theorem 5.3 ---------------

Y = norm(Q_i * F(X_i, β_i, κ_i), Inf)
condition("1. Y ≥ ‖Q F_pd(X̄)‖∞ (Theorem 5.1)", isbounded(Y) && isguaranteed(Y), "Y ∈ $(fmt(Y))")

R = 1e-8
X_R = interval.(X_bar, R; format = :midpoint)
Z₁ = opnorm(exact(I) - Q_i * DF(X_R, β_i, κ_i), Inf)
condition("2. Z ≥ sup over B(X̄, r⋆ = $R) of ‖I - Q DF_pd(X)‖∞ (Theorem 5.1)", isbounded(Z₁) && isguaranteed(Z₁), "Z ∈ $(fmt(Z₁))")

# RadiiPolynomial.interval_of_existence(Y, Z, R) checks 0 ≤ Y, 0 ≤ Z < 1 and
# sup(Y/(1 - Z)) ≤ R, and returns [sup(Y/(1 - Z)), R]; so r₀ := inf(ie) is a
# binary64 upper bound of Y/(1 - Z). Strictness r₀ < r⋆ is checked on top.
ie, proved = interval_of_existence(Y, Z₁, R)
condition("3. Z < 1 and r₀ = Y/(1 - Z) < r⋆ = $R (Theorem 5.1)", proved & isguaranteed(ie) & (inf(ie) < R),
    "r₀ ≤ $(inf(ie)); interval of existence $(fmt(ie))")

# ---- Theorem 5.2 (Gershgorin), step (2) of Theorem 5.3 -----------------------
# A_r encloses A = D_x f(x⋆, λ⋆) over the ball B(X̄, r₀); V approximates the
# eigenvectors of its midpoint (auxiliary, recomputed in floating point), and
# 𝓜 encloses V⁻¹ A V. Since V is real, each disk 𝒟_k is centred on the real
# axis and is represented by its trace on ℝ, the real interval D[k].

x_r, v_r, λ_r, w_r, μ_r, h_r = parts(interval.(X_bar, inf(ie); format = :midpoint))
A_r = Df(x_r, λ_r, β_i, κ_i)

V_i = interval.(eigvecs(mid.(A_r)))
𝓜 = inv(V_i) * A_r * V_i
D = [𝓜[k,k] + interval(-1, 1) * sum(abs(𝓜[k,j]) for j = 1:n if j ≠ k) for k = 1:n]
k₀ = argmin(abs.(mid.(D) .+ 1))
spectrum_verified = all(isguaranteed, D) && in_interval(-1, D[k₀]) && all(isdisjoint_interval(D[k], D[k₀]) && sup(abs(D[k])) < 1 for k = 1:n if k ≠ k₀)

others = [k for k = 1:n if k ≠ k₀]
check_4a = (eltype(V_i) === Interval{Float64}) && all(isguaranteed, D) && in_interval(-1, D[k₀])
check_4b = all(isdisjoint_interval(D[k], D[k₀]) for k ∈ others)
check_4c = all(sup(abs(D[k])) < 1 for k ∈ others)
(check_4a && check_4b && check_4c) == spectrum_verified || error("internal: Gershgorin sub-conditions disagree with the original test")

condition("4a. disk 𝒟_$k₀ contains -1 (Theorem 5.2)", check_4a,
    "D[$k₀] = $(fmt(D[k₀])); largest disk radius $(maximum(radius, D)); V real: $(eltype(V_i) === Interval{Float64})")
condition("4b. 𝒟_$k₀ is disjoint from the other $(n-1) disks (Theorem 5.2)", check_4b,
    "other disks ⊆ $(fmt(reduce(hull, D[others])))")
condition("4c. the other $(n-1) disks lie in {|z| < 1} (Theorem 5.2)", check_4c,
    "max |z| over the other disks ≤ $(maximum(k -> sup(abs(D[k])), others))")

# ---- normal-form coefficient, step (3) of Theorem 5.3 ------------------------
# C ≡ 0 since f is quadratic, so (2.5) reads c = -w⋆ᵀ B(v⋆, h⋆).

c = -w_r' * (B(v_r, λ_r, β_i, κ_i) * h_r)
condition("5. c = -w⋆ᵀB(v⋆,h⋆) < 0 (Theorem 5.3 (3), supercritical)", isguaranteed(c) && sup(c) < 0, "c ∈ $(fmt(c))")

# ---- crossing speed (4.5), with ẋ from (2.3) ---------------------------------

ẋ = -inv(A_r - exact(I)) * f(x_r, λ_r, β_i, κ_i) / λ_r
ν̇ = w_r' * (B(v_r, λ_r, β_i, κ_i) * ẋ + A_r * v_r / λ_r)
condition("6. ν̇_c = w⋆ᵀB(ẋ,v⋆) + w⋆ᵀA_λ v⋆ < 0 (Theorem 4.1, (4.5); Proposition 2.2)", isguaranteed(ν̇) && sup(ν̇) < 0, "ν̇_c ∈ $(fmt(ν̇))")

# ---- the enclosures stated in Theorem 6.1 ------------------------------------
# compared exactly with the decimal endpoints of the theorem

condition("7a. r₀ ≤ 2.1e-14 (Theorem 6.1)", exact_le(inf(ie), 21, 15), "r₀ ≤ $(inf(ie))")
condition("7b. λ⋆ ∈ [2.70423246926594, 2.70423246926599] (Theorem 6.1)",
    isguaranteed(λ_r) && exact_ge(inf(λ_r), 270423246926594, 14) && exact_le(sup(λ_r), 270423246926599, 14), "λ⋆ ∈ $(fmt(λ_r)); a⋆,₀ ∈ $(fmt(x_r[1]))")
condition("7c. c ∈ [-18.6289496098, -18.6289496097] (Theorem 6.1)",
    exact_ge(inf(c), -186289496098, 10) && exact_le(sup(c), -186289496097, 10), "c ∈ $(fmt(c))")
condition("7d. ν̇_c ∈ [-0.9888090931, -0.9888090930] (Theorem 6.1)",
    exact_ge(inf(ν̇), -9888090931, 10) && exact_le(sup(ν̇), -9888090930, 10), "ν̇_c ∈ $(fmt(ν̇))")

# ---- verdict --------------------------------------------------------------

all_pass = all(verdicts)
println()
println(all_pass ? "ALL $(length(verdicts)) CONDITIONS PASS" : "$(count(!, verdicts)) OF $(length(verdicts)) CONDITIONS FAILED")
exit(all_pass ? 0 : 1)
