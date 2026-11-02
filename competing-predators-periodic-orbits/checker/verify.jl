# CHECKER — the proof.
#
# Reads the frozen inputs (floating-point approximations of the periodic orbits
# and of their Floquet normal forms on the Chebyshev nodes) and evaluates the
# conditions of statement.md: the Newton–Kantorovich bounds of Theorem 2.1 for
# the family of periodic orbits (Section 2.2) and for the Floquet normal form
# (Appendix A), the sign and monotonicity checks of Section 2.3, and the
# localisation of the Floquet exponents of Section 3. Nothing is searched for
# or solved here: only evaluated. The approximate inverses A, B are recomputed
# from the inputs in floating point (see statement.md).
#
# Usage (from the package root):
#   julia --project=checker checker/verify.jl
# Prints each condition, in the order of statement.md, prefixed pass or fail.
# Exit code 0 iff every condition holds.

import Pkg
Pkg.instantiate()

using RadiiPolynomial, LinearAlgebra, SHA, JSON
include(joinpath(@__DIR__, "src", "CompetingPredatorsCertificate.jl"))
using .CompetingPredatorsCertificate

BLAS.set_num_threads(min(2, Sys.CPU_THREADS))
IntervalArithmetic.configure(; nthreads = min(2, Sys.CPU_THREADS)) # threads of the interval matrix product
setdisplay(:infsup; sigdigits = 17, decorations = false) # print enclosures with all digits

verdicts = Bool[]
function condition(label, ok, enclosure)
    push!(verdicts, ok)
    println(ok ? "pass " : "fail ", label, ": ", enclosure)
    return ok
end

finite(x::Interval) = !isempty_interval(x) & isbounded(x)

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

# constants of the proof (Theorem 1.1 and Sections 2–3)

K     = 20 # Fourier truncation order, periodic orbits
K_lin = 20 # Fourier truncation order, Floquet normal form
N     = 30 # Chebyshev truncation order
npts  = CompetingPredatorsCertificate._npts(N) # Chebyshev nodes of the frozen grids (33)

X_po  = ScalarSpace()^3 × Fourier(K, 1.0)^3     # (τ, ζ₁, ζ₂, u₁, u₂, u₃)
X_lin = ScalarSpace()^9 × Fourier(K_lin, 1.0)^9 # (vec C, vec V)

d_po, d_lin = dimension(X_po), dimension(X_lin)
filesize(joinpath(datadir, "x_po_grid.bin"))  == npts * d_po  * 16 || error("x_po_grid.bin has the wrong size")
filesize(joinpath(datadir, "x_lin_grid.bin")) == npts * d_lin * 16 || error("x_lin_grid.bin has the wrong size")

function read_grid(f, X)
    raw = Vector{Float64}(undef, npts * dimension(X) * 2)
    read!(joinpath(datadir, f), raw)
    z = reinterpret(ComplexF64, ltoh.(raw))
    return [Sequence(X, z[(j-1)*dimension(X)+1:j*dimension(X)]) for j ∈ 1:npts]
end

x_po_grid  = read_grid("x_po_grid.bin",  X_po)
x_lin_grid = read_grid("x_lin_grid.bin", X_lin)
println("loaded x̄ (periodic orbits) and ῡ (Floquet normal form) on $npts Chebyshev nodes, dimensions $d_po and $d_lin")
println()

# parameters (1.4) of the paper, the range κ ∈ [κ₁, κ₂] = [92, 129], weight ν

γ = m₁ = m₂ = y₁ = y₂ = I"1.0"
d₁, d₂, a₁, a₂ = I"0.8", I"0.5", I"10.0", I"41.0"

κ₁ = interval(92.0)
κ₂ = interval(129.0)

ν = I"1.1"
R = 1e-6

# ---- conditions 1 and 2: Theorem 2.1 for F (periodic orbits) and G (Floquet) --

true_τ, true_ζ₁, true_ζ₂, true_C, bpo, blin =
    prove_contraction(x_po_grid, x_lin_grid, κ₁, κ₂, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂); K = K, K_lin = K_lin, N = N, ν = ν, R = R)

condition("1a. Y (Section 2.2.1, Lemma 2.2) is a finite bound", finite(bpo.Y), "Y ∈ $(bpo.Y), with ‖A‖ ∈ $(bpo.opnorm_A)")
condition("1b. Z₁ (Section 2.2.2) is a finite bound", finite(bpo.Z₁), "Z₁ ∈ $(bpo.Z₁)")
condition("1c. Z₂ (Section 2.2.3, Lemma 2.3 with R = $R) is a finite bound, the Lemma 2.3 denominators positive",
    finite(bpo.Z₂) & finite(bpo.rig_g₁) & finite(bpo.rig_g₂) & (inf(bpo.rig_g₁) > 0) & (inf(bpo.rig_g₂) > 0),
    "Z₂ ∈ $(bpo.Z₂); sup ‖(u₃ + αⱼ)⁻¹‖ ≤ 𝔠₁ ∈ $(bpo.rig_g₁), 𝔠₂ ∈ $(bpo.rig_g₂)")
check_ie_po = !isempty_interval(bpo.ie) & (0 ≤ inf(bpo.ie)) & (sup(bpo.ie) ≤ R) & (inf(bpo.ie) ≤ 1e-10)
condition("1d. Theorem 2.1 for F: the interval of existence in [0, R] is non-empty, and its left end r ≤ 1e-10", check_ie_po,
    "interval of existence $(bpo.ie); proved radius r = $(inf(bpo.r))")

condition("2a. Y for G (Appendix A) is a finite bound", finite(blin.Y), "Y ∈ $(blin.Y), with ‖B‖ ∈ $(blin.opnorm_A)")
condition("2b. Z₁ for G (Appendix A) is a finite bound", finite(blin.Z₁), "Z₁ ∈ $(blin.Z₁)")
condition("2c. Z₂ = 2‖B‖ for G (Appendix A, R = ∞) is a finite bound", finite(blin.Z₂), "Z₂ ∈ $(blin.Z₂)")
check_ie_lin = !isempty_interval(blin.ie) & (0 ≤ inf(blin.ie))
condition("2d. Theorem 2.1 for G: the interval of existence is non-empty", check_ie_lin,
    "interval of existence $(blin.ie); proved radius r = $(inf(blin.r))")

# ---- condition 3: positivity and boundary crossings (Section 2.3) ----------

h₁⁻, h₁⁺ = I"-0.95", I"-0.8"
h₂⁻, h₂⁺ = I"0.8", I"0.95"

# for display only: the enclosure whose sign `sign_derivative` (Lemma 2.6) takes
function derivative_enclosure(f, η, ν)
    N = order(f)[1]
    val = zero(Interval{Float64})
    err_finite = zero(Interval{Float64})
    for n ∈ 1:N
        val += n * real(f[(n,0)]) * sin(n * acos(η))
        err_finite += n * abs(sin(n * acos(η)))
    end
    return interval(val, total_error(f) * (err_finite + (ν - 1)*N + ν)/((ν - 1)^2 * ν^N); format = :midpoint)
end
hull_over(g, a, b, n) = reduce(hull, (g(η) for η ∈ mince(interval(a, b), n)))
fmt(x) = "[$(inf(x)), $(sup(x))]"

ok = prove_crossing(true_τ, true_ζ₁, true_ζ₂; a = h₁⁻, b = h₁⁺, n = 100, ν = ν)
condition("3a. crossing of Q₁ on η ∈ [h₁⁻, h₁⁺] = [-0.95, -0.8] (Lemmas 2.5, 2.6; 100 pieces)", ok,
    "ζ₂(h₁⁻) ∈ $(fmt(real(true_ζ₂(h₁⁻, 0)))), ζ₂(h₁⁺) ∈ $(fmt(real(true_ζ₂(h₁⁺, 0)))); " *
    "τ ∈ $(fmt(hull_over(η -> real(true_τ(η, 0)), h₁⁻, h₁⁺, 100))), ζ₁ ∈ $(fmt(hull_over(η -> real(true_ζ₁(η, 0)), h₁⁻, h₁⁺, 100))), " *
    "Σ n ζ₂,ₙ sin(n acos η) ∈ $(fmt(hull_over(η -> derivative_enclosure(true_ζ₂, η, ν), h₁⁻, h₁⁺, 100)))")

ok = prove_interior(true_τ, true_ζ₁, true_ζ₂; a = h₁⁺, b = h₂⁻, n = 100)
condition("3b. positivity on η ∈ [h₁⁺, h₂⁻] = [-0.8, 0.8] (Lemma 2.5; 100 pieces)", ok,
    "τ ∈ $(fmt(hull_over(η -> real(true_τ(η, 0)), h₁⁺, h₂⁻, 100))), ζ₁ ∈ $(fmt(hull_over(η -> real(true_ζ₁(η, 0)), h₁⁺, h₂⁻, 100))), " *
    "ζ₂ ∈ $(fmt(hull_over(η -> real(true_ζ₂(η, 0)), h₁⁺, h₂⁻, 100)))")

ok = prove_crossing(true_τ, true_ζ₂, true_ζ₁; a = h₂⁻, b = h₂⁺, n = 100, ν = ν)
condition("3c. crossing of Q₂ on η ∈ [h₂⁻, h₂⁺] = [0.8, 0.95] (Lemmas 2.5, 2.6; 100 pieces)", ok,
    "ζ₁(h₂⁻) ∈ $(fmt(real(true_ζ₁(h₂⁻, 0)))), ζ₁(h₂⁺) ∈ $(fmt(real(true_ζ₁(h₂⁺, 0)))); " *
    "τ ∈ $(fmt(hull_over(η -> real(true_τ(η, 0)), h₂⁻, h₂⁺, 100))), ζ₂ ∈ $(fmt(hull_over(η -> real(true_ζ₂(η, 0)), h₂⁻, h₂⁺, 100))), " *
    "Σ n ζ₁,ₙ sin(n acos η) ∈ $(fmt(hull_over(η -> derivative_enclosure(true_ζ₁, η, ν), h₂⁻, h₂⁺, 100)))")

# ---- condition 4: Floquet exponents (Section 3) -----------------------------
# prove_inside_box(C; a, b, at, box) returns (number of Gershgorin enclosures of
# the spectrum of C(at) inside the box, whether (C(η) - zI)⁻¹ is bounded for all
# z on the boundary of the box and all η ∈ [a, b]). The enclosures inside the
# box must also be disjoint from the others, so that the box contains exactly
# that many eigenvalues of C(at).

boxes_disjoint(z, w) = isdisjoint_interval(real(z), real(w)) | isdisjoint_interval(imag(z), imag(w))

# The eigenvalues are counted at η = `at` ∈ [a, b]: η = 0 for the interior interval
# [h₁⁺, h₂⁻] (Section 3, step 2), the left end point a for the two outer intervals.
function stability(label, interval_label, a, b, nexp, box; at = a)
    left_side, right_side, bottom_side, top_side = box
    count, contour = prove_inside_box(true_C; a = a, b = b, at = at, left_side = left_side, right_side = right_side, bottom_side = bottom_side, top_side = top_side)
    Λ = eigvals(evaluate.(true_C, Ref((at, 0))))
    Ω = complex(interval(left_side, right_side), interval(bottom_side, top_side))
    inside = [issubset_interval(λ, Ω) for λ ∈ Λ]
    # the Gershgorin discs inside Ω must be disjoint from those outside, otherwise their
    # union only bounds the total number of eigenvalues, not the number inside Ω
    separated = true
    for i ∈ eachindex(Λ), j ∈ eachindex(Λ)
        if inside[i] && !inside[j]
            separated &= boxes_disjoint(Λ[i], Λ[j])
        end
    end
    boxstr = "Ω = [$left_side, $right_side] × [$bottom_side, $top_side]i"
    condition("$(label).1 $nexp Floquet exponent(s) of C($(interval_label)) in $boxstr (Gershgorin)", (count == nexp) & separated,
        "enclosures of the spectrum of C($(interval_label)): $(join(string.(Λ), ", ")); $count inside Ω")
    condition("$(label).2 no Floquet exponent on ∂Ω for η ∈ [$(inf(a)), $(sup(b))] (adaptive subdivision)", contour,
        contour ? "(C(η) - zI)⁻¹ bounded for all z ∈ ∂Ω" : "could not bound (C(η) - zI)⁻¹ on ∂Ω")
end

stability("4a", "h₁⁻", h₁⁻, h₁⁺, 1, (-0.15, -0.05, -0.1, 0.1))
stability("4b", "0", h₁⁺, h₂⁻, 2, (-0.1, -0.0025, -0.1, 0.1); at = interval(0))
stability("4c", "h₂⁻", h₂⁻, h₂⁺, 1, (-0.04, -0.02, -0.1, 0.1))

# ---- verdict --------------------------------------------------------------

all_pass = all(verdicts)
println()
println(all_pass ? "ALL $(length(verdicts)) CONDITIONS PASS" : "$(count(!, verdicts)) OF $(length(verdicts)) CONDITIONS FAILED")
exit(all_pass ? 0 : 1)
