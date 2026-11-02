# CHECKER — the proof.
#
# For each of the 18 cases (Bianchi type VIII and IX, N₁N₂N₃ = 10^k for
# k = 1, 0, -1, …, -7), reads the frozen approximate zero x̄ of F and the
# constants C₁, C₂, ν, R, and evaluates the conditions of statement.md with the
# original proof code of HLPeriodicOrbits_VIII_IX.jl (checker/src/proof.jl):
# the bounds Y, Z₁, Z₂ of Proposition B.1, the radii polynomial, and the checks
# (ii) giving the hypotheses of Lemma 4.2, plus the reality of the zero and
# the Bianchi type. Nothing is searched for or solved
# here: only evaluated (the approximate inverse A = inv(mid(DF(x̄))) is
# recomputed, as in the original; see statement.md).
#
# Usage (from the package root; Julia 1.13, see checker/Manifest.toml):
#   julia --project=checker checker/verify.jl <case>     (one of the 18 cases)
# This is the `run` command of entry.json, whose `cases` field declares
# case = 1:18: the checker is run once per case, prints the conditions of that
# case in the order of statement.md, prefixed pass or fail, and exits 0 iff
# they all hold. The theorem holds iff all 18 runs exit 0.

import Pkg
Pkg.instantiate()

using RadiiPolynomial, LinearAlgebra, SHA, JSON
include(joinpath(@__DIR__, "src", "HLCertificate.jl"))
using .HLCertificate

BLAS.set_num_threads(min(2, Sys.CPU_THREADS))
# Interval matrix products use IntervalArithmetic 1.0's default :fast mode
# (Rump's midpoint-radius product, rigorous); the original's entrywise products
# (:slow) are 2-3x slower and gave the same verdicts. Integer powers stay
# correctly rounded, as in IntervalArithmetic 0.20 (which used BigFloat).
IntervalArithmetic.configure(; matmul = :fast, power = :slow)
precision(BigFloat) == 256 || error("the proof uses the default BigFloat precision of 256 bits")
rounding(BigFloat) == RoundNearest || error("the proof uses the default BigFloat rounding mode")

verdicts = Bool[]
function condition(label, ok, enclosure)
    push!(verdicts, ok)
    println(ok ? "pass " : "fail ", label, ": ", enclosure)
    flush(stdout)
    return ok
end

# every digit of both endpoints, with the decoration and the IntervalArithmetic
# 1.0 "guaranteed" flag (false = a non-interval number entered the computation;
# here only exact constants, such as the integer coefficients of the original
# formulas; see statement.md)
setdisplay(:full)
full(x::Interval) = sprint(show, x)
full(x) = repr(x)

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

# ---- the 18 cases (constants of the proof: type, k, Fourier order n, number format)

const CASES = [
    # case  type    k    file tag  n    format
    ( 1, "VIII",  1, "10",    32, :mp256),
    ( 2, "VIII",  0, "1",     32, :f64),
    ( 3, "VIII", -1, "1e-1",  32, :f64),
    ( 4, "VIII", -2, "1e-2",  32, :f64),
    ( 5, "VIII", -3, "1e-3",  32, :f64),
    ( 6, "VIII", -4, "1e-4",  64, :f64),
    ( 7, "VIII", -5, "1e-5",  64, :f64),
    ( 8, "VIII", -6, "1e-6", 128, :f64),
    ( 9, "VIII", -7, "1e-7", 256, :mp256),
    (10, "IX",    1, "10",    32, :f64),
    (11, "IX",    0, "1",     64, :f64),
    (12, "IX",   -1, "1e-1",  64, :f64),
    (13, "IX",   -2, "1e-2",  64, :f64),
    (14, "IX",   -3, "1e-3", 128, :f64),
    (15, "IX",   -4, "1e-4", 192, :f64),
    (16, "IX",   -5, "1e-5", 192, :f64),
    (17, "IX",   -6, "1e-6", 192, :mp256),
    (18, "IX",   -7, "1e-7", 256, :mp256),
]

exact(v) = Rational{BigInt}(v)                     # exact value of a binary float
pow10(k) = k ≥ 0 ? big(10)^k // 1 : 1 // big(10)^(-k)
contains_exactly(X, q) = exact(inf(X)) ≤ q ≤ exact(sup(X))

# |x̄ - †x̄|_X, with †x = (conj(a₋ₖ))ₖ for the five Fourier components and
# conj for γ, η₁, η₂ (condition c.5): the symmetrisation of the original code
# updates the coefficients in place, in increasing order of k, so x̄ is only
# approximately conjugate-symmetric.
function asymmetry(x, ν)
    n = order(space(space(x)[1]))
    d = copy(x)
    for j ∈ 1:5
        a, b = component(component(x, 1), j), component(component(d, 1), j)
        for k ∈ -n:n
            b[k] = a[k] - conj(a[-k])
        end
    end
    for i ∈ 1:3
        p = component(x, 2)[i]
        component(d, 2)[i] = p - conj(p)
    end
    X = NormedCartesianSpace((NormedCartesianSpace(ℓ¹(GeometricWeight(ν)), ℓ∞()), NormedCartesianSpace(ℓ∞(), ℓ∞())), ℓ∞())
    return norm(d, X)
end

length(ARGS) == 1 || error("usage: julia --project=checker checker/verify.jl <case>, case = 1, …, 18")
selected = [parse(Int, ARGS[1])]
issubset(selected, [c for (c,) ∈ CASES]) || error("cases are numbered 1 to 18")

for (c, type, k, tag, n, fmt) ∈ CASES
    c ∈ selected || continue
    file = "$(type)_N1N2N3_$(tag).bin"
    lab(j) = "$c.$j $type N₁N₂N₃=10^$k"

    x, C₁, C₂, ν, R = read_case(joinpath(datadir, file), n, fmt)

    # condition c.0 — the stored constants are the stated ones
    ok0 = contains_exactly(C₁, 1 // 1) && contains_exactly(C₂, pow10(k)) && inf(ν) ≥ 1 && R > 0
    condition(lab(0) * " inputs ($file, n = $n, $fmt)", ok0,
        "C₁ = $(full(C₁)) ∋ 1, C₂ = $(full(C₂)) ∋ 10^$k, ν = $(full(ν)), R = $R")

    out = prove_case(x, C₁, C₂, ν, R)   # symmetrises x in place, as the original
    Y, Z₁, Z₂, ie = out[:Y], out[:Z₁], out[:Z₂], out[:ie]

    # conditions c.1–c.4 — Proposition B.1 (radii polynomial)
    condition(lab(1) * " Y (Prop. B.1) is a finite bound", isfinite(sup(Y)), "Y ∈ $(full(Y))")
    condition(lab(2) * " Z₁ (Prop. B.1) is a finite bound", isfinite(sup(Z₁)), "Z₁ ∈ $(full(Z₁))")
    condition(lab(3) * " Z₂ (Prop. B.1) is a finite bound", isfinite(sup(Z₂)), "Z₂ ∈ $(full(Z₂))")
    condition(lab(4) * " radii polynomial: nonempty interval of existence in [0, R]", !isempty_interval(ie),
        "interval of existence $(full(ie)), r₀ = inf = $(isempty_interval(ie) ? "—" : inf(ie))")

    # condition c.5 — the zero is real (Proposition B.1, last part): †x̃ is a
    # zero of F with |†x̃ - x̄| ≤ r₀ + |x̄ - †x̄|, and F has a unique zero in
    # the ball of radius sup(ie) around x̄, so r₀ + |x̄ - †x̄| ≤ sup(ie) gives
    # †x̃ = x̃.
    # conditions c.6–c.11 — the hypotheses of Lemma 4.2 (checks 1–6 of the
    # original code), with r₀ the error bound; c.12 — the Bianchi type
    if isempty_interval(ie)
        for (j, name) ∈ ((5, "x̃ is real"), (6, "γ ∈ (0, +∞)"), (7, "Σ₊ ≢ 0"), (8, "N₁ ≢ 0"), (9, "N₃ ≢ 0"),
                         (10, "N₁ has no zero on [0, 2π]"), (11, "N₃ has no zero on [0, 2π]"),
                         (12, "Bianchi type $type"))
            condition(lab(j) * " " * name, false, "not evaluated: condition $c.4 failed")
        end
    else
        δ = asymmetry(x, ν)
        lhs = interval(inf(ie)) + δ
        condition(lab(5) * " x̃ is real: r₀ + |x̄ - †x̄|_X ≤ sup(interval of existence)", sup(lhs) ≤ sup(ie),
            "|x̄ - †x̄|_X ∈ $(full(δ)), r₀ + |x̄ - †x̄|_X ≤ $(sup(lhs)) ≤ $(sup(ie))")
        r₀ = out[:r₀]
        γ  = out[:γ]
        condition(lab(6) * " γ ∈ (0, +∞) (γ real and γ > R ≥ r₀)", out[:check_1],
            "γ ∈ $(full(real(γ))) + $(full(imag(γ)))i, R = $(full(out[:R_ii]))")
        for (j, key, name) ∈ ((7, :check_2, "Σ₊"), (8, :check_3, "N₁"), (9, :check_4, "N₃"))
            at = out[Symbol(key, :_at)]
            condition(lab(j) * " $name is not identically zero", out[key],
                at === nothing ? "no point t ∈ 0:1e-3:2π with 0 ∉ $(name)(t) ± r₀" :
                    "$(name)(t) ∈ $(full(at[2])) ± r₀ at t = $(full(at[1])), r₀ = $r₀")
        end
        for (j, key, name, bnd) ∈ ((10, :check_5, "N₁", :bound_DN₁), (11, :check_6, "N₃", :bound_DN₃))
            steps, f_0, min_f, t_end = out[Symbol(key, :_info)]
            condition(lab(j) * " $name has no zero on [0, 2π]", out[key],
                "$steps steps up to t = $t_end, $(name)(0) ∈ $(f_0 === nothing ? "—" : full(f_0)) ± r₀, min |$(name)(tⱼ)| ≥ $min_f, sup |$(name)'| ≤ $(sup(out[bnd]))")
        end
        # condition c.12 — the signs of N₁, N₃ (constant by c.10, c.11) give the
        # Bianchi type, N₂ = 10^k / (N₁N₃) having the sign of N₁N₃: type IX iff
        # N₁ > 0 and N₃ > 0; type VIII otherwise.
        N₁_0, N₃_0 = out[:check_5_info][2], out[:check_6_info][2]
        sgn(f) = f === nothing ? 0 : inf(f - r₀) > 0 ? 1 : sup(f + r₀) < 0 ? -1 : 0
        s₁, s₃ = sgn(N₁_0), sgn(N₃_0)
        ok12 = out[:check_5] && out[:check_6] && s₁ != 0 && s₃ != 0 &&
            (type == "IX" ? (s₁ == 1 && s₃ == 1) : !(s₁ == 1 && s₃ == 1))
        signs = s₁ == 0 || s₃ == 0 ? "undetermined" :
            join((s₁ > 0 ? "+" : "-", s₁ * s₃ > 0 ? "+" : "-", s₃ > 0 ? "+" : "-"), ", ")
        condition(lab(12) * " Bianchi type $type", ok12, "signs of (N₁, N₂, N₃) = ($signs)")
    end
end

# ---- verdict --------------------------------------------------------------

all_pass = all(verdicts)
println()
println("case ", selected[1], ": ", all_pass ? "ALL $(length(verdicts)) CONDITIONS PASS" : "$(count(!, verdicts)) OF $(length(verdicts)) CONDITIONS FAILED")
exit(all_pass ? 0 : 1)
