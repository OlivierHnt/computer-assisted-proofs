# CHECKER — the proof.
#
# Reads the frozen floating-point approximations (equilibrium c̄₁, approximate
# eigendecompositions at c₀ and c₁, Taylor coefficients of the local stable
# and unstable manifolds, Chebyshev coefficients of the connecting orbit) and
# evaluates, in interval arithmetic, the conditions of statement.md: for each of
# the six subproblems the Y and Z bounds of the paper and the Radii Polynomial
# Theorem 1.1, each subproblem using the validated enclosures of the previous
# ones. The approximate inverses are recomputed here in floating point from the
# frozen data. Nothing is searched for or solved here: no Newton, no
# eigensolver, no ODE solver.
#
# Usage (from the package root):
#   julia --project=checker checker/verify.jl
# Prints each condition, in the order of statement.md, prefixed pass or fail.
# Exit code 0 iff every condition holds.

import Pkg
Pkg.instantiate()

using RadiiPolynomial, LinearAlgebra, SHA, JSON
include(joinpath(@__DIR__, "src", "ShimizuMoriokaCertificate.jl"))
using .ShimizuMoriokaCertificate

BLAS.set_num_threads(2)
setdisplay(:infsup; sigdigits = 17) # print enclosures with every digit of their endpoints

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

function readbin(name, n) # n little-endian binary64 numbers
    p = joinpath(datadir, name)
    filesize(p) == 8n || error("$name has the wrong size")
    v = Vector{Float64}(undef, n)
    open(io -> read!(io, v), p)
    return ltoh.(v)
end
pairs_to_complex(v) = complex.(v[1:2:end], v[2:2:end])

# constants of the proof (paper, Sections 2–4)

parameters = (I"0.75", I"0.45") # (a, b) = (3/4, 9/20)
K_Q, ν_Q = 15, interval(1.1)    # stable manifold of c₀: Taylor order, weight
K_P, ν_P = 22, interval(1.1)    # unstable manifold of c₁: Taylor order, weight
K_u, μ   = 125, interval(1.05)  # connecting orbit: Chebyshev order, weight
τ        = interval(13)         # integration time

c1_bar = readbin("c1_bar.bin", 3)
eig_c0 = readbin("eig_c0.bin", 3 + 9)
λ0_approx, V0_approx = eig_c0[1:3], reshape(eig_c0[4:12], 3, 3)
eig_c1 = pairs_to_complex(readbin("eig_c1.bin", 2 * (3 + 9)))
λ1_approx, V1_approx = eig_c1[1:3], reshape(eig_c1[4:12], 3, 3)
Q_bar = Sequence((Taylor(K_Q) ⊗ Taylor(K_Q))^3, readbin("Q_bar.bin", 3 * (K_Q + 1)^2))
P_data = Sequence((Taylor(K_P) ⊗ Taylor(K_P))^3, pairs_to_complex(readbin("P_bar.bin", 2 * 3 * (K_P + 1)^2)))

# Symmetry P̄_{(m,n)} = conj(P̄_{(n,m)}) of P̄ (Remark (real_image) of the revised
# manuscript). statement.md deduces P⋆(γ(α)) ∈ ℝ³ from the uniqueness of the
# zero in a ball centred at P̄, which needs P̄ exactly symmetric; the frozen data
# are symmetric up to rounding noise only, so P̄ is symmetrised here, bit for
# bit (condition 5e reports the deviation). RadiiPolynomial's
# `conjugacy_symmetry!` implements a different symmetry and is not used.
function swap_conjugate_symmetrise(P)
    P_sym = copy(P)
    K = order(component(P, 1))[1]
    for i ∈ 1:3, m ∈ 0:K, n ∈ m:K
        a = (component(P, i)[(m, n)] + conj(component(P, i)[(n, m)])) / 2
        m == n && (a = complex(real(a), 0.0))
        component(P_sym, i)[(m, n)] = a
        component(P_sym, i)[(n, m)] = conj(a)
    end
    return P_sym
end
function is_swap_conjugate_symmetric(P)
    K = order(component(P, 1))[1]
    return all(component(P, i)[(m, n)] == conj(component(P, i)[(n, m)]) for i ∈ 1:3, m ∈ 0:K, n ∈ 0:K)
end
P_bar = swap_conjugate_symmetrise(P_data)
P_data_deviation = maximum(abs, coefficients(P_bar - P_data)) # floating point, reported in 5e
P_data_tolerance = 1e-15 # stated in statement.md: the data must be P̄ up to rounding noise
x_bar = Sequence(Chebyshev(K_u) × Chebyshev(K_u) × Chebyshev(K_u) × ScalarSpace() × ScalarSpace() × ScalarSpace(),
                 readbin("x_bar.bin", 3 * (K_u + 1) + 3))
println("loaded c̄₁, the approximate eigenpairs at c₀ and c₁, Q̄ (Taylor order $K_Q), P̄ (Taylor order $K_P, symmetrised), x̄ (Chebyshev order $K_u)")
println()

finite(x) = isbounded(x) && !isempty_interval(x)

# ---- 1. the equilibrium c₁ (Section 2.1, Theorem 1.1) -----------------------

c0_star = interval([0, 0, 0]) # exact equilibrium

eq = validate(EquilibriumProblem(f, Df, parameters), c1_bar)
c1_star = eq.c_star

condition("1a. c₁ equilibrium: Y = ‖A f(c̄₁)‖₁ is a finite bound", finite(eq.Y), "Y ∈ $(eq.Y)")
condition("1b. c₁ equilibrium: Z = ‖I - A Df(Box(c̄₁, R))‖₁ is a finite bound", finite(eq.Z), "Z ∈ $(eq.Z), R = $(eq.R)")
condition("1c. c₁ equilibrium: Theorem 1.1 (Z < 1, Y/(1-Z) ≤ R)", eq.proof_success,
    "interval of existence $(eq.ie); c₁⋆ ∈ $c1_star")

# ---- 2.–3. the eigenpairs at c₀ and c₁ (Section 2.2) ------------------------

eig0 = validated_eigen(Df(c0_star, parameters), λ0_approx, V0_approx)
eig1 = validated_eigen(Df(c1_star, parameters), λ1_approx, V1_approx)
λ0, V0 = eig0.eig
λ1, V1 = eig1.eig
# each eigenpair is validated on ℂ³ × ℂ (ℝ³ × ℝ when real) with the 1-norm; a
# conjugate pair is validated once, the second eigenpair being its conjugate
fmt_certs(e, field) = join(["λ$(j): $(c.conjugate_of == 0 ? getfield(c, field) : "conj. of λ$(c.conjugate_of)")" for (j, c) ∈ enumerate(e.certs)], "; ")
validated(e) = [c for c ∈ e.certs if c.conjugate_of == 0]

# the type of the equilibria, which includes the hypothesis Re(λ₁)Re(λ₂) > 0 of
# Proposition 3.3 for the two eigenvalues used for each manifold
spectrum_c0 = eltype(λ0) <: Interval && sup(λ0[1]) < 0 && sup(λ0[2]) < 0 && inf(λ0[3]) > 0
spectrum_c1 = isthinzero(imag(λ1[1])) && sup(real(λ1[1])) < 0 &&
    inf(real(λ1[2])) > 0 && inf(real(λ1[3])) > 0 && !in_interval(0, imag(λ1[2])) &&
    isequal_interval(real(λ1[2]), real(λ1[3])) && isequal_interval(imag(λ1[2]), -imag(λ1[3]))

# Conjugacy of the eigenpairs of c₁ (needed for the symmetry of P⋆, condition
# 5e): the data must satisfy λ̄₁, v̄₁ real and (λ̄₃, v̄₃) = (conj λ̄₂, conj v̄₂)
# bit for bit, and the enclosures must then be conjugate too (real parts equal,
# imaginary parts opposite).
eigdata_c1_conjugate = imag(λ1_approx[1]) == 0 && all(imag.(V1_approx[:,1]) .== 0) &&
    λ1_approx[3] == conj(λ1_approx[2]) && V1_approx[:,3] == conj(V1_approx[:,2])
eigvecs_c1_conjugate = all(isequal_interval.(real.(V1[:,2]), real.(V1[:,3]))) &&
    all(isequal_interval.(imag.(V1[:,2]), -imag.(V1[:,3]))) && all(isthinzero.(imag.(V1[:,1])))

for (n, name, e, spectrum, type) ∈ ((2, "c₀", eig0, spectrum_c0, "c₀ is a saddle: λ₁, λ₂ < 0 < λ₃, all real"),
                                    (3, "c₁", eig1, spectrum_c1 && eigdata_c1_conjugate && eigvecs_c1_conjugate,
                                     "c₁ is a saddle-focus: λ₁ < 0 real, λ₂ = conj(λ₃), Re λ₂ > 0, Im λ₂ ≠ 0; v₁ real, v₃ = conj(v₂) (data bitwise and enclosures)"))
    condition("$(n)a. $name eigenpairs: Y = ‖A F_eig(x̄)‖₁ is a finite bound for each validated eigenpair",
        all(finite(c.Y) for c ∈ validated(e)), fmt_certs(e, :Y))
    condition("$(n)b. $name eigenpairs: Z = ‖I − A DF_eig(Box(x̄, R))‖₁ is a finite bound for each validated eigenpair",
        all(finite(c.Z) for c ∈ validated(e)), fmt_certs(e, :Z))
    condition("$(n)c. $name eigenpairs: Theorem 1.1 (Z < 1, Y/(1-Z) ≤ R) for each validated eigenpair", all(c.success for c ∈ validated(e)),
        "r (1-norm on (v, λ), hence on each entry) = " * fmt_certs(e, :r))
    condition("$(n)d. $type", spectrum, "λ ∈ $(e.eig.values)")
end
println("c₁ eigenvectors: v₂ ∈ $(V1[:,2])")
println("                 v₃ ∈ $(V1[:,3])")

# ---- 4.–5. the local stable manifold of c₀ and unstable manifold of c₁ ------
# (Section 3, Proposition 3.3, Theorem 1.1)

λ0_stable = λ0[1:2];
V0_stable = exact(0.5) * V0[:,1:2];
mQ = validate(ManifoldProblem(f, Df, parameters, c0_star, λ0_stable, V0_stable), Q_bar, ν_Q)

λ1_unstable = λ1[2:3];
V1_unstable = exact(0.15) * V1[:,2:3];
mP = validate(ManifoldProblem(f, Df, parameters, c1_star, λ1_unstable, V1_unstable), P_bar, ν_P)

for (n, name, m) ∈ ((4, "Wˢ(c₀) manifold Q", mQ), (5, "Wᵘ(c₁) manifold P", mP))
    condition("$(n)a. $name: Y (Section 3, Step 4) is a finite bound", finite(m.Y), "Y ∈ $(m.Y)")
    condition("$(n)b. $name: Z₀ (Proposition 3.3) is a finite bound", finite(m.Z₀), "Z₀ ∈ $(m.Z₀)")
    condition("$(n)c. $name: Z₁(R) (Proposition 3.3) is a finite bound", finite(m.Z₁), "Z₁ ∈ $(m.Z₁), R = $(m.R)")
    condition("$(n)d. $name: Theorem 1.1 (Z = Z₀ + Z₁ < 1, Y/(1-Z) ≤ R)", m.proof_success,
        "Z ∈ $(m.Z), interval of existence $(m.ie); r = $(m.r)")
end
condition("5e. Wᵘ(c₁) manifold P: P̄ is exactly conjugacy-symmetric, P̄₍ₘ,ₙ₎ = conj(P̄₍ₙ,ₘ₎), and the data P_bar.bin are P̄ up to $(P_data_tolerance)",
    is_swap_conjugate_symmetric(P_bar) && P_data_deviation ≤ P_data_tolerance,
    "max |P_bar.bin - P̄| = $(P_data_deviation) (symmetrised: $(is_swap_conjugate_symmetric(P_bar)), data as read: $(is_swap_conjugate_symmetric(P_data)))")

# ---- 6. the connecting orbit (Section 4, Lemma 4.3, Proposition 4.4, Theorem 1.1)

Q, r_Q = mQ.P, mQ.r
P, r_P = mP.P, mP.r

pb_het = HeteroclinicProblem(f, Df, parameters,
    P, [Derivative(1, 0) .* P Derivative(0, 1) .* P], r_P, ν_P,
    Q, [Derivative(1, 0) .* Q Derivative(0, 1) .* Q], r_Q, ν_Q,
    τ);

het = validate(pb_het, x_bar, μ)

condition("6a. orbit: θ̄₁ ± R, θ̄₂ ± R in the domain (-ν_Q, ν_Q) of Lemma 4.3 (applied on the ball)", het.θ_in_domain, "θ₁ ∈ $(het.θ_ball[1]), θ₂ ∈ $(het.θ_ball[2]), ν_Q = $(ν_Q)")
condition("6b. orbit: Y (Section 4, Step 4) is a finite bound", finite(het.Y), "Y ∈ $(het.Y)")
condition("6c. orbit: Z₀(R) (Proposition 4.4, B(α, θ) over the ball) is a finite bound", finite(het.Z₀), "Z₀(R) ∈ $(het.Z₀)")
condition("6d. orbit: Z₁(R) (Proposition 4.4) is a finite bound", finite(het.Z₁), "Z₁(R) ∈ $(het.Z₁), R = $(het.R)")
condition("6e. orbit: Theorem 1.1 (Z = Z₀(R) + Z₁(R) < 1, Y/(1-Z) ≤ R)", het.proof_success,
    "Z ∈ $(het.Z), interval of existence $(het.ie); r = $(het.r), ᾱ = $(x_bar[end-2])")

# ---- verdict --------------------------------------------------------------

all_pass = all(verdicts)
println()
println(all_pass ? "ALL $(length(verdicts)) CONDITIONS PASS" : "$(count(!, verdicts)) OF $(length(verdicts)) CONDITIONS FAILED")
exit(all_pass ? 0 : 1)
