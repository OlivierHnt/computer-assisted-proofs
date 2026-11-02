# Case II — SEIR model (paper, Section 4.2).
#
#   S' = μ − ν S I − μ S,  I' = ν S(t−τ) I(t−τ) e^{−μτ} − (γ + μ) I,
#   F(x, y) = (μ − ν x₁ x₂ − μ x₁,  ν e^{−μτ} y₁ y₂ − (γ + μ) x₂)
#
# with μ = 0.1, ν = 0.5, γ = 0.1, τ = 1, the disease-free equilibrium
# c = (1, 0), one unstable eigenvalue λ, eigenvector ξ = s (ν, −(λ + μ)).
# Taylor coefficients are real-interval sequences on Taylor(N).

module SEIR

using RadiiPolynomial
include("common.jl")

# parameters: the tightest binary64 intervals around the decimals
const μ, ν, γ, τ = I"0.1", I"0.5", I"0.1", interval(1)
const c = (interval(1), interval(0))                 # equilibrium (S, I) = (1, 0)
const s = I"1.8"                                     # eigenvector scaling, ξ = s (ν, −(λ + μ))
const N = 200                                        # Taylor order (Section 4.2)
const r_star = 1e-3                                  # radius r* of Theorem 3.1 (binary64 nearest to 10⁻³)
const r_λ = 1e-10                                    # radius of the interval on which f″ is bounded (eigenvalue proof)

const S = Taylor(N)

# characteristic matrix at c = (1, 0) (Section 4.2), with νe := ν e^{−μτ}:
#   Ψ(z) = [−μ − z, −ν; 0, −μ − γ + νe e^{−zτ} − z],
# so the characteristic function is f(z) = det Ψ(z) = (z + μ) g(z) with
# g(z) = z + γ + μ − νe e^{−zτ}; f′ = g + (z + μ) g′, f″ = 2 g′ + (z + μ) g″
const νe = ν * exp(-μ * τ)
g(z) = z + γ + μ - νe * exp(-z * τ)
g′(z) = exact(1) + τ * νe * exp(-z * τ)
g″(z) = -τ * τ * νe * exp(-z * τ)
f(z) = (z + μ) * g(z)
f′(z) = g(z) + (z + μ) * g′(z)
f″(z) = exact(2) * g′(z) + (z + μ) * g″(z)

eigenvalue_bounds(λ̄) = eigenvalue_bounds(f, f′, f″, λ̄, r_λ)

# the eigenvector ξ = s (ν, −(λ + μ)) of Section 4.2
eigenvector(λ) = (s * ν, -s * (λ + μ))

# condition II.0d: Ψ(λ) ξ ∋ 0, i.e. ξ is an eigenvector for the enclosed λ
function check_eigenvector(λ)
    ξ = eigenvector(λ)
    Ψξ = ((-μ - λ) * ξ[1] - ν * ξ[2], (-μ - γ + νe * exp(-λ * τ) - λ) * ξ[2])
    return all(containszero, Ψξ)
end

# the delayed coefficients b_k = a_k e^{−kλτ} (eq. (19))
delayed(x, λ, τ) = Sequence(space(x), [x[k] * exp(-exact(k) * λ * τ) for k ∈ indices(space(x))])

# condition II.0e and the coefficients: the homological equations (10).
# A₀ = c, A₁ = ξ (Lemma 1.2), and for 2 ≤ k ≤ N the nonlinear part of F gives
# F̂_k = (−ν \hat{(a¹ * a²)}_k,  νe \hat{(b¹ * b²)}_k) (eq. (17)), and with
# Ψ(z)^{-1} = (1/f(z)) [Ψ₂₂, ν; 0, Ψ₁₁], z = kλ, A_k = −Ψ(z)^{-1} F̂_k reads
#   a¹_k = (Ψ₂₂(z) ν â_k − ν νe b̂_k) / f(z),   a²_k = (z + μ) νe b̂_k / f(z),
# order by order on the enclosures of the coefficients of order < k. By
# inclusion isotonicity the intervals contain the exact Π_N(A) for every λ in
# the enclosure, provided every f(kλ) is nonzero, which is checked.
function solve_homological_equations(λ)
    a¹, a², b¹, b² = (zeros(Interval{Float64}, S) for _ ∈ 1:4)
    ξ = eigenvector(λ)
    a¹[0], a²[0] = c
    a¹[1], a²[1] = ξ
    b¹[0], b²[0] = c
    b¹[1], b²[1] = ξ[1] * exp(-λ * τ), ξ[2] * exp(-λ * τ)

    nonresonant = true
    min_abs_f = Inf
    for k ∈ 2:N
        z = exact(k) * λ
        fz = f(z)
        nonresonant &= !containszero(fz)
        min_abs_f = min(min_abs_f, inf(abs(fz)))
        â, b̂ = hat_product(k, a¹, a²), hat_product(k, b¹, b²)
        a¹[k] = ((-μ - γ + νe * exp(-z * τ) - z) * ν * â - ν * νe * b̂) / fz
        a²[k] = (z + μ) * νe * b̂ / fz
        b¹[k], b²[k] = a¹[k] * exp(-z * τ), a²[k] * exp(-z * τ)
    end
    return (a¹, a²), nonresonant, min_abs_f
end

# the bounds of Lemma 3.3 and the resolvent bound, on Π_N(A)
function bounds(A, λ)
    a¹, a² = A
    b¹, b² = delayed(a¹, λ, τ), delayed(a², λ, τ)
    r = interval(r_star)

    # 1. C_N (Proposition A.1): K₁ = D₁F(c,c) = [−μ −ν; 0 −(γ+μ)], K₂ = D₂F(c,c) = [0 0; 0 ν e^{−μτ}],
    #    ∞-norms ‖K₁‖ = max(μ + ν, γ + μ), ‖K₂‖ = ν e^{−μτ}; ζ = (N+1) λ
    C = resolvent_bound(exact(N + 1) * λ, max(μ + ν, γ + μ), νe, τ)

    # 2. Y (eq. (30)): F̂ = (−ν a¹ * a², ν e^{−μτ} b¹ * b²) beyond order N (full
    #    products, on Taylor(2N)); Y = C_N max over the two components of Σ_{k > N} |F̂_k|
    Y = C * max(ν * tail_norm(a¹ * a², N), νe * tail_norm(b¹ * b², N))

    # 3. Z₁ (eq. (31)): ∂_{x₁}F̂₁ = −ν a², ∂_{x₂}F̂₁ = −ν a¹; ∂_{y₁}F̂₂ = ν e^{−μτ} b², ∂_{y₂}F̂₂ = ν e^{−μτ} b¹
    Z₁ = C * max(ν * (norm(a¹, 1) + norm(a², 1)), νe * (norm(b¹, 1) + norm(b², 1)))

    # 4. Z₂ (eq. (32), see hessian_sum): |F₁| ∋ ν x₁ x₂, |F₂| ∋ ν e^{−μτ} y₁ y₂,
    #    exponents (x₁, x₂, y₁, y₂) = (1, 1, 0, 0) and (0, 0, 1, 1): the sums are 2ν and
    #    2ν e^{−μτ}, independent of z; Z₂ = 2 C_N ν
    z = (norm(a¹, 1) + r, norm(a², 1) + r, norm(b¹, 1) + r, norm(b², 1) + r)
    Z₂ = C * max(hessian_sum([(ν, (1, 1, 0, 0))], z), hessian_sum([(νe, (0, 0, 1, 1))], z))

    return C, Y, Z₁, Z₂
end

end # module
