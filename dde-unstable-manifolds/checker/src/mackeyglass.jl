# Case III — Mackey–Glass equation, ρ = 3, in its polynomial form (paper,
# Example 1.8, eq. (12), Section 4.3).
#
#   u' = −γ u + β u(t−τ) v(t−τ),  v' = −ρ u^{ρ−1} v² (−γ u + β u(t−τ) v(t−τ)),  v = 1/(1 + u^ρ),
#   F(x, y) = (−γ x₁ + β y₁ y₂,  ργ x₁^ρ x₂² − ρβ x₁^{ρ−1} x₂² y₁ y₂)
#
# with γ = 1, β = 2, τ = 2, equilibrium c = (0, 1), one unstable eigenvalue λ,
# eigenvector ξ = s (1, −ρ (u*)^{ρ−1}) = s (1, 0). Taylor coefficients are
# real-interval sequences on Taylor(N).

module MackeyGlass

using RadiiPolynomial
include("common.jl")

# parameters (exact binary64 numbers)
const γ, β, τ = interval(1), interval(2), interval(2)
const ρ = 3
const c = (interval(0), interval(1))                 # equilibrium (u, v) = (0, 1)
const s = I"0.7"                                     # eigenvector scaling, ξ = s (1, 0)
const N = 600                                        # Taylor order (Section 4.3)
const r_star = 1e-3                                  # radius r* of Theorem 3.1 (binary64 nearest to 10⁻³)
const r_λ = 1e-10                                    # radius of the interval on which f″ is bounded (eigenvalue proof)

const S = Taylor(N)

# characteristic matrix at c = (0, 1): Ψ(z) = diag(−γ + β e^{−zτ} − z, −z)
# (D₁F(c,c) = diag(−γ, 0), D₂F(c,c) = diag(β, 0)), so the characteristic
# function is f(z) = det Ψ(z) = z g(z) with g(z) = z + γ − β e^{−zτ} (eq. (33));
# f′ = g + z g′, f″ = 2 g′ + z g″
g(z) = z + γ - β * exp(-z * τ)
g′(z) = exact(1) + β * τ * exp(-z * τ)
g″(z) = -β * τ * τ * exp(-z * τ)
f(z) = z * g(z)
f′(z) = g(z) + z * g′(z)
f″(z) = exact(2) * g′(z) + z * g″(z)

eigenvalue_bounds(λ̄) = eigenvalue_bounds(f, f′, f″, λ̄, r_λ)

# the eigenvector ξ = s (1, 0) of Section 4.3
eigenvector(λ) = (s, interval(0))

# condition III.0d: Ψ(λ) ξ ∋ 0, i.e. ξ is an eigenvector for the enclosed λ
function check_eigenvector(λ)
    ξ = eigenvector(λ)
    Ψξ = (-g(λ) * ξ[1], -λ * ξ[2])
    return all(containszero, Ψξ)
end

# the delayed coefficients b_k = a_k e^{−kλτ} (eq. (19))
delayed(x, λ, τ) = Sequence(space(x), [x[k] * exp(-exact(k) * λ * τ) for k ∈ indices(space(x))])

# condition III.0e and the coefficients: the homological equations (10).
# A₀ = c, A₁ = ξ (Lemma 1.2), and for 2 ≤ k ≤ N the nonlinear part of F gives
#   F̂_k = (β \hat{(b¹ * b²)}_k,  ρ (γ \hat{((a¹)^ρ * (a²)²)}_k − β \hat{((a¹)^{ρ−1} * (a²)² * b¹ * b²)}_k)),
# and Ψ(z) is diagonal, so with z = kλ, A_k = −Ψ(z)^{-1} F̂_k reads
#   a¹_k = −F̂¹_k / Ψ₁₁(z) = β b̂_k / g(z),   a²_k = F̂²_k / z,
# order by order on the enclosures of the coefficients of order < k. By
# inclusion isotonicity the intervals contain the exact Π_N(A) for every λ in
# the enclosure, provided every f(kλ) = kλ g(kλ) is nonzero, which is checked.
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
        b̂ = hat_product(k, b¹, b²)
        ĝ = hat_product(k, a¹, a¹, a¹, a², a²)             # (a¹)^ρ (a²)², ρ = 3
        ĥ = hat_product(k, a¹, a¹, a², a², b¹, b²)         # (a¹)^{ρ−1} (a²)² b¹ b²
        a¹[k] = β * b̂ / g(z)
        a²[k] = exact(ρ) * (γ * ĝ - β * ĥ) / z
        b¹[k], b²[k] = a¹[k] * exp(-z * τ), a²[k] * exp(-z * τ)
    end
    return (a¹, a²), nonresonant, min_abs_f
end

# the bounds of Lemma 3.3 and the resolvent bound, on Π_N(A)
function bounds(A, λ)
    a¹, a² = A
    b¹, b² = delayed(a¹, λ, τ), delayed(a², λ, τ)
    r = interval(r_star)

    # 1. C_N (Proposition A.1): ‖K₁‖ = γ, ‖K₂‖ = β (∞-norms of the diagonal matrices above); ζ = (N+1) λ
    C = resolvent_bound(exact(N + 1) * λ, γ, β, τ)

    # the Cauchy powers and products of the coefficients (full products)
    a¹², a¹³, a²² = a¹ * a¹, a¹ * a¹ * a¹, a² * a²
    b¹b² = b¹ * b²

    # 2. Y (eq. (30)): beyond order N the hat products are the full products;
    #    F̂¹ = β b¹ * b² (on Taylor(2N)), F̂² = ρ (γ (a¹)³ * (a²)² − β (a¹)² * (a²)² * b¹ * b²) (on Taylor(6N));
    #    Y = C_N max over the two components of Σ_{k > N} |F̂_k|
    F̂¹ = β * b¹b²
    F̂² = exact(ρ) * (γ * (a¹³ * a²²) - β * (a¹² * a²² * b¹b²))
    Y = C * max(tail_norm(F̂¹, N), tail_norm(F̂², N))

    # 3. Z₁ (eq. (31)): the ℓ¹ norms of the partial derivatives of the nonlinear part,
    #    ∂_{y₁}F₁ = β y₂, ∂_{y₂}F₁ = β y₁;
    #    ∂_{x₁}F₂ = ρ (ργ x₁^{ρ−1} x₂² − (ρ−1)β x₁^{ρ−2} x₂² y₁ y₂),  ∂_{x₂}F₂ = 2ρ (γ x₁^ρ x₂ − β x₁^{ρ−1} x₂ y₁ y₂),
    #    ∂_{y₁}F₂ = −ρβ x₁^{ρ−1} x₂² y₂,  ∂_{y₂}F₂ = −ρβ x₁^{ρ−1} x₂² y₁,   at (x, y) = (a, b)
    ∂x₁F₂ = exact(ρ) * (exact(ρ) * γ * (a¹² * a²²) - exact(ρ - 1) * β * (a¹ * a²² * b¹b²))
    ∂x₂F₂ = exact(2ρ) * (γ * (a¹³ * a²) - β * (a¹² * a² * b¹b²))
    ∂y₁F₂ = exact(ρ) * β * (a¹² * a²² * b²)
    ∂y₂F₂ = exact(ρ) * β * (a¹² * a²² * b¹)
    Z₁ = C * max(β * (norm(b², 1) + norm(b¹, 1)),
                 norm(∂x₁F₂, 1) + norm(∂x₂F₂, 1) + norm(∂y₁F₂, 1) + norm(∂y₂F₂, 1))

    # 4. Z₂ (eq. (32), see hessian_sum): |F₁| ∋ β y₁ y₂, exponents (x₁, x₂, y₁, y₂) = (0, 0, 1, 1);
    #    |F₂| = ργ x₁^ρ x₂² + ρβ x₁^{ρ−1} x₂² y₁ y₂, exponents (ρ, 2, 0, 0) and (ρ−1, 2, 1, 1);
    #    z = (‖a¹‖₁ + r*, ‖a²‖₁ + r*, ‖b¹‖₁ + r*, ‖b²‖₁ + r*). The first component gives 2β.
    z = (norm(a¹, 1) + r, norm(a², 1) + r, norm(b¹, 1) + r, norm(b², 1) + r)
    Z₂ = C * max(hessian_sum([(β, (0, 0, 1, 1))], z),
                 hessian_sum([(exact(ρ) * γ, (ρ, 2, 0, 0)), (exact(ρ) * β, (ρ - 1, 2, 1, 1))], z))

    return C, Y, Z₁, Z₂
end

end # module
