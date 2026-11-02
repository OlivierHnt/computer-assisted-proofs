# Case I — delayed van der Pol system (paper, Sections 2.5, 4.1, A.1).
#
#   u' = v,  v' = −a u − b u(t−τ) + μ (1 − u(t−τ)²) v,
#   F(x, y) = (x₂, −a x₁ − b y₁ + μ (1 − y₁²) x₂)                        (22)
#
# with a = 0.9, b = 0.1, μ = 1, τ = 2, equilibrium c = 0, two unstable
# eigenvalues Λ = (λ, λ̄), eigenvectors ξ = s (1, λ), ξ̄. Taylor coefficients are
# sequences on Taylor(N) ⊗ Taylor(N), complex intervals, with
# A_{(α₂,α₁)} = conj(A_{(α₁,α₂)}); the entries with |α| > N of the square
# storage are zero, so a sequence is Π_N(A).

module VanDerPol

using RadiiPolynomial
include("common.jl")

# parameters: the tightest binary64 intervals around the decimals
const a, b, μ, τ = I"0.9", I"0.1", interval(1), interval(2)
const c = (interval(0), interval(0))                 # equilibrium (u, v) = (0, 0)
const s = interval(2)                                # eigenvector scaling, ξ = s (1, λ); Section 4.1 does not print s
const N = 99                                         # Taylor order (Section 4.1)
const r_star = 1e-3                                  # radius r* of Theorem 3.1 (binary64 nearest to 10⁻³)
const r_λ = 1e-10                                    # radius of the disc on which f″ is bounded (eigenvalue proof)

const S = Taylor(N) ⊗ Taylor(N)

# ⟨α,Λ⟩ = α₁ λ + α₂ λ̄
Λ(α, λ) = exact(α[1]) * λ + exact(α[2]) * conj(λ)

# characteristic matrix (23), Ψ(z) = [−z, 1; −a − b e^{−zτ}, μ − z], and its
# determinant f(z) = det Ψ(z) = z² − μ z + a + b e^{−zτ} (the characteristic
# function), with its first two derivatives
f(z) = z * z - μ * z + a + b * exp(-z * τ)
f′(z) = exact(2) * z - μ - b * τ * exp(-z * τ)
f″(z) = exact(2) + b * τ * τ * exp(-z * τ)

eigenvalue_bounds(λ̄) = eigenvalue_bounds(f, f′, f″, λ̄, r_λ)

# the eigenvector ξ = s (1, λ) of Section 4.1 (with the complex scaling s + 0i)
eigenvector(λ) = (complex(s, interval(0)), s * λ)

# condition I.0d: Ψ(λ) ξ ∋ 0, i.e. ξ is an eigenvector for the enclosed λ
function check_eigenvector(λ)
    ξ = eigenvector(λ)
    Ψξ = (-λ * ξ[1] + ξ[2], -(a + b * exp(-λ * τ)) * ξ[1] + (μ - λ) * ξ[2])
    return all(containszero, Ψξ)
end

# the delayed coefficients b¹_α = a¹_α e^{−⟨α,Λ⟩τ} (eq. (19))
delayed(a¹, λ, τ) = Sequence(space(a¹), [a¹[α] * exp(-Λ(α, λ) * τ) for α ∈ indices(space(a¹))])

# condition I.0e and the coefficients: the homological equations (10), (25).
# A₀ = c, A_{e₁} = ξ, A_{e₂} = ξ̄ (Lemma 1.2), and for 2 ≤ |α| ≤ N, with
# F̂_α = (0, −μ \hat{(a² * b¹ * b¹)}_α) (eq. (24)) and the inverse of Ψ(z),
#   A_α = −Ψ(⟨α,Λ⟩)^{-1} F̂_α = −(μ / f(z)) ĥ_α (1, z),   ĥ = \hat{(b¹ * b¹ * a²)}, z = ⟨α,Λ⟩,
# order by order on the enclosures of the coefficients of order < |α|. By
# inclusion isotonicity the intervals contain the exact Π_N(A) for every λ in
# the enclosure, provided every f(⟨α,Λ⟩) is nonzero, which is checked (its
# interval must exclude zero).
#
# Conjugate symmetry: the parameters are real and Λ = (λ, λ̄), so the equation
# at (α₂, α₁) is the complex conjugate of the one at (α₁, α₂) and
# A_{(α₂,α₁)} = conj(A_{(α₁,α₂)}) exactly. The coefficients are computed for
# α₁ ≤ α₂ and the mirror entry is set to the conjugate enclosure (conjugation
# of a complex interval is exact).
function solve_homological_equations(λ)
    a¹, a², b¹ = (zeros(Complex{Interval{Float64}}, S) for _ ∈ 1:3)
    ξ = eigenvector(λ)
    a¹[(0, 0)], a²[(0, 0)] = complex(c[1], interval(0)), complex(c[2], interval(0))
    a¹[(1, 0)], a²[(1, 0)] = ξ
    a¹[(0, 1)], a²[(0, 1)] = conj.(ξ)
    b¹[(0, 0)] = a¹[(0, 0)]
    b¹[(1, 0)] = a¹[(1, 0)] * exp(-λ * τ)
    b¹[(0, 1)] = conj(b¹[(1, 0)])

    nonresonant = true
    min_abs_f = Inf
    for n ∈ 2:N, α₁ ∈ 0:n÷2
        α = (α₁, n - α₁)
        z = Λ(α, λ)
        fz = f(z)
        nonresonant &= !containszero(fz)
        min_abs_f = min(min_abs_f, inf(abs(fz)))
        ĥ = hat_product(α, b¹, b¹, a²)
        a¹[α] = -μ * inv(fz) * ĥ
        a²[α] = -μ * inv(fz) * z * ĥ
        b¹[α] = a¹[α] * exp(-z * τ)
        mirror = (α[2], α[1])
        a¹[mirror], a²[mirror], b¹[mirror] = conj(a¹[α]), conj(a²[α]), conj(b¹[α])
    end
    return (a¹, a²), nonresonant, min_abs_f
end

# the bounds of Lemma 3.3 and the resolvent bound, on Π_N(A)
function bounds(A, λ)
    a², b¹ = A[2], delayed(A[1], λ, τ)
    r = interval(r_star)

    # 1. C_N (Proposition A.1): K₁ = D₁F(0,0) = [0 1; −a μ], K₂ = D₂F(0,0) = [0 0; −b 0],
    #    ∞-norms ‖K₁‖ = max(1, |a| + |μ|), ‖K₂‖ = |b|; Re λ = Re λ̄ so ζ = (N+1) Re λ
    C = resolvent_bound(exact(N + 1) * real(λ), max(interval(1), abs(a) + abs(μ)), abs(b), τ)

    # 2. Y (eq. (30)): the only nonlinear monomial is −μ y₁² x₂, so
    #    Y = C_N Σ_{|α| > N} |μ (b¹ * b¹ * a²)_α| (for |α| > N the hat product is
    #    the full product, as Π_N(A)_α = 0); the product lives on Taylor(3N)⊗Taylor(3N),
    #    so the sum runs over all N+1 ≤ |α| ≤ 3N
    Y = C * abs(μ) * tail_norm(b¹ * b¹ * a², N)

    # 3. Z₁ (eq. (31)): ∂_{x₂}F̂₂ = −μ b¹ * b¹, ∂_{y₁}F̂₂ = −2μ b¹ * a², the other
    #    partial derivatives of the nonlinear part vanish
    Z₁ = C * abs(μ) * (norm(b¹ * b¹, 1) + exact(2) * norm(b¹ * a², 1))

    # 4. Z₂ (eq. (32), see hessian_sum): |F₂| has the nonlinear monomial μ y₁² x₂,
    #    exponents (x₁, x₂, y₁, y₂) = (0, 1, 2, 0); z = (‖a¹‖₁ + r*, ‖a²‖₁ + r*, ‖b¹‖₁ + r*, ‖b²‖₁ + r*)
    #    (b² does not occur). The sum is 2 x₂ + 4 y₁ at z.
    z = (norm(A[1], 1) + r, norm(a², 1) + r, norm(b¹, 1) + r, norm(delayed(A[2], λ, τ), 1) + r)
    Z₂ = C * hessian_sum([(abs(μ), (0, 1, 2, 0))], z)

    return C, Y, Z₁, Z₂
end

end # module
