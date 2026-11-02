# Shared pieces of the three cases, derived from the paper (Hénot, Lessard,
# Mireles James, JDDE 2022, revised version; numbering of that version).
#
# Taylor coefficients are RadiiPolynomial `Sequence`s on `Taylor(N)` (one
# parameter) or `Taylor(N) ⊗ Taylor(N)` (two parameters), with interval
# entries; `*` is the Cauchy product (eq. (14)), computed on the full product
# space, and `norm(x, 1)` is the unweighted ℓ¹ norm Σ_α |x_α| of Section 3.1.
# Included inside each case module.

# |α| for a one- or two-dimensional multi-index
order(α::Integer) = α
order(α::Tuple) = sum(α)

# Σ_{|α| > N} |x_α|: the ℓ¹ norm of π_∞(x) (Section 3.2), for x on a Taylor
# space or a tensor product of Taylor spaces
function tail_norm(x::Sequence, N::Integer)
    s = abs(zero(eltype(x)))
    for α ∈ indices(space(x))
        order(α) > N && (s += abs(x[α]))
    end
    return s
end

# the hat product (eq. (15)) at the multi-index α: the α-th coefficient of the
# Cauchy product of the factors with their α-th coefficients set to zero, so it
# only involves coefficients of order < |α|. Coefficients outside the box
# {β ≤ α} do not contribute to the α-th coefficient, hence the factors and the
# intermediate products are restricted to the box.
box(α::Integer) = Taylor(α)
box(α::Tuple{Integer,Integer}) = Taylor(α[1]) ⊗ Taylor(α[2])
function hat_product(α, factors::Sequence...)
    s = box(α)
    lower = map(factors) do x
        y = project(x, s)
        y[α] = zero(eltype(y))
        return y
    end
    p = foldl((u, v) -> project(u * v, s), lower)
    return p[α]
end

# ---- eigenvalue enclosure ----------------------------------------------------
#
# The unstable eigenvalue is a zero of the characteristic function
# f(λ) = det Ψ(λ), Ψ(λ) = λI − K₁ − e^{−λτ}K₂ (eq. (2)); each case module gives
# f, f′, f″ in closed form. The enclosure comes from Theorem 3.1 in dimension
# one (the "rigorous Newton method" of the paper, ref. [38]): with λ̄ the frozen
# approximate eigenvalue and A the binary64 number nearest to 1/f′(λ̄) (A ≠ 0),
# the zeros of f are the fixed points of T(λ) = λ − A f(λ). For |λ − λ̄| ≤ r ≤ r_λ,
# h = λ − λ̄,
#   |T(λ̄) − λ̄| = |A f(λ̄)| =: Y,
#   |T′(λ)| = |1 − A f′(λ)| ≤ |1 − A f′(λ̄)| + |A| |f′(λ) − f′(λ̄)| ≤ Z₁ + Z₂ r,
#        Z₁ := |1 − A f′(λ̄)|,   Z₂ := |A| sup_{|λ−λ̄| ≤ r_λ} |f″(λ)|,
# (f′(λ) − f′(λ̄) = ∫₀¹ f″(λ̄ + t h) h dt), hence
#   |T(λ) − λ̄| ≤ Y + Z₁ r + Z₂ r²/2.
# So T maps the closed disc D_r(λ̄) into itself when p(r) = Z₂ r²/2 + (Z₁ − 1) r + Y ≤ 0
# and is a contraction on it when Z₁ + Z₂ r < 1: f has exactly one zero in
# D_r(λ̄). interval_of_existence(Y, Z₁, Z₂, r_λ) certifies both inequalities on
# an interval of radii; ε is its infimum. For real λ (cases II, III) the same
# runs in real arithmetic on the segment [λ̄ − r, λ̄ + r].
#
# The sup of |f″| is bounded by the interval evaluation of f″ on `ball(λ̄, r_λ)`,
# a real interval, or a complex interval (a square) containing the disc.
ball(x::Interval, r::Real) = x + interval(-r, r)
ball(z::Complex{<:Interval}, r::Real) = complex(ball(real(z), r), ball(imag(z), r))

# the binary64 number nearest to the midpoint of an interval, as a thin interval
thin(x::Interval) = interval(mid(x))
thin(z::Complex{<:Interval}) = complex(thin(real(z)), thin(imag(z)))

function eigenvalue_bounds(f, f′, f″, λ̄, r_λ)
    A = thin(inv(f′(λ̄)))
    Y = abs(A * f(λ̄))
    Z₁ = abs(exact(1) - A * f′(λ̄))
    Z₂ = abs(A) * abs(f″(ball(λ̄, r_λ)))        # its sup bounds |A| sup |f″| on the disc
    return Y, Z₁, Z₂
end

# ---- resolvent and second-derivative bounds ----------------------------------

# Proposition A.1: for z with |z| − ‖K₁‖ − e^{−τ Re z}‖K₂‖ > 0, Ψ(z) is
# invertible and ‖Ψ(z)^{-1}‖ ≤ 1/(|z| − ‖K₁‖ − e^{−τ Re z}‖K₂‖). For
# z = ⟨α,Λ⟩ with |α| ≥ N+1 one has |z| ≥ Re z ≥ ζ := (N+1) min_j Re λ_j, and the
# denominator increases with |z| and with Re z, so
#   C_N = sup_{|α| ≥ N+1} ‖Ψ(⟨α,Λ⟩)^{-1}‖ ≤ C := 1/(ζ − ‖K₁‖ − e^{−τζ}‖K₂‖)
# as soon as the denominator at ζ is positive (then hypotheses (34), (35) hold
# for all these z). The caller checks inf(C) > 0 and C bounded, which is
# equivalent to the denominator being positive.
resolvent_bound(ζ, normK₁, normK₂, τ) = inv(ζ - normK₁ - normK₂ * exp(-ζ * τ))

# Z₂ bound (Lemma 3.3, proof): with U = (Π_N A + H, τ(Π_N A + H)), H ∈ B_{r*}(0),
#   D²F̂_j(U)[(V,τV),(W,τW)] = Σ_{p,q} ∂_p∂_q F_j(U) * (V,τV)_p * (W,τW)_q
# summed over all ORDERED pairs (p, q) of the 2d variables (x₁..x_d, y₁..y_d),
# so by (29) and ‖τ‖ = 1
#   ‖D²F̂_j(U)‖ ≤ Σ_{p,q} ∂_p∂_q |F_j|(z),   z = (‖π_N a^i‖₁ + r*, ‖π_N τ(a^i)‖₁ + r*)_i.
# Each mixed pair {x_i, y_j} therefore appears twice, as ∂_{x_i}∂_{y_j} and as
# ∂_{y_j}∂_{x_i}; eq. (32) of the paper prints it once (see statement.md,
# "Notes on the paper"). For a monomial
# |c| x^κ y^γ with exponent vector e = (κ, γ) ∈ N^{2d},
#   ∂_p∂_q z^e = e_p (e_q − δ_{pq}) z^{e − e_p − e_q}.
# `monomials` lists the monomials of |F_j| as (|c|, e), `z` the 2d-vector above.
function hessian_sum(monomials, z)
    s = zero(eltype(z))
    for (c, e) ∈ monomials, p ∈ eachindex(e), q ∈ eachindex(e)
        n = e[p] * (e[q] - (p == q))
        n == 0 && continue
        f = collect(e)
        f[p] -= 1
        f[q] -= 1
        s += c * exact(n) * prod(pown(z[i], f[i]) for i ∈ eachindex(f))   # pown: the guaranteed integer power
    end
    return s
end

# the interval [0, 0] test for the components of a complex interval
containszero(x::Interval) = in_interval(0, x)
containszero(z::Complex{<:Interval}) = containszero(real(z)) && containszero(imag(z))
