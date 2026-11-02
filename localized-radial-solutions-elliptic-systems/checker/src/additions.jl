# Checks ADDED for the proof package; not part of the original code.
# (Ported to IntervalArithmetic 1.0: interval/issubset_interval/in_interval/…)
#
# The original code takes the hypotheses below for granted (they are
# established in the paper, Sections 3–4 and Appendix B). They are checked here
# so that "all conditions hold ⇒ theorem" can be read off the checker output.
# None of them changes the arithmetic of the original conditions.

_full_(x) = _full(x)
_isequal(a::Complex{<:Interval}, b::Complex{<:Interval}) = isequal_interval(real(a), real(b)) && isequal_interval(imag(a), imag(b))
_hull(xs) = reduce(hull, xs)

_big(x::Interval) = interval(big(inf(x)), big(sup(x)))
_contains(stored::Interval, exact::Interval) = issubset_interval(_big(exact), _big(stored))
_contains(stored::Complex{<:Interval}, exact::Complex{<:Interval}) = _contains(real(stored), real(exact)) && _contains(imag(stored), imag(exact))
_contains(stored::Complex{<:Interval}, exact::Interval) = _contains(real(stored), exact) && in_interval(0, imag(stored))

_thin(x::Interval, v) = inf(x) == v && sup(x) == v

# ---- hypotheses of Theorem 3.3 on the constants ---------------------------

function check_parameters(data)
    x = data.x
    n_T = order(space(space(component(x, 3)))) + data.n_T_pad
    n_C = order(space(space(component(x, 4)))) + data.n_C_pad
    r_max = exp(-inv(interval(n_T + 1))) # e^{-1/(n_T+1)}
    ok = n_T ≥ 2 && n_C ≥ 1 &&
        0 < data.r_star ≤ inf(r_max) &&
        data.ν > 1 && data.l > 0 && data.L > 0 && data.ℒ_x > 0 && data.ℒ_y > 0 && data.ϱ > 0
    msg = "n_T = $n_T, n_C = $n_C, r_* = $(data.r_star) ≤ e^{-1/(n_T+1)} ∈ $(_full_(r_max)), ν = $(data.ν), ℓ = $(data.l), L = $(data.L), ℒ_x = $(data.ℒ_x), ℒ_y = $(data.ℒ_y), ϱ = $(data.ϱ)"
    return ok, msg
end

# ---- x̄ = S(x̄) (hypothesis of Lemma 3.4 and Corollary 3.5) -----------------

function check_symmetry(data)
    x, Λ, Γ = data.x, coefficients(data.Λ), coefficients(data.Γ)
    q = nspaces(space(component(x, 1)))
    if eltype(x) <: Real
        ok = eltype(Λ) <: Interval{<:Real} && eltype(Γ) <: Interval{<:Real}
        return ok, "x̄, Λ, Γ real: S is complex conjugation, x̄ = S(x̄) trivially"
    end
    # complex case: here q = 2 with one conjugate pair (Swift–Hohenberg)
    η = coefficients(component(x, 1))
    rest = [coefficients(component(x, 2)); coefficients(component(x, 3)); coefficients(component(x, 4))]
    ok_pair = q == 2 && η[1] == conj(η[2]) &&
        _isequal(Λ[1,1], conj(Λ[2,2])) && all(_isequal(Γ[i,1], conj(Γ[i,2])) for i ∈ 1:q)
    ok_real = all(iszero ∘ imag, rest)
    msg = "η̄₁ = conj(η̄₂): $(η[1] == conj(η[2])), Λ₁₁ = conj(Λ₂₂) and Γ₍:,1₎ = conj(Γ₍:,2₎): $(ok_pair), ϕ̄, v̄, w̄ real: $(ok_real) (max |imag| = $(maximum(abs ∘ imag, rest)))"
    return ok_pair && ok_real, msg
end

# ---- equation parameters, equilibrium and eigen-data ------------------------
# The stored intervals must enclose the parameter values of the theorem, an
# exact zero c of N, and an exact solution (Λ, Γ) of -DN(c) Γ = Γ Λ² with
# Re Λᵢᵢ > 0 and Γ invertible (Section 2).

_I(p, q = 1) = interval(p) / interval(q)

function check_equilibrium(eq::KleinGordon, data)
    c, Λ, Γ = coefficients(data.c), coefficients(data.Λ), coefficients(data.Γ)
    ok_par = _contains(eq.β₁, _I(1)) && _contains(eq.β₂, _I(1))
    ok_eig = _thin(c[1], 0.0) && _thin(Λ[1,1], 1.0) && _thin(Γ[1,1], 1.0)
    return ok_par && ok_eig,
        "β₁ ∋ 1, β₂ ∋ 1: $ok_par; c = 0, Λ = Γ = 1 exactly (N(0) = 0, -DN(0) = 1): $ok_eig"
end

# principal square root of a complex interval a + ib with a < 0 (componentwise formula)
function _csqrt(z::Complex{<:Interval})
    a, b = real(z), imag(z)
    m = sqrt(a^2 + b^2)
    re = sqrt((m + a) / 2)
    im_ = sqrt((m - a) / 2)
    sup(b) < 0 && return complex(re, -im_)
    inf(b) > 0 && return complex(re, im_)
    error("branch not determined")
end

# the closed-form values are enclosed in 256-bit BigFloat interval arithmetic
_Ib(p, q = 1) = interval(big(p)) / interval(big(q))

check_equilibrium(eq::SwiftHohenberg, data) = setprecision(() -> _check_equilibrium_sh(eq, data), BigFloat, 256)

function _check_equilibrium_sh(eq, data)
    c, Λ, Γ = coefficients(data.c), coefficients(data.Λ), coefficients(data.Γ)
    β₁, β₂, β₃, β₄ = _Ib(-3, 5), sqrt(_Ib(6)), _Ib(-1), _Ib(1)
    ok_par = _contains(eq.β₁, β₁) && _contains(eq.β₂, β₂) && _contains(eq.β₃, β₃) && _contains(eq.β₄, β₄)
    # c = 0; -DN(0) = [-β₄ 1; β₁ -β₄]; eigenvector (1, g) with g² = β₁, eigenvalue λ² = g - β₄
    g = complex(zero(β₁), -sqrt(-β₁))  # first column; second column is the conjugate
    λ = _csqrt(complex(-β₄, zero(β₁)) + g)
    ok_c = all(cᵢ -> _thin(cᵢ, 0.0), c)
    ok_eig = _contains(Λ[1,1], λ) && _contains(Λ[2,2], conj(λ)) &&
        _contains(Γ[1,1], _Ib(1)) && _contains(Γ[1,2], _Ib(1)) &&
        _contains(Γ[2,1], g) && _contains(Γ[2,2], conj(g)) &&
        inf(real(λ)) > 0 && sup(imag(λ)) < 0
    return ok_par && ok_c && ok_eig,
        "β ∋ (-3/5, √6, -1, 1): $ok_par; c = 0: $ok_c; Λ₁₁ ∋ √(-β₄ - i√(-β₁)) ∈ $(_full_(interval(Float64, real(λ)))) + $(_full_(interval(Float64, imag(λ))))i, Γ ∋ [1 1; -i√(-β₁) i√(-β₁)], conjugate pair: $ok_eig"
end

function check_equilibrium(eq::FitzHughNagumo, data)
    c, Λ, Γ = coefficients(data.c), coefficients(data.Λ), coefficients(data.Γ)
    ϵ, β₁, β₂, β₃, β₄ = _I(3, 10), _I(1, 2), _I(1, 2), _I(1), _I(3)
    ok_par = _contains(eq.ϵ, ϵ) && _contains(eq.β₁, β₁) && _contains(eq.β₂, β₂) && _contains(eq.β₃, β₃) && _contains(eq.β₄, β₄)
    # c = (c_*, c_*, c_*), c_* = -(5 + √145)/20, a root of υ - υ³ - ϵ(β₁ + (β₂+β₃)υ)
    # = -(υ - 1/2)(υ² + υ/2 - 3/10)  (for these parameter values)
    c_star_big = setprecision(() -> -(_Ib(5) + sqrt(_Ib(145))) / _Ib(20), BigFloat, 256)
    ok_c = all(cᵢ -> _contains(cᵢ, c_star_big), c)
    c_star = interval(Float64, c_star_big) # binary64 enclosure for the eigenproblem below (port: rounded outward from the BigFloat one, as IntervalArithmetic 1.0 gives a 1-ulp wider binary64 enclosure)
    # eigenpairs of -DN(c) by the Newton–Kantorovich argument of Appendix B
    M = [inv(ϵ^2)*(1 - 3c_star^2)  -β₂/ϵ  -β₃/ϵ ;
         interval(1.0)              interval(-1.0)  interval(0.0) ;
         inv(β₄^2)                  interval(0.0)   -inv(β₄^2)]
    oks = Bool[]; msgs = String[]
    for j ∈ 1:3
        ok_j, r_j = _eig_nk(M, mid(Λ[j,j]), mid.(Γ[:,j]), 1e-10)
        λ_encl = interval(mid(Λ[j,j])) + interval(-r_j, r_j)
        g_encl = interval.(mid.(Γ[:,j])) .+ interval(-r_j, r_j)
        ok_j = ok_j && _contains(Λ[j,j], λ_encl) && all(_contains.(Γ[:,j], g_encl)) && inf(Λ[j,j]) > 0
        push!(oks, ok_j); push!(msgs, "λ$j ∈ $(_full_(λ_encl)) (r = $r_j)")
    end
    distinct = all(isdisjoint_interval(Λ[i,i], Λ[j,j]) for i ∈ 1:3 for j ∈ i+1:3)
    ok_eig = all(oks) && distinct
    return ok_par && ok_c && ok_eig,
        "(ϵ, β) ∋ (3/10, 1/2, 1/2, 1, 3): $ok_par; c ∋ -(5+√145)/20 ∈ $(_full_(c_star)): $ok_c; eigenpairs (Appendix B): " * join(msgs, ", ") * ", distinct: $distinct"
end

# Newton–Kantorovich for F_eig(λ, g) = (g₁ - 1, M g + λ² g) (Appendix B), in the
# max norm on ℝ × ℝ³. A_eig = inv(mid(DF_eig(x̄))) is a floating-point auxiliary.
function _eig_nk(M, λ̄::Float64, ḡ::Vector{Float64}, ϱe::Float64)
    # port: products written as explicit sums (the generic matrix product of
    # IntervalArithmetic 0.20 used by the faithful version), not the `:fast`
    # BLAS-based product of IntervalArithmetic 1.0, which is wider on 4×4 data
    mv(A, x) = [sum(A[i,k] * x[k] for k ∈ axes(A, 2)) for i ∈ axes(A, 1)]
    mm(A, B) = [sum(A[i,k] * B[k,j] for k ∈ axes(A, 2)) for i ∈ axes(A, 1), j ∈ axes(B, 2)]
    F_eig(λ, g) = [g[1] - 1; mv(M, g) + λ^2*g]
    DF_eig(λ, g) = [interval(0.0) interval(1.0) interval(0.0) interval(0.0); 2λ*g  M + λ^2*interval.(Matrix{Float64}(I, 3, 3))]
    λ₀, g₀ = interval(λ̄), interval.(ḡ)
    A = inv(mid.(DF_eig(λ₀, g₀)))
    Y = maximum(sup.(abs.(mv(A, F_eig(λ₀, g₀)))))
    B = interval(-ϱe, ϱe)
    E = interval.(Matrix{Float64}(I, 4, 4)) - mm(A, DF_eig(λ₀ + B, g₀ .+ B))
    Z = maximum(sup(sum(abs.(E[i,:]))) for i ∈ 1:4)
    Z < 1 || return false, Inf
    r = sup(interval(Y) / (interval(1) - interval(Z)))
    return r ≤ ϱe, r
end

# ---- Proposition 2.2 with a ψ̂ including |c| and |Γ| -------------------------
# (1/2) |Λ⁻¹Γ⁻¹ (DN(c + Γζ') - DN(c)) Γ| with |ζ'| ≤ (1+ℒ_y)μ, so
# |Γζ'| ≤ |Γ|(1+ℒ_y)μ =: s. In each case DN(c+ζ) - DN(c) has one nonzero
# entry, bounded by B(s) below; ψ̂ = B(s)/2 · |Λ⁻¹Γ⁻¹| |Γ|.

_B(eq::KleinGordon, c, s) = (2abs(eq.β₁) + 3abs(eq.β₂)*(2abs(c[1]) + s))*s
_B(eq::SwiftHohenberg, c, s) = (2abs(eq.β₂) + 3abs(eq.β₃)*(2abs(c[1]) + s))*s
_B(eq::FitzHughNagumo, c, s) = 3inv(abs(pow(eq.ϵ, 2)))*(2abs(c[1]) + s)*s

function ψ̂_strict(eq, c, Λ, Γ, μ, ℒ_y)
    nΓ = opnorm(coefficients(Γ), Inf)
    s = nΓ * (1 + ℒ_y) * μ
    return _B(eq, c, s) / 2 * opnorm(coefficients(inv(Γ * Λ)), Inf) * nΓ
end

function check_LP_strict(eq, c, Λ, Γ, λ̂, δ, μ, ℒ_x, ℒ_y)
    ψ̂_ = ψ̂_strict(eq, c, Λ, Γ, μ, ℒ_y)
    a = ((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y)
    c₁ = sup(a) < inf(λ̂)
    c₂ = inf(ℒ_y * (2λ̂ - a)) ≥ sup(a)
    rhs₃ = (inv(2λ̂ - a) + a / ((2λ̂ - (3(eq.d - 1) * δ / 2 + 2ψ̂_) * (1 + ℒ_y)) * (2λ̂ - a))) *
        ((eq.d - 1) * (1 + ℒ_y) / 2 + ((eq.d - 1) * δ / 2 + ψ̂_) * ℒ_x)
    c₃ = inf(ℒ_x) ≥ sup(rhs₃)
    return c₁ && c₂ && c₃, "ψ̂ ∈ $(_full_(ψ̂_)); inequalities: $c₁, $c₂, $c₃ (a ∈ $(_full_(a)), ℒ_x bound ∈ $(_full_(rhs₃)))"
end
