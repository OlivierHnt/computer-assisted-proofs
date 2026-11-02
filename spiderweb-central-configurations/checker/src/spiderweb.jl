# The spiderweb map f of the paper (Section 6, first written in Section 2.4,
# eq. (fiphi)), its Jacobian Df and its second derivative D²f, evaluated in
# interval arithmetic. The formulas are those printed in Section 6.1
# ("Computation of the bounds"); each was rederived and checked against finite
# differences (see the report of the revision and tests/ in the scratch work).
#
# Conventions. Circles i = 1, …, n, radii r_i; l half-lines at angles
# θ_k = 2πk/l, k = 0, …, l-1; masses m_i on circle i, m_0 at the centre; λ < 0.
#
#   f_i(r) = λ r_i + m_i ζ_l / (2^{3/2} r_i²) + m_0 / r_i²
#            + Σ_{j ≠ i} m_j Σ_{k=0}^{l-1} g(r_i, r_j, θ_k),
#
#   g(a, b, θ) = (a − b cos θ) / (a² + b² − 2ab cos θ)^{3/2},
#   ζ_l = Σ_{k=1}^{l-1} (1 − cos θ_k)^{-1/2}.
#
# f_i is a function of r_i and of the other r_j, so its second derivative has
# only three kinds of non-zero entries: ∂²_{ii} f_i, ∂²_{ij} f_i = ∂²_{ji} f_i
# and ∂²_{jj} f_i (j ≠ i); ∂²_{jl} f_i = 0 for j ≠ l both ≠ i.
#
# All arguments are intervals (IntervalArithmetic 1.x); no bare floating-point
# number enters the arithmetic, so every result carries the "guaranteed" flag.

module Spiderweb

using IntervalArithmetic, LinearAlgebra

export f_Df_D2f

# ---- the pair term g and its derivatives (Section 6.1) -----------------------
# a = r_i, b = r_j; c, c2, c3 = cos θ, cos 2θ, cos 3θ; d = (a² + b² − 2ab cos θ)^{1/2}.
# The numerators are those printed in the paper, which agree with the symbolic
# derivatives of g (checked with a computer algebra system).

const two, three, four = interval(2), interval(3), interval(4)

dist(a, b, c) = sqrt(a^2 + b^2 - two*a*b*c)

g(a, b, c)         = (a - b*c) / dist(a, b, c)^3
g_a(a, b, c, c2)   = -(four*a^2 + b^2 - interval(8)*a*b*c + three*b^2*c2) / (two*dist(a, b, c)^5)
g_b(a, b, c, c2)   = -(-four*(a^2 + b^2)*c + a*b*(interval(7) + c2)) / (two*dist(a, b, c)^5)
g_aa(a, b, c, c2)  = three*(a - b*c)*(four*a^2 - b^2 - interval(8)*a*b*c + interval(5)*b^2*c2) / (two*dist(a, b, c)^7)
g_ab(a, b, c, c2, c3) = -three*(a*(interval(8)*a^2 + interval(23)*b^2)*c
                                 - b*(interval(20)*a^2 + two*b^2 + (four*a^2 + interval(6)*b^2)*c2 - a*b*c3)) / (four*dist(a, b, c)^7)
g_bb(a, b, c, c2, c3) = g_ab(b, a, c, c2, c3)  # d is symmetric in (a, b); the paper's third case is the second with i ↔ j

# ---- f, Df, D²f --------------------------------------------------------------

"""
    f_Df_D2f(x, ball, m0, m, λ, n, l)

`x`: the point x̄ (vector of n intervals), `ball`: the box B̄_{ρ*}(x̄) (vector of
n intervals), `m0`, `m`, `λ`: the masses and the common factor (intervals).
Returns `(F, DF, D2F)`: f(x̄) and Df(x̄) (n-vector and n × n matrix), and the
n × n² matrix B = (B_1 | … | B_n) of Section 6.1 enclosing D²f over `ball`,
with B[k, (i-1)n + j] ∋ ∂²_{ji} f_k(b) for every b ∈ `ball`.
"""
function f_Df_D2f(x::Vector{<:Interval}, ball::Vector{<:Interval}, m0::Interval, m::Vector{<:Interval}, λ::Interval, n::Int, l::Int)
    θ   = [two*interval(k)*interval(π)/interval(l) for k ∈ 0:l-1]   # θ_k = 2πk/l
    cosθ, cos2θ, cos3θ = cos.(θ), cos.(two .* θ), cos.(three .* θ)
    ζ   = sum(inv(sqrt(interval(1) - cosθ[k+1])) for k ∈ 1:l-1)    # ζ_l
    isqrt2 = sqrt(two)

    F   = zeros(eltype(x), n)
    DF  = zeros(eltype(x), n, n)
    H_ii    = zeros(eltype(x), n)      # ∂²_{ii} f_i on the ball
    H_mixed = zeros(eltype(x), n, n)   # [i, j] = ∂²_{ij} f_i, j ≠ i
    H_pure  = zeros(eltype(x), n, n)   # [i, j] = ∂²_{jj} f_i, j ≠ i

    for i ∈ 1:n
        # own circle and central mass (Sections 2.4, 6.1)
        F[i]     = λ*x[i] + m[i]*ζ/(two*isqrt2*x[i]^2) + m0/x[i]^2
        DF[i, i] = λ - m[i]*ζ/(isqrt2*x[i]^3) - two*m0/x[i]^3
        H_ii[i]  = three*m[i]*ζ/(isqrt2*ball[i]^4) + interval(6)*m0/ball[i]^4
        # the other circles: Σ_k over the l bodies of circle j
        for j ∈ 1:n
            j == i && continue
            for k ∈ 1:l
                c, c2, c3 = cosθ[k], cos2θ[k], cos3θ[k]
                F[i]     += m[j]*g(x[i], x[j], c)
                DF[i, i] += m[j]*g_a(x[i], x[j], c, c2)
                DF[i, j] += m[j]*g_b(x[i], x[j], c, c2)
                H_ii[i]        += m[j]*g_aa(ball[i], ball[j], c, c2)
                H_mixed[i, j]  += m[j]*g_ab(ball[i], ball[j], c, c2, c3)
                H_pure[i, j]   += m[j]*g_bb(ball[i], ball[j], c, c2, c3)
            end
        end
    end
    return F, DF, hessian_matrix(H_ii, H_mixed, H_pure, n)
end

# Unfold D²f into the n × n² matrix B = (B_1 | … | B_n) of Section 6.1:
# B[k, (i-1)n + j] = ∂²_{ji} f_k. Row k holds the Hessian of f_k, which is
# non-zero only at (k, k), (i, i), (i, k) and (k, i) for i ≠ k.
function hessian_matrix(H_ii, H_mixed, H_pure, n)
    B = zeros(eltype(H_ii), n, n*n)
    col(i, j) = (i - 1)*n + j
    for k ∈ 1:n
        B[k, col(k, k)] = H_ii[k]                 # ∂²_{kk} f_k
        for i ∈ 1:n
            i == k && continue
            B[k, col(i, i)] = H_pure[k, i]        # ∂²_{ii} f_k
            B[k, col(i, k)] = H_mixed[k, i]       # ∂²_{ki} f_k
            B[k, col(k, i)] = H_mixed[k, i]       # ∂²_{ik} f_k
        end
    end
    return B
end

end # module
