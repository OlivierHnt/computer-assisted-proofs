# Radii-polynomial check of one case (n, l): Theorem 6.1 with Y₀, Z₀ of
# Section 6.1 and Z₂ of Lemma 6.3, in interval arithmetic; the polynomial is
# handled by RadiiPolynomial.interval_of_existence.

"""
    check_case(n, l, x̄; λ, m0, m, r_star)

Checks the conditions of statement.md for the numerical zero `x̄` of the
spiderweb map with `n` circles and `l` half-lines. Returns a named tuple with
- `in_Rn`: the ball B̄_{ρ*}(x̄) lies in R^n = {0 < r_1 < … < r_n} (so f is C² on it),
- `Y0`, `Z0`, `Z2`: the interval enclosures of the bounds of Theorem 6.1 / Lemma 6.3,
- `ie`, `success`: the interval of existence and the verdict of `interval_of_existence`,
- `gap`: min(inf ball₁, min_i(inf ball_{i+1} − sup ball_i)), positive iff `in_Rn`.
"""
function check_case(n::Int, l::Int, x̄::Vector{Float64};
                    λ::Float64 = -1.0, m0::Float64 = 0.0, m::Vector{Float64} = fill(1/l, n),
                    r_star::Float64 = 1e-5)
    length(x̄) == n == length(m) || error("x̄ and m must have length n")
    x    = interval.(x̄)
    ball = x .+ interval(-r_star, r_star)   # B̄_{ρ*}(x̄), outward rounded

    # Theorem 6.1 needs B̄_{ρ*}(x̄) ⊂ U with f ∈ C¹(U): here U = R^n, where f is C^∞
    gap   = min(inf(ball[1]), minimum(inf(ball[i+1]) - sup(ball[i]) for i ∈ 1:n-1; init = Inf))
    in_Rn = gap > 0

    F, DF, D2F = f_Df_D2f(x, ball, interval(m0), interval.(m), interval(λ), n, l)

    # A = numerical inverse of Df(x̄) (Section 6.1), as a thin interval matrix. A
    # degenerate x̄ (singular or non-finite Df) is a failed case, not an error.
    DF_mid = mid.(DF)
    DF_lu  = all(isfinite, DF_mid) ? lu(DF_mid; check = false) : nothing
    if DF_lu === nothing || !issuccess(DF_lu)
        return (; in_Rn, gap, Y0 = nai(), Z0 = nai(), Z2 = nai(), ie = emptyinterval(), success = false)
    end
    A = interval.(inv(DF_lu))

    Y0 = norm(A*F, Inf)                                 # ‖A f(x̄)‖_∞
    Z0 = opnorm(interval.(Matrix(I, n, n)) - A*DF, Inf) # ‖I − A Df(x̄)‖_∞
    Z2 = opnorm(A*D2F, Inf)                             # ‖A B(B̄_{ρ*}(x̄))‖_∞, Lemma 6.3

    # Radii polynomial. Z₂ bounds ‖A D²f‖ on the whole ball, so Taylor's formula
    # gives the polynomial Y₀ + (Z₀ − 1)ρ + Z₂ρ²/2 of RadiiPolynomial (the paper's
    # Theorem 6.1 has no 1/2; see statement.md, condition 2). The suprema of the
    # enclosures are passed. The library returns an interval I ⊂ [0, ρ*] on which
    # the polynomial is ≤ 0 and Z₀ + Z₂ρ < 1, and `true`; or the empty interval
    # and `false`.
    ie, success = interval_of_existence(interval(sup(Y0)), interval(sup(Z0)), interval(sup(Z2)), r_star)
    success &= isguaranteed(ie)

    return (; in_Rn, gap, Y0, Z0, Z2, ie, success)
end
