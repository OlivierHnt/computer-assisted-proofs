# Compatibility layer for the port to RadiiPolynomial 0.11 / IntervalArithmetic 1.0.

# Printing: IntervalArithmetic 1.0 has no `showfull`; print both endpoints in
# full (shortest round-trip decimal), in the format of IntervalArithmetic 0.20's
# `showfull`, so that outputs of the two versions can be compared verbatim.
_full(x::Interval) = isempty_interval(x) ? "∅" : "Interval($(inf(x)), $(sup(x)))"
_full(z::Complex{<:Interval}) = "$(_full(real(z))) + $(_full(imag(z)))im"
_full(x::AbstractArray) = "[" * join(map(_full, x), ", ") * "]"
_full(x) = string(x)
showfull(x) = print(_full(x))

# `interval_of_existence(Y, Z₁, Z₂, R)` of RadiiPolynomial 0.7.1
# (src/rpa/interval_existence.jl), line by line. RadiiPolynomial 0.11 returns a
# different interval (its right end is min(R, (1 - Z₁)/Z₂) and it works on the
# intervals rather than on their suprema), so the original algorithm is kept:
# with Y, Z₁, Z₂, R replaced by their suprema, it returns I ⊂ [0, R] with
# Y + (Z₁ - 1) r + Z₂ r²/2 ≤ 0 and Z₁ + Z₂ r < 1 for all r ∈ I, or ∅.

_supremum(x::Interval) = sup(x)
_supremum(x::Real) = x

function interval_of_existence_071(Y_::Real, Z₁_::Real, R_::Real)
    Y, Z₁, R = _supremum(Y_), _supremum(Z₁_), _supremum(R_)
    if !(Y ≥ 0 && isfinite(Y) && Z₁ ≥ 0 && isfinite(Z₁) && R ≥ 0)
        return throw(DomainError((Y, Z₁, R), "Y and Z₁ must be positive and finite, R must be positive"))
    elseif Z₁ ≥ 1
        return emptyinterval(Float64)
    else
        r = interval(Y) / (one(Interval{Float64}) - interval(Z₁))
        r_sup = sup(r)
        if 0 ≤ r_sup ≤ R
            return interval(r_sup, R)
        else
            return emptyinterval(Float64)
        end
    end
end

function interval_of_existence_071(Y_::Real, Z₁_::Real, Z₂_::Real, R_::Real)
    Y, Z₁, Z₂, R = _supremum(Y_), _supremum(Z₁_), _supremum(Z₂_), _supremum(R_)
    if Z₂ == 0
        return interval_of_existence_071(Y, Z₁, R)
    elseif !(Y ≥ 0 && isfinite(Y) && Z₁ ≥ 0 && isfinite(Z₁) && Z₂ ≥ 0 && isfinite(Z₂) && R ≥ 0)
        return throw(DomainError((Y, Z₁, Z₂, R), "Y, Z₁ and Z₂ must be positive and finite, R must be positive"))
    else
        b = interval(Z₁) - one(Interval{Float64})
        Δ = b*b - interval(2)*(interval(Z₂)*interval(Y))
        if inf(Δ) < 0 # complex roots
            return emptyinterval(Float64)
        else # real roots
            sqrtΔ = sqrt(Δ)
            r₁ = -(sqrtΔ + b)/interval(Z₂)
            r₁_sup = sup(r₁)
            if 0 ≤ r₁_sup ≤ R && sup(interval(Z₁) + interval(Z₂) * interval(r₁_sup)) < 1
                r₂ = (sqrtΔ - b)/interval(Z₂)
                r₂_inf = inf(r₂)
                if r₁_sup > r₂_inf
                    return emptyinterval(Float64)
                elseif 0 ≤ r₂_inf ≤ R && sup(interval(Z₁) + interval(Z₂) * interval(r₂_inf)) < 1
                    return interval(r₁_sup, r₂_inf)
                elseif sup(interval(Z₁) + interval(Z₂) * interval(R)) < 1
                    return interval(r₁_sup, R)
                else
                    x = prevfloat(inf(-b/interval(Z₂)))
                    if r₁_sup ≤ x ≤ r₂_inf && x ≤ R
                        return interval(r₁_sup, x)
                    else
                        return interval(r₁_sup)
                    end
                end
            else
                return emptyinterval(Float64)
            end
        end
    end
end

# coefficients of a Sequence, or the value itself (evaluation of a sequence on a
# CartesianSpace returns a Sequence or a Vector depending on the version)
_coefs(x::Sequence) = coefficients(x)
_coefs(x::AbstractVector) = x
