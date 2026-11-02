# The proof of one localized radial solution: Proposition 2.2 (Lyapunov–Perron
# contraction), the bounds Y, Z₁, Z₂ and the radii polynomial of Theorem 3.3,
# and the C⁰ error of Corollary 3.5 (arXiv numbering).
#
# Copied from EllipticPDELocalizedRadialSol.jl (src/proof.jl, commit c0fda16),
# ported to RadiiPolynomial 0.11 / IntervalArithmetic 1.0: `Interval(x)` →
# `interval(x)`, `ParameterSpace` → `ScalarSpace`, `isempty` → `isempty_interval`,
# and `interval_of_existence` → `interval_of_existence_071` (the algorithm of
# RadiiPolynomial 0.7.1, reimplemented in compat.jl, since the 0.11 function
# returns a different interval), and in-place broadcasts on Sequence and
# LinearOperator (no longer supported) act on their coefficients. Further
# changes, all I/O and printing only:
#   - `prove(filename)` (JLD2 loading) and `display_data` are replaced by
#     `prove(data, report, label)`, which takes the data already read from
#     data/*.bin (see io.jl);
#   - `__prove` reports each condition through `report.condition(...)` instead
#     of `println`, evaluates every bound even when an earlier condition fails,
#     and records timings;
#   - `isLPcontracting` also returns the three inequalities separately.
# The arithmetic (every expression, its operands and their order) is unchanged;
# F and DF are verbatim.

function prove(data::NamedTuple, report, label::String)
    return _prove(data.equation, data.x, data.n_T_pad, data.n_C_pad, data.c, data.Λ, data.Γ,
        data.l, data.r_star, data.L, data.ℒ_x, data.ℒ_y, data.ν, data.ϱ; report = report, label = label)
end

_prove(eq::Union{KleinGordon{Interval{Float64}},FitzHughNagumo{Interval{Float64}}},
    x::Sequence{<:CartesianProduct,Vector{Float64}}, n_T_pad::Int, n_C_pad::Int,
    c::Sequence{CartesianPower{ScalarSpace},Vector{Interval{Float64}}},
    Λ::LinearOperator{CartesianPower{ScalarSpace},CartesianPower{ScalarSpace},Diagonal{Interval{Float64},Vector{Interval{Float64}}}},
    Γ::LinearOperator{CartesianPower{ScalarSpace},CartesianPower{ScalarSpace},Matrix{Interval{Float64}}},
    l::Float64, r_star::Float64, L::Float64, ℒ_x::Float64, ℒ_y::Float64, ν::Float64, ϱ::Float64; kwargs...) =
        __prove(eq, interval.(x), n_T_pad, n_C_pad, c, Λ, Γ, interval(l), interval(r_star), interval(L), interval(ℒ_x), interval(ℒ_y), interval(ν), interval(ϱ); kwargs...)

_prove(eq::SwiftHohenberg{Interval{Float64}},
    x::Sequence{<:CartesianProduct,Vector{Complex{Float64}}}, n_T_pad::Int, n_C_pad::Int,
    c::Sequence{CartesianPower{ScalarSpace},Vector{Interval{Float64}}},
    Λ::LinearOperator{CartesianPower{ScalarSpace},CartesianPower{ScalarSpace},Diagonal{Complex{Interval{Float64}},Vector{Complex{Interval{Float64}}}}},
    Γ::LinearOperator{CartesianPower{ScalarSpace},CartesianPower{ScalarSpace},Matrix{Complex{Interval{Float64}}}},
    l::Float64, r_star::Float64, L::Float64, ℒ_x::Float64, ℒ_y::Float64, ν::Float64, ϱ::Float64; kwargs...) =
        __prove(eq, interval.(x), n_T_pad, n_C_pad, c, Λ, Γ, interval(l), interval(r_star), interval(L), interval(ℒ_x), interval(ℒ_y), interval(ν), interval(ϱ); kwargs...)

function __prove(eq, x, n_T_pad, n_C_pad, c, Λ, Γ, l, r_star, L, ℒ_x, ℒ_y, ν, ϱ; report, label)
    condition = report.condition
    timings = Pair{String,Float64}[]

    λ̂ = minimum(real(diag(coefficients(Λ))))

    η = component(x, 1)
    v, w = component(x, 3), component(x, 4)

    n_T_num = order(space(space(v)))
    n_T = n_T_num + n_T_pad

    n_C_num = order(space(space(w)))
    n_C = n_C_num + n_C_pad

    δ = inv(l*r_star + L)
    norm_η = norm(η, Inf)
    μ = norm_η + ϱ

    # Proposition 2.2 (original: isLPcontracting)

    LP = isLPcontracting_parts(eq, c, Λ, Γ, λ̂, δ, μ, ℒ_x, ℒ_y)
    condition("$label.1a Proposition 2.2, first inequality: (½(d-1)δ + ψ̂)(1+ℒ_y) < λ̂", LP.condition_1,
        "(½(d-1)δ + ψ̂)(1+ℒ_y) ∈ $(_full(LP.a)), λ̂ ∈ $(_full(λ̂)), with ψ̂ ∈ $(_full(LP.ψ̂)), δ ∈ $(_full(δ)), μ ∈ $(_full(μ))")
    condition("$label.1b Proposition 2.2, ℒ_y inequality: ℒ_y(2λ̂ - a) ≥ a", LP.condition_2,
        "ℒ_y(2λ̂ - a) ∈ $(_full(LP.lhs_2)) ≥ a ∈ $(_full(LP.a))")
    condition("$label.1c Proposition 2.2, ℒ_x inequality", LP.condition_3,
        "ℒ_x = $(_full(ℒ_x)) ≥ $(_full(LP.rhs_3))")
    report.after_LP(LP, (λ̂ = λ̂, δ = δ, μ = μ))

    q = nspaces(space(v))

    #

    t = @elapsed F_ = F(eq, x, c, Λ, Γ, l, r_star, L)
    push!(timings, "F" => t)

    t = @elapsed DF_ = DF(eq, x, n_T, n_C, Λ, Γ, l, r_star, L)
    push!(timings, "DF" => t)

    tail_N = copy(component(F_, 3))
    tail_f = copy(component(F_, 4))
    for i ∈ 1:q
        component(tail_N, i)[0:min(n_T-2, order(component(tail_N, i)))] .= 0
        component(tail_f, i)[0:min(n_C-1, order(component(tail_f, i)))] .= 0
        component(tail_f, i+q+1)[0:min(n_C-1, order(component(tail_f, i+q+1)))] .= 0
    end

    # A: floating-point inverse of the midpoint of the square part of DF(x̄; 0),
    # recomputed here exactly as in the original code (not an input).
    t = @elapsed A = inv(mid.(project(DF_, codomain(DF_), codomain(DF_))))
    push!(timings, "A = inv(mid(DF))" => t)

    #

    X_T = ℓ¹()
    X_T_cartesian = NormedCartesianSpace(X_T, ℓ∞())

    X_C = ℓ¹(GeometricWeight(ν))
    X_C_cartesian = NormedCartesianSpace(X_C, ℓ∞())

    X = NormedCartesianSpace((ℓ∞(), ℓ∞(), X_T_cartesian, X_C_cartesian), ℓ∞())

    #

    t = @elapsed Y = norm(A * F_, X) + max(
        l^2 * norm(tail_N, X_T_cartesian) / interval((n_T + 1) * (n_T + eq.d - 1)),
        L * (ν + inv(ν)) * norm(tail_f, X_C_cartesian) / (4(n_C + 1))
        )
    push!(timings, "Y" => t)
    condition("$label.2 Y bound (Theorem 3.3) is finite", isfinite(sup(Y)), "Y ∈ $(_full(Y))")

    n_C_max = Int((order(space(space(F_)[4]))-1)/n_C_num)*n_C + 1
    t = @elapsed opnorm_A = opnorm(interval.(A), X)
    push!(timings, "‖A‖" => t)
    opnorm_Df_ = opnorm_Df(eq, w, X_C)
    t₁ = @elapsed ADF_I = rsub!(A * DF_, I) # port: `A * DF_ - I` is a lazy operator in RadiiPolynomial 0.11
    push!(timings, "A·DF - I" => t₁)
    t = @elapsed Z₁ = opnorm(ADF_I, X) +
            max(
                l^2 * opnorm_DN(eq, v, X_T) / interval((n_T + 1) * (n_T + eq.d - 1)),
                L * (ν + inv(ν)) * opnorm_Df_ / (4(n_C + 1))
            ) + opnorm_A *
            max(
                2 / ν^(n_C_max + 1), # NOTE: 1 / ν^(n_C_max + 1) is also valid
                r_star^(n_T + 1) * max(1, (n_T + 1) / l) + L * opnorm_Df_ / (2 * ν^(n_C + 2) * ((n_C + 2)^2 - 1))
            )
    ADF_I = nothing
    push!(timings, "Z₁ (rest)" => t)
    condition("$label.3 Z₁ bound (Theorem 3.3) is finite", isfinite(sup(Z₁)), "Z₁ ∈ $(_full(Z₁))")

    Z₂ = (opnorm_A + 1) *
        max(
            l^2 * opnorm_D²N_abs(eq, v, X_T, ϱ),
            L * (1 + ν) * opnorm_D²f_abs(eq, w, X_C, ϱ) / 2
            )
    condition("$label.4 Z₂ bound (Theorem 3.3) is finite", isfinite(sup(Z₂)), "Z₂ ∈ $(_full(Z₂))")

    Q = LinearOperator(ScalarSpace()^q, ScalarSpace()^q × ScalarSpace()^q, [coefficients(Γ) ; coefficients(Γ*Λ)])
    opnorm_AQ = opnorm(component(A, :, 1:2) * Q, ℓ∞(), X)

    ρ = interval_of_existence_071(Y + opnorm_AQ*ℒ_y*norm_η, Z₁ + opnorm_AQ*ℒ_y, Z₂, ϱ)

    condition("$label.5 radii polynomial (Theorem 3.3): interval of existence in [0, ϱ] is nonempty", !isempty_interval(ρ),
        "|A(Γ, ΓΛ, 0, 0)| ∈ $(_full(opnorm_AQ)); interval of existence $(_full(ρ)), ϱ = $(sup(ϱ))")

    error = isempty_interval(ρ) ? Inf : max(inf(ρ), inf(opnorm(coefficients(Γ), Inf)*(ρ + (1 + ℒ_y)*μ)))
    error_sup = isempty_interval(ρ) ? Inf : max(inf(ρ), sup(opnorm(coefficients(Γ), Inf)*(interval(inf(ρ)) + (1 + ℒ_y)*μ)))

    return (; Y, Z₁, Z₂, ρ, error, error_sup, opnorm_AQ, timings, n_T, n_C)
end

#

function F(eq, x, c, Λ, Γ, l, r_star, L)
    η, ϕ = component(x, 1), component(x, 2)
    v, w = component(x, 3), component(x, 4)
    q = nspaces(space(v))

    T = eltype(x)

    N_ = N(eq, v)
    f_ = f(eq, w)

    space_F₃ = Taylor(order(space(space(N_))) + 2)^q
    space_F₄ = Chebyshev(order(space(space(f_))) + 1)^(1 + 2q)
    space_F = ScalarSpace()^q × ScalarSpace()^q × space_F₃ × space_F₄
    F_ = zeros(T, space_F)

    n_T = order(space(space(v)))
    n_T_max = order(space(space_F₃))

    #

    # port: Sequence and LinearOperator are no longer broadcastable in RadiiPolynomial 0.11;
    # the same elementwise operations act on their coefficients
    coefficients(component(F_, 1)) .= _coefs(component(w, 2:(1+q))(one(T))) .- coefficients(c) .- coefficients(Γ * η)
    coefficients(component(F_, 2)) .= _coefs(component(w, (2+q):(1+2q))(one(T))) .+ coefficients(Γ * Λ * η)

    #

    F₃ = component(F_, 3)
    l² = pow(l, 2)
    for i ∈ 1:q
        component(F₃, i)[0] = component(v, i)[0] - ϕ[i]
        component(F₃, i)[1] = component(v, i)[1]
        for n ∈ 2:n_T_max
            component(F₃, i)[n] = l² * component(N_, i)[n-2]
            if n ≤ n_T
                component(F₃, i)[n] += n * (n + eq.d - 2) * component(v, i)[n]
            end
        end
    end

    #

    F₄ = component(F_, 4)
    coefficients(integrate!(F₄, f_)) .*= L/2
    rsub!(F₄, w)
    component(F₄, 1)[0] += inv(l*r_star)
    for i ∈ 1:q, n ∈ 0:n_T
        component(F₄, i+1)[0] += component(v, i)[n] * pow(r_star, n)
        component(F₄, i+q+1)[0] += n * component(v, i)[n] * pow(r_star, n) / l
    end

    #

    return F_
end

function DF(eq, x, n_T, n_C, Λ, Γ, l, r_star, L)
    v, w = component(x, 3), component(x, 4)
    q = nspaces(space(v))

    T = eltype(x)

    DN_ = DN(eq, v, n_T)
    Df_ = Df(eq, w, n_C)

    domain_DF = ScalarSpace()^q × ScalarSpace()^q × Taylor(n_T)^q × domain(Df_)
    space_x_padded = ScalarSpace()^q × ScalarSpace()^q × Taylor(n_T)^q × Chebyshev(n_C)^(1+2q)
    DF_ = zeros(T, domain_DF, space_x_padded)

    #

    coefficients(component(DF_, 1, 1)) .= (-).(coefficients(Γ))
    mul!(component(DF_, 2, 1), Γ, Λ)
    D₄F₁ = component(DF_, 1, 4)
    D₄F₂ = component(DF_, 2, 4)
    for i ∈ 1:q
        coefficients(component(D₄F₁, i, i+1)) .= 2
        component(D₄F₁, i, i+1)[1,0] = 1
        coefficients(component(D₄F₂, i, i+q+1)) .= 2
        component(D₄F₂, i, i+q+1)[1,0] = 1
    end

    #

    D₂F₃ = component(DF_, 3, 2)
    D₃F₃ = component(DF_, 3, 3)
    l² = pow(l, 2)
    for j ∈ 1:q, i ∈ 1:q
        if i == j
            component(D₂F₃, i, i)[0,1] = -1
            component(D₃F₃, i, i)[0,0] = component(D₃F₃, i, i)[1,1] = 1
            for m ∈ 0:n_T, n ∈ 2:n_T
                component(D₃F₃, i, i)[n,m] = l² * component(DN_, i, i)[n-2,m]
                if n == m
                    component(D₃F₃, i, i)[n,n] += n * (n + eq.d - 2)
                end
            end
        else
            for m ∈ 0:n_T, n ∈ 2:n_T
                component(D₃F₃, i, j)[n,m] = l² * component(DN_, i, j)[n-2,m]
            end
        end
    end

    #

    rsub!(mul!(component(DF_, 4, 4), Integral(1), Df_, L/2, false), I)
    D₃F₄ = component(DF_, 4, 3)
    for i ∈ 1:q, n ∈ 0:n_T
        component(D₃F₄, i+1, i)[0,n] = pow(r_star, n)
        component(D₃F₄, i+q+1, i)[0,n] = n * pow(r_star, n) / l
    end

    #

    return DF_
end


# Original isLPcontracting, returning the three inequalities and the
# quantities they compare (same expressions, same order).

function isLPcontracting_parts(eq, c, Λ, Γ, λ̂, δ, μ, ℒ_x, ℒ_y)
    ψ̂_ = ψ̂(eq, c, Λ, Γ, μ, ℒ_y)

    condition_1 = sup(((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y)) < inf(λ̂)

    condition_2 = inf(ℒ_y * (2λ̂ - ((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y))) ≥ sup(((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y))

    rhs_3 = (inv(2λ̂ - ((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y)) +
        ((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y) / ((2λ̂ - (3(eq.d - 1) * δ / 2 + 2ψ̂_) * (1 + ℒ_y)) * (2λ̂ - ((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y)))) *
        ((eq.d - 1) * (1 + ℒ_y) / 2 + ((eq.d - 1) * δ / 2 + ψ̂_) * ℒ_x)

    condition_3 = inf(ℒ_x) ≥ sup(rhs_3)

    a = ((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y)
    lhs_2 = ℒ_y * (2λ̂ - ((eq.d - 1) * δ / 2 + ψ̂_) * (1 + ℒ_y))

    return (; condition_1, condition_2, condition_3, ψ̂ = ψ̂_, a, lhs_2, rhs_3)
end
