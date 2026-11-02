# Connecting orbit (Section 4 of the paper; Lemma 4.3 and Proposition 4.4 in
# the form of the revised manuscript, see statement.md). `validate` takes the
# approximate zero x_bar and returns the Y, Z₀(R), Z₁(R) bounds and the radius
# of Theorem 1.1.
# From ShimizuMoriokaTutorial.jl/src/heteroclinic.jl; Step 2 (ODE solve + Newton) is in search/.

struct HeteroclinicProblem
    f
    Df
    parameters
    P # unstable_manifold
    DP
    r_P
    ν_P
    Q # stable_manifold
    DQ
    r_Q
    ν_Q
    integration_time
end

#-

function validate(pb::HeteroclinicProblem, x_bar::Sequence, μ::Real)
    pb_approx = HeteroclinicProblem(pb.f, pb.Df, mid.(pb.parameters),
        mid.(pb.P), mid.(pb.DP), 0.0, 0.0,
        mid.(pb.Q), mid.(pb.DQ), 0.0, 0.0,
        mid(pb.integration_time))

    # Step 2: x_bar is an input (data/x_bar.bin)

    # Step 3
    A = approximate_inverse(pb_approx, x_bar)

    # Step 4
    X_C = Ell1(GeometricWeight(μ))
    X = NormedCartesianSpace((X_C, X_C, X_C, Ell1(), Ell1(), Ell1()), Ell1())
    #-
    Y = Y_bound(pb, interval(x_bar), A, X)
    R = interval(10sup(Y))
    # Lemma 4.3 is applied to Q⋆ at every θ = (θ₁, θ₂) of the ball B(x̄, R): the
    # intervals θ̄_j ± R must stay in the domain 𝔻_{ν_Q}² ∩ ℝ² of the lemma
    θ_ball = ball(interval(x_bar), R)[2:3]
    θ_in_domain = all(θ -> sup(abs(θ)) < inf(pb.ν_Q), θ_ball)
    Z, Z₀, Z₁ = Z_bound(pb, interval(x_bar), A, R, X_C, X)
    #-
    ie, proof_success = interval_of_existence(Y, Z, R)

    #
    return (; θ_in_domain, θ_ball, Y, Z₀, Z₁, Z, R, ie, proof_success, x = interval.(unpack(x_bar)), r = interval(inf(ie)))
end

# the scalar components (α, θ₁, θ₂) of the ball B(x̄, R) of 𝒳_{C,μ}: since
# ‖x − x̄‖ = Σᵢ‖u⁽ⁱ⁾ − ū⁽ⁱ⁾‖ + |α − ᾱ| + |θ₁ − θ̄₁| + |θ₂ − θ̄₂| ≤ R, each scalar lies
# in the interval centred at its approximation with radius R
function ball(x_bar, R)
    α, θ₁, θ₂ = [unpack(x_bar)[j][1] for j ∈ 4:6]
    return [interval(s, R; format = :midpoint) for s ∈ (α, θ₁, θ₂)]
end





# Step 1

# Rigorous evaluation of a validated parameterization (Lemma 4.3). Given
# ‖P⋆ − P̄‖_{X_{T,ν}} ≤ r and θ = (θ₁, θ₂) in the polydisk 𝔻_ν², real or complex,
#     P⋆(θ) ∈ P̄(θ) + [-r, r]   and   ∂_{θ_j}P⋆(θ) ∈ ∂_{θ_j}P̄(θ) + [-r_j, r_j],
#     r_j = r β_ν(|θ_j|)   (β_ν below),
# the caller passing the radius. For a complex point the radius is added to the
# real and the imaginary part. When θ is a box (Z₀(R) below), the radius is an
# interval and its supremum is used.
error_box(r) = interval(-sup(r), sup(r))

# β_ν of Lemma 4.3: for ν > 1 and 0 ≤ δ < ν, sup_{k ≥ 1} k δ^{k-1} ν^{-k} ≤ β_ν(δ),
#     β_ν(δ) = 1/ν                    for 0 ≤ δ ≤ ν/e,
#     β_ν(δ) = 1/(e δ |ln(δ/ν)|)      for ν/e < δ < ν
# (s ↦ s δ^{s-1} ν^{-s} is maximal at s = 1/|ln(δ/ν)|, which is ≤ 1 iff δ ≤ ν/e;
# the first branch is the value at k = 1). β_ν is constant on [0, ν/e] and
# increasing on [ν/e, ν), so for an interval δ its supremum is β_ν(sup δ).
β(ν, δ::AbstractFloat) = δ ≤ ν / ℯ ? inv(ν) : inv(ℯ * δ * abs(log(δ / ν))) # floating point (approximate inverse only)
β(ν::AbstractFloat, δ::Interval) = β(ν, sup(δ))                                # idem (approximate problem, radius 0)
function β(ν::Interval, δ::Interval)
    threshold = ν / exact(ℯ)
    sup(δ) ≤ inf(threshold) && return inv(ν)
    δ_sup = interval(sup(δ))
    return max(inv(ν), inv(exact(ℯ) * δ_sup * abs(log(δ_sup / ν))))
end
rig_eval(P, θ::NTuple{2,Interval}, r) = Evaluation(θ) .* P .+ error_box(r)
rig_eval(P, θ::NTuple{2,Complex{<:Interval}}, r) = Evaluation(θ) .* P .+ complex(error_box(r), error_box(r))
rig_eval(P, θ::NTuple{2,Union{Float64,ComplexF64}}, r) = Evaluation(θ) .* P # floating point: approximate inverse only (r = 0)
add_radius(v, r::Interval) = v .+ error_box(r)     # real box of radius sup(r) around each entry
add_radius(v, r::AbstractFloat) = v                # floating point: approximate inverse only (r = 0)

# P⋆(γ(α)) and (d/dα) P⋆(γ(α)) are real: P⋆ satisfies P⋆(θ₂*, θ₁*) = P⋆(θ₁, θ₂)*
# (statement.md, "Real-valuedness of P⋆(γ(α))") and γ(α) = (e^{iα}, e^{-iα}) has
# the form (θ, θ*). So the real part of the complex enclosure is a rigorous
# enclosure; the imaginary part is used nowhere.
real_enclosure(z) = real(z)

function F(pb::HeteroclinicProblem, x)
    u = x[1:3]
    α, θ₁, θ₂ = [x[j][1] for j ∈ 4:6]

    γ = (cis(α), cis(-α))
    Pγ = real_enclosure( rig_eval(pb.P, γ, pb.r_P) )
    F₁ = u - Pγ - pb.integration_time / exact(2) * (Integral(1) .* pb.f(u, pb.parameters))

    F₂ = Evaluation(exact(1)) .* u - rig_eval(pb.Q, (θ₁, θ₂), pb.r_Q)

    return [F₁ ; F₂]
end

function DF(pb::HeteroclinicProblem, x)
    u = x[1:3]
    α, θ₁, θ₂ = [x[j][1] for j ∈ 4:6]
    return DF(pb, u, α, θ₁, θ₂)
end

# The derivative at (u, α, θ₁, θ₂): the u-block depends on u only, the α column
# −DP⋆(γ(α))γ'(α) on α only and the θ columns −∂_{θ_j}Q⋆(θ) on θ only. Called
# with u = ū and interval α, θ it is the operator 𝓙(α, θ) of Proposition 4.4.
function DF(pb::HeteroclinicProblem, u, α, θ₁, θ₂)
    DuF₁ = Diagonal([exact(I), exact(I), exact(I)]) - pb.integration_time / exact(2) * (Integral(1) .* Multiplication.(pb.Df(u, pb.parameters)))
    DuF₂ = Diagonal([Evaluation(exact(1)), Evaluation(exact(1)), Evaluation(exact(1))])

    # Lemma 4.3 radii: DP⋆(γ(α))γ'(α) = i e^{iα} ∂_{σ₁}P⋆ − i e^{-iα} ∂_{σ₂}P⋆, with
    # |γ_j(α)| = 1, lies in the real box of radius 2 r_P β_{ν_P}(1) around the real
    # part of DP̄(γ(α))γ'(α) (the exact value is real; the radius is added after
    # taking the real part); r_Q β_{ν_Q}(|θ_j|) for ∂_{θ_j}Q⋆.
    γ = (cis(α), cis(-α))
    Dγ = exact(im) .* (cis(α), -cis(-α))
    r_DP = exact(2) * pb.r_P * β(pb.ν_P, interval(1))
    DP_Dγ_values = add_radius(real_enclosure(sum(i -> (Evaluation(γ) .* pb.DP[:,i]) .* Dγ[i], 1:2)), r_DP)
    DP_Dγ = LinearOperator.(ScalarSpace(), Chebyshev(0), hcat.(DP_Dγ_values))
    r_DQ = [pb.r_Q * β(pb.ν_Q, abs(θ)) for θ ∈ (θ₁, θ₂)]
    DQ_θ = [LinearOperator.(ScalarSpace(), ScalarSpace(), hcat.(rig_eval(pb.DQ[:,1], (θ₁, θ₂), r_DQ[1]))) ;;
            LinearOperator.(ScalarSpace(), ScalarSpace(), hcat.(rig_eval(pb.DQ[:,2], (θ₁, θ₂), r_DQ[2])))]

    return [DuF₁  -DP_Dγ                  fill(zeros(ExactReal{Int}, ScalarSpace(), Chebyshev(0)), 3, 2)
            DuF₂   fill(zeros(ExactReal{Int}, Chebyshev(0), ScalarSpace()), 3, 1)  -DQ_θ]
end





# Step 3

function approximate_inverse(pb_approx::HeteroclinicProblem, x_bar)
    Π = Projection(space(x_bar))
    A_finite = inv(mid(Π * DF(pb_approx, unpack(x_bar)) * Π))
    return interval(A_finite) + (interval(I) - interval(Π))
end





# Step 4

# Y = ‖A F(x̄)‖ with A = A_{≤K}Π_{≤K} + Π_{>K}. The operator product A F(x̄) forms
# the coefficients k ≤ K as A_{≤K}Π_{≤K}F + F_k − F_k, which widens them by
# 2 rad(F_k) in interval arithmetic; they are replaced by A_{≤K}Π_{≤K}F computed
# alone (the coefficients k > K, the identity part, are exact).
function Y_bound(pb::HeteroclinicProblem, x_bar, A, X)
    Fx = F(pb, unpack(x_bar))
    Π = Projection(space(x_bar))
    AF = A * Fx          # (A_{≤K}Π_{≤K} + I − Π_{≤K}) F(x̄): exact on k > K, widened on k ≤ K
    low = (A * Π) * Fx   # A_{≤K} Π_{≤K} F(x̄): the coefficients k ≤ K, exactly
    for i ∈ 1:3, k ∈ 0:order(component(x_bar, i))
        component(AF, i)[k] = component(low, i)[k]
    end
    for i ∈ 4:6
        component(AF, i)[1] = component(low, i)[1]
    end
    return norm(AF, X)
end

# returns Z(R) = Z₀(R) + Z₁(R) and its two terms (Proposition 4.4)
function Z_bound(pb::HeteroclinicProblem, x_bar, A, R, X_C, X)
    K = order(component(x_bar, 1)) # since each component of u_bar have the same truncation order
    Π = Projection(space(x_bar))
    Π_2Kp1 = interval(Projection(Chebyshev(2K+1) × Chebyshev(2K+1) × Chebyshev(2K+1) ×
                                 ScalarSpace() × ScalarSpace() × ScalarSpace()))

    opnorm_A = max(opnorm(A * Π, X), interval(1)) # ‖A‖ = max(‖A_{≤K}‖, 1), eq. (bounds_het)
    Π_0 = interval(Projection(Chebyshev(0) × Chebyshev(0) × Chebyshev(0) × ScalarSpace() × ScalarSpace() × ScalarSpace()))
    opnorm_AΠ0 = opnorm(A * Π_0, X) # ‖A Π_{≤0}‖ of M
    μ = rate(X_C)
    opnorm_L = interval(1) + μ
    first_row_tail_bound_L = inv(μ^interval(K+2) * (interval(K+2)^interval(2) - interval(1)))
    tail_bound_L = inv(μ * interval(2*(K+1))) + μ / (interval(2*(K+3)))

    # Proposition 4.4:
    #   Z₀(R) = max( sup_{x ∈ B(x̄,R)} ‖Π_{≤2K+1} − Π_{≤3K+2} A Π_{≤3K+2} 𝓙(α, θ) Π_{≤2K+1}‖,
    #                M (τ/2) ‖Df(ū)‖ + ‖A‖ μ^{-(2K+2)} ),
    # where 𝓙(α, θ) is DF(x̄) with its α and θ columns evaluated over the ball
    # (the u block stays at ū; its variation is the Z₁(R) term). The supremum is
    # enclosed by evaluating these columns at the intervals ᾱ ± R, θ̄_j ± R. Since
    # 𝓙 Π_{≤2K+1} has range in Π_{≤3K+2} (f quadratic, ū of order K), the finite
    # matrix below is exactly the operator of the formula. On the columns > 2K+1
    # of I − A 𝓙, the u block gives M (τ/2) ‖Df(ū)‖ and the evaluation functionals
    # 𝓔_i give, through A, ‖A‖ μ^{-(2K+2)}; both act on the same columns, so the
    # two bounds add.
    u_bar = unpack(x_bar)[1:3]
    α_ball, θ₁_ball, θ₂_ball = ball(x_bar, R)
    B = DF(pb, u_bar, α_ball, θ₁_ball, θ₂_ball)
    W = pb.Df(u_bar, pb.parameters)
    M = opnorm_AΠ0 * first_row_tail_bound_L + tail_bound_L
    tail_E = opnorm_A * inv(μ^interval(2K+2))
    Z₀ = max(opnorm(Π_2Kp1 - A * (B * Π_2Kp1), X),
             M * pb.integration_time / interval(2) * opnorm(norm.(W, X_C), 1) + tail_E)

    Z₁ = exact(2) * R * pb.integration_time / interval(2) * opnorm_A * opnorm_L

    return Z₀ + Z₁, Z₀, Z₁
end
