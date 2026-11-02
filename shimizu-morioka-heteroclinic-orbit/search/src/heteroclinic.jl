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

function validated_solve(pb::HeteroclinicProblem, α_init::Float64, order::Int, μ::Real)
    pb_approx = HeteroclinicProblem(pb.f, pb.Df, mid.(pb.parameters),
        mid.(pb.P), mid.(pb.DP), 0.0, 0.0,
        mid.(pb.Q), mid.(pb.DQ), 0.0, 0.0,
        mid(pb.integration_time))

    # Step 2
    x_bar, newton_success = approximate_zero(pb_approx, α_init, order)
    @assert newton_success
    @assert abs(x_bar[end-1]) < 1
    @assert abs(x_bar[end])   < 1

    # Step 3
    A = approximate_inverse(pb_approx, x_bar)

    # Step 4
    X_C = Ell1(GeometricWeight(μ))
    X = NormedCartesianSpace((X_C, X_C, X_C, Ell1(), Ell1(), Ell1()), Ell1())
    #-
    Y = Y_bound(pb, interval(x_bar), A, X)
    R = interval(10sup(Y))
    Z = Z_bound(pb, interval(x_bar), A, R, X_C, X)
    #-
    ie, proof_success = interval_of_existence(Y, Z, R; verbose = true)
    @assert proof_success

    #
    return interval.(unpack(x_bar)), interval(inf(ie))
end





# Step 1

rig_eval(P, θ, r) = Evaluation(θ) .* P
rig_eval(P, θ::NTuple{2,Interval}, r) = interval.(Evaluation(θ) .* P, r; format = :midpoint)

function F(pb::HeteroclinicProblem, x)
    u = x[1:3]
    α, θ₁, θ₂ = [x[j][1] for j ∈ 4:6]

    γ = (cis(α), cis(-α))
    Pγ = real( rig_eval(pb.P, γ, pb.r_P) )
    F₁ = u - Pγ - pb.integration_time / exact(2) * (Integral(1) .* pb.f(u, pb.parameters))

    F₂ = Evaluation(exact(1)) .* u - rig_eval(pb.Q, (θ₁, θ₂), pb.r_Q)

    return [F₁ ; F₂]
end

function DF(pb::HeteroclinicProblem, x)
    u = x[1:3]
    α, θ₁, θ₂ = [x[j][1] for j ∈ 4:6]

    DuF₁ = Diagonal([exact(I), exact(I), exact(I)]) - pb.integration_time / exact(2) * (Integral(1) .* Multiplication.(pb.Df(u, pb.parameters)))
    DuF₂ = Diagonal([Evaluation(exact(1)), Evaluation(exact(1)), Evaluation(exact(1))])

    γ = (cis(α), cis(-α))
    Dγ = exact(im) .* (cis(α), -cis(-α))
    DP_Dγ = real(sum(i -> LinearOperator.(ScalarSpace(), Chebyshev(0), hcat.(rig_eval(pb.DP[:,i], γ, pb.r_P/(exact(ℯ) * log(pb.ν_P))) .* Dγ[i])), 1:2))
    DQ_θ = [LinearOperator.(ScalarSpace(), ScalarSpace(), hcat.(rig_eval(pb.DQ[:,1], (θ₁, θ₂), pb.r_Q/(exact(ℯ) * abs(θ₁) * abs(log(θ₁/pb.ν_Q)))))) ;;
            LinearOperator.(ScalarSpace(), ScalarSpace(), hcat.(rig_eval(pb.DQ[:,2], (θ₁, θ₂), pb.r_Q/(exact(ℯ) * abs(θ₂) * abs(log(θ₂/pb.ν_Q))))))]

    return [DuF₁  -DP_Dγ                  fill(zeros(ExactReal{Int}, ScalarSpace(), Chebyshev(0)), 3, 2)
            DuF₂   fill(zeros(ExactReal{Int}, Chebyshev(0), ScalarSpace()), 3, 1)  -DQ_θ]
end





# Step 2

function approximate_zero(pb_approx::HeteroclinicProblem, α_init::Float64, order::Int)
    u_ic = real( Evaluation(cis(α_init), cis(-α_init)) .* pb_approx.P )
    pb_ode = ODEProblem((du, u, p, t) -> du .= pb_approx.f(u, p), u_ic, (0.0, pb_approx.integration_time), pb_approx.parameters)
    sol = solve(pb_ode; alg = Tsit5())

    x_init = zeros(Chebyshev(order) × Chebyshev(order) × Chebyshev(order) × ScalarSpace() × ScalarSpace() × ScalarSpace())
    for j ∈ 1:3
        project!(component(x_init, j),
                real(to_coef(s -> sol((s + 1) * pb_approx.integration_time / 2; idxs = j), Chebyshev(order))))
    end
    x_init[end-2] = α_init
    # approximate coordinates on the stable manifold
    θ₁_init, θ₂_init = (Evaluation(0, 0) .* pb_approx.DQ) \ (sol(pb_approx.integration_time) - Evaluation(0, 0) .* pb_approx.Q)
    x_init[end-1] = θ₁_init
    x_init[end]   = θ₂_init

    F_DF_newton(x) = F(pb_approx, unpack(x)), DF(pb_approx, unpack(x))
    return newton(F_DF_newton, x_init; verbose = true)
end





# Step 3

function approximate_inverse(pb_approx::HeteroclinicProblem, x_bar)
    Π = Projection(space(x_bar))
    A_finite = inv(mid(Π * DF(pb_approx, unpack(x_bar)) * Π))
    return interval(A_finite) + (interval(I) - interval(Π))
end





# Step 4

Y_bound(pb::HeteroclinicProblem, x_bar, A, X) = norm(A * F(pb, unpack(x_bar)), X)

function Z_bound(pb::HeteroclinicProblem, x_bar, A, R, X_C, X)
    K = order(component(x_bar, 1)) # since each component of u_bar have the same truncation order
    Π = Projection(space(x_bar))
    Π_2Kp1 = interval(Projection(Chebyshev(2K+1) × Chebyshev(2K+1) × Chebyshev(2K+1) ×
                                 ScalarSpace() × ScalarSpace() × ScalarSpace()))

    opnorm_A = max(opnorm(A * Π, X), interval(1))
    μ = rate(X_C)
    opnorm_L = interval(1) + μ
    first_row_tail_bound_L = inv(μ^interval(K+2) * (interval(K+2)^interval(2) - interval(1)))
    tail_bound_L = inv(μ * interval(2*(K+1))) + μ / (interval(2*(K+3)))

    B = DF(pb, unpack(x_bar))
    W = pb.Df(unpack(x_bar)[1:3], pb.parameters)
    M = opnorm_A * first_row_tail_bound_L + tail_bound_L
    Z₀ = max(opnorm(Π_2Kp1 - A * (B * Π_2Kp1), X),
             M * pb.integration_time / interval(2) * opnorm(norm.(W, X_C), 1) )

    Z₁ = exact(2) * R * pb.integration_time / interval(2) * opnorm_A * opnorm_L

    return Z₀ + Z₁
end
