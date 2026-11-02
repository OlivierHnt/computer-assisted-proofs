# Local manifolds (Section 3 of the paper, Proposition 3.3). `validate` takes
# the approximate Taylor coefficients P_bar and returns the Y, Z₀, Z₁ bounds and
# the radius of Theorem 1.1.
# From ShimizuMoriokaTutorial.jl/src/manifold.jl; Step 2 (Newton) is in search/.

struct ManifoldProblem
    f
    Df
    parameters
    equilibrium
    eigvals
    eigvecs
end

#-

function validate(pb::ManifoldProblem, P_bar::Sequence, ν::Real)
    pb_approx = ManifoldProblem(pb.f, pb.Df, mid.(pb.parameters), mid.(pb.equilibrium), mid.(pb.eigvals), mid.(pb.eigvecs))

    # Step 2: P_bar is an input (data/Q_bar.bin, data/P_bar.bin)

    # Step 3
    A = approximate_inverse(pb_approx, P_bar)

    # Step 4
    X_T = Ell1(GeometricWeight(ν), GeometricWeight(ν))
    X = NormedCartesianSpace(X_T, Ell1())
    #-
    Y = Y_bound(pb, interval(P_bar), A, X)
    R = interval(10sup(Y))
    Z, Z₀, Z₁ = Z_bound(pb, interval(P_bar), A, R, X_T, X)
    #-
    ie, proof_success = interval_of_existence(Y, Z, R)

    #
    return (; Y, Z₀, Z₁, Z, R, ie, proof_success, P = interval.(unpack(P_bar)), r = interval(inf(ie)))
end





# Step 1

struct L_manifold <: AbstractDiagonalOperator
    λ₁
    λ₂
end

function RadiiPolynomial.getcoefficient(L::L_manifold, (codom, k)::T, (dom, l)::T) where {T<:Tuple{TensorSpace{NTuple{2,Taylor}},NTuple{2,Int}}}
    x = inv(exact(k[1]) * L.λ₁ + exact(k[2]) * L.λ₂)
    k == l && sum(k) ≥ 2 && return x
    return zero(x)
end

function F(pb::ManifoldProblem, u)
    L = L_manifold(pb.eigvals[1], pb.eigvals[2])

    θ₁ = Sequence(Taylor(1) ⊗ Taylor(0), [exact(0), exact(1)])
    θ₂ = Sequence(Taylor(0) ⊗ Taylor(1), [exact(0), exact(1)])
    ϕ = pb.equilibrium + pb.eigvecs[:,1] .* θ₁ + pb.eigvecs[:,2] .* θ₂

    return u - ϕ - L .* pb.f(u, pb.parameters)
end

function DF(pb::ManifoldProblem, u)
    L = L_manifold(pb.eigvals[1], pb.eigvals[2])

    return Diagonal([exact(I), exact(I), exact(I)]) - L .* Multiplication.(pb.Df(u, pb.parameters))
end





# Step 3

function approximate_inverse(pb_approx::ManifoldProblem, P_bar)
    Π = Projection(space(P_bar))
    A_finite = inv(mid(Π * DF(pb_approx, unpack(P_bar)) * Π))
    return interval(A_finite) + (interval(I) - interval(Π))
end





# Step 4

Y_bound(pb::ManifoldProblem, P_bar, A, X) = norm(A * F(pb, unpack(P_bar)), X)

function Z_bound(pb::ManifoldProblem, P_bar, A, R, X_T, X)
    K = order(component(P_bar, 1))[1] # since each component of P_bar have the same truncation order, in both variables
    Π = interval(Projection(space(P_bar)))

    L = Diagonal([L_manifold(pb.eigvals[1], pb.eigvals[2]) for j ∈ 1:3])
    tail_bound_L = inv(interval(K + 1) * minimum(abs ∘ real, pb.eigvals))
    opnorm_AL = max(opnorm(A * (L * Π), X), tail_bound_L)

    B = DF(pb, unpack(P_bar))
    W = pb.Df(unpack(P_bar), pb.parameters)
    Z₀ = max(opnorm(Π - A * (B * Π), X),
             tail_bound_L * opnorm(norm.(W, X_T), 1))

    Z₁ = exact(2) * R * opnorm_AL

    return Z₀ + Z₁, Z₀, Z₁
end
