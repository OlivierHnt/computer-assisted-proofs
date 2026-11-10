# The lift of a finite-dimensional problem to a parameter interval, mocking the validate pipeline.
# For a fixed parameter the unknown is a vector of numbers; over the interval it is a vector of
# Chebyshev series in t ∈ [-1, 1]. The problem and the bounds are the same code: only the
# coefficient ring changes (numbers → series, abs → ℓ¹ norm) and the inverse becomes a type
# whose products split a (sequence + error) into a grid product and a norm product.

using RadiiPolynomial, LinearAlgebra

### Problems

struct FiniteDimensionalProblem{FF,FDF,P}
    F  :: FF
    DF :: FDF
    p  :: P
end
F(pb::FiniteDimensionalProblem, x)  = pb.F(x, pb.p)
DF(pb::FiniteDimensionalProblem, x) = pb.DF(x, pb.p)

# the same problem with the parameter running over [p_lo, p_hi], as the series p(t) = mid + half t
struct ContinuationProblem{P,X}
    pointwise :: P
    X         :: X    # Banach space in t
end
function ContinuationProblem(pb::FiniteDimensionalProblem, (p_lo, p_hi), X)
    p = Sequence(Chebyshev(1), [(p_lo + p_hi) / exact(2), (p_hi - p_lo) / exact(4)])   # T₁ carries the weight 2
    return ContinuationProblem(FiniteDimensionalProblem(pb.F, pb.DF, p), X)
end
F(pb::ContinuationProblem, x)  = F(pb.pointwise, x)
DF(pb::ContinuationProblem, x) = DF(pb.pointwise, x)

### Inverses

struct FiniteInverse{T}
    matrix :: T
end

# A(t) = Σₙ Aₙ Tₙ(t), a matrix of Chebyshev series. Its product with a vector or matrix of
# (sequence + error): the sequence parts through a grid resolving the product degree, exactly;
# the errors through the block norms of A.
struct ChebyshevInverse{T,X}
    coefficients :: T
    X            :: X
end
opnorms(A::ChebyshevInverse) = norm.(A.coefficients, A.X)

function Base.:*(A::ChebyshevInverse, f::AbstractVecOrMat)
    d = maximum(order, A.coefficients) + maximum(order ∘ sequence, f)
    m = nextpow(2, d) + 1
    A_grid = to_grid.(A.coefficients, m)
    f_grid = to_grid.(sequence.(f), m)
    Af = [real(to_coef(sum(A_grid[i,k] .* f_grid[k,j] for k ∈ axes(f, 1)), Chebyshev(d))) for i ∈ axes(A_grid, 1), j ∈ axes(f, 2)]
    e = interval.(sup.(opnorms(A) * total_error.(f)))   # the fast interval matrix product may give a negative lower bound
    return reshape([InfiniteSequence(Af[i], A.X; total_error = e[i]) for i ∈ eachindex(Af)], size(Af))
end
Base.:*(A::ChebyshevInverse, f::AbstractVector) = vec(A * reshape(f, :, 1))

### Bounds: the same two formulas, on the two coefficient rings

Y_bound(pb::FiniteDimensionalProblem, x, A::FiniteInverse) = abs.(interval(A.matrix) * F(pb, interval(x)))
Z_bound(pb::FiniteDimensionalProblem, x, A::FiniteInverse, R) = abs.(interval(I) - interval(A.matrix) * DF(pb, ball(x, R)))
ball(x::AbstractVector{<:Number}, R) = interval.(x, R; format = :midpoint)

# F and DF are evaluated on balls (of radius 0 for Y): a non-polynomial f then returns (sequence + error) itself
Y_bound(pb::ContinuationProblem, x, A::ChebyshevInverse) = norm.(A * F(pb, ball(x, zero(eltype(x[1])), pb.X)), pb.X)
Z_bound(pb::ContinuationProblem, x, A::ChebyshevInverse, R) = norm.(_one_minus(A * DF(pb, ball(x, R, pb.X))), pb.X)
ball(x::AbstractVector, R, X) = [InfiniteSequence(interval(xᵢ), X; total_error = interval(r)) for (xᵢ, r) ∈ zip(x, R)]
ball(x::AbstractVector, R::Number, X) = ball(x, fill(R, length(x)), X)
_one_minus(M::AbstractMatrix) = [i == j ? exact(1) - M[i,j] : -M[i,j] for i ∈ axes(M, 1), j ∈ axes(M, 2)]

function validate(pb, x, A; R = nothing)
    Y = Y_bound(pb, x, A)
    R_ = isnothing(R) ? fill(exact(10maximum(sup.(Y))), length(Y)) : R
    Z = Z_bound(pb, x, A, R_)
    W = zeros(eltype(Z), length(Y), length(Y), length(Y))
    radii, η, success = set_of_radii(Y, Z, W, R_)
    return (; Y, Z, R = R_, radii, success)
end

### Approximation (floating-point arithmetic): Newton at the Chebyshev nodes, then interpolation

nodes(N) = [-cospi(j / N) for j ∈ 0:N]
interpolate(values, N) = interval(real(to_coef(reverse(values), Chebyshev(N))))   # values at nodes(N)

function approximate(pb::FiniteDimensionalProblem, x0, (p_lo, p_hi), N)
    p_nodes = (p_lo + p_hi) / 2 .+ (p_hi - p_lo) / 2 .* nodes(N)
    x_nodes = [first(newton(x -> (pb.F(x, p), pb.DF(x, p)), x0)) for p ∈ p_nodes]
    A_nodes = [inv(pb.DF(x, p)) for (x, p) ∈ zip(x_nodes, p_nodes)]
    x̄ = [interpolate([x[i] for x ∈ x_nodes], N) for i ∈ eachindex(x0)]
    Ā = [interpolate([A[i,j] for A ∈ A_nodes], N) for i ∈ eachindex(x0), j ∈ eachindex(x0)]
    return x̄, Ā
end

### 1. F(x, p) = x² - 2 + p

F₁(x, p)  = [x[1]^2 - exact(2) + p]
DF₁(x, p) = [exact(2) * x[1];;]

#- pointwise, p = 0.3
pb = FiniteDimensionalProblem(F₁, DF₁, 0.3)
x̄, _ = newton(x -> (F(pb, x), DF(pb, x)), [1.3])
cert = validate(FiniteDimensionalProblem(F₁, DF₁, interval(0.3)), x̄, FiniteInverse(inv(DF(pb, x̄))))
println("pointwise p = 0.3:   |x - x̄| ≤ ", inf(cert.radii[1]), "   (exact: ", sqrt(2 - 0.3) - x̄[1], ")")

#- lifted, p ∈ [-1, 1]
N = 16
X = Ell1()
x̄_t, Ā_t = approximate(pb, [1.3], (-1.0, 1.0), N)
cpb  = ContinuationProblem(FiniteDimensionalProblem(F₁, DF₁, nothing), (interval(-1), interval(1)), X)
cert_t = validate(cpb, x̄_t, ChebyshevInverse(Ā_t, X))
println("lifted p ∈ [-1, 1]:  sup |x(p) - x̄(p)| ≤ ", inf(cert_t.radii[1]), "   proved = ", cert_t.success)
println("check at p = 0.3:    |√1.7 - x̄(0.3)| = ", abs(sqrt(2 - 0.3) - mid(x̄_t[1](0.3))))

### 2. a non-polynomial nonlinearity: F(x, p) = exp(x) - 2 - p, exact solution log(2 + p)

# the rigorous exp of a (sequence + error) is RadiiPolynomial's; the inverse never sees exp, it only splits what comes out of it
F₂(x, p)  = [exp(x[1]) - exact(2) - p]
DF₂(x, p) = [exp(x[1]);;]

pb₂ = FiniteDimensionalProblem(F₂, DF₂, nothing)
X₂  = Ell1(GeometricWeight(interval(1.2)))
x̄₂, Ā₂ = approximate(FiniteDimensionalProblem(F₂, DF₂, 0.0), [0.7], (-0.5, 0.5), 10)
cpb₂ = ContinuationProblem(pb₂, (interval(-0.5), interval(0.5)), X₂)
cert₂ = validate(cpb₂, x̄₂, ChebyshevInverse(Ā₂, X₂))
println("exp, p ∈ [-1/2, 1/2]: sup |x(p) - x̄(p)| ≤ ", inf(cert₂.radii[1]), "   proved = ", cert₂.success,
        "   Y = ", sup(cert₂.Y[1]), "   Z = ", sup(cert₂.Z[1]))
println("check at p = 0.3:    |log(2.3) - x̄(t = 0.6)| = ", abs(log(2.3) - mid(x̄₂[1](0.6))))   # p = t/2
