# Periodic orbit of the Lorenz system: first-order Newton-Kantorovich, with DF
# evaluated on the ball through InfiniteSequence / InfiniteLinearOperator.
# Same problem, data and approximate inverse as the RadiiPolynomial example
# docs/src/examples/periodic_orbits/lorenz_po.md; only Step 3 differs.

using RadiiPolynomial, LinearAlgebra
using Capsule: Capsule, InfiniteLinearOperator, ProjectedInverse, opnorms

### Step 1: Formulation

function f(u, params)
    σ, ρ, β = params
    u₁, u₂, u₃ = u
    return [σ*(u₂ - u₁)
            u₁*(ρ - u₃) - u₂
            u₁*u₂ - β*u₃]
end

function Df(u, params)
    σ, ρ, β = params
    u₁, u₂, u₃ = u
    return [-σ*one(u₁)    σ*one(u₂)    zero(u₃)
            ρ-u₃          -one(u₂)     -u₁
            u₂            u₁           -β*one(u₃)]
end

# the unknown x = (u₁, u₂, u₃, τ); the blocks are the primitives

F(u, τ, params, ξ) = [Derivative(1) .* u .- τ .* f(u, params) ; sum(adjoint(ξ[i]) * u[i] for i ∈ 1:3)]

function DF(u, τ, params, ξ)
    O = zero(exact(I))
    L = [Derivative(1)  O              O
         O              Derivative(1)  O
         O              O              Derivative(1)]
    return [L - τ * Multiplication.(Df(u, params))    LinearOperator.(-f(u, params))
            permutedims(adjoint.(ξ))                   LinearOperator(exact(0))]
end

F(x, params, ξ)  = F(unpack(x)[1:3], unpack(x)[4][1], params, ξ)
DF(x, params, ξ) = DF(unpack(x)[1:3], unpack(x)[4][1], params, ξ)

### Step 2: Approximation (floating-point arithmetic)

σ, ρ, β = 10.0, 28.0, 8/3
K = 40

u_guess = [zeros(ComplexF64, Fourier(K, 1.0)) for _ ∈ 1:3]
u_guess[1][1:2:5] = [-2.9 - 4.3im, 1.6 - 1.1im,  0.3 + 0.4im]
u_guess[2][1:2:5] = [-1.2 - 5.4im, 3.0 + 0.8im, -0.4 + 1.1im]
u_guess[3][0:2:4] = [ 23,          3.8 + 4.7im, -1.8 + 0.9im]
u_guess[1][-5:2:-1] .= conj.(u_guess[1][5:-2:1])
u_guess[2][-5:2:-1] .= conj.(u_guess[2][5:-2:1])
u_guess[3][-4:2:0]  .= conj.(u_guess[3][4:-2:0])

ξ = differentiate.(u_guess)
τ_guess = 1.5/(2π)

space_x = Fourier(K, 1.0) × Fourier(K, 1.0) × Fourier(K, 1.0) × ScalarSpace()
x_guess = Sequence(space_x, [coefficients(u_guess[1]) ; coefficients(u_guess[2]) ; coefficients(u_guess[3]) ; τ_guess])

x_bar, converged = newton(x -> (F(x, (σ, ρ, β), ξ), DF(x, (σ, ρ, β), ξ)), x_guess; verbose = true)
conjugacy_symmetry!(x_bar)

#- the approximate inverse A = A_K Π + (I - Π) ∫ on the u-blocks

params_i = (interval(10), interval(28), interval(8)/interval(3))
ξ_i = interval.(ξ)
x_i = interval(x_bar)
u_i, τ_i = unpack(x_i)[1:3], unpack(x_i)[4][1]

Π = Projection(space_x)
A_K = interval(inv(mid(Π * DF(x_i, params_i, ξ_i) * Π)))

X_F = Ell1(GeometricWeight(interval(1)))
Xs  = [X_F, X_F, X_F, Ell1()]

tail = (interval(I) - interval(Projection(Fourier(K, interval(1))))) * Integral(1)
T = [i == j ≤ 3 ? tail : nothing for i ∈ 1:4, j ∈ 1:4]
A = ProjectedInverse(A_K, T, inv(interval(K+1)), Xs)

### Step 3: Bounds (interval arithmetic), first order, one radius per component

M = 4

#- Y bound, one per component
AF = A * F(u_i, τ_i, params_i, ξ_i)
Y  = norm.(AF, Xs)

#- Z bound: the matrix of block norms of I - A DF(x) over the product of balls B(ū_i, R_i) × B(τ̄, R_4),
#  columns split at q = 2K+1
R   = fill(10maximum(sup.(Y)), 4)
u_R = [InfiniteSequence(u_i[i], X_F; total_error = interval(R[i])) for i ∈ 1:3]
τ_R = interval(real(τ_i), R[4]; format=:midpoint)

Πs = Diagonal([[interval(Projection(Fourier(2K+1, interval(1)))) for _ ∈ 1:3] ; interval(Projection(ScalarSpace()))])
D  = Πs - A * (DF(u_R, τ_R, params_i, ξ_i) * Πs)
Z_finite = opnorm.(D, permutedims(Xs), Xs)   # domain of column i, codomain of row m

# columns j > 2K+1 of the u-blocks: (I - A DF)_{m,i} e_j = τ A M(Df(u))_{m,i} e_j, since A ∂ e_j = e_j and ξ' e_j = 0.
# The sequence part of Df(u) sends e_j above K, where A = (I - Π) ∫ has norm 1/(K+1); its error R meets all of A
W = Df(u_R, params_i)
Z_tail = (abs(τ_i) + exact(R[4])) .* (opnorms(A)[:, 1:3] * total_error.(W))
Z_tail = [Z_tail zeros(eltype(Z_tail), 4)]
Z_tail[1:3, 1:3] .+= (abs(τ_i) + exact(R[4])) .* A.tail_bound .* norm.(sequence.(W), X_F)

Z = max.(Z_finite, Z_tail)         # exact per block: every domain is a weighted ℓ¹
W = zeros(eltype(Z), 4, 4, 4)      # first order

r, η, proved = set_of_radii(Y, Z, W, R; verbose = true)

### Step 4: Conclusion

println("Y = ", Y)
println("Z = "); display(Z)
println("R = ", R, "
proved = ", proved, "
radii r = ", r)
