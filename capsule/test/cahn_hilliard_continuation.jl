# A branch of steady states of Cahn-Hilliard, Δu + β(u - u³ - c) = 0 on cosine series, for c ∈ [-0.1, 0.1] and β = 41:
# the edge γ₄ of the region of the RadiiPolynomial example docs/src/examples/continuation/cahn_hilliard_cont.md.
# The pointwise problem of logistic_ivp.jl, lifted along the path c(t) = -0.1 t.

using RadiiPolynomial, LinearAlgebra
using Capsule: Capsule, SteadyStateProblem, ContinuationProblem, ContinuationInverse, NewtonKantorovich, segment, validate, radius

struct PseudoInverseLaplacian <: RadiiPolynomial.AbstractDiagonalOperator end
RadiiPolynomial.getcoefficient(::PseudoInverseLaplacian, (codom, i)::Tuple{SymmetricSpace{<:Fourier},Integer}, (dom, j)::Tuple{SymmetricSpace{<:Fourier},Integer}) =
    (i == j) & !(i == j == 0) ? inv(-(frequency(dom) * exact(i))^2) : zero(frequency(dom))

N_ch(u, (c, β))  = [β * (u[1] - u[1]^3 - c)]
DN_ch(u, (c, β)) = [β * (exact(1) - exact(3) * u[1]^2);;]
problem(c, β) = SteadyStateProblem(Laplacian(), PseudoInverseLaplacian(), N_ch, DN_ch, (c, β), [Ell1()])

### Approximation (floating-point arithmetic): Newton at the Chebyshev nodes t_j = cos(jπ/N), then interpolation

K, N = 20, 16
β = 41.0

u_grid = Vector{Sequence}(undef, N+1)
u_grid[1], _ = newton(u -> (Capsule.F(problem(-0.1, β), [u])[1], Capsule.DF(problem(-0.1, β), [u])[1,1]),
                      Sequence(evensym(Fourier(K, π)), [-0.1022666446473428 ; 0 ; 0.044213651292169816 ; zeros(K-2)]))
for j ∈ 1:N
    pb_j = problem(-0.1cospi(j/N), β)
    u_grid[j+1], converged = newton(u -> (Capsule.F(pb_j, [u])[1], Capsule.DF(pb_j, [u])[1,1]), u_grid[j])
    converged || error("Newton failed at the node $j")
end
x = [interval(to_coef(u_grid, Chebyshev(N)))]

### Proof

cpb  = ContinuationProblem(problem(interval(0), interval(β)), (segment(Chebyshev, interval(0.1), interval(-0.1)), interval(β)), Ell1())
cert = validate(cpb, x, NewtonKantorovich{1}(ContinuationInverse(cpb, x)))
println(cert, "   Z = ", sup(cert.Z[1]))

#- check at c = -0.03, that is t = 0.3, against Newton there: T_n(t) = cos(n acos t)
u_c, _ = newton(u -> (Capsule.F(problem(-0.03, β), [u])[1], Capsule.DF(problem(-0.03, β), [u])[1,1]), u_grid[1])
x_c = [sum(mid(real(x[1][(n, k)])) * (n == 0 ? 1 : 2cos(n * acos(0.3))) for n ∈ 0:N) for k ∈ 0:K]
println("check at c = -0.03:   ‖u - x̄(0.3)‖ = ", sum(k -> (k == 0 ? 1 : 2) * abs(u_c[k] - x_c[k+1]), 0:K))
