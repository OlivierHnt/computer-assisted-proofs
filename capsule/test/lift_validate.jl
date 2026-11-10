# Lifts through the framework: the two problems of lift_prototype.jl, the first again on Taylor series, and an
# initial value problem.

using RadiiPolynomial, LinearAlgebra
using Capsule: Capsule, InitialValueProblem, FiniteDimensionalProblem, ContinuationProblem, ContinuationInverse, FiniteInverse, segment, NewtonKantorovich, validate, radius, solution

### Approximation (floating-point arithmetic): Newton at the Chebyshev nodes, then interpolation

nodes(N) = [-cospi(j / N) for j ∈ 0:N]
interpolate(values, N) = interval(real(to_coef(reverse(values), Chebyshev(N))))   # values at nodes(N)

function approximate(pb::FiniteDimensionalProblem, x0, (p_lo, p_hi), N)
    p_nodes = (p_lo + p_hi) / 2 .+ (p_hi - p_lo) / 2 .* nodes(N)
    x_nodes = [first(newton(x -> (pb.F(x, p), pb.DF(x, p)), x0)) for p ∈ p_nodes]
    A_nodes = [inv(pb.DF(x, p)) for (x, p) ∈ zip(x_nodes, p_nodes)]
    x_bar = [interpolate([x[i] for x ∈ x_nodes], N) for i ∈ eachindex(x0)]
    A_bar = [interpolate([A[i,j] for A ∈ A_nodes], N) for i ∈ eachindex(x0), j ∈ eachindex(x0)]
    return x_bar, A_bar
end

### 1. F(x, p) = x² - 2 + p, p ∈ [-1, 1]

F₁(x, p)  = [x[1]^2 - exact(2) + p]
DF₁(x, p) = [exact(2) * x[1];;]

pb₁ = FiniteDimensionalProblem(F₁, DF₁, nothing)
X₁  = Ell1()
x₁, A₁ = approximate(pb₁, [1.3], (-1.0, 1.0), 16)
cert₁ = validate(ContinuationProblem(pb₁, segment(Chebyshev, interval(-1), interval(1)), X₁), x₁, NewtonKantorovich{1}(ContinuationInverse(FiniteInverse(A₁), X₁)))
println(cert₁)

### 2. F(x, p) = exp(x) - 2 - p, p ∈ [-1/2, 1/2]

F₂(x, p)  = [exp(x[1]) - exact(2) - p]
DF₂(x, p) = [exp(x[1]);;]

pb₂ = FiniteDimensionalProblem(F₂, DF₂, nothing)
X₂  = Ell1(GeometricWeight(interval(1.2)))
x₂, A₂ = approximate(pb₂, [0.7], (-0.5, 0.5), 10)
cert₂ = validate(ContinuationProblem(pb₂, segment(Chebyshev, interval(-0.5), interval(0.5)), X₂), x₂, NewtonKantorovich{1}(ContinuationInverse(FiniteInverse(A₂), X₂)))
println(cert₂, "   Z = ", sup(cert₂.Z[1]))
println("x(p) ∈ ", typeof(solution(cert₂)[1]).name.name, " of radius ", radius(cert₂)[1])

### 3. the first problem again, with Taylor series in t: Newton on the coefficients

s₃ = Taylor(40)
p₃ = Sequence(Taylor(1), [0.0, 1.0])
x₃, _ = newton(x -> (project(x * x - 2 + p₃, s₃), project(Multiplication(2x), s₃, s₃)), Sequence(s₃, [1.4 ; zeros(40)]))
A₃ = project(Multiplication(2x₃), s₃, s₃) \ Sequence(s₃, [1.0 ; zeros(40)])
cert₃ = validate(ContinuationProblem(pb₁, segment(Taylor, interval(-1), interval(1)), X₁), [interval(x₃)], NewtonKantorovich{1}(ContinuationInverse(FiniteInverse([interval(A₃);;]), X₁)))
println(cert₃)
println("check at p = 0.3:    |√1.7 - x̄(0.3)| = ", abs(sqrt(2 - 0.3) - x₃(0.3)))

### 4. a sequence problem: the logistic equation u' = λ u (1 - u), u(0) = 1/2 on [-1, 1], for λ ∈ [0.9, 1.1]

f₄(u, λ)  = [λ * u[1] * (exact(1) - u[1])]
Df₄(u, λ) = [λ * (exact(1) - exact(2) * u[1]);;]
problem₄(λ) = InitialValueProblem(f₄, Df₄, λ, [exact(0.5)], [Ell1(GeometricWeight(interval(1)))])

K₄, N₄ = 27, 10
u_grid = map(0:N₄) do j
    pb_j = problem₄(1 + 0.1cospi(j / N₄))
    first(newton(u -> (Capsule.F(pb_j, [u])[1], Capsule.DF(pb_j, [u])[1,1]), zeros(Taylor(K₄))))
end
x₄ = [interval(to_coef(u_grid, Chebyshev(N₄)))]

cpb₄  = ContinuationProblem(problem₄(interval(1)), segment(Chebyshev, interval(0.9), interval(1.1)), Ell1())
cert₄ = validate(cpb₄, x₄, NewtonKantorovich{1}(ContinuationInverse(cpb₄, x₄)))
println(cert₄, "   Z = ", sup(cert₄.Z[1]))
# λ = 1 is t = 0, where T_n(0) = cos(nπ/2), and the time 1 sums the Taylor coefficients
u₄ = sum(mid(real(x₄[1][(n, k)])) * (n == 0 ? 1 : 2cospi(n / 2)) for n ∈ 0:N₄, k ∈ 0:K₄)
println("check at λ = 1, time 1:   |1/(1 + e⁻¹) - x̄| = ", abs(inv(1 + exp(-1)) - u₄))
