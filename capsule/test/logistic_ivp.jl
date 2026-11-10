# The logistic equation u' = u(1 - u), u(0) = 1/2 on [-2, 2], as in the RadiiPolynomial docs
# (getting_started/logistic_equation.md), through the framework: first order, componentwise.

using RadiiPolynomial, LinearAlgebra
using Capsule: Capsule, InitialValueProblem, SteadyStateProblem, ProjectedInverse, NewtonKantorovich, validate, radius, solution

### Step 1: Formulation, on the blocks u = (u₁)

f(u, p)  = [u[1] * (exact(1) - u[1])]
Df(u, p) = [exact(1) - exact(2) * u[1];;]

pb = InitialValueProblem(f, Df, nothing, [exact(0.5)], [Ell1(GeometricWeight(interval(2)))])

### Step 2: Approximation (floating-point arithmetic)

K = 27
u_bar, converged = newton(u -> (Capsule.F(pb, [u])[1], Capsule.DF(pb, [u])[1,1]), zeros(Taylor(K)); verbose = true)

### Step 3-4: Bounds and conclusion

cert = validate(pb, [u_bar], NewtonKantorovich{1}(ProjectedInverse(pb, [u_bar])))
println(cert)
println("u ∈ ", solution(cert)[1])

#- the same equation on Chebyshev, u(-1) = 1/2 on [-1, 1]: the exact solution is 1/(1 + e^{-(t+1)})
pb_c = InitialValueProblem(f, Df, nothing, [exact(0.5)], [Ell1(GeometricWeight(interval(1.5)))])
u_c, converged_c = newton(u -> (Capsule.F(pb_c, [u])[1], Capsule.DF(pb_c, [u])[1,1]), zeros(Chebyshev(40)); verbose = false)
cert_c = validate(pb_c, [u_c], NewtonKantorovich{1}(ProjectedInverse(pb_c, [u_c])))
println(cert_c)
sol_c = solution(cert_c)[1]
u1 = sequence(sol_c)(interval(1)) + interval(-1, 1) * total_error(sol_c)
println("Chebyshev u(1) ∈ ", u1, "  exact ", inv(interval(1) + exp(interval(-2))), "  contained: ", issubset_interval(inv(interval(1) + exp(interval(-2))), u1))
sol_t = solution(cert)[1]
u2 = sequence(sol_t)(interval(2)) + interval(-1, 1) * total_error(sol_t)
println("Taylor    u(2) ∈ ", u2, "  contained: ", issubset_interval(inv(interval(1) + exp(interval(-2))), u2))

#- the docs' second-order proof, for comparison
Π = Projection(Taylor(K)); A_i = interval(inv(Π * Capsule.DF(pb, [u_bar])[1,1] * Π)) + (interval(I) - interval(Π))
u_i = interval(u_bar); ν = interval(2); X_T = Ell1(GeometricWeight(ν))
Y  = norm(A_i * Capsule.F(pb, [u_i])[1], X_T)
Z₁ = opnorm(interval(Projection(Taylor(K+1))) - A_i * Capsule.DF(pb, [u_i])[1,1] * Projection(Taylor(K+1)), X_T)
Z₂ = max(opnorm(A_i * Π, X_T), interval(1)) * ν * interval(2)
ie, proved = interval_of_existence(Y, Z₁, Z₂, Inf)
println("docs, second order: radius ", inf(ie), ", proved = ", proved)


### A planar vector field: Lotka-Volterra, u₁' = u₁(1 - u₂), u₂' = u₂(u₁ - 1), u(0) = (2, 1), on [-1, 1]

f_lv(u, p)  = [u[1] * (exact(1) - u[2]), u[2] * (u[1] - exact(1))]
Df_lv(u, p) = [exact(1) - u[2]    -u[1]
               u[2]               u[1] - exact(1)]

ν_lv  = interval(1)
pb_lv = InitialValueProblem(f_lv, Df_lv, nothing, [exact(2), exact(1)], [Ell1(GeometricWeight(ν_lv)), Ell1(GeometricWeight(ν_lv))])

K_lv = 40
x_lv, converged_lv = newton(x -> (Capsule.F(pb_lv, unpack(x)), Capsule.DF(pb_lv, unpack(x))), zeros(Taylor(K_lv)^2); verbose = true)
u_lv = unpack(x_lv)

cert_lv = validate(pb_lv, u_lv, NewtonKantorovich{1}(ProjectedInverse(pb_lv, u_lv)))
println(cert_lv)
println("Z = "); display(cert_lv.Z)


### A steady state: Cahn-Hilliard, Δu + β(u - u³ - c) = 0 on cosine series (docs, examples/steady_states/cahn_hilliard.md)

struct PseudoInverseLaplacian <: RadiiPolynomial.AbstractDiagonalOperator end
RadiiPolynomial.getcoefficient(::PseudoInverseLaplacian, (codom, i)::Tuple{SymmetricSpace{<:Fourier},Integer}, (dom, j)::Tuple{SymmetricSpace{<:Fourier},Integer}) =
    (i == j) & !(i == j == 0) ? inv(-(frequency(dom) * exact(i))^2) : zero(frequency(dom))

N_ch(u, (c, β))  = [β * (u[1] - u[1]^3 - c)]
DN_ch(u, (c, β)) = [β * (exact(1) - exact(3) * u[1]^2);;]

K_ch = 32
c, β = 0.0, 40.15625
pb_ch_float = SteadyStateProblem(Laplacian(), PseudoInverseLaplacian(), N_ch, DN_ch, (c, β), [Ell1()])              # for Newton
pb_ch       = SteadyStateProblem(Laplacian(), PseudoInverseLaplacian(), N_ch, DN_ch, (interval(c), interval(β)), [Ell1()])  # for the proof

u_init = Sequence(evensym(Fourier(K_ch, π)), [c ; 0.45 ; zeros(K_ch-1)])
u_ch, converged_ch = newton(u -> (Capsule.F(pb_ch_float, [u])[1], Capsule.DF(pb_ch_float, [u])[1,1]), u_init; verbose = true)

cert_ch = validate(pb_ch, [u_ch], NewtonKantorovich{1}(ProjectedInverse(pb_ch, [u_ch])))
println(cert_ch)

#- the docs' second-order proof, for comparison
Π = Projection(space(u_ch)); u_i = interval(u_ch); c_i, β_i = interval(c), interval(β)
F_ch(u) = Laplacian() * u + β_i * (u - u^3 - c_i); DF_ch(u) = Laplacian() + Multiplication(β_i * (exact(1) - exact(3) * u^2))
A_ch = interval(inv(Π * DF_ch(u_ch) * Π)) + PseudoInverseLaplacian() * (interval(I) - interval(Π))
Y_ch = norm(A_ch * F_ch(u_i), Ell1())
Π_3K = interval(Projection(evensym(Fourier(3K_ch+1, π))))
Z₁_ch = max(opnorm(Π_3K - A_ch * (DF_ch(u_i) * Π_3K), Ell1()), β_i * (interval(1) + interval(3) * norm(u_i, Ell1())^2) / (interval(π) * interval(K_ch+2))^2)
R_ch = exact(10sup(Y_ch))
Z₂_ch = exact(3) * β_i * max(opnorm(interval(inv(Π * DF_ch(u_ch) * Π)), Ell1()), interval(1)/interval(π)^2) * (exact(2) * norm(u_i, Ell1()) + R_ch)
ie_ch, proved_ch = interval_of_existence(Y_ch, Z₁_ch, Z₂_ch, R_ch)
println("docs, second order: radius ", inf(ie_ch), ", proved = ", proved_ch)
