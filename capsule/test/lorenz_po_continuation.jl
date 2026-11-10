# The branch of periodic orbits of Lorenz through the orbit of lorenz_po_first_order.jl, for ρ ∈ [28, 28 + δ]:
# a problem with a scalar unknown, the period, lifted along a path of parameters.

include(joinpath(@__DIR__, "lorenz_po_first_order.jl"))
using Capsule: Capsule, PeriodicOrbitProblem, ContinuationProblem, ContinuationInverse, NewtonKantorovich, segment, validate, radius

### Approximation (floating-point arithmetic): Newton at the Chebyshev nodes t_j = cos(jπ/N), then interpolation

N_t, δ = 6, 0.02
x_grid = map(0:N_t) do j
    params_j = (σ, ρ + δ * (1 + cospi(j / N_t)) / 2, β)
    x_j, converged_j = newton(x -> (F(x, params_j, ξ), DF(x, params_j, ξ)), x_bar)
    converged_j || error("Newton failed at the node $j")
    conjugacy_symmetry!(x_j)
end
x_t = [[interval(to_coef([component(x_j, i) for x_j ∈ x_grid], Chebyshev(N_t))) for i ∈ 1:3] ;
       interval(real(to_coef([real(component(x_j, 4)[1]) for x_j ∈ x_grid], Chebyshev(N_t))))]

### Proof

path   = (interval(σ), segment(Chebyshev, interval(ρ), interval(ρ) + interval(δ)), interval(8)/interval(3))
cpb    = ContinuationProblem(PeriodicOrbitProblem(f, Df, params_i, ξ_i, Xs), path, Ell1())
cert_t = validate(cpb, x_t, NewtonKantorovich{1}(ContinuationInverse(cpb, x_t)))
println(cert_t)
println("Z = "); display(cert_t.Z)
