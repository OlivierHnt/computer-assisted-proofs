# A steady state of the Shigesada-Kawasaki-Teramoto system, Δ Φ(v, λ) + R(v) = 0 on cosine series, at the last node of
# the arc u_bar1a of chebyshev-rigorous-continuation, then the branch through it over an interval of λ.
# The diffusion is nonlinear: the tail of the inverse is M(W) Δ⁻¹.

include(joinpath(@__DIR__, "skt_problem.jl"))
using Capsule: Capsule, QuasilinearProblem, ProjectedInverse, ContinuationProblem, ContinuationInverse, NewtonKantorovich, segment, validate, radius

#- the last of the 17 nodes of the arc: (v₁, v₂, λ), cosine series of order 80
data = open(io -> read!(io, Vector{Float64}(undef, 17(2K + 3))), joinpath(@__DIR__, "..", "..", "chebyshev-rigorous-continuation", "data", "skt_u_bar1a.bin"))[16(2K+3)+1:end]
v, λ = [Sequence(S, data[1:K+1]), Sequence(S, data[K+2:2K+2])], data[2K+3]

pb   = problem(interval(λ))
A    = ProjectedInverse(pb, v, [evensym(Fourier(K_A, π)), evensym(Fourier(K_A, π))])
cert = validate(pb, v, NewtonKantorovich{1}(A))
println("λ = ", λ)
println(cert)
println("Z = "); display(cert.Z)

### The branch for λ ∈ [λ, λ + δ]: Newton at the Chebyshev nodes, then interpolation

N, δ = 8, 1e-4
S² = S × S
u_grid = Vector{Sequence}(undef, N+1)
u = Sequence(S², [coefficients(v[1]) ; coefficients(v[2])])
for j ∈ N:-1:0
    pb_j = problem(λ + δ * (1 + cospi(j / N)) / 2)
    global u, converged = newton(u -> (Capsule.F(pb_j, unpack(u)), Capsule.DF(pb_j, unpack(u))), u)
    converged || error("Newton failed at the node $j")
    u_grid[j+1] = u
end
x = [interval(to_coef([component(u, i) for u ∈ u_grid], Chebyshev(N))) for i ∈ 1:2]

cpb    = ContinuationProblem(problem(interval(λ)), segment(Chebyshev, interval(λ), interval(λ) + interval(δ)), Ell1())
A_t    = ContinuationInverse(cpb, x, [evensym(Fourier(K_A, π)), evensym(Fourier(K_A, π))])
cert_t = validate(cpb, x, NewtonKantorovich{1}(A_t))
println(cert_t)
println("Z = "); display(cert_t.Z)
