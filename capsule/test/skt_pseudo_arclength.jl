# The arc u_bar1a of steady states of the Shigesada-Kawasaki-Teramoto system (chebyshev-rigorous-continuation), by
# pseudo-arclength: the parameter λ is an unknown, and (v₁, v₂, λ) lies on the hyperplane through a point of the arc
# orthogonal to its tangent. First at one node, then along the whole arc, the arclength being the path.

include(joinpath(@__DIR__, "skt_problem.jl"))
using Capsule: Capsule, PseudoArclengthProblem, ProjectedInverse, ContinuationProblem, ContinuationInverse, NewtonKantorovich, validate, radius

pseudo_arclength(tangent, point) = PseudoArclengthProblem(problem(nothing), (v, λ) -> v, nothing, (tangent = tangent, point = point))

#- the 17 nodes of the arc and their tangents: (v₁, v₂, λ), cosine series of order 80
N = 16
data = reshape(open(io -> read!(io, Vector{Float64}(undef, 2(N + 1) * (2K + 3))), joinpath(@__DIR__, "..", "..", "chebyshev-rigorous-continuation", "data", "skt_u_bar1a.bin")), 2K + 3, 2(N + 1))
blocks(c) = [Sequence(S, c[1:K+1]), Sequence(S, c[K+2:2K+2]), c[2K+3]]
u_grid, u̇_grid = [blocks(data[:,j]) for j ∈ 1:N+1], [blocks(data[:,N+1+j]) for j ∈ 1:N+1]
spaces_A = [evensym(Fourier(K_A, π)), evensym(Fourier(K_A, π)), ScalarSpace()]

### At the first node, next to the bifurcation point

pa   = pseudo_arclength(interval.(u̇_grid[1]), interval.(u_grid[1]))
cert = validate(pa, u_grid[1], NewtonKantorovich{1}(ProjectedInverse(pa, u_grid[1], spaces_A)))
println(cert)

### Along the arc: the Chebyshev interpolants of the nodes and of the tangents

series(grid) = [[interval(to_coef([u[i] for u ∈ grid], Chebyshev(N))) for i ∈ 1:2] ; interval(real(to_coef([u[3] for u ∈ grid], Chebyshev(N))))]
x_t, ẋ_t = series(u_grid), series(u̇_grid)

cpb    = ContinuationProblem(pseudo_arclength(nothing, nothing), (tangent = ẋ_t, point = x_t), Ell1(); dims = 1)
cert_t = validate(cpb, x_t, NewtonKantorovich{1}(ContinuationInverse(cpb, x_t, spaces_A)))
println(cert_t)
println("Z = "); display(cert_t.Z)
