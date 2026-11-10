# The Lorenz periodic orbit through the framework: the same data and Steps 1-2 as the template,
# then one `validate` call in place of its Step 3.

include(joinpath(@__DIR__, "lorenz_po_first_order.jl"))
using Capsule: Capsule, validate, radius, PeriodicOrbitProblem, NewtonKantorovich, solution

pb  = PeriodicOrbitProblem(f, Df, params_i, ξ_i, Xs)
x   = [unpack(x_bar)[1:3] ; real(unpack(x_bar)[4][1])]      # the blocks, τ real
alg = NewtonKantorovich{1}(ProjectedInverse(pb, x))
cert = validate(pb, x, alg)

println(cert)
println("template radii : ", inf.(r))
println("framework radii: ", radius(cert))
println("solution: ", typeof.(solution(cert)))
