module Capsule

using RadiiPolynomial
using LinearAlgebra

export validate, NewtonKantorovich, Certificate,
    FiniteDimensionalProblem, PeriodicOrbitProblem, InitialValueProblem, SteadyStateProblem, QuasilinearProblem,
    PseudoArclengthProblem,
    ContinuationProblem, segment,
    FiniteInverse, ProjectedInverse, ContinuationInverse,
    InfiniteLinearOperator

# `F`, `DF` and `radius` are common names in a proof script: they are public and not exported
public F, DF, radius, solution, ReturnCode, Y_bound, Z_bound, opnorms, ball,
    multipliers, prefactor, leading, tail, tail_bound

include("radiipolynomial.jl")

include("contract.jl")
include("certificate.jl")

include("finite.jl")

include("sequence/projected_inverse.jl")
include("sequence/bounds.jl")
include("sequence/periodic_orbit.jl")
include("sequence/initial_value.jl")
include("sequence/steady_state.jl")
include("sequence/quasilinear.jl")
include("sequence/pseudo_arclength.jl")
include("sequence/not_written.jl")

include("continuation/series.jl")
include("continuation/series_operator.jl")
include("continuation/continuation.jl")
include("continuation/inverse.jl")

end
