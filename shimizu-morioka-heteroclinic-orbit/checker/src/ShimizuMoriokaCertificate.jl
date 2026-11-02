# The verification part of ShimizuMoriokaTutorial.jl (O. Hénot), with the
# searches (Newton, LAPACK eigendecomposition, ODE solve) removed: see search/.

module ShimizuMoriokaCertificate

using RadiiPolynomial
using LinearAlgebra

#

include("system.jl")
    export f, Df

include("equilibrium.jl")
    export EquilibriumProblem

include("eigen.jl")
    export validated_eigen

include("manifold.jl")
    export ManifoldProblem

include("heteroclinic.jl")
    export HeteroclinicProblem

#

    export validate

end
