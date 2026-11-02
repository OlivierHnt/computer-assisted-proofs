# The proof functions of the original repository (BHWPeriodicCoexistence.jl),
# without the numerical continuation (numerics.jl, OrdinaryDiffEq), which lives
# in search/.
module CompetingPredatorsCertificate

    using RadiiPolynomial
    using LinearAlgebra

include("zero_finding_problems.jl") # F, DF for the periodic orbits and the Floquet normal form
include("checks.jl")                # grid2cheb, cheb2grid, prove_interior, prove_crossing, prove_inside_box
include("contraction.jl")           # prove_contraction = prove_po + prove_lin

    export prove_contraction, prove_interior, prove_crossing, prove_inside_box, sign_derivative

end
