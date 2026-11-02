# Proof functions of the checker: the zero-finding maps, the Newton–Kantorovich
# bounds and the a posteriori checks. From OlivierHnt/BHWPeriodicCoexistence.jl
# (commit 6effd7d), without the numerical continuation, which lives in search/.
module CompetingPredatorsCertificate

    using RadiiPolynomial
    using LinearAlgebra

include("zero_finding_problems.jl") # F, DF for the periodic orbits and the Floquet normal form
include("checks.jl")                # grid2cheb, cheb2grid, prove_interior, prove_crossing, prove_inside_box
include("contraction.jl")           # prove_contraction = prove_po + prove_lin

    export prove_contraction, prove_interior, prove_crossing, prove_inside_box, sign_derivative

end
