module EllipticPDELocalizedRadialSol

# The equations, the zero-finding map F and its derivative DF, the proof
# (Proposition 2.2, Theorem 3.3, Corollary 3.5 of the paper), the reader of
# data/*.bin (io.jl) and the checks of the hypotheses (additions.jl).
# From OlivierHnt/EllipticPDELocalizedRadialSol.jl, commit c0fda16.

    using RadiiPolynomial
    import LinearAlgebra: Diagonal, diag

include("compat.jl")

include("klein_gordon.jl")
    export KleinGordon

include("swift_hohenberg.jl")
    export SwiftHohenberg

include("fitzhugh_nagumo.jl")
    export FitzHughNagumo

include("proof.jl")
    export prove

include("io.jl")
    export CASES, read_case

include("additions.jl")
    export check_parameters, check_symmetry, check_equilibrium, ψ̂_strict, check_LP_strict

end
