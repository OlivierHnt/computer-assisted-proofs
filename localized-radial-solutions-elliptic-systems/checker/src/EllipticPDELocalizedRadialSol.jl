module EllipticPDELocalizedRadialSol

# The equations, the zero-finding map F and its derivative DF, and the proof
# (Proposition 2.2, Theorem 3.3, Corollary 3.5 of the paper), copied from
# EllipticPDELocalizedRadialSol.jl (commit c0fda16) with the JLD2 loading
# replaced by raw binary I/O (io.jl), ported to RadiiPolynomial 0.11 /
# IntervalArithmetic 1.0 (compat.jl). The module keeps its original name so
# that the original JLD2 files deserialise into these types (search/).

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
    export check_parameters, check_symmetry, check_equilibrium, check_eigendata, ψ̂_strict, check_LP_strict

end
