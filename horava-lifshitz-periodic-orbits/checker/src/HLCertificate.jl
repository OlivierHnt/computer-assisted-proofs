module HLCertificate

# The map F, its derivative, the bounds Y, Z₁, Z₂, the radii polynomial and
# the checks (ii) of the original HLPeriodicOrbits_VIII_IX.jl (src/proof.jl),
# plus the raw-binary reader of the frozen inputs. Nothing here searches:
# the only computed auxiliary object is the approximate inverse
# A = inv(mid(π²ⁿ DF(x̄) π²ⁿ)), as in the original.

using RadiiPolynomial

include("proof.jl")
    export prove_case

include("witness_io.jl")
    export read_case, write_case, expected_filesize

end
