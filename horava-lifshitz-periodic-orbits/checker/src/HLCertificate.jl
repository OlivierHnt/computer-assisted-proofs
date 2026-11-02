module HLCertificate

# The map F, its derivative, the bounds Y, Z₁, Z₂, the radii polynomial and
# the checks of Lemma 4.2 (proof.jl), plus the reader and writer of the
# frozen inputs (witness_io.jl). Nothing here searches: the only computed
# auxiliary object is the approximate inverse A = inv(mid(π²ⁿ DF(x̄) π²ⁿ)).

using RadiiPolynomial

include("proof.jl")
    export prove_case

include("witness_io.jl")
    export read_case, write_case, expected_filesize

end
