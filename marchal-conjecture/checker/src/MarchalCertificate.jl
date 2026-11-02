module MarchalCertificate

# The P₁₂ symmetric spaces, the map F and its derivative, and the witness I/O.
# Everything else the proof needs (Chebyshev interpolation, products, norms,
# the radii polynomial) is RadiiPolynomial. No plotting dependencies.

using RadiiPolynomial

include("spaces.jl")
    export evencos, oddcos, evensin, oddsin, evenfourier, marchal_space,
        impose_reality!, norm_space, tail, tail!

include("marchal.jl")
    export F_marchal!, DF_marchal!

# ------------------------------------------------------------------
# Input I/O.
#   x_nodes.bin : for each node, the d complex interval coefficients of x̄ in
#                 the symmetric representation, stored as
#                 (inf re, sup re, inf im, sup im) — four binary64 each.
# Raw little-endian IEEE binary64 throughout; no text round-trip.
# The approximate inverses A on the nodes are not stored: they are cheap to
# recompute from x̄ (see approximate_inverse).
# ------------------------------------------------------------------

function write_witness(dir, x_nodes::Vector{<:Sequence})
    mkpath(dir)
    open(joinpath(dir, "x_nodes.bin"), "w") do io
        for x ∈ x_nodes, c ∈ coefficients(x)
            write(io, inf(real(c)), sup(real(c)), inf(imag(c)), sup(imag(c)))
        end
    end
    return nothing
end

function read_witness(dir, npts::Int, d::Int)
    x_nodes = open(joinpath(dir, "x_nodes.bin"), "r") do io
        map(1:npts) do _
            raw = read!(io, Matrix{Float64}(undef, 4, d))
            map(1:d) do i
                complex(interval(raw[1,i], raw[2,i]), interval(raw[3,i], raw[4,i]))
            end
        end
    end
    return x_nodes
end

# Floating-point inverse of DF(x, Ω) at the midpoint of an interval node value
# x; any matrix would do, the proof only needs it to be the same operator
# wherever it is used.
function approximate_inverse(x::Sequence, Ω::Float64, K::Int)
    X = marchal_space(K, 1.0)
    xf = Sequence(X, mid.(coefficients(x)))
    return Matrix{ComplexF64}(coefficients(inv(DF_marchal!(zeros(ComplexF64, X, X), xf, Ω))))
end

    export write_witness, read_witness, approximate_inverse

end
