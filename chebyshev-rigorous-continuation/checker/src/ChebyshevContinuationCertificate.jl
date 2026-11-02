module ChebyshevContinuationCertificate

# The four examples of Breden–Hénot, "Efficient Rigorous Continuation via
# Chebyshev Series Expansion I" (arXiv:2609.22689), and the I/O of their frozen
# inputs. Everything else (Chebyshev interpolation, products, norms, the radii
# polynomial) is RadiiPolynomial. No plotting dependencies.

using RadiiPolynomial

include("cube_root.jl")
include("square_root.jl")
include("cahn_hilliard.jl")
include("skt.jl")

export CubeRoot, SquareRoot, CahnHilliard, SKT,
    write_floats, read_floats,
    write_cube_root, read_cube_root,
    write_square_root, read_square_root,
    write_cahn_hilliard, read_cahn_hilliard,
    write_skt, read_skt

# ------------------------------------------------------------------
# Inputs: raw little-endian IEEE binary64, no header, no text round-trip.
#   cube_root.bin              x_j, j = 0..8                              (9 values)
#   square_root.bin            (x_j, λ_j), j = 0..32, then (ẋ_j, λ̇_j)       (2·33 + 2·33 values)
#   cahn_hilliard_<case>.bin   v_{j₁,j₂,k}, k fastest, then j₁, then j₂   (21·17·17 values)
#   skt_<arc>.bin              u_j (163 values each), j = 0..N, then u̇_j   (2·163·(N+1) values)
# ------------------------------------------------------------------

ENDIAN_BOM == 0x04030201 || error("big-endian hosts are not supported by this reader")

write_floats(path, x::AbstractVector{Float64}) = open(io -> write(io, x), path, "w")

function read_floats(path, n::Int)
    filesize(path) == 8n || error("$(basename(path)) has $(filesize(path)) bytes, expected $(8n)")
    return open(io -> read!(io, Vector{Float64}(undef, n)), path, "r")
end

# cube root
write_cube_root(path, x_grid::Vector{Float64}) = write_floats(path, x_grid)
read_cube_root(path) = read_floats(path, CubeRoot.N + 1)

# square root
write_square_root(path, u_grid, u̇_grid) = write_floats(path, [reduce(vcat, u_grid) ; reduce(vcat, u̇_grid)])
function read_square_root(path)
    n = SquareRoot.N + 1
    x = read_floats(path, 4n)
    u_grid = [x[2j-1:2j] for j ∈ 1:n]
    u̇_grid = [x[2n+2j-1:2n+2j] for j ∈ 1:n]
    return u_grid, u̇_grid
end

# Cahn–Hilliard
write_cahn_hilliard(path, v_grid) = write_floats(path, reduce(vcat, [collect(Float64, coefficients(v)) for v ∈ v_grid]))
function read_cahn_hilliard(path)
    s = evensym(Fourier(CahnHilliard.K, π))
    d = CahnHilliard.K + 1
    n₁, n₂ = CahnHilliard.N₁ + 1, CahnHilliard.N₂ + 1
    x = reshape(read_floats(path, d * n₁ * n₂), d, n₁, n₂)
    v_grid = Matrix{Sequence}(undef, n₁, n₂)
    for j₂ ∈ 1:n₂, j₁ ∈ 1:n₁
        v_grid[j₁,j₂] = Sequence(s, x[:,j₁,j₂])
    end
    return v_grid
end

# SKT
write_skt(path, u_grid, u̇_grid) = write_floats(path,
    [reduce(vcat, [collect(Float64, coefficients(u)) for u ∈ u_grid]) ;
     reduce(vcat, [collect(Float64, coefficients(u)) for u ∈ u̇_grid])])
function read_skt(path, name)
    X = SKT.X_u()
    d = dimension(X)
    n = SKT.arc_N(name) + 1
    x = reshape(read_floats(path, 2d * n), d, 2n)
    u_grid = Vector{Sequence}(undef, n)
    u̇_grid = Vector{Sequence}(undef, n)
    for j ∈ 1:n
        u_grid[j] = Sequence(X, x[:,j])
        u̇_grid[j] = Sequence(X, x[:,n+j])
    end
    return u_grid, u̇_grid
end

end
