# SEARCH — not part of the proof.
#
# Runs the original pipeline of ShimizuMoriokaTutorial.jl (script/main.jl,
# without the plot): Newton for the equilibrium c₁, the LAPACK eigendecompositions
# at c₀ and c₁, Newton for the Taylor coefficients of the local stable manifold
# of c₀ and of the local unstable manifold of c₁, and an ODE solve followed by
# Newton for the Chebyshev coefficients of the connecting orbit. The
# validations of the original code are run too, because each approximation is
# computed from the validated objects of the previous step (mid points of the
# validated equilibrium, eigenpairs, manifolds), exactly as in the original.
#
# Writes the floating-point approximations to data/ (raw little-endian
# binary64, see statement.md). The checker does not care how they were
# obtained.
#
# Usage (from the package root; this is the `produce` command of entry.json):
#   julia --project=search search/produce.jl
#   julia --project=checker tools/hash_data.jl

import Pkg
Pkg.instantiate()

using RadiiPolynomial, LinearAlgebra
include(joinpath(@__DIR__, "src", "ShimizuMoriokaTutorial.jl")) # original code (verbatim)
using .ShimizuMoriokaTutorial
const SMT = ShimizuMoriokaTutorial

BLAS.set_num_threads(2)

datadir = joinpath(@__DIR__, "..", "data")

writebin(name, v) = open(io -> write(io, htol.(v)), joinpath(datadir, name), "w")
interleave(z::AbstractArray{<:Complex}) = vec(permutedims(hcat(real.(vec(z)), imag.(vec(z)))))

# Parameters
parameters = (I"0.75", I"0.45") # (a, b)

# Section 2: Equilibria and eigenspaces

c0_star = interval([0, 0, 0])

pb_c1 = EquilibriumProblem(f, Df, parameters)
c1_bar, ok = SMT.approximate_zero(EquilibriumProblem(f, Df, mid.(parameters)), [1.0, 0.0, 1.0]) # the Newton step of validated_solve
@assert ok
c1_star = validated_solve(pb_c1, [1.0, 0.0, 1.0])

# approximate eigenpairs: what `eigen` of IntervalArithmetic computes from the
# mid point matrix before validating it
λ0_approx, V0_approx = eigen!(mid.(Df(c0_star, parameters)))
λ1_approx, V1_approx = eigen!(mid.(Df(c1_star, parameters)))

λ0, V0 = eigen(Df(c0_star, parameters))
λ1, V1 = eigen(Df(c1_star, parameters))

# Section 3

λ0_stable = λ0[1:2];
V0_stable = exact(0.5) * V0[:,1:2];
pb_stable_manifold = ManifoldProblem(f, Df, parameters, c0_star, λ0_stable, V0_stable);
ν_Q = interval(1.1)
Q_bar, r_Q = validated_solve(pb_stable_manifold, 15, ν_Q);

λ1_unstable = λ1[2:3];
V1_unstable = exact(0.15) * V1[:,2:3];
pb_unstable_manifold = ManifoldProblem(f, Df, parameters, c1_star, λ1_unstable, V1_unstable);
ν_P = interval(1.1)
P_bar, r_P = validated_solve(pb_unstable_manifold, 22, ν_P);

# Section 4

pb_het = HeteroclinicProblem(f, Df, parameters,
    P_bar, [Derivative(1, 0) .* P_bar Derivative(0, 1) .* P_bar], r_P, ν_P,
    Q_bar, [Derivative(1, 0) .* Q_bar Derivative(0, 1) .* Q_bar], r_Q, ν_Q,
    interval(13));

x_bar, r = validated_solve(pb_het, 0.6, 125, interval(1.05));

# validated_solve returns the Newton outputs of the manifolds and of the orbit
# as thin intervals; their mid points are the floating-point approximations.

Q_float = reduce(vcat, [coefficients(mid.(c)) for c ∈ Q_bar]) # 3 × 16² Float64
P_float = reduce(vcat, [coefficients(mid.(c)) for c ∈ P_bar]) # 3 × 23² ComplexF64
x_float = reduce(vcat, [coefficients(mid.(c)) for c ∈ x_bar]) # 3 × 126 + 3 Float64
@assert all(isthin, reduce(vcat, coefficients.(Q_bar)))
@assert all(isthin, reduce(vcat, coefficients.(P_bar)))
@assert all(isthin, reduce(vcat, coefficients.(x_bar)))
@assert eltype(λ0_approx) == eltype(V0_approx) == Float64
@assert eltype(λ1_approx) == eltype(V1_approx) == ComplexF64
@assert eltype(Q_float) == eltype(x_float) == Float64 && eltype(P_float) == ComplexF64

mkpath(datadir)
writebin("c1_bar.bin", c1_bar)
writebin("eig_c0.bin", [λ0_approx; vec(V0_approx)])
writebin("eig_c1.bin", interleave([λ1_approx; vec(V1_approx)]))
writebin("Q_bar.bin", Q_float)
writebin("P_bar.bin", interleave(P_float))
writebin("x_bar.bin", x_float)
println("wrote ", join(sort(filter(endswith(".bin"), readdir(datadir))), ", "), " to ", normpath(datadir))
