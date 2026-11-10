using RadiiPolynomial, LinearAlgebra
using Capsule: Capsule, validate, solution, FiniteInverse, NewtonKantorovich, FiniteDimensionalProblem, DF, InfiniteLinearOperator
setdisplay(:infsup; sigdigits = 6)

# ---- 1. finite-dimensional: the Shimizu-Morioka equilibrium ------------------

include(joinpath(@__DIR__, "..", "..", "shimizu-morioka-heteroclinic-orbit", "checker", "src", "system.jl"))
c1_bar = let v = Vector{Float64}(undef, 3)
    open(io -> read!(io, v), joinpath(@__DIR__, "..", "..", "shimizu-morioka-heteroclinic-orbit", "data", "c1_bar.bin"))
    ltoh.(v)
end

pb   = FiniteDimensionalProblem(f, Df, (I"0.75", I"0.45"))
alg  = NewtonKantorovich{1}(FiniteInverse(inv(mid.(DF(pb, c1_bar)))))
cert = validate(pb, c1_bar, alg)           # try: validate(pb, c1_bar, alg; R = 1e-10)
display(cert)
println("c1 ∈ ", solution(cert))

# ---- 2. operators with an error ---------------------------------------------

X = Ell1(GeometricWeight(interval(1.1)))
a = InfiniteSequence(Sequence(Taylor(2), interval.([1.0, 0.5, 0.25])), X; total_error = interval(0.1))

M = Multiplication(a)                     # Multiplication(ā) with error 0.1
P = interval(I) - interval(2) * M         # lazy, error 0.2
Q = P * Projection(Taylor(4))             # dense on Taylor(4), error kept
println("Q: ", typeof(Q.operator).name.name, ", error ", Q.error)
println("opnorm(Q)     = ", opnorm(Q, X))
println("opnorm(known) = ", opnorm(Q.operator, X))

b = Sequence(Taylor(1), interval.([1.0, -1.0]))
println("P * b has total error ", total_error(P * b))

# the first-order Z of a toy problem F(u) = u - u*u/4 on the ball ‖u - ū‖ ≤ R:
# DF(u) = I - M(u)/2 evaluated on the ball, A = A_K Π + (I - Π) as in the manifold step
R  = 0.01
u_ball = InfiniteSequence(Sequence(Taylor(2), interval.([0.3, 0.1, 0.05])), X; total_error = interval(R))
B  = interval(I) - interval(0.5) * Multiplication(u_ball)      # DF on the ball, error R/2
Π  = Projection(Taylor(2))
A_K = inv(mid.(coefficients(Π * B.operator * Π)))
A  = LinearOperator(Taylor(2), Taylor(2), interval.(A_K)) + (interval(I) - interval(Π))
Z_finite = opnorm(interval(Π) - A * (B * Π), X)                 # columns ≤ 2, exact for Taylor
Z_tail   = norm(sequence(u_ball), X) / interval(2) + B.error    # columns > 2: ‖M(u)/2‖ on the ball
println("Z = max(", Z_finite, ", ", Z_tail, ")")
