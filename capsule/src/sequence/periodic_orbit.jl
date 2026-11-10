# PeriodicOrbitProblem
# x = (u..., τ): the components of the orbit and the period scaling

struct PeriodicOrbitProblem{Ff,FDf,P,Ph,BX} <: AbstractSequenceProblem
    f     :: Ff
    Df    :: FDf
    p     :: P
    phase :: Ph
    X     :: BX
end

function F(pb::PeriodicOrbitProblem, x)
    u, τ = x[1:end-1], x[end]
    return [Derivative(1) .* u .- τ .* pb.f(u, pb.p) ; sum(adjoint(pb.phase[i]) * u[i] for i ∈ eachindex(u))]
end

function DF(pb::PeriodicOrbitProblem, x)
    u, τ = x[1:end-1], x[end]
    L = [i == j ? Derivative(1) : zero(exact(I)) for i ∈ eachindex(u), j ∈ eachindex(u)]
    return [L - τ * Multiplication.(pb.Df(u, pb.p))    LinearOperator.(-pb.f(u, pb.p))
            permutedims(adjoint.(pb.phase))            LinearOperator(exact(0))]
end

F(pb::PeriodicOrbitProblem, x::Sequence) = F(pb, _blocks(x))
DF(pb::PeriodicOrbitProblem, x::Sequence) = DF(pb, _blocks(x))
_blocks(x::Sequence) = (b = unpack(x); [b[1:end-1] ; b[end][1]])

# L = ∂, B = τ, N = f; the phase functionals vanish beyond their order
function multipliers(pb::PeriodicOrbitProblem, x)
    n = length(x) - 1
    return [pb.Df(x[1:n], pb.p)  fill(exact(0), n) ; permutedims(pb.phase)  exact(0)]
end
prefactor(::PeriodicOrbitProblem, x) = x[end]
_with_parameters(pb::PeriodicOrbitProblem, p) = PeriodicOrbitProblem(pb.f, pb.Df, p, pb.phase, pb.X)
tail(::PeriodicOrbitProblem, s::Fourier) = (interval(I) - interval(Projection(s))) * Integral(1)
tail_bound(::PeriodicOrbitProblem, s::Fourier) = inv(interval(order(s) + 1) * frequency(s))
