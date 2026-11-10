# SteadyStateProblem
# x = (u...): F(u) = L u + N(u), DF(u) = L + M(DN(u))

struct SteadyStateProblem{TL,TLi,FN,FDN,P,BX} <: AbstractSequenceProblem
    L    :: TL    # linear part, a diagonal operator (e.g. Laplacian())
    Linv :: TLi   # its pseudo-inverse, nonincreasing in the mode beyond the projection
    N    :: FN    # nonlinearity on the blocks, and its derivative
    DN   :: FDN
    p    :: P
    X    :: BX
end

F(pb::SteadyStateProblem, x) = Ref(pb.L) .* x .+ pb.N(x, pb.p)

function DF(pb::SteadyStateProblem, x)
    Lblocks = [i == j ? pb.L : zero(exact(I)) for i ∈ eachindex(x), j ∈ eachindex(x)]
    return Lblocks + Multiplication.(pb.DN(x, pb.p))
end

# B = I
multipliers(pb::SteadyStateProblem, x) = pb.DN(x, pb.p)
prefactor(::SteadyStateProblem, x) = interval(1)
_with_parameters(pb::SteadyStateProblem, p) = SteadyStateProblem(pb.L, pb.Linv, pb.N, pb.DN, p, pb.X)
tail(pb::SteadyStateProblem, s::SequenceSpace) = pb.Linv * (interval(I) - interval(Projection(s)))
tail_bound(pb::SteadyStateProblem, s::SequenceSpace) = _linv_bound(pb, s, order(s))
