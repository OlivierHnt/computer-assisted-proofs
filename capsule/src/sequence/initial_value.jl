# InitialValueProblem
# x = (u...): the Taylor components of the solution; F(u) = u - u0 - ∫ f(u)

struct InitialValueProblem{Ff,FDf,P,U,BX} <: AbstractSequenceProblem
    f  :: Ff
    Df :: FDf
    p  :: P
    u0 :: U   # initial condition, one number per component
    X  :: BX
end

F(pb::InitialValueProblem, x) = x .- pb.u0 .- Integral(1) .* pb.f(x, pb.p)

function DF(pb::InitialValueProblem, x)
    Id = [i == j ? exact(I) : zero(exact(I)) for i ∈ eachindex(x), j ∈ eachindex(x)]
    return Id - Integral(1) .* Multiplication.(pb.Df(x, pb.p))
end

# L = I, B = ∫, N = f
multipliers(pb::InitialValueProblem, x) = pb.Df(x, pb.p)
prefactor(::InitialValueProblem, x) = Integral(1)
_with_parameters(pb::InitialValueProblem, p) = InitialValueProblem(pb.f, pb.Df, p, pb.u0, pb.X)
tail(::InitialValueProblem, s::SequenceSpace) = interval(I) - interval(Projection(s))
tail_bound(::InitialValueProblem, s::SequenceSpace) = interval(1)
