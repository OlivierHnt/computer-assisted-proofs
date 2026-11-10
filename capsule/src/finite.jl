# Finite-dimensional tier: numbers, abs, a matrix inverse

struct FiniteInverse{T} <: AbstractInverse
    matrix :: T
end

Base.:*(A::FiniteInverse, f::AbstractVecOrMat) = interval(A.matrix) * f
opnorms(A::FiniteInverse) = abs.(interval(A.matrix))

#- FiniteDimensionalProblem

struct FiniteDimensionalProblem{FF,FDF,P} <: AbstractZeroFindingProblem
    F  :: FF
    DF :: FDF
    p  :: P   # parameters
end

F(pb::FiniteDimensionalProblem, x) = pb.F(x, pb.p)
DF(pb::FiniteDimensionalProblem, x) = pb.DF(x, pb.p)
_X(::FiniteDimensionalProblem, x) = nothing
_with_parameters(pb::FiniteDimensionalProblem, p) = FiniteDimensionalProblem(pb.F, pb.DF, p)
