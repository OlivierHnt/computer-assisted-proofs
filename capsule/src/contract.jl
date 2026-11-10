# The contract. A problem gives F(pb, x) and DF(pb, x), which return (known part + error) entries. An inverse gives
# A * f and A * C on known parts, and opnorms(A), its block norms, which the errors go through.
# The two bounds below are written once and dispatch on the inverse.

abstract type AbstractZeroFindingProblem end

# x is the vector of blocks (sequences, and scalars as numbers); pb.X the vector of their Banach spaces.
# F(u) = L u + B N(u), with L invertible and B N relatively compact with respect to L: the tail of A is L⁻¹ (I − Π)
abstract type AbstractSequenceProblem <: AbstractZeroFindingProblem end

abstract type AbstractInverse end

abstract type AbstractCAPAlgorithm end

struct NewtonKantorovich{N,Inv<:AbstractInverse} <: AbstractCAPAlgorithm
    inverse :: Inv
    function NewtonKantorovich{N}(inverse::Inv) where {N,Inv}
        N ∈ (1, 2) || throw(ArgumentError("order must be 1 or 2"))
        return new{N,Inv}(inverse)
    end
end

#- what a problem states

function F end
function DF end

# for a sequence problem: DN, B and L⁻¹ beyond the projection
function multipliers end   # the entries of DN(u), multiplication operators (numbers where there is none)
function prefactor end     # B
function leading end       # the entries of DΦ(u), for a nonlinearity Φ under L
tail(::AbstractSequenceProblem, ::ScalarSpace) = nothing        # L⁻¹ (I − Π) on a block
tail_bound(::AbstractSequenceProblem, ::ScalarSpace) = interval(0)   # and its norm

#- what an inverse gives, besides its products

function opnorms end

#- the proof

function validate(pb::AbstractZeroFindingProblem, x, alg::NewtonKantorovich{1}; R=nothing)
    Y = Y_bound(pb, x, alg.inverse)
    R_ = isnothing(R) ? fill(exact(10maximum(sup.(Y))), length(Y)) : R
    Z = Z_bound(pb, x, alg.inverse, R_)
    W = zeros(eltype(Z), length(Y), length(Y), length(Y))
    return Certificate(pb, alg, x, Y, Z, W, R_)
end

# one bound per component: the known part through A, the error through the block norms of A
function Y_bound(pb::AbstractZeroFindingProblem, x, A::AbstractInverse)
    Fx = F(pb, point(pb, x))
    return _norm.(A * _known.(Fx), _X(pb, x)) + opnorms(A) * _error.(Fx)
end

# one bound per pair of components, on the product of the component balls
# The method is chosen by the kind of the problem: a problem along a path of parameters is of the kind of its
# pointwise problem, and runs the same method on series
Z_bound(pb::AbstractZeroFindingProblem, x, A::AbstractInverse, R) = _Z_bound(_kind(pb), pb, x, A, R)
_kind(pb::AbstractZeroFindingProblem) = pb
_pointwise(pb::AbstractZeroFindingProblem) = pb

function _Z_bound(::AbstractZeroFindingProblem, pb, x, A, R)
    C = DF(pb, ball(pb, x, R))
    return _norm.(_identity_minus(A * _known.(C)), _X(pb, x)) + opnorms(A) * _error.(C)
end
_identity_minus(M::AbstractMatrix) = [i == j ? exact(1) - M[i,j] : -M[i,j] for i ∈ axes(M, 1), j ∈ axes(M, 2)]

point(::AbstractZeroFindingProblem, x) = interval.(x)
ball(pb::AbstractZeroFindingProblem, x, R) = _ball.(x, R, _X(pb, x))
_ball(a::Number, r, _) = interval(a, r; format=:midpoint)
_ball(a::Sequence, r, X) = InfiniteSequence(interval(a), X; total_error = interval(r))
_X(pb::AbstractZeroFindingProblem, x) = pb.X   # the Banach spaces of the blocks
_norm(a::Number, _) = abs(a)
_norm(a, X) = norm(a, X)
