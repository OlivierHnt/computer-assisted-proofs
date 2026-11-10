# Continuation tier: any pointwise problem along a path of parameters p(t), t the continuation variable(s).
# A number becomes a series in t, a sequence in a space S a sequence in T ⊗ S: t comes first. The problem is of the
# kind of its pointwise problem: F, DF, the statements and the bounds are the pointwise ones, run on `Series`.

struct ContinuationProblem{P,BX} <: AbstractZeroFindingProblem
    pointwise :: P     # the pointwise problem, its parameters series in t
    X         :: BX    # Banach space in t
    dims      :: Int   # number of continuation variables
end

ContinuationProblem(pb::AbstractZeroFindingProblem, p, X::BanachSpace; dims = _path_dims(p)) =
    ContinuationProblem(_with_parameters(pb, p), X, dims)
_path_dims(p::Sequence) = length(grid_size(space(p)))
_path_dims(p::Tuple) = _path_dims(first(Iterators.filter(q -> q isa Sequence, p)))

# the path t ∈ [-1, 1] ↦ [a, b]
segment(::Type{Chebyshev}, a, b) = Sequence(Chebyshev(1), [(a + b) / exact(2), (b - a) / exact(4)])   # T₁ carries the weight 2
segment(::Type{Taylor}, a, b) = Sequence(Taylor(1), [(a + b) / exact(2), (b - a) / exact(2)])

_kind(cpb::ContinuationProblem) = _kind(cpb.pointwise)
_pointwise(cpb::ContinuationProblem) = cpb.pointwise
_ring(cpb::ContinuationProblem) = _with_parameters(cpb.pointwise, _series(cpb.pointwise.p, cpb.dims))

F(cpb::ContinuationProblem, x) = F(_ring(cpb), _series.(x, cpb.dims))
DF(cpb::ContinuationProblem, x) = DF(_ring(cpb), _series.(x, cpb.dims))
multipliers(cpb::ContinuationProblem, x) = multipliers(_ring(cpb), _series.(x, cpb.dims))
prefactor(cpb::ContinuationProblem, x) = prefactor(_ring(cpb), _series.(x, cpb.dims))
leading(cpb::ContinuationProblem, x) = leading(_ring(cpb), _series.(x, cpb.dims))

_sequence(a::Series) = a.a
_sequence(a::Union{Sequence,InfiniteSequence}) = a

point(cpb::ContinuationProblem, x) = _series.(interval.(_sequence.(x)), cpb.dims)
# F is evaluated on balls of radius 0: a non-polynomial f then returns (sequence + error) itself
point(cpb::ContinuationProblem{<:FiniteDimensionalProblem}, x) = ball(cpb, x, zeros(Int, length(x)))
ball(cpb::ContinuationProblem, x, R) = _series.(_ball.(_sequence.(x), R, _X(cpb, x)), cpb.dims)

# the norm of a block: that of the path for a series of numbers, with the weights of the block after it otherwise
function _X(cpb::ContinuationProblem, x)
    X = _X(cpb.pointwise, x)
    return [_lift(cpb.X, isnothing(X) ? nothing : X[i], _inner(space(_sequence(x[i])), cpb.dims), cpb.dims) for i ∈ eachindex(x)]
end
_lift(Xt, Xs, ::ScalarSpace, D) = Xt
_lift(Xt::Ell1{IdentityWeight}, Xs::Ell1{IdentityWeight}, ::SequenceSpace, D) = Ell1()
_lift(Xt::Ell1{<:Weight}, Xs::Ell1{<:Weight}, ::SequenceSpace, D) = Ell1((ntuple(_ -> weight(Xt), D)..., weight(Xs)))

# the parameters at the nodes of a grid in t
_nodes(p::Sequence, m) = _inner(space(p), length(m)) isa ScalarSpace ? real.(to_grid(p, m)) : to_grid(p, m)
_nodes(p::Number, m) = fill(p, m)
_nodes(p::Tuple, m) = map(tuple, _nodes.(p, Ref(m))...)
_nodes(p::AbstractVector, m) = (g = _nodes.(p, Ref(m)); [[gᵢ[j] for gᵢ ∈ g] for j ∈ CartesianIndices(m)])
_nodes(p::NamedTuple, m) = (g = map(q -> _nodes(q, m), p); [map(gᵢ -> gᵢ[j], g) for j ∈ CartesianIndices(m)])

# what the prefactor B is when it is a series of numbers
_size(B::Series) = norm(B.a, banachspace(B.a))
_beyond(B::Series, X, s) = _size(B)
_opnorm(B::Series, X, ::Any) = _size(B)
