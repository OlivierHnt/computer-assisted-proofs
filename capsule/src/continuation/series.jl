# A block of the unknown along a path of parameters: a series in the continuation variable(s) t, whose coefficients
# are numbers (a sequence on T) or sequences of a space S (a sequence on T ⊗ S, t first). It is what a number or a
# sequence becomes: the pointwise F, DF and bounds run on it unchanged.

struct Series{A<:Union{Sequence,InfiniteSequence}}
    a    :: A
    dims :: Int   # number of continuation variables
end

Base.broadcastable(s::Series) = Ref(s)

_lead(s, D) = first(RadiiPolynomial._lead_inner(s, Val(D)))
_inner(s, D) = last(RadiiPolynomial._lead_inner(s, Val(D)))

_full_space(s::Series) = space(s.a)
_path_space(s::Series) = _lead(space(s.a), s.dims)
RadiiPolynomial.space(s::Series) = _inner(space(s.a), s.dims)   # ScalarSpace for a series of numbers
_space(s::Series) = space(s)
_isscalar(s::Series) = space(s) isa ScalarSpace

RadiiPolynomial.sequence(s::Series) = Series(sequence(s.a), s.dims)
RadiiPolynomial.total_error(s::Series) = total_error(s.a)
LinearAlgebra.norm(s::Series, X::BanachSpace) = norm(s.a, X)
_pad(s::Series, space_) = Series(project(s.a, (_isscalar(s) ? _lead(space_, s.dims) : space_) ∪ space(s.a)), s.dims)
_series(a::Series, D) = a
_series(a::Union{Sequence,InfiniteSequence}, D) = Series(a, D)
_series(a::Number, D) = a
_series(a::Tuple, D) = _series.(a, D)
_series(a::AbstractVector, D) = _series.(a, D)
_series(a::NamedTuple, D) = map(q -> _series(q, D), a)

#- the ring: a series of numbers is constant in the variables of S

_embed(a::Sequence, b, D) = Sequence(space(a) ⊗ RadiiPolynomial._zero_space(_inner(space(b), D)), coefficients(a))
_embed(a::InfiniteSequence, b::InfiniteSequence, D) =
    InfiniteSequence(_embed(sequence(a), b, D), banachspace(b); total_error = total_error(a))
_embed(a::InfiniteSequence, b::Sequence, D) =
    throw(ArgumentError("a ball of series of numbers cannot multiply a sequence without a Banach space"))
_pair(x::Series, y::Series) =
    _isscalar(x) == _isscalar(y) ? (x.a, y.a) : _isscalar(x) ? (_embed(x.a, y.a, x.dims), y.a) : (x.a, _embed(y.a, x.a, x.dims))

for f ∈ (:+, :-, :*)
    @eval Base.$f(x::Series, y::Series) = Series($f(_pair(x, y)...), x.dims)
    @eval Base.$f(x::Series, y::Number) = Series($f(x.a, y), x.dims)
    @eval Base.$f(x::Number, y::Series) = Series($f(x, y.a), y.dims)
end
Base.:+(x::Series) = x
Base.:-(x::Series) = Series(-x.a, x.dims)
Base.:/(x::Series, y::Number) = Series(x.a / y, x.dims)
Base.:^(x::Series, n::Integer) = Series(x.a ^ n, x.dims)
Base.zero(x::Series) = Series(zero(x.a), x.dims)
Base.one(x::Series) = Series(one(x.a), x.dims)
for f ∈ (:exp, :log, :sqrt, :inv, :sin, :cos)
    @eval Base.$f(x::Series) = Series($f(x.a), x.dims)
end
Base.:*(x::Series, M::AbstractArray) = x .* M
Base.:*(M::AbstractArray, x::Series) = M .* x

#- an operator of the problem acts along the variables of S, at the nodes of a grid in t

_node_values(a::Sequence, m) = _inner(space(a), length(m)) isa ScalarSpace ? [Sequence(ScalarSpace(), [v]) for v ∈ to_grid(a, m)] : to_grid(a, m)

function Base.:*(L::AbstractLinearOperator, s::Series)
    s.a isa Sequence || throw(ArgumentError("an operator of the problem acts on a series, not on a ball of series"))
    T = _path_space(s)
    return Series(to_coef(Ref(L) .* _node_values(s.a, grid_size(T)), T), s.dims)
end
