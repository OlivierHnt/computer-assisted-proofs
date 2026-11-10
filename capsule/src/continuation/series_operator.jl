# An operator along a path of parameters, with an error: what Multiplication of a ball becomes. Its known part is a
# polynomial in t with operator coefficients, given by its values at the nodes of any grid in t resolving `space`.

# (one concrete type, whatever the operations it comes from, so that it sits in a matrix)
struct SeriesOperator <: AbstractLinearOperator
    at    :: Function   # m -> the operators at the nodes of a grid of size m
    space               # space in t of the known part
    error               # a real number
    X                   # Banach spaces of the blocks it acts between (nothing: no error to carry)
    Y
    bound               # bound of the norm of the known part (nothing: unknown)
    basis               # the space of the multiplier behind it, for the norm of an operator in front
end

_at(P::SeriesOperator, m) = P.at(m)
_at(L::AbstractLinearOperator, m) = fill(L, m)
_operator(P::SeriesOperator) = SeriesOperator(P.at, P.space, zero(P.error), P.X, P.Y, P.bound, P.basis)
_error(P::SeriesOperator) = P.error
_series_error(P::SeriesOperator) = P.error
_series_error(::AbstractLinearOperator) = interval(0)
_series_space(P::SeriesOperator) = P.space

_banachspace(a::InfiniteSequence) = banachspace(a)
_banachspace(::Sequence) = nothing
_bound(a::InfiniteSequence) = norm(sequence(a), banachspace(a))
_bound(::Sequence) = nothing
# the Banach space of a series of numbers next to a block, and of a block along the variables of S only
_path_banachspace(X::Ell1{<:Tuple}, D) = D == 1 ? Ell1(weight(X)[1]) : Ell1(weight(X)[1:D])
_path_banachspace(X, D) = X
_block_banachspace(X::Ell1{<:Tuple}) = Ell1(last(weight(X)))
_block_banachspace(X) = X

function RadiiPolynomial.Multiplication(w::Series)
    w̄, X = sequence(w.a), _banachspace(w.a)
    return SeriesOperator(m -> Multiplication.(to_grid(w̄, m)), _lead(space(w̄), w.dims), _error(w.a), X, X, _bound(w.a), space(w))
end
# the block as an operator on a series of numbers, and as a functional
function RadiiPolynomial.LinearOperator(a::Series)
    ā, X = sequence(a.a), _banachspace(a.a)
    return SeriesOperator(m -> LinearOperator.(to_grid(ā, m)), _lead(space(ā), a.dims), _error(a.a), _path_banachspace(X, a.dims), X, _bound(a.a), nothing)
end
function Base.adjoint(a::Series)
    ā, X = sequence(a.a), _banachspace(a.a)
    return SeriesOperator(m -> adjoint.(to_grid(ā, m)), _lead(space(ā), a.dims), _error(a.a), X, _path_banachspace(X, a.dims), _bound(a.a), nothing)
end

_sum(a, b) = isnothing(a) || isnothing(b) ? nothing : a + b
_err(P::SeriesOperator, factor) = isnothing(P.X) ? P.error : factor() * P.error

Base.:-(P::SeriesOperator) = SeriesOperator(m -> .-(P.at(m)), P.space, P.error, P.X, P.Y, P.bound, P.basis)
for f ∈ (:+, :-)
    @eval Base.$f(P::SeriesOperator, Q::SeriesOperator) =
        SeriesOperator(m -> broadcast($f, P.at(m), Q.at(m)), P.space ∪ Q.space, P.error + Q.error, P.X, P.Y, _sum(P.bound, Q.bound), nothing)
    for TL ∈ (:AbstractLinearOperator, :LinearOperator, :Projection, :UniformScalingOperator)
        @eval Base.$f(L::$TL, P::SeriesOperator) = SeriesOperator(m -> broadcast($f, Ref(L), P.at(m)), P.space, P.error, P.X, P.Y, nothing, nothing)
        @eval Base.$f(P::SeriesOperator, L::$TL) = SeriesOperator(m -> broadcast($f, P.at(m), Ref(L)), P.space, P.error, P.X, P.Y, nothing, nothing)
    end
end
for TL ∈ (:AbstractLinearOperator, :LinearOperator, :Projection, :UniformScalingOperator)
    @eval Base.:*(L::$TL, P::SeriesOperator) =
        SeriesOperator(m -> Ref(L) .* P.at(m), P.space, _err(P, () -> _opnorm(L, _block_banachspace(P.Y), P.basis)), P.X, P.Y, nothing, nothing)
    @eval Base.:*(P::SeriesOperator, L::$TL) =
        SeriesOperator(m -> P.at(m) .* Ref(L), P.space, _err(P, () -> opnorm(L, _block_banachspace(P.X))), P.X, P.Y, nothing, nothing)
end
Base.:*(c::Number, P::SeriesOperator) = SeriesOperator(m -> c .* P.at(m), P.space, abs(c) * P.error, P.X, P.Y, isnothing(P.bound) ? nothing : abs(c) * P.bound, P.basis)
Base.:*(P::SeriesOperator, c::Number) = c * P

# a series of numbers in front: the known parts multiply at the nodes, the errors as in a Banach algebra
function Base.:*(τ::Series, P::SeriesOperator)
    _isscalar(τ) || throw(ArgumentError("only a series of numbers multiplies an operator"))
    τ̄ = sequence(τ.a)
    at, T = (m -> to_grid(τ̄, m) .* P.at(m)), codomain(*, space(τ̄), P.space)
    τ.a isa Sequence && return SeriesOperator(at, T, P.error, P.X, P.Y, nothing, P.basis)
    norm_τ = norm(τ̄, banachspace(τ.a))
    error = norm_τ * P.error + total_error(τ.a) * (P.bound + P.error)
    return SeriesOperator(at, T, error, P.X, P.Y, norm_τ * P.bound, P.basis)
end

# an operator and a block, both series in t
function Base.:*(P::SeriesOperator, s::Series)
    s.a isa Sequence || throw(ArgumentError("an operator of the problem acts on a series, not on a ball of series"))
    T = codomain(*, P.space, _path_space(s))
    m = grid_size(T)
    return Series(to_coef(P.at(m) .* _node_values(s.a, m), T), s.dims)
end

_coefficients(P::SeriesOperator) = to_coef(P.at(grid_size(P.space)), P.space)
LinearAlgebra.opnorm(P::SeriesOperator, X::BanachSpace, Y::BanachSpace) = opnorm(_coefficients(P), X, Y) + P.error
