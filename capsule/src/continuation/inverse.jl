# The inverse along a path of parameters: a pointwise inverse whose entries are series in t. Its products go node by
# node on a grid in t resolving the product, and its block norms are those of series.

struct ContinuationInverse{A<:AbstractInverse,BX} <: AbstractInverse
    pointwise :: A     # the pointwise inverse, its entries series in t
    X         :: BX    # Banach space in t
    dims      :: Int   # number of continuation variables
end

#- a matrix of series

ContinuationInverse(A::FiniteInverse, X) = ContinuationInverse(A, X, length(grid_size(space(first(A.matrix)))))

function Base.:*(A::ContinuationInverse{<:FiniteInverse}, f::AbstractVecOrMat{<:Series})
    g = _sequence.(f)
    s = reduce(∪, [codomain(*, space(a), space(b)) for a ∈ A.pointwise.matrix, b ∈ g])
    m = grid_size(s)
    A_grid, g_grid = to_grid.(A.pointwise.matrix, Ref(m)), to_grid.(g, Ref(m))
    return [Series(_real(to_coef(sum(A_grid[i,k] .* g_grid[k,j] for k ∈ axes(g, 1)), s), eltype(eltype(g))), A.dims) for i ∈ axes(A_grid, 1), j ∈ axes(g, 2)]
end
Base.:*(A::ContinuationInverse{<:FiniteInverse}, f::AbstractVector{<:Series}) = vec(A * reshape(f, :, 1))
_real(a, ::Type{<:Real}) = real(a)
_real(a, ::Type) = a

opnorms(A::ContinuationInverse{<:FiniteInverse}) = norm.(A.pointwise.matrix, A.X)

#- A(t) = A_K(t) Π + T(t)

# the pointwise inverses of the nodes, interpolated; `s...` goes to the pointwise constructor
function ContinuationInverse(cpb::ContinuationProblem{<:AbstractSequenceProblem}, x, s...)
    pb, D = cpb.pointwise, cpb.dims
    T = reduce(∪, _lead.(space.(x), D))
    m = grid_size(T)
    x_grid, p_grid = _block_values.(x, Ref(m)), _nodes(pb.p, m)
    A_grid = [ProjectedInverse(_with_parameters(pb, p_grid[j]), [x_grid[i][j] for i ∈ eachindex(x)], s...) for j ∈ CartesianIndices(m)]
    return ContinuationInverse(_interpolate(pb, A_grid, T, _X(cpb, x)), cpb.X, D)
end
_block_values(a::Sequence, m) = _inner(space(a), length(m)) isa ScalarSpace ? real.(to_grid(a, m)) : to_grid(a, m)

function _interpolate(pb, A_grid, T, X)
    A = first(A_grid)
    finite = _to_coef([A.finite for A ∈ A_grid], T)
    isnothing(A.W) && return ProjectedInverse(finite, A.T, A.tail_bound, pb.X)
    # a quasilinear tail M(W(t)) L⁻¹ (I - Π): W is interpolated, the tails L⁻¹ (I - Π) of the blocks are kept
    s, kind = block_spaces(A), _kind(pb)
    q = _sequence_blocks(s)
    W = [to_coef([A.W[i,j] for A ∈ A_grid], T) for i ∈ axes(A.W, 1), j ∈ axes(A.W, 2)]
    tails = [sₗ isa ScalarSpace ? nothing : kind.Linv * (interval(I) - interval(Projection(sₗ))) for sₗ ∈ s]
    tail_bound = fill(interval(0), size(A.tail_bound))
    tail_bound[q,q] = norm.(W, X[q]) .* permutedims([_linv_bound(kind, sₗ, order(sₗ)) for sₗ ∈ s[q]])
    return ProjectedInverse(finite, tails, tail_bound, pb.X, W)
end
_to_coef(grid::AbstractArray{<:LinearOperator}, T) = to_coef(grid, T)
_to_coef(grid::AbstractArray{<:AbstractMatrix}, T) = [to_coef([F[i,j] for F ∈ grid], T) for i ∈ axes(first(grid), 1), j ∈ axes(first(grid), 2)]
_to_grid(F::LinearOperator, m) = to_grid(F, m)
_to_grid(F::AbstractMatrix{<:LinearOperator}, m) = (G = to_grid.(F, Ref(m)); [[_as_operator(G[i,j][l]) for i ∈ axes(F, 1), j ∈ axes(F, 2)] for l ∈ CartesianIndices(m)])
_to_grid(F::AbstractMatrix{<:Sequence}, m) = (G = to_grid.(F, Ref(m)); [[G[i,j][l] for i ∈ axes(F, 1), j ∈ axes(F, 2)] for l ∈ CartesianIndices(m)])

# the pointwise inverses at the nodes of a grid in t
function _nodes(P::ProjectedInverse, m)
    finite = _to_grid(P.finite, m)
    isnothing(P.W) && return [ProjectedInverse(finite[j], P.T, P.tail_bound, P.X) for j ∈ CartesianIndices(m)]
    W, q = _to_grid(P.W, m), findall(!isnothing, P.T)
    return [ProjectedInverse(finite[j], _tails(W[j], P.T, q), P.tail_bound, P.X, W[j]) for j ∈ CartesianIndices(m)]
end
function _tails(W, tails, q)
    T = Matrix{Any}(nothing, length(tails), length(tails))
    for (a, k) ∈ enumerate(q), (b, l) ∈ enumerate(q)
        T[k,l] = Multiplication(W[a,b]) * tails[l]
    end
    return T
end

_blocks(A::ContinuationInverse{<:ProjectedInverse}) = _unpack(A.pointwise.finite)
_path_space(A::ContinuationInverse{<:ProjectedInverse}) = _lead(codomain(first(_blocks(A))), A.dims)
block_spaces(A::ContinuationInverse{<:ProjectedInverse}) = [_inner(domain(B), A.dims) for B ∈ _blocks(A)[1,:]]
block_banachspaces(A::ContinuationInverse{<:ProjectedInverse}) = _lift.(Ref(A.X), A.pointwise.X, block_spaces(A), A.dims)
tail_bounds(A::ContinuationInverse{<:ProjectedInverse}) = A.pointwise.tail_bound
tail_multipliers(A::ContinuationInverse{<:ProjectedInverse}) = Series.(A.pointwise.W, A.dims)

function Base.:*(A::ContinuationInverse{<:ProjectedInverse}, f::AbstractVector{<:Series})
    T = codomain(*, _path_space(A), reduce(∪, _path_space.(f)))
    m = grid_size(T)
    A_grid, f_grid = _nodes(A.pointwise, m), [_node_values(fᵢ.a, m) for fᵢ ∈ f]
    Af_grid = [A_grid[j] * [f_grid[i][j] for i ∈ eachindex(f)] for j ∈ CartesianIndices(m)]
    return [Series(to_coef([Af[i] for Af ∈ Af_grid], T), A.dims) for i ∈ eachindex(f)]
end

# the known parts at the nodes, the errors through the block norms of A
function Base.:*(A::ContinuationInverse{<:ProjectedInverse}, C::AbstractMatrix)
    T = codomain(*, _path_space(A), reduce(∪, [_series_space(c) for c ∈ C if c isa SeriesOperator]))
    m = grid_size(T)
    A_grid, C_grid = _nodes(A.pointwise, m), _at.(C, Ref(m))
    AC_grid = [_apply(A_grid[j], [C_grid[k,i][j] for k ∈ axes(C, 1), i ∈ axes(C, 2)]) for j ∈ CartesianIndices(m)]
    errors, X = opnorms(A) * _series_error.(C), block_banachspaces(A)
    return [_dense(to_coef([AC[k,i] for AC ∈ AC_grid], T), T, errors[k,i], X[i], X[k]) for k ∈ axes(AC_grid[1], 1), i ∈ axes(C, 2)]
end
_dense(coefficients, T, error, X, Y) = SeriesOperator(m -> _as_operator.(to_grid(coefficients, m)), T, error, X, Y, nothing, nothing)
# the values at the nodes of an operator on series of numbers come back as sequences, or numbers
_as_operator(a::LinearOperator) = a
_as_operator(a::Sequence) = LinearOperator(a)
_as_operator(a::Number) = LinearOperator(ScalarSpace(), ScalarSpace(), [a;;])

function opnorms(A::ContinuationInverse{<:ProjectedInverse})
    X = block_banachspaces(A)
    return max.(opnorm.(_blocks(A), permutedims(X), X), _tail_norms(A.pointwise))
end
head_norms(A::ContinuationInverse{<:ProjectedInverse}) = opnorms(A)

# block norms of A L
function opnorms(A::ContinuationInverse{<:ProjectedInverse}, L)
    X, T, s = block_banachspaces(A), _path_space(A), block_spaces(A)
    q = _sequence_blocks(s)
    AL_grid = [[(F[k,l] * L) * Projection(domain(F[k,l])) for k ∈ axes(F, 1), l ∈ q] for F ∈ _to_grid(A.pointwise.finite, grid_size(T))]
    N = fill(interval(0), size(A.pointwise.finite))
    N[:,q] = opnorm.(_to_coef(AL_grid, T), permutedims(X[q]), X)
    N[q,q] = max.(N[q,q], norm.(A.pointwise.W, X[q]))
    return N
end
