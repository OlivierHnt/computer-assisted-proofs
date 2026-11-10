# The inverse of the sequence tier: blocks of sequences, weighted ℓ¹ norms, A = A_K Π + T

struct ProjectedInverse{F,T,B,Xs,TW} <: AbstractInverse
    finite     :: F   # A_K on its projection space Π (an operator on the product space, or its matrix of blocks)
    T          :: T   # matrix of blocks (or nothing), the tail part: A = A_K Π + T
    tail_bound :: B   # bounds of ‖T Π'‖, Π' = I − Π, one per block (or per pair of blocks)
    X          :: Xs  # Banach spaces of the blocks
    W          :: TW  # the multipliers of a quasilinear tail M(W) L⁻¹ (I − Π), or nothing
end
ProjectedInverse(finite, T, tail_bound, X) = ProjectedInverse(finite, T, tail_bound, X, nothing)
_unpack(A::LinearOperator) = unpack(A)
_unpack(A::AbstractMatrix) = A

# A = A_K Π + T, with one tail block per sequence block: tail(pb, space) and tail_bound(pb, space) come from the problem
function ProjectedInverse(pb::AbstractSequenceProblem, x)
    x = interval.(x)
    Π = Projection(_product(_space.(x)))
    A_K = _approximate_inverse(Π * DF(pb, x) * Π)
    T = [i == j ? tail(pb, _space(x[i])) : nothing for i ∈ eachindex(x), j ∈ eachindex(x)]
    return ProjectedInverse(A_K, T, [tail_bound(pb, _space(xᵢ)) for xᵢ ∈ x], pb.X)
end
_product(spaces) = length(spaces) == 1 ? spaces[1]^1 : reduce(×, spaces)
# inverse of the midpoint, on the interval spaces themselves: `mid` would turn a non-thin frequency into a float
_approximate_inverse(M::LinearOperator) = LinearOperator(domain(M), codomain(M), interval(inv(mid.(coefficients(M)))))

_add(A, ::Nothing) = A
_add(A, T) = A + T
_lazy(A::ProjectedInverse) = _add.(_unpack(A.finite), A.T)   # A_K + T blockwise, nothing = no tail

# the error of a block is not banded: it meets all of A, through the block norms of A
Base.:*(A::ProjectedInverse, f::AbstractVector) = (_lazy(A) * Diagonal(Projection.(space.(f)))) * f
Base.:*(A::ProjectedInverse, C::AbstractMatrix) = _with_error.(_apply(A, _operator.(C)), opnorms(A) * _error.(C), permutedims(A.X), A.X)
_apply(A::ProjectedInverse, C::AbstractMatrix) = (_lazy(A) * Diagonal(Projection.(_codomain(C)))) * C
_codomain(C::AbstractMatrix) = [reduce(∪, codomain.(view(C, i, :))) for i ∈ axes(C, 1)]   # row by row
_with_error(L, e, X, Y) = InfiniteLinearOperator(L, e, X, Y)

# block norms of A = A_K Π + T: max of the two parts, exact for weighted ℓ¹ domains
opnorms(A::ProjectedInverse) = max.(opnorm.(_unpack(A.finite), permutedims(A.X), A.X), _tail_norms(A))
_tail_norms(A::ProjectedInverse) = isnothing(A.W) ? _tail_norm.(A.T, A.tail_bound) : A.tail_bound
_tail_norm(::Nothing, b) = zero(b)
_tail_norm(T, b) = b
# block norms of A Π_0, Π_0 the projection on the mode 0 (the tail of A does not see it)
function head_norms(A::ProjectedInverse)
    B = _unpack(A.finite)
    return [opnorm(B[m,l] * Projection(_head(domain(B[m,l]))), A.X[l], A.X[m]) for m ∈ axes(B, 1), l ∈ axes(B, 2)]
end
_head(s::ScalarSpace) = s
_head(::Taylor) = Taylor(0)
_head(::Chebyshev) = Chebyshev(0)
_head(s::Fourier) = Fourier(0, frequency(s))
