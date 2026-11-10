# The Z bound of a semilinear sequence problem: the columns up to a cut are computed, the others bounded through the tail of A

# the columns up to the cut are computed; beyond it, (I - A DF) e_j = -A B M(w) e_j with w the entries of DN
function _Z_bound(kind::AbstractSequenceProblem, pb, x, A, R)
    X = _X(pb, x)
    x_R = ball(pb, x, R)
    W, B, s = multipliers(pb, x_R), prefactor(pb, x_R), _space.(x_R)
    Πs = cut(B, s, W)
    D = Πs - A * (DF(pb, x_R) * Πs)
    return max.(opnorm.(D, permutedims(X), X), Z_tail(B, s, W, A, _pointwise(pb).X, X))   # exact per block: weighted ℓ¹ domains
end

# the sequence part of a multiplier only meets the tail of A, its error meets all of A;
# s and X are the spaces of the blocks and their Banach spaces, Xn the Banach spaces the multipliers are measured in
function Z_tail(B, s, W, A, X, Xn)
    S = [_sequence_norm(W[m,i], Xn[m]) for m ∈ axes(W, 1), i ∈ axes(W, 2)]
    return _tail_term(B, A, X, s, S) + opnorms(A) * ([_opnorm(B, X[l], s[l]) for l ∈ eachindex(s)] .* _error.(W))
end
_beyond_term(B, A, X, s, S) = [_beyond(B, X[m], s[i]) for m ∈ eachindex(s), i ∈ eachindex(s)] .* tail_bounds(A) .* S
_tail_term(B, A, X, s, S) = _beyond_term(B, A, X, s, S)
# what B writes below the projection meets A_K
_tail_term(B::Integral, A, X, s, S) =
    _beyond_term(B, A, X, s, S) + head_norms(A) * ([_leak(B, X[l], s[i]) for l ∈ eachindex(s), i ∈ eachindex(s)] .* S)
tail_bounds(A::ProjectedInverse) = A.tail_bound

# the columns beyond K_i + bandwidth of the multipliers (+ the downward shift of B) only see the tail of A
cut(B, s, W) =
    Diagonal([interval(Projection(_shift(B, reduce(∪, [_cut(s[i], W[m,i]) for m ∈ axes(W, 1)])))) for i ∈ eachindex(s)])
_cut(s::ScalarSpace, w) = s
_cut(s::SequenceSpace, ::Number) = s
_cut(s::SequenceSpace, w) = _cut(s, space(w))
_cut(s::SequenceSpace, w::VectorSpace) = codomain(*, _bandwidth(w), s)
_bandwidth(::Taylor) = Taylor(0)   # lower triangular products
_bandwidth(s::SequenceSpace) = s
_bandwidth(s::TensorSpace) = TensorSpace(map(_bandwidth, spaces(s)))
_bandwidth(s::SymmetricSpace) = SymmetricSpace(_bandwidth(desymmetrize(s)), symmetry(s))

# the diagonal of the pseudo-inverse at the first mode beyond k, a bound of ‖Linv Π_{>k}‖
_linv_bound(pb, s, k) = abs(RadiiPolynomial.getcoefficient(pb.Linv, (s, k + 1), (s, k + 1)))

#- what B does beyond the order of a space: its downward shift, the norm of what stays above, of what leaks to the mode 0

_shift(B, s) = s
_shift(::Integral, s::Chebyshev) = Chebyshev(order(s) + 1)   # ∫ e_j touches the mode j - 1
_beyond(B::Number, X, s) = abs(B)
_beyond(::Integral, X, s::Taylor) = rate(X) / interval(order(s) + 2)
_beyond(::Integral, X, s::Chebyshev) = (rate(X) / interval(order(s) + 3) + inv(rate(X) * interval(order(s) + 1))) / exact(2)
_leak(::Integral, X, ::Taylor) = zero(rate(X))
_leak(::Integral, X, s::Chebyshev) = inv(pown(rate(X), order(s) + 2) * interval((order(s) + 2)^2 - 1))   # the value at -1
