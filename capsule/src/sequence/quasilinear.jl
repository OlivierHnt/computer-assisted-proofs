# QuasilinearProblem
# x = (u...): F(u) = L Φ(u) + N(u), DF(u) = L M(DΦ(u)) + M(DN(u)). The leading part is not L: beyond the projection
# the inverse is M(W) L⁻¹ with W an approximate inverse of the matrix DΦ(u), and I - M(W) M(DΦ(u)) does not decay

struct QuasilinearProblem{TL,TLi,FP,FDP,FN,FDN,P,BX} <: AbstractSequenceProblem
    L    :: TL    # a diagonal operator
    Linv :: TLi   # its pseudo-inverse, nonincreasing in the mode beyond the projection
    Φ    :: FP    # nonlinearity under L, and its derivative
    DΦ   :: FDP
    N    :: FN    # nonlinearity, and its derivative
    DN   :: FDN
    p    :: P
    X    :: BX
end

F(pb::QuasilinearProblem, x) = Ref(pb.L) .* pb.Φ(x, pb.p) .+ pb.N(x, pb.p)
DF(pb::QuasilinearProblem, x) = Ref(pb.L) .* Multiplication.(pb.DΦ(x, pb.p)) + Multiplication.(pb.DN(x, pb.p))

leading(pb::QuasilinearProblem, x) = pb.DΦ(x, pb.p)
multipliers(pb::QuasilinearProblem, x) = pb.DN(x, pb.p)
_with_parameters(pb::QuasilinearProblem, p) = QuasilinearProblem(pb.L, pb.Linv, pb.Φ, pb.DΦ, pb.N, pb.DN, p, pb.X)

# A = A_K Π + M(W) L⁻¹ - Π M(W) L⁻¹ Π on the spaces s: A_K the compression of the inverse on the doubled spaces,
# W the inverse of DΦ(u) at the nodes of a grid, interpolated. A scalar block has no tail
ProjectedInverse(pb::QuasilinearProblem, x, s = _space.(x)) = _quasilinear_inverse(pb, pb, x, s)

function _quasilinear_inverse(kind, pb, x, s)
    x, s = interval.(x), _interval.(s)
    q = _sequence_blocks(s)
    Π, Π₂ = Projection(_product(s)), Projection(_product(_doubled.(s)))
    A_K = unpack(Π * _approximate_inverse(Π₂ * DF(pb, x) * Π₂) * Π)
    W = _pointwise_inverse(leading(pb, x)[q, q], s[q])
    finite, T, tail_bound = Matrix{LinearOperator}(A_K), Matrix{Any}(nothing, size(A_K)), fill(interval(0), size(A_K))
    for (a, m) ∈ enumerate(q), (b, l) ∈ enumerate(q)
        G = (Multiplication(W[a,b]) * kind.Linv) * Projection(s[l])
        finite[m,l] = A_K[m,l] + (G - Projection(s[m]) * G)
        T[m,l] = Multiplication(W[a,b]) * (kind.Linv * (interval(I) - interval(Projection(s[l]))))
        tail_bound[m,l] = norm(W[a,b], pb.X[m]) * _linv_bound(kind, s[l], order(s[l]))
    end
    return ProjectedInverse(finite, T, tail_bound, pb.X, W)
end
_sequence_blocks(s) = findall(sᵢ -> sᵢ isa SequenceSpace, s)
_interval(s::ScalarSpace) = s
_interval(s::SequenceSpace) = interval(s)
_doubled(s::ScalarSpace) = s
_doubled(s::SequenceSpace) = codomain(*, s, s)
function _pointwise_inverse(U, s)
    m = grid_size(s[1])
    U_grid = [to_grid(mid(U[i,j]), m) for i ∈ axes(U, 1), j ∈ axes(U, 2)]
    W_grid = [inv([U_grid[i,j][l] for i ∈ axes(U, 1), j ∈ axes(U, 2)]) for l ∈ CartesianIndices(m)]
    return [interval(real(to_coef([W[i,j] for W ∈ W_grid], s[i]))) for i ∈ axes(U, 1), j ∈ axes(U, 2)]
end

# block norms of A L, L acting on the sequence blocks: on the columns of the projection, and M(W) beyond
function opnorms(A::ProjectedInverse, L)
    s = block_spaces(A)
    q = _sequence_blocks(s)
    N = fill(interval(0), size(A.finite))
    for m ∈ axes(N, 1), l ∈ q
        N[m,l] = opnorm((A.finite[m,l] * L) * Projection(s[l]), A.X[l], A.X[m])
    end
    N[q,q] = max.(N[q,q], norm.(A.W, A.X[q]))
    return N
end

# the errors of DΦ on the ball sit under L: they go through the block norms of A L, those of DN through the ones of A
function _Z_bound(kind::QuasilinearProblem, pb, x, A, R)
    X, x = _X(pb, x), point(pb, x)
    V, W = leading(pb, x), multipliers(pb, x)
    # the centre is padded to the space of its multipliers: the product of two balls is truncated to their spaces
    Ws = reduce(∪, [_full_space(w) for w ∈ [V ; W] if !(w isa Number)])
    x_R = ball(pb, _pad.(x, Ref(Ws)), R)
    s = block_spaces(A)
    q = _sequence_blocks(s)
    Πs = cut(interval(1), s, [V ; W])
    known = opnorm.(_operator.(Πs - A * (DF(pb, x) * Πs)), permutedims(X), X)
    # beyond the cut, (I - A DF) e_j = (I - M(W_A) M(DΦ)) e_j - M(W_A) L⁻¹ M(DN) e_j on the sequence blocks
    W_A = tail_multipliers(A)
    tail = tail_bounds(A) * [_sequence_norm(W[m,i], X[m]) for m ∈ axes(W, 1), i ∈ axes(W, 2)]
    tail[q,q] += norm.(_identity_minus([sum(W_A[a,c] * V[q[c],q[b]] for c ∈ eachindex(q)) for a ∈ eachindex(q), b ∈ eachindex(q)]), X[q])
    errors = opnorms(A, kind.L) * _error.(leading(pb, x_R)) + opnorms(A) * _error.(multipliers(pb, x_R))
    return max.(known, tail) + errors
end
block_spaces(A::ProjectedInverse) = [domain(A.finite[1,i]) for i ∈ axes(A.finite, 2)]
tail_multipliers(A::ProjectedInverse) = A.W
_full_space(a::Sequence) = space(a)
_pad(a::Sequence, s) = project(a, s ∪ space(a))
_pad(a::Number, s) = a
