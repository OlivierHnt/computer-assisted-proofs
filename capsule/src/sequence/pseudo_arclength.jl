# PseudoArclengthProblem
# x = (u..., λ): the problem F(u, λ) = L Φ(u, λ) + N(u, λ) = 0 with its parameter λ as an unknown, and (u, λ) on the
# hyperplane through `point` orthogonal to `tangent`. Along a branch the path is the arclength: the parameters of
# this problem are the tangent and the point.

struct PseudoArclengthProblem{P,FP,FN,TP,BX} <: AbstractSequenceProblem
    problem :: P     # its parameter is the scalar λ
    DλΦ     :: FP    # derivatives with respect to λ, one block per component (nothing: none)
    DλN     :: FN
    p       :: TP    # (tangent = (u̇..., λ̇), point = (w..., λ_w))
    X       :: BX
end
PseudoArclengthProblem(pb::AbstractSequenceProblem, DλΦ, DλN, p::NamedTuple) = PseudoArclengthProblem(pb, DλΦ, DλN, p, [pb.X ; Ell1()])

_with_parameters(pa::PseudoArclengthProblem, p) = PseudoArclengthProblem(pa.problem, pa.DλΦ, pa.DλN, p, pa.X)
_kind(pa::PseudoArclengthProblem) = _kind(pa.problem)

function F(pa::PseudoArclengthProblem, x)
    u, λ = x[1:end-1], x[end]
    t, w = pa.p.tangent, pa.p.point
    return [F(_with_parameters(pa.problem, λ), u) ; sum(adjoint(t[i]) * (u[i] - w[i]) for i ∈ eachindex(u)) + t[end] * (λ - w[end])]
end

function DF(pa::PseudoArclengthProblem, x)
    u, λ = x[1:end-1], x[end]
    t = pa.p.tangent
    return [DF(_with_parameters(pa.problem, λ), u)    LinearOperator.(_DλF(pa, u, λ))
            permutedims(adjoint.(t[1:end-1]))          LinearOperator(t[end])]
end
_DλF(pa, u, λ) = isnothing(pa.DλΦ) ? pa.DλN(u, λ) : isnothing(pa.DλN) ? Ref(pa.problem.L) .* pa.DλΦ(u, λ) : Ref(pa.problem.L) .* pa.DλΦ(u, λ) .+ pa.DλN(u, λ)

# the statements of the problem, bordered by the derivatives in λ and by the tangent
leading(pa::PseudoArclengthProblem, x) = _bordered(leading(_with_parameters(pa.problem, x[end]), x[1:end-1]), pa.DλΦ, x, fill(interval(0), length(x) - 1))
multipliers(pa::PseudoArclengthProblem, x) = _bordered(multipliers(_with_parameters(pa.problem, x[end]), x[1:end-1]), pa.DλN, x, pa.p.tangent[1:end-1])
_bordered(M, Dλ, x, row) = [M  (isnothing(Dλ) ? fill(interval(0), length(row)) : Dλ(x[1:end-1], x[end])) ; permutedims(row)  interval(0)]
prefactor(pa::PseudoArclengthProblem, x) = prefactor(_with_parameters(pa.problem, x[end]), x[1:end-1])
tail(pa::PseudoArclengthProblem, s::SequenceSpace) = tail(pa.problem, s)
tail_bound(pa::PseudoArclengthProblem, s::SequenceSpace) = tail_bound(pa.problem, s)

ProjectedInverse(pa::PseudoArclengthProblem{<:QuasilinearProblem}, x, s = _space.(x)) = _quasilinear_inverse(pa.problem, pa, x, s)
