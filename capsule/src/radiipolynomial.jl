# Additions to RadiiPolynomial.jl, to move there: these methods extend types this package does not own

#- an operator with an error

struct InfiniteLinearOperator{T<:AbstractLinearOperator,R<:Real,Xx<:BanachSpace,Yy<:BanachSpace} <: AbstractLinearOperator
    operator    :: T
    error       :: R
    X :: Xx
    Y :: Yy
end

RadiiPolynomial.Multiplication(a::InfiniteSequence) =
    InfiniteLinearOperator(Multiplication(sequence(a)), total_error(a), banachspace(a), banachspace(a))
RadiiPolynomial.LinearOperator(a::InfiniteSequence) =
    InfiniteLinearOperator(LinearOperator(sequence(a)), total_error(a), Ell1(), banachspace(a))

#- known part and error: what every tier splits on

_known(a::Number) = a
_known(a) = sequence(a)
_operator(L) = L
_operator(P::InfiniteLinearOperator) = P.operator
_error(::AbstractLinearOperator) = exact(0)
_error(P::InfiniteLinearOperator) = P.error
_error(a) = total_error(a)
_error(a::Number) = zero(a)
_sequence_norm(a, X) = norm(sequence(a), X)
_sequence_norm(a::Number, X) = abs(a)
_space(a::Union{Sequence,InfiniteSequence}) = space(a)
_space(::Number) = ScalarSpace()
_basis(ℳ::Multiplication) = space(sequence(ℳ))
_basis(::AbstractLinearOperator) = nothing

#- arithmetic

RadiiPolynomial.domain(P::InfiniteLinearOperator) = domain(P.operator)
RadiiPolynomial.codomain(P::InfiniteLinearOperator) = codomain(P.operator)

RadiiPolynomial.project(P::InfiniteLinearOperator, dom::VectorSpace, codom::VectorSpace) =
    InfiniteLinearOperator(project(P.operator, dom, codom), P.error, P.X, P.Y)
RadiiPolynomial.project(P::InfiniteLinearOperator, dom::VectorSpace, codom::VectorSpace, ::Type{T}) where {T} =
    InfiniteLinearOperator(project(P.operator, dom, codom, T), P.error, P.X, P.Y)

Base.zero(P::InfiniteLinearOperator) = InfiniteLinearOperator(zero(P.operator), zero(P.error), P.X, P.Y)

Base.:-(P::InfiniteLinearOperator) = InfiniteLinearOperator(-P.operator, P.error, P.X, P.Y)
Base.:+(P::InfiniteLinearOperator, Q::InfiniteLinearOperator) = InfiniteLinearOperator(P.operator + Q.operator, P.error + Q.error, P.X, P.Y)
Base.:+(L::AbstractLinearOperator, P::InfiniteLinearOperator) = InfiniteLinearOperator(L + P.operator, P.error, P.X, P.Y)
Base.:+(P::InfiniteLinearOperator, L::AbstractLinearOperator) = InfiniteLinearOperator(P.operator + L, P.error, P.X, P.Y)

Base.:*(P::InfiniteLinearOperator, Q::InfiniteLinearOperator) =
    InfiniteLinearOperator(P.operator * Q.operator, opnorm(P.operator, P.X, P.Y) * Q.error + P.error * opnorm(Q.operator, Q.X, Q.Y) + P.error * Q.error, Q.X, P.Y)
Base.:*(L::AbstractLinearOperator, P::InfiniteLinearOperator) = InfiniteLinearOperator(L * P.operator, _opnorm(L, P.Y, _basis(P.operator)) * P.error, P.X, P.Y)
Base.:*(P::InfiniteLinearOperator, L::AbstractLinearOperator) = InfiniteLinearOperator(P.operator * L, P.error * opnorm(L, P.X), P.X, P.Y)
Base.:*(A::LinearOperator, P::InfiniteLinearOperator) = InfiniteLinearOperator(A * P.operator, opnorm(A, P.Y) * P.error, P.X, P.Y)
Base.:*(P::InfiniteLinearOperator, A::LinearOperator) = InfiniteLinearOperator(P.operator * A, P.error * opnorm(A, P.X), P.X, P.Y)
Base.:*(Π::Projection, P::InfiniteLinearOperator) = InfiniteLinearOperator(Π * P.operator, P.error, P.X, P.Y)
Base.:*(P::InfiniteLinearOperator, Π::Projection) = InfiniteLinearOperator(P.operator * Π, P.error, P.X, P.Y)

Base.:*(P::InfiniteLinearOperator, b::Sequence) =
    InfiniteSequence(P.operator * b, P.Y; total_error = P.error * norm(b, P.X))
Base.:*(P::InfiniteLinearOperator, b::InfiniteSequence) =
    InfiniteSequence(P.operator * sequence(b), P.Y;
        total_error = opnorm(P.operator, P.X, P.Y) * total_error(b) + P.error * (sequence_norm(b) + total_error(b)))

#- norms

LinearAlgebra.opnorm(P::InfiniteLinearOperator, X::BanachSpace, Y::BanachSpace) = opnorm(P.operator, X, Y) + P.error
LinearAlgebra.opnorm(P::InfiniteLinearOperator, X::BanachSpace) = opnorm(P, X, X)
LinearAlgebra.opnorm(P::InfiniteLinearOperator) = opnorm(P, P.X, P.Y)

LinearAlgebra.opnorm(S::AbstractLinearOperator, X::BanachSpace, Y::BanachSpace) =
    X == Y ? opnorm(S, X) : throw(ArgumentError("the operator norm of a lazy operator between different spaces is not defined"))
LinearAlgebra.opnorm(J::UniformScalingOperator, ::BanachSpace) = abs(J.λ)
LinearAlgebra.opnorm(::Projection, ::BanachSpace) = interval(1)
LinearAlgebra.opnorm(S::Negate, X::BanachSpace) = opnorm(S.A, X)
LinearAlgebra.opnorm(S::Add, X::BanachSpace) = opnorm(S.A, X) + opnorm(S.B, X)
LinearAlgebra.opnorm(S::ComposedOperator, X::BanachSpace) = opnorm(S.outer, X) * opnorm(S.inner, X)

# the norm of ∫ depends on the basis, read from the multiplier behind it
_opnorm(L, X, ::Any) = opnorm(L, X)
_opnorm(B::Number, X, ::Any) = abs(B)
_opnorm(L::Integral{Int}, X::Ell1{<:GeometricWeight}, ::Taylor) = pown(rate(X), order(L))
_opnorm(L::Integral{Int}, X::Ell1{<:GeometricWeight}, s::Fourier) = pown(inv(interval(frequency(s))), order(L))
_opnorm(L::Integral{Int}, X::Ell1{<:GeometricWeight}, ::Chebyshev) = pown(interval(1) + rate(X), order(L))
