# The Shigesada-Kawasaki-Teramoto system on cosine series, Δ Φ(v, λ) + R(v) = 0: the diffusion is nonlinear

using RadiiPolynomial, LinearAlgebra
using Capsule: QuasilinearProblem

struct PseudoInverseLaplacian <: RadiiPolynomial.AbstractDiagonalOperator end
RadiiPolynomial.getcoefficient(::PseudoInverseLaplacian, (codom, i)::Tuple{SymmetricSpace{<:Fourier},Integer}, (dom, j)::Tuple{SymmetricSpace{<:Fourier},Integer}) =
    (i == j) & !(i == j == 0) ? inv(-(frequency(dom) * exact(i))^2) : zero(frequency(dom))

const d₁₁, d₁₂, d₂₁, d₂₂ = exact(0), exact(3), exact(0), exact(0)
const r₁, a₁, b₁ = exact(5), exact(3), exact(1)
const r₂, b₂, a₂ = exact(2), exact(1), exact(3)

Φ(v, λ)  = [(λ + d₁₁*v[1] + d₁₂*v[2]) * v[1],
            (λ + d₂₁*v[1] + d₂₂*v[2]) * v[2]]
DΦ(v, λ) = [exact(2)*d₁₁*v[1]+d₁₂*v[2]+λ                      d₁₂*v[1]
                              d₂₁*v[2]      exact(2)*d₂₂*v[2]+d₂₁*v[1]+λ]
R(v, λ)  = [(r₁ - a₁*v[1] - b₁*v[2]) * v[1],
            (r₂ - b₂*v[1] - a₂*v[2]) * v[2]]
DR(v, λ) = [r₁-exact(2)*a₁*v[1]-b₁*v[2]                       -b₁*v[1]
                               -b₂*v[2]    r₂-exact(2)*a₂*v[2]-b₂*v[1]]

problem(λ) = QuasilinearProblem(Laplacian(), PseudoInverseLaplacian(), Φ, DΦ, R, DR, λ, [Ell1(), Ell1()])

K, K_A = 80, 165
S = evensym(Fourier(K, π))
