# ------------------------------------------------------------------
# The symmetric Fourier spaces of the P₁₂ family (paper, Section 2).
#
# RadiiPolynomial's SymmetricSpace is the subspace invariant under a group
# acting on the coefficients by (g ⋅ a)_k = α_g(k) a_{β_g(k)}, with β_g a map
# of the index lattice and α_g(k) = e^{iπ(θ + φk)}. Two generators are needed:
#   flip(θ)   : a_k = e^{iπθ} a_{-k}        θ = 0: cosine series, θ = 1: sine series
#   parity(θ) : a_k = e^{iπ(θ + k)} a_k     θ = 0: even modes only, θ = 1: odd modes only
# One coefficient per orbit is stored (k ≥ 0 for the cosine and sine series,
# in increasing order); norms, products, derivatives and projections account
# for the full orbit.
# ------------------------------------------------------------------

_flip(θ::Rational{Int})   = GroupElement(LatticeAut([-1;;]), Cocycle(θ, Rational{Int}[0//1]))
_parity(θ::Rational{Int}) = GroupElement(LatticeAut([ 1;;]), Cocycle(θ, Rational{Int}[1//1]))

evencos(s::Fourier)     = SymmetricSpace(s, Group(_flip(0//1), _parity(0//1))) # Σ_{k even} a_k cos(kt)
oddcos(s::Fourier)      = SymmetricSpace(s, Group(_flip(0//1), _parity(1//1))) # Σ_{k odd}  a_k cos(kt)
evensin(s::Fourier)     = SymmetricSpace(s, Group(_flip(1//1), _parity(0//1))) # Σ_{k even} a_k sin(kt)
oddsin(s::Fourier)      = SymmetricSpace(s, Group(_flip(1//1), _parity(1//1))) # Σ_{k odd}  a_k sin(kt)
evenfourier(s::Fourier) = SymmetricSpace(s, Group(_parity(0//1)))              # Σ_{k even} a_k e^{ikt}

"""
    marchal_space(Ku, Kv, Kw, ω = 1)

The space of the unknown x = (a, β, α, u, v, w): three scalars, u ∈ EvenCos ×
EvenSin × OddCos, v ∈ EvenSin × EvenCos × OddSin and w ∈ EvenFourier, with the
Fourier truncation orders `Ku`, `Kv`, `Kw` and frequency `ω`.
"""
function marchal_space(Ku::Int, Kv::Int, Kw::Int, ω = 1)
    U = evencos(Fourier(Ku, ω)) × evensin(Fourier(Ku, ω)) × oddcos(Fourier(Ku, ω))
    V = evensin(Fourier(Kv, ω)) × evencos(Fourier(Kv, ω)) × oddsin(Fourier(Kv, ω))
    W = evenfourier(Fourier(Kw, ω))
    return CartesianProduct(ScalarSpace(), ScalarSpace(), ScalarSpace(), U, V, W)
end
marchal_space(K::Int, ω = 1) = marchal_space(K, K, K, ω)

# ------------------------------------------------------------------
# Real-valuedness. A real function has a_{-k} = conj(a_k), so the stored
# coefficients are real on the cosine spaces, purely imaginary on the sine
# spaces, and conjugate in pairs (k, -k) on EvenFourier. The Chebyshev
# coefficients in Ω of a real function are real too, so the same rule applies
# to Chebyshev ⊗ (symmetric Fourier). For interval data, the average of two
# enclosures of the same value is an enclosure.
# ------------------------------------------------------------------

# the angle θ (in units of π) of the relation a_k = e^{iπθ} a_{-k} at the index
# `k`, or `nothing` without reflection. The angle is evaluated at `k`: a
# composite such as flip ∘ parity carries its sign in the phase φ.
_dot(φ, k::Int) = φ[1] * k
_dot(φ, k::Tuple) = sum(φ .* k)
function _flip_angle(G::Group, k)
    for g ∈ elements(G)
        if g.lattice_aut.matrix[end,end] == -1
            return mod(g.cocycle.angle + _dot(g.cocycle.phase, k), 2)
        end
    end
    return nothing
end

_reflect(k::Int) = -k
_reflect(k::Tuple) = Base.setindex(k, -last(k), length(k))

function impose_reality!(x::Sequence{<:CartesianSpace})
    for xᵢ ∈ eachcomponent(x)
        impose_reality!(xᵢ)
    end
    return x
end

function impose_reality!(x::Sequence{<:Union{ScalarSpace,Chebyshev}})
    coefficients(x) .= real.(coefficients(x))
    return x
end

function impose_reality!(x::Sequence{<:SymmetricSpace})
    s = space(x)
    isempty(indices(s)) && return x
    θ = _flip_angle(symmetry(s), first(indices(s)))
    if θ == 0           # cosine series
        coefficients(x) .= real.(coefficients(x))
    elseif θ == 1       # sine series
        coefficients(x) .= complex.(zero(real(eltype(x))), imag.(coefficients(x)))
    else                # no reflection: conjugate pairs
        for k ∈ indices(space(x))
            if last(k) > 0
                v = (x[k] + conj(x[_reflect(k)])) / exact(2)
                x[k] = v
                x[_reflect(k)] = conj(v)
            elseif last(k) == 0
                x[k] = real(x[k])
            end
        end
    end
    return x
end

# ------------------------------------------------------------------
# The Banach space of the paper: ℓ¹ with geometric weight ν in the Fourier
# variable, unweighted in the Chebyshev variable, summed over the components.
# ------------------------------------------------------------------

norm_space(::ScalarSpace, ν) = Ell1(IdentityWeight())
norm_space(::Chebyshev, ν) = Ell1(IdentityWeight())
norm_space(::Fourier, ν) = Ell1(GeometricWeight(ν))
norm_space(::TensorSpace{<:Tuple{Chebyshev,Fourier}}, ν) = Ell1(IdentityWeight(), GeometricWeight(ν))
norm_space(s::SymmetricSpace, ν) = norm_space(desymmetrize(s), ν)
norm_space(s::CartesianProduct, ν) = NormedCartesianSpace(map(sᵢ -> norm_space(sᵢ, ν), spaces(s)), Ell1())
norm_space(s::CartesianPower, ν) = NormedCartesianSpace(norm_space(space(s), ν), Ell1())

# ------------------------------------------------------------------
# Tails: zero the Fourier modes |k| ≤ K, keeping the rest.
# ------------------------------------------------------------------

function tail!(x::Sequence{<:CartesianSpace}, K::Int)
    for xᵢ ∈ eachcomponent(x)
        tail!(xᵢ, K)
    end
    return x
end

function tail!(x::Sequence{<:SequenceSpace}, K::Int)
    for k ∈ indices(space(x))
        if abs(last(k)) ≤ K
            x[k] = zero(eltype(x))
        end
    end
    return x
end

# a component constant in t (a scalar, or its Chebyshev expansion in Ω) has no Fourier tail
function tail!(x::Sequence{<:Union{ScalarSpace,Chebyshev}}, ::Int)
    coefficients(x) .= zero(eltype(x))
    return x
end

tail(x, K::Int) = tail!(copy(x), K)
