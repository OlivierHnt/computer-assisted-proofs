# The map F(x, Ω) of Section 3 and its derivative, acting directly on the
# symmetric spaces of `marchal_space`. The unknown is x = (a, β, α, u, v, w).
#
# The shifts by 2ζ and 4ζ, ζ = 2π/3, are 2/3 and 4/3 of the period: as rational
# numbers of periods they are exact roots of unity, so no π enters F or DF.
#
# Every product, shift and derivative is carried out on the symmetric
# representation; `project!` onto a component of F then restricts the result
# to the symmetry of that component (it averages over orbits, which is the
# identity on sequences that already carry the symmetry).

const 𝒮  = Shift(2//3) # u ↦ u(⋅ + 2ζ)
const 𝒮² = Shift(4//3) # u ↦ u(⋅ + 4ζ)
const ℰ₀ = Evaluation(0)

# reflection t ↦ -t on EvenFourier, i.e. a_k ↦ a_{-k}
_reflect(a::Sequence{<:SymmetricSpace{<:Fourier}}) = Sequence(space(a), reverse(coefficients(a)))

# in-place scaling of an operator block (LinearAlgebra.rmul! compares the scalar with one, which intervals forbid)
function _scale!(A::LinearOperator, s::Number)
    coefficients(A) .*= s
    return A
end

function _reflection(s::SymmetricSpace{<:Fourier}, ::Type{T}) where {T}
    R = zeros(T, s, s)
    for k ∈ indices(s)
        R[-k,k] = one(T)
    end
    return R
end

function F_marchal!(F, x, Ω)
    a = x[1]
    β = x[2]
    α = x[3]
    u, v, w = component(x, 4), component(x, 5), component(x, 6)
    ux, uy, uz = eachcomponent(u)
    vx, vy, vz = eachcomponent(v)

    # amplitude equation

    F[1] = ℰ₀ * uz - exact(1)

    # zero average

    F[2] = ux[0]

    # initial condition for the polynomial embedding

    M₁ux = ux - 𝒮 * ux
    M₁uy = uy - 𝒮 * uy
    M₁uz = uz - 𝒮 * uz

    ℰ₀M₁ux = ℰ₀ * M₁ux
    ℰ₀M₁uy = ℰ₀ * M₁uy
    ℰ₀M₁uz = ℰ₀ * M₁uz

    ℰ₀w = ℰ₀ * w
    ℰ₀w² = ℰ₀w * ℰ₀w

    F[3] = ℰ₀w² * (ℰ₀M₁ux * ℰ₀M₁ux + ℰ₀M₁uy * ℰ₀M₁uy + a * ℰ₀M₁uz * ℰ₀M₁uz) - exact(1)

    # equations of motion

    M₂ux = ux - 𝒮² * ux
    M₂uy = uy - 𝒮² * uy
    M₂uz = uz - 𝒮² * uz

    w³ = w * w * w
    Rw³ = _reflect(w³)

    Fu = component(F, 4)
    project!(component(Fu, 1), β - Ω^2 * ux + exact(2) * Ω * vy + differentiate(vx) + w³ * M₁ux + Rw³ * M₂ux)
    project!(component(Fu, 2),   - Ω^2 * uy - exact(2) * Ω * vx + differentiate(vy) + w³ * M₁uy + Rw³ * M₂uy)
    project!(component(Fu, 3),                                    differentiate(vz) + w³ * M₁uz + Rw³ * M₂uz)

    # v = u'

    project!(component(F, 5), differentiate(u) - v)

    # polynomial embedding

    M₁vx = vx - 𝒮 * vx
    M₁vy = vy - 𝒮 * vy
    M₁vz = vz - 𝒮 * vz

    project!(component(F, 6), α + differentiate(w) + w³ * (M₁ux * M₁vx + M₁uy * M₁vy + a * M₁uz * M₁vz))

    return F
end

function DF_marchal!(DF, x, Ω)
    a = x[1]
    u, v, w = component(x, 4), component(x, 5), component(x, 6)
    ux, uy, uz = eachcomponent(u)
    vx, vy, vz = eachcomponent(v)

    CoefType = eltype(DF)
    one_ = exact(true)
    Id = UniformScaling(one(real(CoefType))) # identity in the arithmetic of DF (LinearAlgebra.I would mix in Bool)

    coefficients(DF) .= zero(CoefType)

    ℰ₀𝒮 = ℰ₀ * 𝒮 # u ↦ u(2ζ)

    # derivative of the amplitude equation

    project!(component(component(DF, 1, 4), 3), ℰ₀)

    # derivative of the zero average

    component(component(DF, 2, 4), 1)[1,0] = one(CoefType)

    # derivative of the initial condition for the polynomial embedding

    M₁ux = ux - 𝒮 * ux
    M₁uy = uy - 𝒮 * uy
    M₁uz = uz - 𝒮 * uz

    ℰ₀M₁ux = ℰ₀ * M₁ux
    ℰ₀M₁uy = ℰ₀ * M₁uy
    ℰ₀M₁uz = ℰ₀ * M₁uz

    ℰ₀w = ℰ₀ * w
    ℰ₀w² = ℰ₀w * ℰ₀w

    D₃u = component(DF, 3, 4)
    _scale!(sub!(component(D₃u, 1), ℰ₀, ℰ₀𝒮),     exact(2) * ℰ₀w² * ℰ₀M₁ux)
    _scale!(sub!(component(D₃u, 2), ℰ₀, ℰ₀𝒮),     exact(2) * ℰ₀w² * ℰ₀M₁uy)
    _scale!(sub!(component(D₃u, 3), ℰ₀, ℰ₀𝒮), a * exact(2) * ℰ₀w² * ℰ₀M₁uz)

    _scale!(project!(component(DF, 3, 6), ℰ₀), exact(2) * ℰ₀w * (ℰ₀M₁ux * ℰ₀M₁ux + ℰ₀M₁uy * ℰ₀M₁uy + a * ℰ₀M₁uz * ℰ₀M₁uz))

    component(DF, 3, 1)[1,1] = ℰ₀w² * ℰ₀M₁uz * ℰ₀M₁uz

    # derivative of the equations of motion

    M₂ux = ux - 𝒮² * ux
    M₂uy = uy - 𝒮² * uy
    M₂uz = uz - 𝒮² * uz

    w² = w * w
    w³ = w² * w
    Rw² = _reflect(w²)
    Rw³ = _reflect(w³)

    Duu = component(DF, 4, 4)
    Duv = component(DF, 4, 5)
    Duw = component(DF, 4, 6)

    project!(component(Duu, 1, 1), UniformScaling(-Ω^2))
    project!(component(Duv, 1, 1), Derivative(1))
    project!(component(Duv, 1, 2), UniformScaling(exact(2) * Ω))

    project!(component(Duu, 2, 2), UniformScaling(-Ω^2))
    project!(component(Duv, 2, 1), UniformScaling(exact(-2) * Ω))
    project!(component(Duv, 2, 2), Derivative(1))

    project!(component(Duv, 3, 3), Derivative(1))

    # w³ (I - 𝒮) + Rw³ (I - 𝒮²): the shifted spaces are kept exactly in the
    # intermediate step and the result is restricted afterwards
    for i ∈ 1:3
        mul!(component(Duu, i, i), Multiplication(w³),  Id - 𝒮,  one_, one_)
        mul!(component(Duu, i, i), Multiplication(Rw³), Id - 𝒮², one_, one_)
    end

    R = _reflection(domain(component(Duw, 1)), CoefType)
    for (i, M₁uᵢ, M₂uᵢ) ∈ ((1, M₁ux, M₂ux), (2, M₁uy, M₂uy), (3, M₁uz, M₂uz))
        project!(component(Duw, i), Multiplication(exact(3) * w² * M₁uᵢ))
        mul!(component(Duw, i), Multiplication(exact(3) * Rw² * M₂uᵢ), R, one_, one_)
    end

    component(component(DF, 4, 2), 1)[0,1] = one(CoefType)

    # derivative of v = u'

    project!(component(DF, 5, 4), Derivative(1))
    project!(component(DF, 5, 5), UniformScaling(exact(-1)))

    # derivative of the polynomial embedding

    M₁vx = vx - 𝒮 * vx
    M₁vy = vy - 𝒮 * vy
    M₁vz = vz - 𝒮 * vz

    Dwu = component(DF, 6, 4)
    mul!(component(Dwu, 1), Multiplication(    w³ * M₁vx), Id - 𝒮, one_, exact(false))
    mul!(component(Dwu, 2), Multiplication(    w³ * M₁vy), Id - 𝒮, one_, exact(false))
    mul!(component(Dwu, 3), Multiplication(a * w³ * M₁vz), Id - 𝒮, one_, exact(false))

    Dwv = component(DF, 6, 5)
    mul!(component(Dwv, 1), Multiplication(    w³ * M₁ux), Id - 𝒮, one_, exact(false))
    mul!(component(Dwv, 2), Multiplication(    w³ * M₁uy), Id - 𝒮, one_, exact(false))
    mul!(component(Dwv, 3), Multiplication(a * w³ * M₁uz), Id - 𝒮, one_, exact(false))

    add!(component(DF, 6, 6), Derivative(1), Multiplication(exact(3) * w² * (M₁ux * M₁vx + M₁uy * M₁vy + a * M₁uz * M₁vz)))

    component(DF, 6, 3)[0,1] = one(CoefType)

    project!(component(DF, 6, 1), w³ * M₁uz * M₁vz)

    return DF
end
