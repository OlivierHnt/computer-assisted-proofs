# Validated eigenpairs (Section 2.2 of the paper, in the form of the manuscript of
# 28 October 2026, "paper-R2").
#
# For each approximate eigenpair (λ̄, v̄) of M = Df(c⋆), the zero-finding problem is
#     F_eig(v, λ) = (M v − λ v, v^{(l*)} − 1)        on X = ℂ³ × ℂ (or ℝ³ × ℝ),
# with the 1-norm ‖(v, λ)‖_X = Σ_i |v^{(i)}| + |λ| and l* the index of the
# largest-modulus entry of v̄ (so that v̄^{(l*)} = 1 after scaling v̄). The
# procedure is the one of the equilibrium (Section 2.1): A = inv(mid DF(x̄)),
# Y = ‖A F(x̄)‖₁, R = 10 sup(Y), Z = ‖I − A DF(Box(x̄, R))‖₁, and Theorem 1.1
# through `interval_of_existence(Y, Z, R)`: Z < 1 and Y/(1 − Z) ≤ R give a zero
# within r = Y/(1 − Z) of x̄ in the 1-norm, hence each entry within r.
# The 1-norm ball B(x̄, R) is contained in the box of (complex) intervals of
# half-side R around x̄ in every real and imaginary part, so Z over that box is an
# upper bound of sup_{B(x̄, R)} ‖I − A DF(x)‖₁.
#
# A real eigenpair (data with zero imaginary parts, bit for bit) is validated in
# ℝ³ × ℝ. For a complex-conjugate pair one eigenpair is validated and the other
# is defined as its conjugate (paper-R2, Section 2.2): M being real, the
# conjugate of a zero of F_eig is a zero of F_eig with the conjugate data and the
# same l*. The original code validated all three eigenpairs at once on ℂ¹² with
# the sup-norm (IntervalArithmetic's `eigen` for interval matrices).

F_eig(M, v, λ, l) = [M * v - λ * v ; v[l] - one(λ)]

function DF_eig(M, v, λ, l)
    n = length(v)
    T = promote_type(eltype(M), eltype(v), typeof(λ))
    D = zeros(T, n+1, n+1)
    D[1:n, 1:n] .= M - λ * interval(LinearAlgebra.I) # `λ * I` would multiply the interval by a Bool and lose the guarantee flag
    D[1:n, n+1] .= -v
    D[n+1, l] = one(T)
    return D
end

_realdata(λ::Real, v::AbstractVector{<:Real}) = true
_realdata(λ::Complex, v::AbstractVector{<:Complex}) = imag(λ) == 0 && all(imag.(v) .== 0)

_box(x̄::AbstractVector{<:Interval}, R) = x̄ .+ interval(-R, R)
_box(x̄::AbstractVector{<:Complex{<:Interval}}, R) = x̄ .+ complex(interval(-R, R), interval(-R, R))

# Certificate of one eigenpair; λ̄, v̄ are the approximate (floating-point) data.
function validated_eigenpair(M::AbstractMatrix{<:Interval}, λ̄::Number, v̄::AbstractVector)
    if _realdata(λ̄, v̄)
        λ̄, v̄ = real(λ̄), real.(v̄)
    end
    l = argmax(abs.(v̄))
    scale = v̄[l]
    v̄ = v̄ ./ scale
    v̄[l] = one(eltype(v̄))        # v̄^{(l*)} = 1 exactly (the division need not return exactly 1)
    x̄ = interval.([v̄ ; λ̄])        # thin intervals
    v, λ = x̄[1:3], x̄[4]

    A = interval(inv(mid.(DF_eig(M, v, λ, l))))
    Y = LinearAlgebra.norm(A * F_eig(M, v, λ, l), 1)
    R = 10sup(Y)
    box = _box(x̄, R)
    Z = LinearAlgebra.opnorm(interval(LinearAlgebra.I) - A * DF_eig(M, box[1:3], box[4], l), 1)
    ie, success = interval_of_existence(Y, Z, R)
    r = inf(ie)

    # enclosures of the eigenpair, the eigenvector scaled back to the data's scaling
    if success
        true_λ = _box([λ], r)[1]
        true_v = _box(v, r) .* interval(scale)
    else
        nai_x = eltype(x̄) <: Complex ? complex(nai(Interval{Float64}), nai(Interval{Float64})) : nai(Interval{Float64})
        true_λ = nai_x
        true_v = fill(nai_x, length(v))
    end
    return (; λ = true_λ, v = true_v, Y, Z, R, ie, success, r, l, conjugate_of = 0)
end

# Certificates of the n approximate eigenpairs (λ̄_j, v̄_j), j = 1, …, n. An
# eigenpair whose data are, bit for bit, the conjugate of an already validated
# one is defined as its conjugate instead of being validated again.
function validated_eigen(M::AbstractMatrix{<:Interval}, λ̄::Vector, V̄::Matrix)
    n = length(λ̄)
    certs = Vector{Any}(undef, n)
    for j ∈ 1:n
        i = findfirst(i -> λ̄[j] == conj(λ̄[i]) && V̄[:,j] == conj(V̄[:,i]) && λ̄[j] != λ̄[i], 1:j-1)
        if i === nothing
            certs[j] = validated_eigenpair(M, λ̄[j], V̄[:,j])
        else
            c = certs[i]
            certs[j] = (; λ = conj(c.λ), v = conj.(c.v), c.Y, c.Z, c.R, c.ie, c.success, c.r, c.l, conjugate_of = i)
        end
    end
    T = eltype(λ̄) <: Complex ? Complex{Interval{Float64}} : Interval{Float64}
    λ = T[convert(T, c.λ) for c ∈ certs]
    V = T[convert(T, c.v[i]) for i ∈ 1:n, c ∈ certs]
    return (; eig = LinearAlgebra.Eigen(λ, V), certs)
end
