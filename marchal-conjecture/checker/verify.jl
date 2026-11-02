# CHECKER — the proof.
#
# Reads the frozen inputs (x̄ as intervals and A as floats on the Chebyshev
# nodes), encloses everything in interval arithmetic, and evaluates the
# conditions of statement.md: the input enclosures, the bounds of Lemmas 3.2,
# 3.3, 3.4 of the paper, the radii-polynomial inequality, and the two sign
# conditions of Lemma 2.2 (iv). Nothing is searched for or solved here: only
# evaluated.
#
# Usage (from the package root):
#   julia --project=checker checker/verify.jl
# Prints each condition, in the order of statement.md, prefixed pass or fail.
# Exit code 0 iff every condition holds.

import Pkg
Pkg.instantiate()

using RadiiPolynomial, LinearAlgebra, SHA, JSON
include(joinpath(@__DIR__, "src", "MarchalCertificate.jl"))
using .MarchalCertificate

BLAS.set_num_threads(Sys.CPU_THREADS)

verdicts = Bool[]
function condition(label, ok, enclosure)
    push!(verdicts, ok)
    println(ok ? "pass " : "fail ", label, ": ", enclosure)
    return ok
end

# ---- load and authenticate the frozen inputs --------------------------------

root    = joinpath(@__DIR__, "..")
datadir = joinpath(root, "data")
hashes  = JSON.parsefile(joinpath(datadir, "hashes.json"))

for f ∈ sort(collect(keys(hashes)))
    p = joinpath(datadir, f)
    isfile(p) || error("$f is missing: run `julia --project=checker tools/fetch_data.jl` first")
    bytes2hex(open(sha256, p)) == hashes[f]["sha256"] || error("sha256 mismatch for $f: not the recorded input")
    println("pass sha256 ", f)
end

# constants of the proof (paper, Section 3)
K = 70        # Fourier truncation order
N = 20        # Chebyshev truncation order
npts = N + 1  # Chebyshev nodes
ν = I"1.1"    # geometric weight
r = I"1e-6"   # radius

# spaces (symmetric representation, see src/spaces.jl)

X  = marchal_space(K, interval(1))          # the unknown x = (a, β, α, u, v, w)
XF = marchal_space(4K, K, 5K, interval(1))  # F(x̄), untruncated
XD = marchal_space(2K, interval(1))         # domain of DF(x̄) for the Z₁ bound
d  = dimension(X)

filesize(joinpath(datadir, "x_nodes.bin")) == npts * d * 4 * 8 || error("x_nodes.bin has the wrong size")

x_raw   = read_witness(datadir, npts, d)
x_nodes = [Sequence(X, v) for v ∈ x_raw] # interval-valued node values of x̄
println("loaded x̄ (intervals) on $npts nodes, d = $d")

# Auxiliary, not an input: the floating-point inverse of DF at the midpoint of
# each node value, as search/ computes it; the conditions hold for whatever
# matrix is used, as long as the same one is used throughout.
Ω_nodes = [0.5 + 0.5cospi(j/N) for j ∈ 0:N]
A_nodes = [LinearOperator(X, X, approximate_inverse(xⱼ, Ωⱼ, K)) for (xⱼ, Ωⱼ) ∈ zip(x_nodes, Ω_nodes)]
println("computed the approximate inverses A on the $npts nodes")
println()

# ---- condition 0: the inputs enclose the stated data -------------------------
# The proof relies on the interpolant passing through the Lagrange triangle at
# Ω = 1 and through a point with vanishing x-components at Ω = 0, and on x̄
# being real-valued with β = α = 0. The inputs store intervals; here we verify
# they enclose those exact values and carry those properties. Downstream the
# STORED intervals are used, not these reconstructions.

widest(cs) = maximum(c -> max(radius(real(c)), radius(imag(c))), cs)

encloses(stored, exact) = issubset_interval(real(exact), real(stored)) && issubset_interval(imag(exact), imag(stored))

# at Ω = 1, the periodic orbit is the Lagrange triangle

triangle = zeros(Complex{Interval{Float64}}, X)
let u = component(triangle, 4), v = component(triangle, 5), w = component(triangle, 6)
    component(u, 1)[2] =                 inv(sqrt(cbrt(interval(3)))) / interval(2)
    component(u, 2)[2] = interval(-im) * inv(sqrt(cbrt(interval(3)))) / interval(2)
    component(u, 3)[1] = interval(1) / interval(2)

    component(v, 1)[2] = interval( im) * inv(sqrt(cbrt(interval(3))))
    component(v, 2)[2] =                 inv(sqrt(cbrt(interval(3))))
    component(v, 3)[1] = interval( im) / interval(2)

    w[0] = inv(cbrt(interval(3)))
end

check_tri = all(encloses.(coefficients(x_nodes[1]), coefficients(triangle)))
condition("0a. node Ω = 1 encloses the Lagrange triangle", check_tri,
    "w₀ ∈ $(real(component(x_nodes[1], 6)[0])) ∋ 3^(-1/3); widest stored coefficient radius $(widest(coefficients(x_nodes[1])))")

# at Ω = 0, the orbit lies in the yz-plane

Ux₀ = coefficients(component(component(x_nodes[npts], 4), 1))
Vx₀ = coefficients(component(component(x_nodes[npts], 5), 1))
check_planar = all(encloses(c, zero(c)) for c ∈ Ux₀) &&
               all(encloses(c, zero(c)) for c ∈ Vx₀)
condition("0b. node Ω = 0 encloses U_x = V_x = 0", check_planar,
    "max |stored coefficient| of U_x, V_x: $(maximum(c -> max(mag(real(c)), mag(imag(c))), [Ux₀; Vx₀]))")

# the node values are real-valued functions (cosine coefficients real, sine
# coefficients imaginary, conjugate pairs on w) with β = α = 0; the Chebyshev
# interpolant of such data has the same properties, which the DFT below imposes

check_real = all(x_nodes) do xⱼ
    all(isequal_interval.(coefficients(impose_reality!(copy(xⱼ))), coefficients(xⱼ))) &&
        encloses(xⱼ[2], zero(xⱼ[2])) && isthinzero(real(xⱼ[2])) && isthinzero(real(xⱼ[3]))
end
condition("0c. every node is real-valued with β = α = 0", check_real,
    "β, α ∈ $(reduce(hull, [f(x[i]) for x ∈ x_nodes for i ∈ (2, 3) for f ∈ (real, imag)])) (real and imaginary parts) at all $npts nodes")

# ---- Chebyshev interpolants of x̄ and A ---------------------------------------
# x̄ is the Chebyshev interpolant of the interval node values; its coefficients
# are enclosed by the DFT in 256-bit interval arithmetic, and so are its values
# on the finer grids used below. A is the (floating-point) interpolant of the
# node matrices; it is any operator, so no enclosure is needed for its
# construction, only for its use.

_big(z) = complex(interval(BigFloat, real(z)), interval(BigFloat, imag(z)))
_f64(z) = complex(interval(Float64, real(z)), interval(Float64, imag(z)))

m = 7N + 1 # Chebyshev–Lobatto nodes resolving the degree 7N polynomials AF(x̄)
Ω_grid = [inv(interval(2)) + cospi(interval(j) / interval(m - 1)) / interval(2) for j ∈ 0:m-1]

x̄_big, x_grid_big = setprecision(BigFloat, 256) do
    x̄_big = to_coef(map(xⱼ -> Sequence(X, _big.(coefficients(xⱼ))), x_nodes), Chebyshev(N))
    impose_reality!(x̄_big)
    coefficients(component(x̄_big, 2)) .= interval(BigFloat, 0) # β = 0
    coefficients(component(x̄_big, 3)) .= interval(BigFloat, 0) # α = 0
    x̄_big, to_grid(x̄_big, m)
end

x̄     = Sequence(space(x̄_big), _f64.(coefficients(x̄_big)))
x_grid = map(xⱼ -> Sequence(X, _f64.(coefficients(xⱼ))), x_grid_big)

A = let A_float = to_coef(A_nodes, Chebyshev(N))
    LinearOperator(domain(A_float), codomain(A_float), interval.(coefficients(A_float)))
end
A_grid = to_grid(A, m)

opnorm_A = opnorm(A, norm_space(domain(A), ν), norm_space(codomain(A), ν))

# ------------------------------------------------------------------------
# Lemma 3.2 — Y bound

F_grid  = map((xⱼ, Ωⱼ) -> F_marchal!(zeros(Complex{Interval{Float64}}, XF), xⱼ, Ωⱼ), x_grid, Ω_grid)
AF_grid = map((Aⱼ, Fⱼ) -> Aⱼ * project(Fⱼ, X), A_grid, F_grid)

AF     = to_coef(AF_grid, Chebyshev(7N))
tail_F = tail(to_coef(F_grid, Chebyshev(6N)), K)

Y = norm(AF, norm_space(space(AF), ν)) + norm(tail_F, norm_space(space(tail_F), ν)) / interval(K + 1)
condition("1. Y (Lemma 3.2) is a finite bound", isbounded(Y), "Y ∈ $Y")

# ------------------------------------------------------------------------
# Lemma 3.3 — Z₁ bound

DF_grid    = map((xⱼ, Ωⱼ) -> DF_marchal!(zeros(Complex{Interval{Float64}}, XD, X), xⱼ, Ωⱼ), x_grid, Ω_grid)
ADF_I_grid = map((Aⱼ, DFⱼ) -> project(Aⱼ * DFⱼ - UniformScaling(interval(1)), XD, X), A_grid, DF_grid)

ADF_I = to_coef(ADF_I_grid, Chebyshev(6N))

opnorm_ADF_I = opnorm(ADF_I, norm_space(domain(ADF_I), ν), norm_space(codomain(ADF_I), ν))

#--- the explicit bounds M₁, M₂, computed on the full Fourier representation
# (Chebyshev(N) ⊗ Fourier) of the components of x̄

TN(k) = Chebyshev(N) ⊗ Fourier(k, interval(1))
𝒳 = Ell1(IdentityWeight(), GeometricWeight(ν))

Ω = Sequence(Chebyshev(1) ⊗ Fourier(0, interval(1)), [inv(interval(2)), inv(interval(4))]) # Ω = 1/2 + s/2

a = Sequence(TN(0), coefficients(component(x̄, 1)))
ux, uy, uz = (project(c, TN(K)) for c ∈ eachcomponent(component(x̄, 4)))
vx, vy, vz = (project(c, TN(K)) for c ∈ eachcomponent(component(x̄, 5)))
w = project(component(x̄, 6), TN(K))

Δux1 = ux - shift(ux, (0, 2//3))
Δuy1 = uy - shift(uy, (0, 2//3))
Δuz1 = uz - shift(uz, (0, 2//3))

Δux2 = ux - shift(ux, (0, 4//3))
Δuy2 = uy - shift(uy, (0, 4//3))
Δuz2 = uz - shift(uz, (0, 4//3))

Δvx1 = vx - shift(vx, (0, 2//3))
Δvy1 = vy - shift(vy, (0, 2//3))
Δvz1 = vz - shift(vz, (0, 2//3))

w² = w * w
w³ = w² * w

w³Δvx1  =     w³ * Δvx1
w³Δvy1  =     w³ * Δvy1
aw³Δvz1 = a * w³ * Δvz1

w³Δux1  =     w³ * Δux1
w³Δuy1  =     w³ * Δuy1
aw³Δuz1 = a * w³ * Δuz1

w²3Δux1 = interval(3) * w² * Δux1
w²3Δuy1 = interval(3) * w² * Δuy1
w²3Δuz1 = interval(3) * w² * Δuz1

w²3Δux2 = interval(3) * w² * Δux2
w²3Δuy2 = interval(3) * w² * Δuy2
w²3Δuz2 = interval(3) * w² * Δuz2

w²3LaΔuΔv = interval(3) * w² * (Δux1 * Δvx1 + Δuy1 * Δvy1 + a * Δuz1 * Δvz1)

M_1 = max(
    norm(Ω^2, 𝒳) + sqrt(interval(3)) * (interval(2) * norm(w³, 𝒳) + norm( w³Δvx1, 𝒳)),
    norm(Ω^2, 𝒳) + sqrt(interval(3)) * (interval(2) * norm(w³, 𝒳) + norm( w³Δvy1, 𝒳)),
                   sqrt(interval(3)) * (interval(2) * norm(w³, 𝒳) + norm(aw³Δvz1, 𝒳)),
    interval(1) + norm(interval(2) * Ω, 𝒳) + sqrt(interval(3)) * norm( w³Δux1, 𝒳),
    interval(1) + norm(interval(2) * Ω, 𝒳) + sqrt(interval(3)) * norm( w³Δuy1, 𝒳),
    interval(1) +                            sqrt(interval(3)) * norm(aw³Δuz1, 𝒳),
    norm(w²3Δux1, 𝒳) + norm(w²3Δuy1, 𝒳) + norm(w²3Δuz1, 𝒳) +
        norm(w²3Δux2, 𝒳) + norm(w²3Δuy2, 𝒳) + norm(w²3Δuz2, 𝒳) +
        norm(w²3LaΔuΔv, 𝒳)
    )

w_0    = w(nothing, interval(0))
Δux1_0 = Δux1(nothing, interval(0))
Δuy1_0 = Δuy1(nothing, interval(0))
Δuz1_0 = Δuz1(nothing, interval(0))

M_2 = max(
    sqrt(interval(3)) * (               interval(2) * norm(w_0^2 * Δux1_0, 𝒳)  / ν^interval(2K+1) + interval(2) * norm(tail(w³, K), 𝒳) + norm(tail( w³Δvx1, K), 𝒳)),
    sqrt(interval(3)) * (               interval(2) * norm(w_0^2 * Δuy1_0, 𝒳)  / ν^interval(2K+1) + interval(2) * norm(tail(w³, K), 𝒳) + norm(tail( w³Δvy1, K), 𝒳)),
    sqrt(interval(3)) * ((interval(1) + interval(2) * norm(w_0^2 * Δuz1_0, 𝒳)) / ν^interval(2K+1) + interval(2) * norm(tail(w³, K), 𝒳) + norm(tail(aw³Δvz1, K), 𝒳)),
    sqrt(interval(3)) * norm(tail( w³Δux1, K), 𝒳),
    sqrt(interval(3)) * norm(tail( w³Δuy1, K), 𝒳),
    sqrt(interval(3)) * norm(tail(aw³Δuz1, K), 𝒳),
    interval(2) * norm(w_0 * (Δux1_0 * Δux1_0 + Δuy1_0 * Δuy1_0 + a * Δuz1_0 * Δuz1_0), 𝒳) / ν^interval(2K+1) +
        norm(tail(w²3Δux1, K), 𝒳) + norm(tail(w²3Δuy1, K), 𝒳) + norm(tail(w²3Δuz1, K), 𝒳) +
        norm(tail(w²3Δux2, K), 𝒳) + norm(tail(w²3Δuy2, K), 𝒳) + norm(tail(w²3Δuz2, K), 𝒳) +
        norm(tail(w²3LaΔuΔv, K), 𝒳)
    )

#---

Z₁ = max(
    opnorm_ADF_I + max(norm(tail(w³ * Δuz1 * Δvz1, K), 𝒳), M_1) / interval(K + 1),
    opnorm_A * M_2 + M_1 / interval(K + 1)
    )
condition("2. Z₁ (Lemma 3.3) is a finite bound", isbounded(Z₁), "Z₁ ∈ $Z₁")

# ------------------------------------------------------------------------
# Lemma 3.4 — Z₂ bound

â  = norm(a,  𝒳) + r
ûx = norm(ux, 𝒳) + r
ûy = norm(uy, 𝒳) + r
ûz = norm(uz, 𝒳) + r
v̂x = norm(vx, 𝒳) + r
v̂y = norm(vy, 𝒳) + r
v̂z = norm(vz, 𝒳) + r
ŵ  = norm(w,  𝒳) + r

opnorm_D²F =
    max(
        interval(6) * max(ŵ^2 * ûz, ŵ * ûz^2) + interval(3) * max(ŵ^3 * v̂z, ŵ^3 * ûz, interval(3) * ŵ^2 * ûz * v̂z),
        interval(6) * max(              ŵ^2, interval(2) * ŵ * ûx)     + interval(6) * sqrt(interval(3)) * ŵ^2 + interval(3) * max(              ŵ^3, interval(3) * ŵ^2 * v̂x),
        interval(6) * max(              ŵ^2, interval(2) * ŵ * ûy)     + interval(6) * sqrt(interval(3)) * ŵ^2 + interval(3) * max(              ŵ^3, interval(3) * ŵ^2 * v̂y),
        interval(6) * max(ŵ^2 * ûz, â * ŵ^2, interval(2) * â * ŵ * ûz) + interval(6) * sqrt(interval(3)) * ŵ^2 + interval(3) * max(ŵ^3 * v̂z, â * ŵ^3, interval(3) * â * ŵ^2 * v̂z),
                                                                                                    interval(3) * max(              ŵ^3, interval(3) * ŵ^2 * ûx),
                                                                                                    interval(3) * max(              ŵ^3, interval(3) * ŵ^2 * ûy),
                                                                                                    interval(3) * max(ŵ^3 * ûz, â * ŵ^3, interval(3) * â * ŵ^2 * ûz),
        interval(6) * max(ŵ * ûz^2, interval(2) * ŵ * ûx, interval(2) * ŵ * ûy, interval(2) * ŵ * â * ûz, ûx^2 + ûy^2 + â * ûz^2) +
            interval(6) * sqrt(interval(3)) * (max(ŵ^2, interval(2) * ŵ * ûx) + max(ŵ^2, interval(2) * ŵ * ûy) + max(ŵ^2, interval(2) * ŵ * ûz)) +
            interval(9) * max(ŵ^2 * ûz * v̂z,
                              ŵ^2 * v̂x, ŵ^2 * v̂y, ŵ^2 * â * v̂z,
                              ŵ^2 * ûx, ŵ^2 * ûy, ŵ^2 * â * ûz,
                              interval(2) * ŵ * (ûx * v̂x + ûy * v̂y + ûz * v̂z * â),)
        )

Z₂ = max(opnorm_A, inv(interval(K + 1))) * opnorm_D²F
condition("3. Z₂ (Lemma 3.4) is a finite bound", isbounded(Z₂), "Z₂ ∈ $Z₂")

# ---- radii polynomial -----------------------------------------------------
#
# Uses RadiiPolynomial.interval_of_existence, i.e. the library's own
# convention for p(r) = Y + (Z₁ - 1) r + Z₂ r²/2 (check against Theorem
# [radii polynomial] in the paper). The radius r must lie in the
# returned interval of existence.

ie, success_rp = interval_of_existence(Y, Z₁, Z₂, sup(r))
check_rp = success_rp && inf(ie) ≤ inf(r) && sup(r) ≤ sup(ie)
condition("4. radii polynomial: r = $(inf(r)) lies in the interval of existence", check_rp, "interval of existence $ie")

# ------------------------------------------------------------------------
# figure-eight shape, Lemma 2.2 (iv)
# The orbit at Ω = 0 is the node value x_nodes[npts]; the validated functions
# are represented on Fourier(2K) so that their products are computed without
# truncation.

x₀ = x_nodes[npts]
F2K = Fourier(2K, interval(1))
𝒳ν = Ell1(GeometricWeight(ν))

Uy = project(component(component(x₀, 4), 2), F2K)
Uz = sqrt(real(x₀[1])) * project(component(component(x₀, 4), 3), F2K)

Vy = project(component(component(x₀, 5), 2), F2K)
Vz = sqrt(real(x₀[1])) * project(component(component(x₀, 5), 3), F2K)

validated(f) = InfiniteSequence(f, 𝒳ν; total_error = r)

μ = validated(Uy) * validated(Vz) - validated(Uz) * validated(Vy)

# prove that μ is negative on [0, 3/2]

μ_vals = map(t -> real(μ(t)), mince(interval(0, 1.5), 2_000))
check_neg = all(v -> sup(v) < 0, μ_vals)
condition("5a. μ < 0 on [0, 3/2] (2000 pieces)", check_neg,
    "μ([0, 3/2]) ⊆ $(interval(minimum(inf, μ_vals), maximum(sup, μ_vals)))")

# prove that μ''' is positive on [3/2, π/2]
# For ‖f - f̄‖_ν ≤ ε the k-th coefficient of f - f̄ is at most ε ν^{-|k|}, so the
# unweighted ℓ¹ error of the third derivative is at most ε sup_{k ≥ 0} k³ ν^{-k} ≤
# ε (3 / (e log ν))³, which bounds the pointwise error of μ''' on the real line.

C₃ = (interval(3) / (exp(interval(1)) * log(ν)))^3
∂³μ = InfiniteSequence(differentiate(sequence(μ), 3), Ell1(); total_error = C₃ * total_error(μ))

∂³μ_vals = map(t -> real(∂³μ(t)), mince(interval(1.5, sup(interval(π) / interval(2))), 100))
check_pos = all(v -> inf(v) > 0, ∂³μ_vals)
condition("5b. μ''' > 0 on [3/2, π/2] (100 pieces)", check_pos,
    "μ'''([3/2, π/2]) ⊆ $(interval(minimum(inf, ∂³μ_vals), maximum(sup, ∂³μ_vals)))")

# ---- verdict --------------------------------------------------------------

all_pass = all(verdicts)
println()
println(all_pass ? "ALL $(length(verdicts)) CONDITIONS PASS" : "$(count(!, verdicts)) OF $(length(verdicts)) CONDITIONS FAILED")
exit(all_pass ? 0 : 1)
