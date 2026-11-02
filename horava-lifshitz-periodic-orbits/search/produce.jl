# SEARCH — not part of the proof.
#
# Regenerates the 18 approximate zeros x̄ of F (Bianchi type VIII and IX,
# N₁N₂N₃ = 10^k, k = 1, 0, …, -7) and writes them, with the constants
# C₁, C₂, ν, R of each case, to data/*.bin in the format read by the checker.
#
# The original repository ships only the data, not the code that produced
# them; this script is a reconstruction: natural continuation in log₁₀(N₁N₂N₃)
# from the seeds of search/seeds.jl (the N₁N₂N₃ = 10 orbits to 6 digits), with
# Newton's method on the truncated F in binary64, followed for the 256-bit
# cases by a chord-Newton refinement in 256-bit BigFloat arithmetic.
# Any other method of producing the candidates would do: the checker does not
# care how the data were obtained.
#
# Usage (from the package root; this is the `produce` command of entry.json):
#   julia --project=search search/produce.jl [outdir]      (default outdir: data)
# then `julia --project=checker tools/hash_data.jl` to record the new hashes.

import Pkg
Pkg.instantiate()

using RadiiPolynomial, LinearAlgebra
include(joinpath(@__DIR__, "..", "checker", "src", "HLCertificate.jl"))
using .HLCertificate
const F! = HLCertificate.F!
const DF! = HLCertificate.DF!
include(joinpath(@__DIR__, "seeds.jl"))

BLAS.set_num_threads(min(2, Sys.CPU_THREADS))

outdir = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "..", "data")
mkpath(outdir)

# type, k, file tag, Fourier order n, number format, weight ν, radius R
# (the constants of the original data)
const CASES = Dict(
    "VIII" => [( 1, "10",    32, :mp256, "1.3",  "1e-3"),
               ( 0, "1",     32, :f64,   "1.3",  "1e-3"),
               (-1, "1e-1",  32, :f64,   "1.2",  "1e-3"),
               (-2, "1e-2",  32, :f64,   "1.1",  "1e-3"),
               (-3, "1e-3",  32, :f64,   "1.1",  "1e-4"),
               (-4, "1e-4",  64, :f64,   "1.1",  "1e-4"),
               (-5, "1e-5",  64, :f64,   "1.1",  "1e-4"),
               (-6, "1e-6", 128, :f64,   "1.07", "5e-7"),
               (-7, "1e-7", 256, :mp256, "1.05", "1e-9")],
    "IX"   => [( 1, "10",    32, :f64,   "1.1",  "1e-3"),
               ( 0, "1",     64, :f64,   "1.15", "1e-3"),
               (-1, "1e-1",  64, :f64,   "1.1",  "1e-3"),
               (-2, "1e-2",  64, :f64,   "1.07", "1e-3"),
               (-3, "1e-3", 128, :f64,   "1.03", "1e-4"),
               (-4, "1e-4", 192, :f64,   "1.03", "1e-4"),
               (-5, "1e-5", 192, :f64,   "1.03", "1e-4"),
               (-6, "1e-6", 192, :mp256, "1.03", "1e-9"),
               (-7, "1e-7", 256, :mp256, "1.03", "1e-10")])

const STEPS_PER_DECADE = 10

space_x(n) = Fourier(n, 1.0)^5 × ScalarSpace()^3

# reality of the Fourier series and of the parameters (as imposed by the checker)
function symmetrize!(x)
    n = order(space(space(x)[1]))
    for x₁ⱼ ∈ eachcomponent(component(x, 1))
        c = coefficients(x₁ⱼ)
        for i ∈ eachindex(c) # in place, k = -n, …, n, as the RadiiPolynomial 0.7.2 broadcast
            c[i] = (c[i] + conj(c[end+1-i])) / 2
        end
        x₁ⱼ[0] = real(x₁ⱼ[0])
    end
    coefficients(component(x, 2)) .= real.(coefficients(component(x, 2)))
    return x
end

# Newton's method in binary64; stops when the step is at rounding level
function newton64!(x, c; tol = 1e-14, maxiter = 30)
    s = space(x)
    F  = zeros(ComplexF64, s)
    DF = zeros(ComplexF64, s, s)
    nδ_old = Inf
    for i ∈ 1:maxiter
        F!(F, x, 1.0, c)
        nF = norm(coefficients(F), Inf)
        isfinite(nF) || return false
        scale = max(1, norm(coefficients(x), Inf))
        DF!(DF, x, 1.0, c)
        δ = coefficients(DF) \ coefficients(F)
        nδ = norm(δ, Inf)
        # stagnation of the step at rounding level, with a small residual
        i ≥ 3 && nδ ≥ nδ_old / 2 && nF ≤ 1e-10 * scale && return true
        coefficients(x) .-= δ
        symmetrize!(x)
        nδ ≤ tol * scale && return true
        nδ_old = nδ
    end
    println("      Newton: no convergence, last |F|∞ = ", norm(coefficients(F), Inf))
    return false
end

# chord-Newton in 256-bit arithmetic: x ← x - A F(x), A = DF(x)⁻¹ in binary64
function refine256(x64, c::BigFloat; maxiter = 30)
    s = space(x64)
    DF = zeros(ComplexF64, s, s)
    DF!(DF, x64, 1.0, Float64(c))
    A = inv(coefficients(DF))
    x = Sequence(s, Complex{BigFloat}.(coefficients(x64)))
    F = zeros(Complex{BigFloat}, s)
    nF_old = big(Inf)
    for _ ∈ 1:maxiter
        F!(F, x, one(BigFloat), c)
        nF = norm(coefficients(F), Inf)
        nF ≥ nF_old / 2 && break             # no more progress: rounding level reached
        nF_old = nF
        coefficients(x) .-= A * coefficients(F)
        symmetrize!(x)
    end
    println("      256-bit refinement: |F|∞ = ", Float64(nF_old), " before the last step")
    return x
end

# the tightest binary64 interval containing 10^k (the C₂ of the original data)
function enclosure10(k)
    q = k ≥ 0 ? big(10)^k // 1 : 1 // big(10)^(-k)
    f = Float64(q)
    lo = Rational{BigInt}(f) ≤ q ? f : prevfloat(f)
    hi = Rational{BigInt}(f) ≥ q ? f : nextfloat(f)
    return interval(lo, hi)
end

thin(x) = Sequence(space(x), complex.(interval.(real.(coefficients(x))), interval.(imag.(coefficients(x)))))

function write_output(type, k, tag, fmt, νstr, Rstr, x)
    file = joinpath(outdir, "$(type)_N1N2N3_$(tag).bin")
    C₂ = enclosure10(k)
    if fmt === :f64
        write_case(file, thin(x), interval(1.0), C₂, interval(parse(Float64, νstr)), parse(Float64, Rstr))
    else
        x256 = refine256(x, k ≥ 0 ? BigFloat(10)^k : inv(BigFloat(10)^(-k)))
        ν = interval(BigFloat(νstr, RoundDown), BigFloat(νstr, RoundUp))
        write_case(file, thin(x256), interval(BigFloat(1.0)), interval(BigFloat(inf(C₂)), BigFloat(sup(C₂))),
            ν, BigFloat(parse(Float64, Rstr)))
    end
    println("    wrote ", file)
end

for type ∈ ("VIII", "IX")
    γ₀, comps = SEEDS[type]
    m = (length(comps[1]) - 1) ÷ 2
    cases = CASES[type]

    x = zeros(ComplexF64, space_x(cases[1][3]))
    for j ∈ 1:5, (i, kk) ∈ enumerate(-m:m)
        component(component(x, 1), j)[kk] = comps[j][i]
    end
    component(x, 2)[1] = γ₀

    log_c = Float64(cases[1][1])
    for (k, tag, n, fmt, νstr, Rstr) ∈ cases
        x = project(x, space_x(n))   # n is nondecreasing along the continuation
        targets = k == cases[1][1] ? [log_c] : [log_c - (log_c - k) * s / STEPS_PER_DECADE for s ∈ 1:STEPS_PER_DECADE]
        for t ∈ targets
            newton64!(x, exp10(t)) || error("Newton failed for type $type at log₁₀(N₁N₂N₃) = $t")
        end
        log_c = Float64(k)
        println("$type, N₁N₂N₃ = 10^$k (n = $n): γ = ", real(component(x, 2)[1]))
        write_output(type, k, tag, fmt, νstr, Rstr, x)
    end
end
