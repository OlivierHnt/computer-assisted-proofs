# PRODUCER — not part of the proof.
#
# Computes the only inputs of the checker: for each of the three examples, the
# approximate unstable eigenvalue λ̄, a binary64 zero of the characteristic
# function f(λ) = det Ψ(λ) found by a plain floating-point Newton iteration
# started from the value printed in the paper (Sections 4.1–4.3). The checker
# proves an enclosure of the eigenvalue around λ̄ and recomputes everything else
# (eigenvector, Taylor coefficients, bounds) in interval arithmetic, so nothing
# here needs to be rigorous.
#
# Output, in data/ (raw little-endian binary64, no header):
#   vdp_lambda.bin   16 bytes: Re λ̄, Im λ̄ (case I, complex eigenvalue)
#   seir_lambda.bin   8 bytes: λ̄ (case II)
#   mg_lambda.bin     8 bytes: λ̄ (case III)
#
# Usage (from the package root):
#   julia --project=search search/produce.jl [output directory, default data/]
#   julia --project=checker tools/hash_data.jl

datadir = length(ARGS) ≥ 1 ? ARGS[1] : joinpath(@__DIR__, "..", "data")   # optional: another output directory
mkpath(datadir)

# Newton iteration on f, stopped when the iterate no longer changes (or after
# 50 steps); f′ is the derivative of f
function newton(f, f′, λ)
    for _ ∈ 1:50
        λ_new = λ - f(λ) / f′(λ)
        λ_new == λ && break
        λ = λ_new
    end
    return λ
end

function save(name, values...)
    path = joinpath(datadir, name)
    open(path, "w") do io
        foreach(x -> write(io, htol(Float64(x))), values)
    end
    println("wrote ", path, " (", filesize(path), " bytes)")
end

# I. delayed van der Pol (Section 4.1): Ψ(z) = [−z, 1; −a − b e^{−zτ}, μ − z],
#    f(z) = z² − μ z + a + b e^{−zτ}
let a = 0.9, b = 0.1, μ = 1.0, τ = 2.0
    f(z) = z^2 - μ * z + a + b * exp(-z * τ)
    f′(z) = 2z - μ - b * τ * exp(-z * τ)
    λ̄ = newton(f, f′, 0.5218379440444361 + 0.8056375092559337im)
    println("I.   λ̄ = ", λ̄, ",  |f(λ̄)| = ", abs(f(λ̄)))
    save("vdp_lambda.bin", real(λ̄), imag(λ̄))
end

# II. SEIR (Section 4.2): f(z) = (z + μ)(z + γ + μ − ν e^{−μτ} e^{−zτ})
let μ = 0.1, ν = 0.5, γ = 0.1, τ = 1.0
    f(z) = (z + μ) * (z + γ + μ - ν * exp(-μ * τ) * exp(-z * τ))
    f′(z) = (z + γ + μ - ν * exp(-μ * τ) * exp(-z * τ)) + (z + μ) * (1 + τ * ν * exp(-μ * τ) * exp(-z * τ))
    λ̄ = newton(f, f′, 0.1784703538844722)
    println("II.  λ̄ = ", λ̄, ",  |f(λ̄)| = ", abs(f(λ̄)))
    save("seir_lambda.bin", λ̄)
end

# III. Mackey–Glass, ρ = 3 (Section 4.3): f(z) = z (z + γ − β e^{−zτ})
let γ = 1.0, β = 2.0, τ = 2.0
    f(z) = z * (z + γ - β * exp(-z * τ))
    f′(z) = (z + γ - β * exp(-z * τ)) + z * (1 + β * τ * exp(-z * τ))
    λ̄ = newton(f, f′, 0.2393001697495648)
    println("III. λ̄ = ", λ̄, ",  |f(λ̄)| = ", abs(f(λ̄)))
    save("mg_lambda.bin", λ̄)
end
