# Used by checker/verify.jl: the constants of the proof, the frozen inputs,
# and the printing of the conditions of one case. Requires RadiiPolynomial, LinearAlgebra, SHA, JSON and the
# files spiderweb.jl, radii_check.jl of this directory to be loaded.

# constants of the proof (paper, Sections 6.1–6.2; original main.jl)
const λ  = -1.0           # common factor of the central configuration
const m0 = 0.0            # central mass
const N_RANGE = 5:100     # number of circles n
const L_RANGE = 18:200    # number of half-lines l
mass(n, l) = fill(1/l, n) # m_i = 1/l (binary64), i = 1, …, n
const r_star = 1e-5       # ρ*, radius of the ball on which D²f is enclosed

# One thread each for LAPACK (the inverse of mid Df) and for the interval
# matrix products (IntervalArithmetic's default `matmul = :fast`, Rump's
# algorithm on a directed-rounding OpenBLAS, which would otherwise use every
# core), so that cases can run in parallel processes. The results do not
# depend on these settings, only the run time does.
BLAS.set_num_threads(1)
IntervalArithmetic.configure(; nthreads = 1)

const verdicts = Bool[]
function condition(label, ok, enclosure)
    push!(verdicts, ok)
    println(ok ? "pass " : "fail ", label, ": ", enclosure)
    flush(stdout)
    return ok
end

# Verify the SHA-256 of every frozen file in data/ against data/hashes.json
# and return the function x̄(n, l) reading data/radii.bin (layout: for
# l ∈ L_RANGE (outer), n ∈ N_RANGE (inner), the n radii x̄(n, l)[1:n]).
function load_inputs(datadir)
    hashes = JSON.parsefile(joinpath(datadir, "hashes.json"))
    for f ∈ sort(collect(keys(hashes)))
        p = joinpath(datadir, f)
        isfile(p) || error("$f is missing: run `julia --project=checker tools/fetch_data.jl` first")
        bytes2hex(open(sha256, p)) == hashes[f]["sha256"] || error("sha256 mismatch for $f: not the recorded input")
        println("pass sha256 ", f)
    end
    nfloats = length(L_RANGE) * sum(N_RANGE)
    filesize(joinpath(datadir, "radii.bin")) == 8nfloats || error("radii.bin has the wrong size")
    raw = ltoh.(reinterpret(Float64, read(joinpath(datadir, "radii.bin"))))
    offset = Dict{Tuple{Int,Int},Int}()
    let o = 0
        for l ∈ L_RANGE, n ∈ N_RANGE
            offset[(n, l)] = o
            o += n
        end
    end
    println("loaded x̄(n, l) for n ∈ $N_RANGE, l ∈ $L_RANGE ($nfloats binary64)")
    return (n, l) -> raw[offset[(n, l)] .+ (1:n)]
end

# Check one case and print its two conditions (statement.md); returns true iff both hold.
function report_case(n, l, x̄)
    res = check_case(n, l, x̄; λ, m0, m = mass(n, l), r_star)
    ok1 = condition("1 (n,l) = ($n,$l): ball B̄_ρ*(x̄) ⊂ R^n", res.in_Rn,
                    "min(inf ball₁, gaps between consecutive balls) = $(res.gap)")
    ok2 = condition("2 (n,l) = ($n,$l): radii polynomial (Theorem 6.1)", res.success,
                    "Y₀ ≤ $(sup(res.Y0)), Z₀ ≤ $(sup(res.Z0)), Z₂ ≤ $(sup(res.Z2)); interval of existence [$(inf(res.ie)), $(sup(res.ie))]")
    return ok1 && ok2
end
