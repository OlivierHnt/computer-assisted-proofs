# Writes data/*.bin; not part of the proof. The approximate solutions x̄, the
# constants and the enclosures were computed by the authors and are stored as
# JLD2 files in OlivierHnt/EllipticPDELocalizedRadialSol.jl (commit c0fda16).
# This script loads them into field-for-field mirror types (module Legacy),
# writes their binary64 fields as they are (signed zeros included), and checks
# the result against data/hashes.json and the checker's reader.
#
# Usage (from the package root; the `produce` command of entry.json):
#   julia --project=search search/convert_jld2.jl [dir]
# where dir holds the four .jld2 files; without it they are downloaded and
# their SHA-256 checked against the values below.

import Pkg
Pkg.instantiate()

using Downloads, SHA
import JLD2
using RadiiPolynomial
include(joinpath(@__DIR__, "..", "checker", "src", "EllipticPDELocalizedRadialSol.jl"))
using .EllipticPDELocalizedRadialSol

module Legacy # the types of the JLD2 files, field for field
    struct Interval{T}; lo::T; hi::T; end
    struct ParameterSpace end
    struct Taylor; order::Int; end
    struct Chebyshev; order::Int; end
    struct CartesianPower{S}; space::S; n::Int; end
    struct CartesianProduct{T}; spaces::T; end
    struct Sequence{S,V}; space::S; coefficients::V; end
    struct LinearOperator{D,C,M}; domain::D; codomain::C; coefficients::M; end
    struct KleinGordon{T}; d::Int; β₁::T; β₂::T; end
    struct SwiftHohenberg{T}; d::Int; β₁::T; β₂::T; β₃::T; β₄::T; end
    struct FitzHughNagumo{T}; d::Int; ϵ::T; β₁::T; β₂::T; β₃::T; β₄::T; end
end

const TYPEMAP = Dict{String,Any}(
    "IntervalArithmetic.Interval" => Legacy.Interval,
    "RadiiPolynomial.ParameterSpace" => Legacy.ParameterSpace,
    "RadiiPolynomial.Taylor" => Legacy.Taylor,
    "RadiiPolynomial.Chebyshev" => Legacy.Chebyshev,
    "RadiiPolynomial.CartesianPower" => Legacy.CartesianPower,
    "RadiiPolynomial.CartesianProduct" => Legacy.CartesianProduct,
    "RadiiPolynomial.Sequence" => Legacy.Sequence,
    "RadiiPolynomial.LinearOperator" => Legacy.LinearOperator,
    "EllipticPDELocalizedRadialSol.KleinGordon" => Legacy.KleinGordon,
    "EllipticPDELocalizedRadialSol.SwiftHohenberg" => Legacy.SwiftHohenberg,
    "EllipticPDELocalizedRadialSol.FitzHughNagumo" => Legacy.FitzHughNagumo,
)

const COMMIT = "c0fda16ef4d5a7110d8b74ced2da7a17a5db78ab"
const URL = "https://raw.githubusercontent.com/OlivierHnt/EllipticPDELocalizedRadialSol.jl/$COMMIT/data/"
const JLD2_SHA256 = Dict(
    "klein_gordon_m_0" => "7847f4dcd47e7abb2ef2377a82e32e55f90a5e6f4a791d829f39b469a5994c88",
    "klein_gordon_m_1" => "e23e5f55bbac21588b0f1ca97a25bcd9225a34692c6a03344fd52a16708f8779",
    "swift_hohenberg"  => "fa36f4b54cc1ea512d842df867a832f283b3ed63be588a240da33b54dd753689",
    "fitzhugh_nagumo"  => "74908d0951a494ffdfc32278e20f5d4bbc801acc2f8f2ba5d77b71b01dfb9b75",
)
const EQ_FIELDS = Dict(KleinGordon => (:β₁, :β₂), SwiftHohenberg => (:β₁, :β₂, :β₃, :β₄), FitzHughNagumo => (:ϵ, :β₁, :β₂, :β₃, :β₄))

srcdir = length(ARGS) ≥ 1 ? ARGS[1] : mktempdir()
datadir = joinpath(@__DIR__, "..", "data")
hashes = isfile(joinpath(datadir, "hashes.json")) ? read(joinpath(datadir, "hashes.json"), String) : ""

# raw binary64 fields of the loaded objects
raw_interval(x) = (x.lo, x.hi)
raw_cinterval(z) = (z.re.lo, z.re.hi, z.im.lo, z.im.hi)
raw_entry(x, cplx) = cplx ? raw_cinterval(x) : raw_interval(x)
raw_float(x::Float64) = (x,)
raw_float(z::ComplexF64) = (real(z), imag(z))
_diag(D) = hasproperty(D, :diag) ? D.diag : [D[i,i] for i ∈ 1:size(D, 1)]

for case ∈ CASES
    src = joinpath(srcdir, case.name * ".jld2")
    if !isfile(src)
        println("downloading ", case.name, ".jld2")
        Downloads.download(URL * case.name * ".jld2", src)
    end
    bytes2hex(open(sha256, src)) == JLD2_SHA256[case.name] || error("sha256 mismatch for $src")

    d = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
        JLD2.load(src; typemap = TYPEMAP)
    end
    eq = d["equation"]
    eq.d == case.d && d["n_T_pad"] == case.n_T_pad && d["n_C_pad"] == case.n_C_pad || error("integer constants differ for $(case.name)")

    vals = Float64[]
    append!(vals, Float64[d["l"], d["r_star"], d["L"], d["ℒ_x"], d["ℒ_y"], d["ν"], d["ϱ"]])
    for f ∈ EQ_FIELDS[case.eqtype]; append!(vals, raw_interval(getproperty(eq, f))); end
    foreach(cᵢ -> append!(vals, raw_interval(cᵢ)), d["c"].coefficients)
    foreach(λᵢ -> append!(vals, raw_entry(λᵢ, case.complex)), _diag(d["Λ"].coefficients))
    foreach(γᵢⱼ -> append!(vals, raw_entry(γᵢⱼ, case.complex)), d["Γ"].coefficients) # column-major
    foreach(xᵢ -> append!(vals, raw_float(xᵢ)), d["x"].coefficients)

    out = joinpath(datadir, case.name * ".bin")
    open(io -> write(io, htol.(vals)), out, "w")

    # checks: size, reader, recorded hash
    back = read_case(out, case)
    same = coefficients(back.x) == d["x"].coefficients &&
        all(f -> (inf(getproperty(back.equation, f)), sup(getproperty(back.equation, f))) == raw_interval(getproperty(eq, f)), EQ_FIELDS[case.eqtype])
    same || error("$(case.name): the checker does not read back the JLD2 values")
    h = bytes2hex(open(sha256, out))
    println("wrote ", relpath(out), " (", filesize(out), " bytes), sha256 ", h,
        occursin(h, hashes) ? " = data/hashes.json" : " (NOT the recorded hash)")
end
