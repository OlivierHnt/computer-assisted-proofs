# ------------------------------------------------------------------
# Reading and writing the frozen inputs.
#
# One file per case, data/<TYPE>_N1N2N3_<c>.bin, holding 7 + 4d numbers, in
# this order:
#
#     inf C₁, sup C₁, inf C₂, sup C₂, inf ν, sup ν, R,
#     then for each of the d coefficients of x̄ (in the order of
#     coefficients(x̄)): inf re, sup re, inf im, sup im
#
# with d = 5(2n+1) + 3 for x̄ ∈ Fourier(n, 1)^5 × ScalarSpace()^3.
#
# Number formats (no header, little-endian):
#   :f64   — IEEE binary64, 8 bytes per number;
#   :mp256 — a 256-bit binary floating-point number (MPFR / BigFloat with
#            precision 256), 48 bytes per number:
#              UInt64 s, Int64 E, UInt64 M₀, M₁, M₂, M₃
#            value = (-1)^s · M · 2^(E - 256), M = M₀ + M₁ 2^64 + M₂ 2^128 + M₃ 2^192,
#            with 2^255 ≤ M < 2^256 for a nonzero value and M = 0, E = 0 for ±0.
# ------------------------------------------------------------------

const MP_PREC = 256

nwords(n::Int) = 7 + 4 * (5 * (2n + 1) + 3)

bytes_per_number(fmt::Symbol) = fmt === :f64 ? 8 : fmt === :mp256 ? 48 : error("unknown format $fmt")

expected_filesize(n::Int, fmt::Symbol) = nwords(n) * bytes_per_number(fmt)

# --- one number ------------------------------------------------------

write_number(io, v::Float64) = write(io, htol(v))

function write_number(io, v::BigFloat)
    precision(v) == MP_PREC || error("BigFloat of precision $(precision(v)) ≠ $MP_PREC")
    isfinite(v) || error("non-finite BigFloat")
    s = UInt64(signbit(v))
    if iszero(v)
        E, M = Int64(0), big(0)
    else
        E = Int64(exponent(v) + 1)                       # |v| = M 2^(E-256), M ∈ [2^255, 2^256)
        M = BigInt(setprecision(() -> ldexp(abs(v), MP_PREC - E), BigFloat, MP_PREC))
        big(2)^(MP_PREC-1) ≤ M < big(2)^MP_PREC || error("unnormalised mantissa")
    end
    write(io, htol(s), htol(E))
    mask = big(typemax(UInt64))
    for i ∈ 0:3
        write(io, htol(UInt64((M >> (64i)) & mask)))
    end
    return 48
end

read_number(io, ::Val{:f64}) = ltoh(read(io, Float64))

function read_number(io, ::Val{:mp256})
    s = ltoh(read(io, UInt64))
    E = ltoh(read(io, Int64))
    M = big(0)
    for i ∈ 0:3
        M += big(ltoh(read(io, UInt64))) << (64i)
    end
    (s == 0 || s == 1) || error("bad sign word")
    M == 0 || big(2)^(MP_PREC-1) ≤ M < big(2)^MP_PREC || error("unnormalised mantissa")
    return setprecision(BigFloat, MP_PREC) do
        v = iszero(M) ? BigFloat(0) : ldexp(BigFloat(M), Int(E) - MP_PREC) # exact: M has 256 bits
        s == 1 ? -v : v
    end
end

# --- one case ----------------------------------------------------------

function write_case(path, x::Sequence, C₁::Interval, C₂::Interval, ν::Interval, R::Real)
    open(path, "w") do io
        for v ∈ (inf(C₁), sup(C₁), inf(C₂), sup(C₂), inf(ν), sup(ν), R)
            write_number(io, v)
        end
        for c ∈ coefficients(x)
            write_number(io, inf(real(c))); write_number(io, sup(real(c)))
            write_number(io, inf(imag(c))); write_number(io, sup(imag(c)))
        end
    end
    return path
end

"""
    read_case(path, n, fmt) -> (x, C₁, C₂, ν, R)

Read one case: `x` a `Sequence` on
`Fourier(n, 1.0)^5 × ScalarSpace()^3` with `Complex{Interval{T}}`
coefficients, `C₁, C₂, ν :: Interval{T}` and `R :: T`, where `T = Float64`
for `fmt = :f64` and `T = BigFloat` (precision 256) for `fmt = :mp256`.
"""
function read_case(path, n::Int, fmt::Symbol)
    filesize(path) == expected_filesize(n, fmt) || error("$path has the wrong size for n = $n, format $fmt")
    space_x = Fourier(n, 1.0)^5 × ScalarSpace()^3
    d = dimension(space_x)
    open(path, "r") do io
        rd() = read_number(io, Val(fmt))
        C₁ = interval(rd(), rd())
        C₂ = interval(rd(), rd())
        ν  = interval(rd(), rd())
        R  = rd()
        v = map(1:d) do _
            a, b, c, e = rd(), rd(), rd(), rd()
            complex(interval(a, b), interval(c, e))
        end
        eof(io) || error("trailing bytes in $path")
        return Sequence(space_x, v), C₁, C₂, ν, R
    end
end
