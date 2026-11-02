# Reader of the frozen inputs, written by search/convert_jld2.jl. One file per
# case, little-endian IEEE binary64, no header:
#
#   1. scalars        l, r_star, L, ℒ_x, ℒ_y, ν, ϱ                      7 × binary64
#   2. equation       the interval parameters of the equation, in the
#                     order of the struct fields (KG: β₁, β₂; SH: β₁…β₄;
#                     FHN: ϵ, β₁…β₄), each (inf, sup)                     2 per parameter
#   3. c              q intervals (inf, sup)                              2q
#   4. Λ              the q diagonal entries                              2q (real) / 4q (complex)
#   5. Γ              q × q, column-major                                 2q² (real) / 4q² (complex)
#   6. x̄              coefficients of x̄ = (η̄, ϕ̄, v̄, w̄), in this order:
#                     η̄ (q), ϕ̄ (q), v̄₁…v̄_q (each n_T_num+1 Taylor coefficients,
#                     order 0 first), w̄₁…w̄_{1+2q} (each n_C_num+1 Chebyshev
#                     coefficients, order 0 first)                       1 (real) / 2 (complex) per coefficient
#
# Real interval: (inf, sup). Complex interval: (inf re, sup re, inf im, sup im).
# Complex float: (re, im).
# The integers (d, q, truncation orders, paddings) are constants of the proof
# and are written below, not stored.

const CASES = [
    (name = "klein_gordon_m_0", eqtype = KleinGordon,    d = 3, q = 1, n_T_num = 60,  n_C_num = 85,  n_T_pad = 100, n_C_pad = 200,   complex = false),
    (name = "klein_gordon_m_1", eqtype = KleinGordon,    d = 3, q = 1, n_T_num = 150, n_C_num = 200, n_T_pad = 100, n_C_pad = 15000, complex = false),
    (name = "swift_hohenberg",  eqtype = SwiftHohenberg, d = 2, q = 2, n_T_num = 30,  n_C_num = 120, n_T_pad = 100, n_C_pad = 230,   complex = true),
    (name = "fitzhugh_nagumo",  eqtype = FitzHughNagumo, d = 2, q = 3, n_T_num = 60,  n_C_num = 250, n_T_pad = 100, n_C_pad = 10000, complex = false),
]

_x_space(case) = ScalarSpace()^case.q × ScalarSpace()^case.q × Taylor(case.n_T_num)^case.q × Chebyshev(case.n_C_num)^(1 + 2case.q)

_nparams(case) = fieldcount(case.eqtype) - 1 # all fields but d

function expected_filesize(case)
    q = case.q
    s = case.complex ? 2 : 1
    n = 7 + 2_nparams(case) + 2q + 2s*q + 2s*q^2 + s*dimension(_x_space(case))
    return 8n
end

_read_f64(io) = ltoh(read(io, Float64))
_read_interval(io, ::Val{false}) = (a = _read_f64(io); b = _read_f64(io); interval(a, b))
_read_interval(io, ::Val{true}) = (a = _read_interval(io, Val(false)); b = _read_interval(io, Val(false)); complex(a, b))
_read_float(io, ::Val{false}) = _read_f64(io)
_read_float(io, ::Val{true}) = (a = _read_f64(io); b = _read_f64(io); complex(a, b))

function read_case(path, case)
    filesize(path) == expected_filesize(case) || error("$path has the wrong size")
    q, cplx = case.q, Val(case.complex)
    open(path, "r") do io
        l, r_star, L, ℒ_x, ℒ_y, ν, ϱ = ntuple(_ -> _read_f64(io), 7)
        params = [_read_interval(io, Val(false)) for _ ∈ 1:_nparams(case)]
        equation = case.eqtype(case.d, params...)
        c = Sequence(ScalarSpace()^q, [_read_interval(io, Val(false)) for _ ∈ 1:q])
        Λ = LinearOperator(ScalarSpace()^q, ScalarSpace()^q, Diagonal([_read_interval(io, cplx) for _ ∈ 1:q]))
        Γ = LinearOperator(ScalarSpace()^q, ScalarSpace()^q, reshape([_read_interval(io, cplx) for _ ∈ 1:q^2], q, q))
        X = _x_space(case)
        x = Sequence(X, [_read_float(io, cplx) for _ ∈ 1:dimension(X)])
        eof(io) || error("$path has trailing bytes")
        return (; equation, x, n_T_pad = case.n_T_pad, n_C_pad = case.n_C_pad, c, Λ, Γ, l, r_star, L, ℒ_x, ℒ_y, ν, ϱ)
    end
end
