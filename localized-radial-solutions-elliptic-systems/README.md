# Localized radial solutions of semilinear elliptic systems — proof package

Localized radial solutions of the cubic Klein–Gordon equation on ℝ³ (two of
them), the Swift–Hohenberg equation on ℝ² and a three-component
FitzHugh–Nagumo system on ℝ², each within an explicit C⁰ distance of a
Taylor–Chebyshev approximation. This package records the computer-assisted
proof as statement + inputs + checker; see `statement.md` for the theorem and
the conditions.

Check the proof (Julia 1.13):

```
julia --project=checker checker/verify.jl
```

This runs the four cases. The Klein–Gordon m = 0 and Swift–Hohenberg
cases take minutes on a laptop; the Klein–Gordon m = 1 and FitzHugh–Nagumo
cases build dense interval matrices of size about 45 900 × 47 100 and
72 200 × 74 000 and need a large server (roughly 0.4 TB and 1 TB of memory,
extrapolated; the paper reports under 2 h and under 6 h on 32 threads with the
original versions). Single cases:
`checker/verify.jl kg0 sh` (names `kg0`, `kg1`, `sh`, `fhn`; a partial run
exits with code 3 when everything it evaluated passes).

Interval matrix products are BLAS products under directed rounding
(IntervalArithmetic's bundled OpenBLAS); `CHECKER_BLAS_THREADS=n` sets the
number of BLAS threads, and condition 0 tests that they honour the rounding mode.

Regenerate the inputs from the original JLD2 files (optional, not part of the
proof; downloads them from the original repository):

```
julia --project=search search/convert_jld2.jl
julia --project=checker tools/hash_data.jl
```

J. B. van den Berg, O. Hénot, J.-P. Lessard, *Constructive proofs for
localised radial solutions of semilinear elliptic systems on ℝ^d*,
Nonlinearity 36 (2023) 6476, [doi:10.1088/1361-6544/ad0278](https://doi.org/10.1088/1361-6544/ad0278);
[arXiv:2211.16445](https://arxiv.org/abs/2211.16445).
Original code: [OlivierHnt/EllipticPDELocalizedRadialSol.jl](https://github.com/OlivierHnt/EllipticPDELocalizedRadialSol.jl)
(Julia 1.9.3, RadiiPolynomial 0.7.1); this package runs it, ported, on Julia 1.13
and RadiiPolynomial 0.11.
