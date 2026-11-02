# Localized radial solutions of semilinear elliptic systems — proof package

Localized radial solutions of the cubic Klein–Gordon equation on ℝ³ (two of
them), the Swift–Hohenberg equation on ℝ² and a three-component
FitzHugh–Nagumo system on ℝ², each within an explicit C⁰ distance of a
Taylor–Chebyshev approximation; this package records the computer-assisted
proof as statement + inputs + checker, see `statement.md`.

Check the proof:

```
julia --project=checker checker/verify.jl {case}
```

The checker is run per case, `julia --project=checker checker/verify.jl {case}`
with `case` one of `kg0`, `kg1`, `sh`, `fhn` (Klein–Gordon m = 0 and m = 1,
Swift–Hohenberg, FitzHugh–Nagumo); it exits 0 iff the conditions of that case
hold, and the theorem needs every case. An automated check runs the single
case `test_case` of `entry.json` (`"test_case": {"case": "kg0"}`), or the first
declared value when there is none. See
[capsule-math/spec](https://github.com/capsule-math/spec), §4–§5. The cases
`kg1` and `fhn` need a large server (about 0.4 TB and 1 TB of memory).

J. B. van den Berg, O. Hénot, J.-P. Lessard, *Constructive proofs for
localised radial solutions of semilinear elliptic systems on ℝ^d*,
Nonlinearity 36 (2023) 6476, [doi:10.1088/1361-6544/ad0278](https://doi.org/10.1088/1361-6544/ad0278);
[arXiv:2211.16445](https://arxiv.org/abs/2211.16445).
