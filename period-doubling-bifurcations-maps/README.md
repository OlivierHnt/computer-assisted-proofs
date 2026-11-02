# Period-doubling bifurcation in the Kot–Schaffer model — proof package

The 10-mode Galerkin projection of the Kot–Schaffer growth-dispersal model
undergoes a supercritical period-doubling bifurcation at
$\lambda_\star \approx 2.7042$, with an attracting 2-cycle for $\lambda > \lambda_\star$.
This package records the computer-assisted proof as statement + inputs +
checker; see `statement.md` for the theorem and the conditions.

Check the proof (Julia 1.13.3; with juliaup, `julia +1.13.3 ...` or the
`release` channel if it is 1.13.3):

```
julia --project=checker checker/verify.jl
```

Regenerate the input (optional, not part of the proof):

```
julia --project=search search/find_candidate.jl
julia --project=checker tools/hash_data.jl
```

O. Hénot, J.-P. Lessard, *Rigorous verification of period-doubling bifurcations
in maps*, [arXiv:2610.22986](https://arxiv.org/abs/2610.22986) (2026), Theorem 6.1.
Original code: [OlivierHnt/PeriodDoublingMaps.jl](https://github.com/OlivierHnt/PeriodDoublingMaps.jl).
