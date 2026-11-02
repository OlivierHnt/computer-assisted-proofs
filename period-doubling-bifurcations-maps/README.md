# Period-doubling bifurcation in the Kot–Schaffer model — proof package

The 10-mode Galerkin projection of the Kot–Schaffer growth-dispersal model
undergoes a supercritical period-doubling bifurcation at
$\lambda_\star \approx 2.7042$, with an attracting 2-cycle for $\lambda > \lambda_\star$;
this package records the computer-assisted proof as statement + inputs +
checker, see `statement.md`.

Check the proof:

```
julia --project=checker checker/verify.jl
```

O. Hénot, J.-P. Lessard, *Rigorous verification of period-doubling bifurcations
in maps*, [arXiv:2610.22986](https://arxiv.org/abs/2610.22986) (2026), Theorem 6.1.
