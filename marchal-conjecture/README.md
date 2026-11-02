# Marchal's conjecture — proof package

The $P_{12}$ family of the three-body problem connects the Lagrange triangle to
the figure-eight choreography; this package records the computer-assisted proof
as statement + inputs + checker, see `statement.md`.

Check the proof:

```
julia --project=checker checker/verify.jl
```

R. Calleja, C. García-Azpeitia, O. Hénot, J.-P. Lessard, J. D. Mireles James,
*From the Lagrange triangle to the figure eight choreography: Proof of Marchal's
conjecture*, Trans. Amer. Math. Soc. (2026), [doi:10.1090/tran/9779](https://doi.org/10.1090/tran/9779);
[arXiv:2406.17564](https://arxiv.org/abs/2406.17564).
