# Transverse heteroclinic orbit in the Shimizu–Morioka system — proof package

For $a = 3/4$, $b = 9/20$, the Shimizu–Morioka system has a transverse
heteroclinic orbit from the saddle-focus $c_1^\star \approx (0.67082, 0, 1)$ to the
saddle at the origin; this package records the computer-assisted proof as
statement + inputs + checker, see `statement.md`.

Check the proof:

```
julia --project=checker checker/verify.jl
```

O. Hénot, A. Takayasu, *Computer-Assisted Proofs in Dynamical Systems: A Case
Study of a Heteroclinic Orbit in the Shimizu–Morioka System*,
[arXiv:2605.07500](https://arxiv.org/abs/2605.07500).
