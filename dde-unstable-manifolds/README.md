# Validated unstable manifolds of DDEs — proof package

Validated error bounds for the Taylor parameterization of the unstable manifold
of an equilibrium of three delay differential equations (a delayed van der Pol
system, an SEIR model, the Mackey–Glass equation with ρ = 3); this package
records the computer-assisted proof as statement + inputs + checker, see
`statement.md`.

Check the proof:

```
julia --project=checker checker/verify.jl
```

O. Hénot, J.-P. Lessard, J. D. Mireles James, *Parameterization of unstable
manifolds for DDEs: formal series solutions and validated error bounds*,
J. Dyn. Differ. Equ. 34 (2022), 1285–1324,
[doi:10.1007/s10884-021-10002-8](https://doi.org/10.1007/s10884-021-10002-8).
