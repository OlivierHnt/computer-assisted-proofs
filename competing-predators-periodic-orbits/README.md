# Stable periodic orbits of two competing predators — proof package

For the predator–prey system of Butler and Waltman with two predators competing
for one prey (Holling type II, $a_1 = 10$, $a_2 = 41$, $d_1 = 0.8$, $d_2 = 0.5$,
$m_j = y_j = \gamma = 1$), a global family of stable positive periodic orbits,
parametrised by the carrying capacity $\kappa$, connects the boundary limit
cycle in $Q_1 = \{X_2 = 0\}$ to the one in $Q_2 = \{X_1 = 0\}$. This package
records the computer-assisted proof as statement + inputs + checker; see
`statement.md` for the theorem and the conditions.

Check the proof (Julia 1.13; the checker pins RadiiPolynomial 0.11.13 and
IntervalArithmetic 1.0.12, see `statement.md` for how it relates to the
original run on Julia 1.11.4 / RadiiPolynomial 0.8.20):

```
julia --project=checker checker/verify.jl
```

Regenerate the inputs (optional, not part of the proof):

```
julia --project=search search/produce.jl
julia --project=checker tools/hash_data.jl
```

K. E. M. Church, J.-Y. Dai, O. Hénot, P. Lappicy, N. Vassena,
*Global continuation of stable periodic orbits in systems of competing predators*,
[arXiv:2504.03058](https://arxiv.org/abs/2504.03058).
Original code: [OlivierHnt/BHWPeriodicCoexistence.jl](https://github.com/OlivierHnt/BHWPeriodicCoexistence.jl).
