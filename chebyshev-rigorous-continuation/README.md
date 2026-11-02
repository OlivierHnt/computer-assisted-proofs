# Rigorous Chebyshev continuation of the Cahn–Hilliard steady-state manifold and SKT branches — proof package

Existence, local uniqueness and continuity of solution branches via Chebyshev
series in the parameter: a two-parameter manifold of Cahn–Hilliard steady states,
four branches of Shigesada–Kawasaki–Teramoto steady states (fourteen arcs that
join), and the two toy examples of the paper (cube root, square root). This
package records the computer-assisted proof as statement + inputs + checker; see
`statement.md`.

Check the proof:

```
julia --project=checker checker/verify.jl {example}
```

The checker is run per case, `julia --project=checker checker/verify.jl {example}`
with `example` one of `cube_root`, `square_root`, `cahn_hilliard`, `skt`; it exits
0 iff the conditions of that case hold, and the theorem needs every case. An
automated check runs the single case `test_case` of entry.json
(`"test_case": {"example": "cahn_hilliard"}`), or the first declared value when
there is none. See [capsule-math/spec](https://github.com/capsule-math/spec), §4–§5.

M. Breden, O. Hénot, *Efficient Rigorous Continuation via Chebyshev Series
Expansion I*, preprint (2026), [arXiv:2609.22689](https://arxiv.org/abs/2609.22689).
