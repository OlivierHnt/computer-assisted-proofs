# Bianchi VIII and IX periodic orbits in Hořava–Lifshitz cosmologies — proof package

For $v = 0$, the Bianchi models of Hořava–Lifshitz gravity have periodic
orbits of Bianchi type VIII and IX along which $N_1N_2N_3 = 10^k$, for each
$k = 1, 0, -1, \dots, -7$ (Theorem 4.1 of the paper); this package records the
computer-assisted proof as statement + inputs + checker, see `statement.md`.

Check the proof:

```
julia --project=checker checker/verify.jl {case}
```

The checker is run per case, `julia --project=checker checker/verify.jl {case}`
with `case = 1, …, 18` (`cases` of `entry.json`); it exits 0 iff the
conditions of that case hold, and the theorem needs every case. An automated
check runs the single case `test_case` of `entry.json`, or the first declared
value when there is none (here case 1). See
[capsule-math/spec](https://github.com/capsule-math/spec), §4–§5. Cases 9, 17
and 18 take several hours each and 2–3.5 GB of memory.

K. E. M. Church, O. Hénot, P. Lappicy, J.-P. Lessard, H. Sprink, *Periodic
orbits in Hořava-Lifshitz cosmologies*, Gen. Relativ. Gravit. 55 (2023), 2,
[doi:10.1007/s10714-022-03054-8](https://doi.org/10.1007/s10714-022-03054-8);
[arXiv:2203.03763](https://arxiv.org/abs/2203.03763).
