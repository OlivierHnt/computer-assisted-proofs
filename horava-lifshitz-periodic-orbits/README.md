# Bianchi VIII and IX periodic orbits in Hořava–Lifshitz cosmologies — proof package

For $v = 0$, the vacuum Bianchi models of Hořava–Lifshitz gravity have periodic
orbits of Bianchi type VIII and IX along which $N_1N_2N_3 = 10^k$, for each
$k = 1, 0, -1, \dots, -7$ (Theorem 4.1 of the paper; 18 cases). This package
records the computer-assisted proof as statement + inputs + checker; see
`statement.md` for the theorem and the conditions.

Check the proof (Julia 1.13, as recorded in `checker/Manifest.toml`):

```
julia --project=checker checker/verify.jl <case>      # case = 1, …, 18
```

`entry.json` declares `cases = 1:18`: the checker is run once per case and
exits with code 0 iff the conditions of that case hold; the theorem is proved
when all 18 runs exit 0. The runs are independent and can be executed in
parallel.

The checker runs on the current releases of its dependencies
(RadiiPolynomial.jl 0.11.13, IntervalArithmetic.jl 1.0.12). The original
proof ran on Julia 1.9.3 with RadiiPolynomial.jl 0.7.2 and
IntervalArithmetic.jl 0.20.9; this package's checker is a line-by-line
translation of that code to the new APIs (see `statement.md`).

Each run verifies the SHA-256 of the frozen inputs in `data/`, then prints
the conditions `c.0`–`c.12` of its case `c`. It recomputes the approximate
inverse (up to 5128 × 5128) and evaluates interval matrix products of that
size: the 15 cases with `n ≤ 192` take about 30 minutes in total with two
threads (IntervalArithmetic's midpoint–radius product). Cases 9, 17 and 18
(`n = 256`, 256-bit inputs, where the library falls back to entrywise
products) take several hours each, single-threaded, and 2–3.5 GB of memory.

Status of the runs with this checker (Julia 1.13, current library versions):
cases 1–8 and 10–16 pass all their conditions; cases 9, 17 and 18 have not
been rerun to completion with it. The previous version of this checker, a
line-by-line port on the original library versions (Julia 1.9.3,
RadiiPolynomial 0.7.2, IntervalArithmetic 0.20.9), passed all 18 cases with
output identical to the original code's.

Regenerate the inputs (optional, not part of the proof; does not reproduce the
frozen files bit for bit):

```
julia --project=search search/produce.jl
julia --project=checker tools/hash_data.jl
```

K. E. M. Church, O. Hénot, P. Lappicy, J.-P. Lessard, H. Sprink, *Periodic
orbits in Hořava-Lifshitz cosmologies*, Gen. Relativ. Gravit. 55 (2023), 2,
[doi:10.1007/s10714-022-03054-8](https://doi.org/10.1007/s10714-022-03054-8);
[arXiv:2203.03763](https://arxiv.org/abs/2203.03763).
Original code: [OlivierHnt/HLPeriodicOrbits_VIII_IX.jl](https://github.com/OlivierHnt/HLPeriodicOrbits_VIII_IX.jl).
