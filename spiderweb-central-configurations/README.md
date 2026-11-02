# Spiderweb central configurations — proof package

For every $5 \le n \le 100$ circles and $18 \le \ell \le 200$ half-lines there is
an $n \times \ell$ spiderweb central configuration with equal masses; this
package records the computer-assisted proof as statement + inputs + checker,
see `statement.md`.

Check the proof:

```
julia --project=checker checker/verify.jl 20 109
```

The checker is run per case, `julia --project=checker checker/verify.jl {n} {l}`
with `"cases": {"n": "5:100", "l": "18:200"}` in `entry.json`; it exits 0 iff
the conditions of that case hold, and the theorem needs every case. An
automated check runs the single case `test_case` of `entry.json`,
$(n, \ell) = (20, 109)$. See
[capsule-math/spec](https://github.com/capsule-math/spec), §4–§5. All
$17\,568$ cases take about 8 hours of CPU.

O. Hénot, C. Rousseau, *Spiderweb central configurations*, Qual. Theory Dyn.
Syst. 18 (2019), 1135–1160, [doi:10.1007/s12346-019-00330-y](https://doi.org/10.1007/s12346-019-00330-y);
[arXiv:1810.09915](https://arxiv.org/abs/1810.09915).
