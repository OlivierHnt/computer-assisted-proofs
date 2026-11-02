# Spiderweb central configurations — proof package

For every $5 \le n \le 100$ circles and $18 \le \ell \le 200$ half-lines there is
an $n \times \ell$ spiderweb central configuration with equal masses, proved
by the radii-polynomial method around a numerical zero, one case at a time.
This package records the computer-assisted proof as statement + inputs +
checker; see `statement.md` for the theorem and the conditions.

Check the proof (Julia 1.13). The checker takes one case $(n, \ell)$:

```
julia --project=checker checker/verify.jl 20 109       # the case (n, l) = (20, 109); exit 0 iff it passes
```

The theorem is the conjunction over all $96 \times 183 = 17\,568$ cases
$n \in 5{:}100$, $\ell \in 18{:}200$: `entry.json` declares
`"run": "julia --project=checker checker/verify.jl {n} {l}"` with
`"cases": {"n": "5:100", "l": "18:200"}`, i.e. one run per element of the
Cartesian product, and the proof is complete iff every run exits 0. The cases
are independent and can run in parallel (each process uses one thread). One
case is the function `check_case(n, l, x̄)` of `checker/src/radii_check.jl`;
the map, its Jacobian and its second derivative are `f_Df_D2f` in
`checker/src/spiderweb.jl`; the constants, the input loading (with the hash
check) and the printing are in `checker/src/common.jl`.

Run time: about 5 s per process start (package loading, hash check of the
7 MB input) plus 0.1 s for the smallest cases and 8 s for
$(n, \ell) = (100, 200)$ on one core (Apple M-series), dominated by the
$O(n^2 \ell)$ interval evaluation of $f, Df, D^2 f$ and the $n \times n^2$
interval product $A \cdot D^2 f$; the computation of all cases is estimated at
about 8 h on one core (from the 240 cases with $n \ge 21$ below: 1.9 s each on
average), plus the start-up of each process. The subset $n \in 5{:}20$,
$\ell \in 18{:}200$ plus $n \in 21{:}100$, $\ell \in \{18, 109, 200\}$
(3 168 cases) was run for this package by a scratch script looping over the
same `check_case` in two single-threaded processes; every condition passes,
in 8 min 9 s of wall time.

Regenerate the inputs (optional, not part of the proof; about 1 h):

```
julia --project=search search/produce.jl
julia --project=checker tools/hash_data.jl
```

`produce` on Julia 1.13 does not reproduce `data/radii.bin` bit for bit (the
LAPACK solves inside Newton's method differ in the last bits from the
OpenBLAS of Julia 1.7, relative differences below $10^{-13}$). The frozen file
comes from the same script with the authors' 2019 package versions
(IntervalArithmetic 0.15.2, Polynomials 0.5.2) on Julia 1.7.2.

## Differences from the 2019 code

The checker was first a line-for-line port of `spiderweb.radii` (2019,
`spiderweb.zip`). It now departs from it in three places, each rederived
from the paper and checked numerically (see `statement.md`, "Corrections to
the 2019 code"):

1. **Diagonal of $D^2 f$.** In the 2019 `iF_iDF_iDDF`, `iDDFixjxj[i,j] = 6*im0/ixri4`
   (for $i = j$) sits inside the loop over $j$ and resets $\partial^2_{ii} f_i$,
   discarding the terms of the circles $j < i$. The published $Z_2$ was
   therefore not a proven bound. With the complete $D^2 f$, $\sup Z_2$ moves
   by at most a few per cent on the cases compared ($+0.4\,\%$ at $(5, 18)$,
   $-5.7\,\%$ at $(10, 50)$, $-2.2\,\%$ at $(20, 200)$, $-0.7\,\%$ at
   $(50, 109)$, $-0.2\,\%$ at $(100, 18)$) and every case run still passes.
2. **Mixed derivatives.** The 2019 code stored $\partial^2_{ik} f_k$ as
   $\partial^2_{kk} f_i$, which is exact only for equal masses
   ($\partial^2_{ik} f_k / m_i = \partial^2_{kk} f_i / m_k$). They are now
   computed from their own formula (the paper's second case), so the code is
   correct for any masses. Checked against 256-bit finite differences of $f$
   for equal and unequal masses (relative agreement $10^{-13}$); with the
   2019 storage and unequal masses exactly the $2n(n-1)$ mixed entries are
   wrong.
3. **Radii polynomial.** The floating-point root test (`Polynomials.roots`,
   smallest root $\le \rho_*$) is replaced by
   `RadiiPolynomial.interval_of_existence(Y₀, Z₀, Z₂, ρ*)`, which certifies
   $q(\rho) = Y_0 + (Z_0 - 1)\rho + Z_2\rho^2/2 \le 0$ and $Z_0 + Z_2\rho < 1$
   on the returned interval. $Z_2$ is passed as it is: it bounds
   $\|A\,D^2 f\|$ on the whole ball, so it bounds the second-order Taylor
   remainder by $Z_2\rho^2/2$; the paper's polynomial
   $Z_2\rho^2 - (1 - Z_0)\rho + Y_0$ (Theorem 6.1, no $1/2$) comes from the
   coarser mean-value estimate and is a valid but weaker sufficient condition.
   The upper end of the interval is now $\min(\rho_*, (1 - Z_0)/Z_2)$ instead
   of the second root of the paper's $p$ (the same quantity to leading
   order); the lower end $\rho_1$ (the radius of the theorem) differs from
   the paper's smallest root by $Z_2\rho_1^2/2$, i.e. relatively by
   $Z_2\rho_1/2$: $1\,\%$ at $n = 100$, and below the outward rounding of the
   quadratic formula ($\approx 10^{-8}$ relative) for small $n$.

Besides, the checker also verifies that the ball of radius $\rho_*$ lies in
$\{0 < r_1 < \dots < r_n\}$ (condition 1), the whole computation is done
with intervals (no bare floating-point operand, the library's "guaranteed"
flag is required), the ball is built by an outward-rounded interval sum, and
the interval matrix products use the library default (Rump's algorithm,
`matmul = :fast`): compared with the entrywise product (`:slow`) the bounds
differ in the 15th digit only. The identifiers $𝑛, ℓ$ of the 2019 code are
written `n`, `l`.

With the 2019 packages, `norm(v, Inf)` of an interval vector (Julia's generic
code, which compares intervals with `>`, and IntervalArithmetic 0.15, where
`>` is false for overlapping intervals) could return an entry that is not the
largest, so the 2019 $Y_0$ was not always an upper bound (96 of the 3 168
cases compared, by up to 23 %). IntervalArithmetic 1.0.12 returns the true
maximum; all those cases pass with it.

O. Hénot, C. Rousseau, *Spiderweb central configurations*, Qual. Theory Dyn.
Syst. 18 (2019), 1135–1160, [doi:10.1007/s12346-019-00330-y](https://doi.org/10.1007/s12346-019-00330-y);
[arXiv:1810.09915](https://arxiv.org/abs/1810.09915).
Original code: `spiderweb.zip` (`main.jl`, `spiderweb.jl`, 2019-03-26) from
http://www.dms.umontreal.ca/~rousseac/spiderweb.zip.
