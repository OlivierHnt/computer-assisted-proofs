# Statement and verification conditions

**Theorem (spiderweb central configurations with equal masses).** Let
$\lambda = -1$, $m_0 = 0$, and, for $\ell \ge 2$, let $\mu_\ell$ be the binary64
number nearest to $1/\ell$. For every number of circles $n \in \{5, \dots, 100\}$
and every number of half-lines $\ell \in \{18, \dots, 200\}$, the map
$f : \mathcal R^n \to \mathbb{R}^n$ of Section 6 with masses
$m_1 = \dots = m_n = \mu_\ell$ has a zero $\tilde r(n,\ell)$ with
$\|\tilde r(n,\ell) - \bar x(n,\ell)\|_\infty \le \rho_1(n,\ell)$, where
$\bar x(n,\ell)$ is the stored numerical zero (`data/radii.bin`) and
$[\rho_1(n,\ell), \rho_2(n,\ell)] \subset (0, 10^{-5}]$ is the interval of
existence printed by the checker; this zero is unique in the closed ball of
radius $\rho$ around $\bar x(n,\ell)$ for every $\rho \in [\rho_1, \rho_2]$.
Since the ball of radius $\rho_* = 10^{-5}$ around $\bar x(n,\ell)$ lies in
$\mathcal R^n = \{0 < r_1 < \dots < r_n\}$, there is an $n \times \ell$
spiderweb central configuration (no central mass) with circles of radii
$\tilde r_1 < \dots < \tilde r_n$ and equal masses. By the scaling (2.1) of
Section 2.1 (equivalence of $(\mathbf r, m, \lambda)$ under homotheties and
changes of units), the same holds for any common positive mass and any
$\lambda < 0$.

Reference: O. Hénot, C. Rousseau, *Spiderweb central configurations*,
Qual. Theory Dyn. Syst. 18 (2019), 1135–1160, doi:10.1007/s12346-019-00330-y;
arXiv:1810.09915. Section, theorem and lemma numbers below are those of the
arXiv source (`article_HR_final.tex`); the journal numbering may differ. The
statement above is the computer-assisted part of Section 6.2 ("Numerical
experimentations with circles of equal mass"), for the range of the authors'
script (`main.jl`: $n \in 5{:}100$, $\ell \in 18{:}200$, $m_i = 1/\ell$). The
cases $\ell \le 9$ (all $n$), $\ell \in \{10, \dots, 18\}$ with $n$ small
(Theorem 4.1) and $n \le 4$ (Theorem 5.1) are proved analytically in the paper
and are not part of this package.

## Setting

- Unknowns: the radii $r = (r_1, \dots, r_n)$ of the circles; the $\ell$ bodies
  on circle $i$ sit at angles $\theta_k = 2\pi k/\ell$, $k = 0, \dots, \ell-1$.
- Map (Section 6, $f$; first written in Section 2.4, eq. (2.6)): for
  $i = 1, \dots, n$,
  $f_i(r) = \lambda r_i + \dfrac{m_i\,\zeta_\ell}{2^{3/2} r_i^2} + \dfrac{m_0}{r_i^2}
  + \sum_{j \ne i} m_j \sum_{k=0}^{\ell-1} g(r_i, r_j, \theta_k)$, with
  $g(a, b, \theta) = \dfrac{a - b\cos\theta}{(a^2 + b^2 - 2ab\cos\theta)^{3/2}}$
  and $\zeta_\ell = \sum_{k=1}^{\ell-1} (1 - \cos\theta_k)^{-1/2}$. It is
  $C^\infty$ on $\{r_i > 0,\ r_i \ne r_j\} \supset \mathcal R^n$.
- $Df$ and $D^2 f$: the explicit formulas of Section 6.1 ("Computation of the
  bounds"), rederived for this package and checked against finite differences
  of $f$ in 256-bit arithmetic (relative agreement $10^{-13}$, for equal and
  for unequal masses). $f_i$ depends on $r_i$ and on the other $r_j$, so
  $D^2 f_i$ has only the entries $\partial^2_{ii} f_i$,
  $\partial^2_{ij} f_i = \partial^2_{ji} f_i$ and $\partial^2_{jj} f_i$
  ($j \ne i$); they are the three cases of the paper, the third being the second
  with $i \leftrightarrow j$. $D^2 f$ is stored as the $n \times n^2$ matrix
  $B = (B_1 \mid \dots \mid B_n)$ of Section 6.1,
  $B[k, (i-1)n + j] = \partial^2_{ji} f_k$. Code: `f_Df_D2f` in
  `checker/src/spiderweb.jl`.
- Norm: $\|\cdot\|_\infty$ on $\mathbb{R}^n$ and the induced operator norms
  (maximal row sum, also for the $n \times n^2$ matrix $AB$).
- Constants: $\lambda = -1$, $m_0 = 0$, $m_i = \mu_\ell$ (`fill(1/l, n)`),
  $\rho_* = 10^{-5}$ (`r_star`).
- Radii polynomial. The paper's Theorem 6.1 uses
  $p(\rho) = Z_2 \rho^2 - (1 - Z_0)\rho + Y_0$ with $Y_0 \ge \|A f(\bar x)\|$,
  $Z_0 \ge \|I - A\,Df(\bar x)\|$ and $Z_2 \rho \ge \|A[Df(c) - Df(\bar x)]\|$
  for all $c$ in the closed ball of radius $\rho \le \rho_*$. The checker
  certifies the Newton–Kantorovich polynomial
  $q(\rho) = Y_0 + (Z_0 - 1)\rho + Z_2 \rho^2/2$ instead (together with
  $Z_0 + Z_2 \rho < 1$), with the same $Z_2$: since $Z_2$ bounds
  $\|A\,D^2 f\|$ on the whole ball, Taylor's formula with integral remainder
  gives $\|A(f(x) - f(\bar x) - Df(\bar x)(x - \bar x))\| \le Z_2 \rho^2/2$,
  while the mean-value bound of Theorem 6.1 still gives the contraction constant
  $Z_0 + Z_2\rho$. As $p = q + Z_2\rho^2/2 \ge q$, the paper's condition
  $p(\rho_0) < 0$ is a valid but coarser sufficient condition.
- $Z_2$ bound (Lemma 6.3): $Z_2 = \sup_{b \in \bar B_{\rho_*}(\bar x)} \|A\,B(b)\|_\infty$.
  Derivation used here: $Df(c) - Df(\bar x) = \int_0^1 D^2 f(\bar x + t(c - \bar x))[c - \bar x]\,dt$,
  so $|(A[Df(c) - Df(\bar x)])_{aj}| \le \rho \sum_i \sup_b |\sum_k A_{ak}\,\partial^2_{ji} f_k(b)|$
  and the row sums give the lemma. The checker evaluates $A B(\bar B)$ with
  $\bar B$ the interval box, which encloses $\sup_b$.
- Arithmetic: IEEE binary64 interval arithmetic with outward rounding,
  IntervalArithmetic.jl 1.0.12 (pinned by RadiiPolynomial.jl 0.11.13;
  correctly rounded `+ - * / sqrt` and integer powers; `cos` via CRlibm;
  interval matrix products with the library default, Rump's midpoint–radius
  product `matmul = :fast`). Every operand is an interval, so the results
  carry the library's "guaranteed" flag, which the checker requires. The
  radii polynomial is handled by `RadiiPolynomial.interval_of_existence`.
  Julia 1.13.

## Inputs (in `data/`)

- Frozen: `radii.bin`, the numerical zeros $\bar x(n, \ell)$ of $f$, for
  $\ell = 18, \dots, 200$ and $n = 5, \dots, 100$ (17 568 cases). Producing
  them takes the authors' Newton continuation in $n$ and $\ell$
  (`proof_𝑛ℓ` of the original code), replayed by
  `julia --project=search search/produce.jl` (the `produce` command of
  `entry.json`), which is not part of the proof: any other data passing the
  conditions below would prove the theorem as well.
- Recomputed by the checker: $A = $ `inv(mid.(DF))`, the floating-point
  (LAPACK, binary64) inverse of the midpoint of the interval enclosure of
  $Df(\bar x)$, then converted to a thin interval matrix. Theorem 6.1 holds
  for any matrix $A$, so the conditions do not depend on its exact value, only
  on the same $A$ being used in $Y_0$, $Z_0$ and $Z_2$; its last bits may vary
  with the platform and BLAS, and so may the last digits of the printed bounds.
- Recomputed by the checker: the ball $\bar B = \bar x + [-\rho_*, \rho_*]^n$
  (interval sum, outward rounded) on which $D^2 f$ is enclosed.

File (raw little-endian IEEE binary64, no header, 922 320 numbers =
7 378 560 bytes): `radii.bin` holds, for $\ell = 18, 19, \dots, 200$ (outer
loop) and $n = 5, 6, \dots, 100$ (inner loop), the $n$ radii
$\bar x_1(n,\ell) < \dots < \bar x_n(n,\ell)$. The block of $(n, \ell)$ starts
at number $(\ell - 18)\cdot 5040 + \sum_{n'=5}^{n-1} n'$. Its SHA-256 digest is
in `data/hashes.json`; the checker verifies it before reading the file.

## Conditions to verify

The checker takes one case: `julia --project=checker checker/verify.jl <n> <l>`
with the parameters $n \in \{5, \dots, 100\}$ (circles) and
$\ell \in \{18, \dots, 200\}$ (half-lines), declared in `entry.json` as
`"cases": {"n": "5:100", "l": "18:200"}`. For the case $(n, \ell)$ it
evaluates in interval arithmetic (`check_case` in `checker/src/radii_check.jl`),
with $\bar x = \bar x(n,\ell)$, $A$ and $\bar B$ as above, the enclosures

- $Y_0 \ni \|A\, f(\bar x)\|_\infty$,
- $Z_0 \ni \|I - A\, Df(\bar x)\|_\infty$,
- $Z_2 \ni \|A\, B(\bar B)\|_\infty$ (Lemma 6.3),

prints their suprema, and checks two conditions, labelled
`1 (n,l) = (n,l)` and `2 (n,l) = (n,l)`; it exits 0 iff both hold.

1. **Ball in the domain.** $\bar B \subset \mathcal R^n$:
   $\inf \bar B_1 > 0$ and $\sup \bar B_i < \inf \bar B_{i+1}$ for
   $i = 1, \dots, n-1$. (Theorem 6.1 needs $\bar B_{\rho_*}(\bar x) \subset U$
   with $f \in C^1(U)$; here $U = \mathcal R^n$, and the zero inherits the
   ordering $0 < \tilde r_1 < \dots < \tilde r_n$.) The printed number is
   $\min(\inf \bar B_1, \min_i (\inf \bar B_{i+1} - \sup \bar B_i))$.
2. **Radii polynomial (Theorem 6.1, Newton–Kantorovich form).**
   `interval_of_existence(sup Y₀, sup Z₀, sup Z₂, ρ*)` succeeds and returns
   a guaranteed interval $I = [\rho_1, \rho_2] \subset [0, \rho_*]$ on which
   $q(\rho) = Y_0 + (Z_0 - 1)\rho + Z_2 \rho^2/2 \le 0$ and $Z_0 + Z_2 \rho < 1$
   (the library's documented guarantee; $\rho_2 = \min(\rho_*, (1 - Z_0)/Z_2)$
   rounded down). Success also entails $Y_0 \ge 0$, $0 \le Z_0 < 1$ and
   $0 \le Z_2 < \infty$ (the library rejects unbounded or negative bounds, and
   the empty interval), and $\rho_1 > 0$ when $Y_0 > 0$. The printed enclosure
   is $I$.

All conditions hold for every case $\Rightarrow$ Theorem. The theorem is the
conjunction of the $96 \times 183 = 17\,568$ statements "f has a unique
zero in $\bar B_\rho(\bar x(n,\ell))$ for $\rho \in I(n,\ell)$", one per
case, and the proof is complete iff every run of the checker exits 0. For one
case, let $\rho \in I$ and $T(x) = x - A f(x)$: for $x, y$ in the closed ball of radius $\rho$,
$\|T(x) - T(y)\| \le \sup_c \|I - A\,Df(c)\|\,\|x - y\| \le (Z_0 + Z_2\rho)\|x - y\|$
with $Z_0 + Z_2 \rho < 1$ (mean value bound, Lemma 6.3); and
$\|T(x) - \bar x\| \le \|A f(\bar x)\| + \|(I - A\,Df(\bar x))(x - \bar x)\| + \|A(f(x) - f(\bar x) - Df(\bar x)(x - \bar x))\| \le Y_0 + Z_0\rho + Z_2\rho^2/2 = q(\rho) + \rho \le \rho$
(Taylor remainder, same $Z_2$). Banach's fixed point theorem gives a unique
fixed point $\tilde r$ in the ball; $Z_0 < 1$ makes $A\,Df(\bar x)$, hence $A$,
invertible, so $f(\tilde r) = 0$. This is the conclusion of Theorem 6.1 under
a weaker hypothesis ($q(\rho) \le 0$ and $Z_0 + Z_2\rho < 1$ instead of
$p(\rho) < 0$). Condition 1 places the zero in $\mathcal R^n$.

**Corrections to the 2019 code (this package).** The checker no longer
reproduces `spiderweb.radii` of the original code: (i) the diagonal entries
$\partial^2_{ii} f_i$ of the 2019 `iF_iDF_iDDF` were assigned (`=`) instead of
accumulated (`+=`) inside a loop, which dropped the terms of the inner circles
$j < i$ (relative error up to $7.8$ on those entries at $(n,\ell) = (5,18)$,
checked against finite differences); (ii) the mixed derivatives
$\partial^2_{ik} f_k$ were filled from $\partial^2_{kk} f_i$, by an identity
exact only for equal masses, and are now computed from their own formula, so
`f_Df_D2f` is correct for any masses; (iii) the floating-point root test on
$p$ (`Polynomials.roots`, smallest root $\le \rho_*$) is replaced by the
certified `interval_of_existence` on the Newton–Kantorovich polynomial $q$. With the published $D^2 f$ the printed
$Z_2$ was not a proven bound of Lemma 6.3; with the corrected one it is, and
all cases run here still pass (the suprema of $Z_2$ change by a few per cent
at most; see README.md).

**Caveat (unchanged).** The full range of 17 568 cases (about 8 h of
computation on one core, see README.md) has not been run in this packaging. The runs recorded in README.md cover
$n \in 5{:}20$, $\ell \in 18{:}200$ and $n \in 21{:}100$,
$\ell \in \{18, 109, 200\}$ (3 168 cases) and all pass. The theorem as stated
follows from the conditions only once `checker/verify.jl <n> <l>` has exited 0
for every $(n, \ell)$ of `cases`.
