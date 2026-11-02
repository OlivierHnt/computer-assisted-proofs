# Statement and verification conditions

**Theorem (Bianchi type VIII and IX periodic orbits).** The Hořava–Lifshitz
Bianchi model with $v = 0$,

$$
\begin{aligned}
\Sigma_+' &= 2\big[(N_3-N_2)^2 - N_1(2N_1-N_2-N_3)\big], &
N_1' &= 4\Sigma_+ N_1,\\
\Sigma_-' &= 2\sqrt3\,(N_3-N_2)(N_1-N_2-N_3), &
N_2' &= -2(\Sigma_+ + \sqrt3\,\Sigma_-)N_2,\\
& & N_3' &= -2(\Sigma_+ - \sqrt3\,\Sigma_-)N_3,
\end{aligned}
$$

possesses nontrivial periodic orbits along which the constraint
$\Sigma_+^2 + \Sigma_-^2 + N_1^2 + N_2^2 + N_3^2 - 2N_1N_2 - 2N_2N_3 - 2N_3N_1 = 1$
is satisfied. In particular, for each $k \in \{1, 0, -1, \dots, -7\}$ there are
periodic solutions along which the product $N_1N_2N_3$ is constant and equal to
$10^k$: one of Bianchi type VIII ($N_1, N_2, N_3$ nonzero, not all of the same
sign) and one of Bianchi type IX ($N_1, N_2, N_3$ nonzero, all of the same
sign).

Reference: K. E. M. Church, O. Hénot, P. Lappicy, J.-P. Lessard, H. Sprink,
*Periodic orbits in Hořava-Lifshitz cosmologies*, Gen. Relativ. Gravit. 55
(2023), 2, doi:10.1007/s10714-022-03054-8; arXiv:2203.03763. This is Theorem
4.1 there. Numbers below (Theorem 4.1, Lemma 4.2, Section 4.1, Proposition B.1,
Appendix B) are those of the arXiv version (v2 source), which defines every
object named here.

The package covers 18 cases, numbered as the checker prints them:

| case | type | $N_1N_2N_3$ | file | $n$ | format | $\nu$ | $R$ |
|---|---|---|---|---|---|---|---|
| 1 | VIII | $10^{1}$  | `VIII_N1N2N3_10.bin`   | 32  | mp256 | 1.3  | 1e-3 |
| 2 | VIII | $10^{0}$  | `VIII_N1N2N3_1.bin`    | 32  | f64   | 1.3  | 1e-3 |
| 3 | VIII | $10^{-1}$ | `VIII_N1N2N3_1e-1.bin` | 32  | f64   | 1.2  | 1e-3 |
| 4 | VIII | $10^{-2}$ | `VIII_N1N2N3_1e-2.bin` | 32  | f64   | 1.1  | 1e-3 |
| 5 | VIII | $10^{-3}$ | `VIII_N1N2N3_1e-3.bin` | 32  | f64   | 1.1  | 1e-4 |
| 6 | VIII | $10^{-4}$ | `VIII_N1N2N3_1e-4.bin` | 64  | f64   | 1.1  | 1e-4 |
| 7 | VIII | $10^{-5}$ | `VIII_N1N2N3_1e-5.bin` | 64  | f64   | 1.1  | 1e-4 |
| 8 | VIII | $10^{-6}$ | `VIII_N1N2N3_1e-6.bin` | 128 | f64   | 1.07 | 5e-7 |
| 9 | VIII | $10^{-7}$ | `VIII_N1N2N3_1e-7.bin` | 256 | mp256 | 1.05 | 1e-9 |
| 10 | IX  | $10^{1}$  | `IX_N1N2N3_10.bin`     | 32  | f64   | 1.1  | 1e-3 |
| 11 | IX  | $10^{0}$  | `IX_N1N2N3_1.bin`      | 64  | f64   | 1.15 | 1e-3 |
| 12 | IX  | $10^{-1}$ | `IX_N1N2N3_1e-1.bin`   | 64  | f64   | 1.1  | 1e-3 |
| 13 | IX  | $10^{-2}$ | `IX_N1N2N3_1e-2.bin`   | 64  | f64   | 1.07 | 1e-3 |
| 14 | IX  | $10^{-3}$ | `IX_N1N2N3_1e-3.bin`   | 128 | f64   | 1.03 | 1e-4 |
| 15 | IX  | $10^{-4}$ | `IX_N1N2N3_1e-4.bin`   | 192 | f64   | 1.03 | 1e-4 |
| 16 | IX  | $10^{-5}$ | `IX_N1N2N3_1e-5.bin`   | 192 | f64   | 1.03 | 1e-4 |
| 17 | IX  | $10^{-6}$ | `IX_N1N2N3_1e-6.bin`   | 192 | mp256 | 1.03 | 1e-9 |
| 18 | IX  | $10^{-7}$ | `IX_N1N2N3_1e-7.bin`   | 256 | mp256 | 1.03 | 1e-10 |

($\nu$ and $R$ are as stored in the files, see Inputs; for the mp256 cases
$\nu$ is the tightest 256-bit enclosure of the decimal value and $R$ the
binary64 value of the decimal, converted to 256 bits.)

## Setting (Section 4.1 and Appendix B)

- Auxiliary system: $2\pi$-periodic solutions of
  $\tfrac{d}{dt}\Sigma_+ = 2\gamma[(N_3-N_2)^2 - N_1(2N_1-N_2-N_3)] + \eta_1\Sigma_+$,
  $\tfrac{d}{dt}\Sigma_- = 2\gamma\sqrt3(N_3-N_2)(N_1-N_2-N_3)$,
  $\tfrac{d}{dt}N_1 = 4\gamma\Sigma_+N_1$,
  $\tfrac{d}{dt}N_2 = -2\gamma(\Sigma_+ + \sqrt3\Sigma_-)N_2 + \eta_2$,
  $\tfrac{d}{dt}N_3 = -2\gamma(\Sigma_+ - \sqrt3\Sigma_-)N_3$,
  with unknowns $\gamma, \eta_1, \eta_2$ (Section 4.1).
- Unknown: $x = (a_1, \dots, a_5, \gamma, \eta_1, \eta_2)$, $a_j$ the Fourier
  coefficients ($e^{ikt}$, $k \in \mathbb Z$) of
  $\Sigma_+, \Sigma_-, N_1, N_2, N_3$.
- Spaces and norms: $X_{\text{Fourier}} = \ell^1_\nu$,
  $|a| = \sum_k |a_k|\nu^{|k|}$; $X = X_{\text{Fourier}}^5 \times \mathbb C^3$
  with the max norm of the eight components (code:
  `NormedCartesianSpace((NormedCartesianSpace(ℓ¹(GeometricWeight(ν)), ℓ∞()), NormedCartesianSpace(ℓ∞(), ℓ∞())), ℓ∞())`).
  $\nu$ per case as in the table.
- Map: $F : X \to X$ of Section 4.1 (code: `F!`, `DF!` in
  `checker/src/proof.jl`): the five equations above with
  $\mathcal D(a)_k = ik\,a_k$, and the three scalar equations
  $\mathcal E(a_1) = 0$ (phase condition $\Sigma_+(0) = 0$),
  $(\Sigma^2 + \Omega_k)(\mathcal E(a_1), \dots, \mathcal E(a_5)) - C_1 = 0$,
  $\mathcal E(a_3)\mathcal E(a_4)\mathcal E(a_5) - C_2 = 0$,
  with $\mathcal E(a) = \sum_k a_k$, $C_1 = 1$, $C_2 = 10^k$.
- Truncation: $\bar x$ has Fourier order $n$ (table above), so
  $F(\bar x)$ has order $2n$; $DF(\bar x)$ is assembled from order $3n$ to
  order $2n$. The operator of Proposition B.1 is
  $\mathcal A = A\pi^{2n} + L|_{\pi^{\infty(2n)}}^{-1}\pi^{\infty(2n)}$, with
  $\|L|_{\pi^{\infty(2n)}}^{-1}\| \le 1/(2n+1)$.
- $A$: the floating-point inverse of the midpoint of $\pi^{2n}DF(\bar x)\pi^{2n}$
  (`inv(mid.(project(DF, space_F, space_F)))`, LAPACK in binary64). It is
  **recomputed by the checker** in binary64 from the frozen $\bar x$, exactly
  as in the original code, and not stored: it takes one LU factorisation to
  obtain (but up to 420 MB to store: a $5128 \times 5128$ complex matrix for
  $n = 256$), Proposition B.1 holds for any $A$, and the conditions below only
  use the $A$ actually computed, not its exact value.
- Bounds (Proposition B.1 and its proof, code `_prove_periodic_orbit`):
  - $Y = \|A F(\bar x)\|_X$;
  - $Z_1 = \|A\,DF(\bar x) - I\|_X + \tfrac{1}{2n+1}\,\zeta_1
    + \|A\|_X\,\nu^{-(3n+1)} \max\big(1,\ 2(|s_1|+|s_2|+3(|s_3|+|s_4|+|s_5|)),\ |s_4||s_5|+|s_3||s_5|+|s_3||s_4|\big)$,
    with $s_j = \mathcal E(\bar a_j)$ and $\zeta_1$ the bound on the
    nonlinear part of $DF(\bar x)$ (`opnorm_Df`);
  - $Z_2 = \big(\|A\|_X + \tfrac{1}{2n+1}\big)\big(\zeta_{\Sigma_+} + \zeta_{\Sigma_-} + \zeta_{N_1} + \zeta_{N_2} + \zeta_{N_3} + \zeta_\gamma + 1\big)$,
    the $\zeta$'s bounding the partial derivatives of $DF$ on the ball of
    radius $R$ around $\bar x$ (`opnorm_∂Σ₊DF`, …, `opnorm_∂γDF`).
- Symmetry: before anything else the checker applies, as the original code,
  the in-place update `x₁ⱼ .= (x₁ⱼ .+ conj.(view(x₁ⱼ, n:-1:-n))) ./ 2` to each
  of the five Fourier components, then replaces $(\bar a_j)_0$ and
  $\gamma, \eta_1, \eta_2$ by their real parts (interval arithmetic). The
  update runs over $k = -n, \dots, n$ in increasing order and reads entries
  already overwritten, i.e. $(\bar a_j)_k \leftarrow ((\bar a_j)_k + \overline{(\bar a_j)_{-k}})/2$
  with the *new* value of $(\bar a_j)_{-k}$ for $k > 0$ (checked against
  RadiiPolynomial 0.7.2 on the frozen data). Since the frozen data are not
  exactly conjugate-symmetric (asymmetry up to $10^{-12}$), the resulting
  $\bar x$ is not exactly the coefficient sequence of real functions, which
  the last part of Proposition B.1 assumes; condition 5 below closes this gap.
  $\bar x$ denotes this updated sequence from now on.
- Arithmetic: interval arithmetic with outward rounding, IntervalArithmetic.jl
  1.0.12 (binary64, correctly rounded basic operations, set-based flavour;
  interval matrix products in the default `matmul = :fast` mode, Rump's
  midpoint-radius product, where the original used entrywise products (for
  the five 256-bit cases the library falls back to entrywise products, since
  the fast mode exists only for binary16/32/64); and
  `power = :slow`, correctly rounded integer powers as in IntervalArithmetic
  0.20) under RadiiPolynomial.jl
  0.11.13, Julia 1.13.3. The original proof (Julia 1.9.3, RadiiPolynomial
  0.7.2, IntervalArithmetic 0.20.9) is the faithful version of this package;
  the code was translated API for API (see the header of
  `checker/src/proof.jl`). The checker prints each enclosure with its
  decoration and IntervalArithmetic 1.0 "guaranteed" flag; the flag is
  `false` whenever a non-interval number entered the computation, which here
  only happens with exact binary constants (the integer coefficients of $F$
  and of the bounds, the Fourier frequency $1$, the grid points $t$ and the
  radius $r_0$ used as exact floats), so it does not affect the enclosures.
  For the five mp256 cases, $F(\bar x)$
  and $Y$ are evaluated in 256-bit MPFR (BigFloat) interval arithmetic; $DF$,
  $A$, $Z_1$, $Z_2$ and the checks 6–11 are evaluated in binary64 on the
  outward-rounded binary64 enclosure of $\bar x, C_1, C_2, \nu, R$, and
  $Z_1, Z_2$ are converted back (outward) to 256-bit intervals. The radii
  polynomial uses `RadiiPolynomial.interval_of_existence` (0.11.13): on the
  upper bounds $y, z_1, z_2$ of $Y, Z_1, Z_2$ it encloses the smaller root
  in the cancellation-free form $r_1 = 2y/(1 - z_1 + \sqrt{(1-z_1)^2 - 2z_2y})$,
  and returns $[\sup r_1, \min(R, (1 - z_1)/z_2)]$ (right endpoint rounded
  down), after checking $z_1 + z_2\sup r_1 < 1$, the computation of
  RadiiPolynomial 0.7.2 (the original) up to rounding (case 1:
  $r_0 = 2.05\cdot10^{-21}$, original $2.1\cdot10^{-21}$). (RadiiPolynomial 0.11.0
  to 0.11.10 evaluated $(1 - Z_1 - \sqrt{\cdot})/Z_2$ on the enclosures, where the
  width of $Z_1$ enters with cancellation and gave $r_0 \approx 3.3\cdot10^{-15}$ in
  case 1; both are valid radii for Proposition B.1.)

## Inputs (in `data/`)

All 18 files are **frozen**: they are the original data of
HLPeriodicOrbits_VIII_IX.jl (`data/*.jld2`, commit 746d8e7), on which the
paper's proof was run, converted bit-exactly to raw binary (the reconstructed
objects were checked bitwise equal to the JLD2 ones). SHA-256 digests are in
`data/hashes.json`; the checker verifies them before reading the files.

Each file `<TYPE>_N1N2N3_<c>.bin` holds $7 + 4d$ numbers, no header, in this
order:

1. $\inf C_1, \sup C_1, \inf C_2, \sup C_2, \inf\nu, \sup\nu, R$;
2. for each of the $d = 5(2n+1) + 3$ coefficients of $\bar x$ (in the order
   $\Sigma_+$ modes $k = -n, \dots, n$, then $\Sigma_-$, $N_1$, $N_2$, $N_3$
   likewise, then $\gamma, \eta_1, \eta_2$): $\inf\operatorname{Re}$,
   $\sup\operatorname{Re}$, $\inf\operatorname{Im}$, $\sup\operatorname{Im}$.

$d = 328, 648, 1288, 1928, 2568$ for $n = 32, 64, 128, 192, 256$. Number
formats (little-endian):

- **f64**: IEEE binary64, 8 bytes per number.
- **mp256**: a binary floating-point number with a 256-bit significand (the
  BigFloat/MPFR numbers of the original data, which have exponents far outside
  the binary64 range), 48 bytes per number: `UInt64 s`, `Int64 E`, then the
  four `UInt64` limbs $M_0, M_1, M_2, M_3$ of
  $M = M_0 + 2^{64}M_1 + 2^{128}M_2 + 2^{192}M_3$; the value is
  $(-1)^s\,M\,2^{E-256}$, with $2^{255} \le M < 2^{256}$ for a nonzero value
  and $M = 0$, $E = 0$ for $\pm 0$ (sign $s$).

What the numbers are: $\bar x$ is a numerical approximate zero of $F$ (stored
as thin intervals; $\eta_1 = \eta_2 = 0$); $C_1 \ni 1$ and $C_2 \ni 10^k$ are
intervals (the tightest binary64 enclosures; for the mp256 cases converted to
256 bits); $\nu$ is the weight and $R$ the radius of the ball on which $Z_2$
is valid.

`search/produce.jl` (the `produce` command of `entry.json`) regenerates
inputs of the same form by Newton continuation in $N_1N_2N_3$; it is not part
of the proof, and its output is not bit-identical to the frozen data (the
original search code is not part of the original repository).

## Conditions to verify

The checker is run once per case $c = 1, \dots, 18$ (`cases` in
`entry.json`); each run prints the conditions $c.0, \dots, c.12$ of its case
in this order, each prefixed `pass` or `fail`, with the enclosure obtained,
and exits 0 iff they all hold. The theorem is the conjunction over the 18
cases. Each is a finite computation in outward-rounded interval
arithmetic on the data above (and on the matrix $A$ computed from it).

0. **Inputs.** $1 \in C_1$ and $10^k \in C_2$ (exact rational comparison of
   the stored endpoints), $\nu \ge 1$, $R > 0$.
1. **$Y$ (Proposition B.1).** $Y = \|A F(\bar x)\|_X$ is a finite bound.
2. **$Z_1$ (Proposition B.1).** The bound $Z_1$ above is finite.
3. **$Z_2$ (Proposition B.1).** The bound $Z_2$ above, valid on the closed ball
   of radius $R$, is finite.
4. **Radii polynomial (Proposition B.1).** The interval of existence returned
   by `interval_of_existence(Y, Z₁, Z₂, R)` is nonempty: every $r$ in it
   satisfies $r \le R$, $Y + (Z_1 - 1)r + \tfrac{Z_2}{2}r^2 \le 0$ and
   $Z_1 + Z_2 r < 1$ (upper bounds of $Y, Z_1, Z_2$ used). With
   $r_0$ its left endpoint, $F$ has a unique zero $\tilde x$ with
   $|\tilde x - \bar x|_X \le r_0$, unique in the ball of radius $\sup I$ for
   the interval of existence $I$ (claim (i) of Section 4.1). Checks 6–11 use
   $r_0$ (rounded up to binary64 for the mp256 cases) as
   the $C^0$ error bound of $\Sigma_+, N_1, N_3$, since
   $\sup_t|f(t)| \le |f|_{X_{\text{Fourier}}}$.
5. **The zero is real (Proposition B.1, last part).** With
   $\dagger x = \big((\overline{(a_j)_{-k}})_k, \overline{\gamma}, \overline{\eta_1}, \overline{\eta_2}\big)$,
   $r_0 + |\bar x - \dagger\bar x|_X \le \sup I$, where $I$ is the interval
   of existence of condition 4. Then $\dagger\tilde x$, which is a zero of $F$
   ($C_1 = 1$ and $C_2 = 10^k$ are real) at distance at most
   $r_0 + |\bar x - \dagger\bar x|_X$ from $\bar x$ ($\dagger$ is an
   isometry of $X$), equals $\tilde x$ by uniqueness of the zero in the ball
   of radius $\sup I$; so $\tilde x$ is real. (Not a check of the original
   code; see Symmetry above.)
6. **$\gamma \in (0, +\infty)$ (Lemma 4.2).** $\operatorname{Im}\bar\gamma = 0$
   and $\operatorname{Re}\bar\gamma > R \ge r_0$, so $\tilde\gamma > 0$.
7. **$\Sigma_+ \not\equiv 0$ (Lemma 4.2).** For some $t$ in the grid
   $0, 0.001, 0.002, \dots \le 2\pi$, the enclosure of the truncated series
   $\sum_{|k|\le n}(\bar a_1)_k e^{ikt}$, widened by $\pm r_0$, excludes $0$.
8. **$N_1 \not\equiv 0$ (Lemma 4.2).** Same for $\bar a_3$.
9. **$N_3 \not\equiv 0$ (Lemma 4.2).** Same for $\bar a_5$.
10. **$N_1$ has no zero on $[0, 2\pi]$ (Lemma 4.2, Appendix B).** With
   $\|N_1'\|_\infty \le 4(|\bar\gamma| + r_0)(|\bar a_1| + r_0)(|\bar a_3| + r_0)$,
   starting at $t_0 = 0$ and stepping
   $t_{j+1} = t_j + (|N_1(t_j)| - r_0)/\|N_1'\|_\infty$ (rounded down), the
   enclosure of $N_1(t_j) \pm r_0$ excludes $0$ at every step until
   $t_j \ge 2\pi$.
11. **$N_3$ has no zero on $[0, 2\pi]$ (Lemma 4.2, Appendix B).** Same, with
    $\|N_3'\|_\infty \le 2(|\bar\gamma| + r_0)(|\bar a_5| + r_0)(|\bar a_1| + r_0 + \sqrt3(|\bar a_2| + r_0))$.
12. **Bianchi type.** The enclosures of $N_1(0) \pm r_0$ and $N_3(0) \pm r_0$
    computed in 10 and 11 have a definite sign, and: for type IX,
    $N_1(0) > 0$ and $N_3(0) > 0$; for type VIII, not both positive. (Not a
    check of the original code, which only labels the cases; it reads off
    the enclosures already computed there.)

Condition $c.0$ ties the stored intervals to the stated constants (not a
check of the original code). Conditions $c.1$–$c.4$ are claim (i) of
Section 4.1, with $c.5$ for the reality of $\tilde x$; $c.6$–$c.11$ (checks
1–6 of the original code) are claim (ii). By Lemma 4.2, $\tilde\eta_1 = \tilde\eta_2 = 0$ and
$\tau \mapsto (\Sigma_+, \Sigma_-, N_1, N_2, N_3)(\tau/\tilde\gamma)$ is a
$2\pi\tilde\gamma$-periodic solution of the system of the theorem; the
constraint and $N_1N_2N_3 = 10^k$ hold at $t = 0$ (last two components of
$F$, for the exact values $C_1 = 1$, $C_2 = 10^k$ contained in the stored
intervals) and are conserved. Since $N_1$, $N_3$ have no zero and
$N_2 = 10^k/(N_1N_3)$, the signs of $N_1, N_2, N_3$ are constant, $N_2$ having
the sign of $N_1N_3$; condition $c.12$ then gives the Bianchi type (IX: all
positive; VIII: not all of the same sign).

All conditions hold ⇒ Theorem.
