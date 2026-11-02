# Statement and verification conditions

**Theorem.** Let $\Delta$ be the Laplacian on $\mathbb{R}^d$.

1. *(Klein–Gordon, Theorem 1.1.)* For $\beta_1 = \beta_2 = 1$, the equation
   $\Delta U - U + \beta_1 U^2 + \beta_2 U^3 = 0$ on $\mathbb{R}^3$ has two
   distinct localized radial solutions, within $C^0$-distance
   $2.9 \times 10^{-7}$ and $5.5 \times 10^{-6}$ of the numerical approximations
   $\bar u$ built from `data/klein_gordon_m_0.bin` and `data/klein_gordon_m_1.bin`.
2. *(Swift–Hohenberg, Theorem 1.2.)* For
   $(\beta_1, \beta_2, \beta_3, \beta_4) = (-\tfrac35, \sqrt6, -1, 1)$, the equation
   $-(\beta_4 + \Delta)^2 U + \beta_1 U + \beta_2 U^2 + \beta_3 U^3 = 0$ on
   $\mathbb{R}^2$ has a localized radial solution within $C^0$-distance
   $2.4 \times 10^{-5}$ of the approximation built from `data/swift_hohenberg.bin`.
   (The paper prints $\beta_3 = -\tfrac1{10}$; the data and the original code use
   $\beta_3 = -1$, which is what is proved here.)
3. *(FitzHugh–Nagumo, Theorem 1.3.)* For $\epsilon = \tfrac3{10}$ and
   $(\beta_1, \beta_2, \beta_3, \beta_4) = (\tfrac12, \tfrac12, 1, 3)$, the system
   $\epsilon^2 \Delta U_1 + U_1 - U_1^3 - \epsilon(\beta_1 + \beta_2 U_2 + \beta_3 U_3) = 0$,
   $\Delta U_2 + U_1 - U_2 = 0$, $\beta_4^2 \Delta U_3 + U_1 - U_3 = 0$ on
   $\mathbb{R}^2$ has a localized radial solution within $C^0$-distance
   $9.8 \times 10^{-7}$ of the approximation built from `data/fitzhugh_nagumo.bin`.

"Localized radial" means $U(x) = u(|x|)$ with $u(r) \to c$, $u'(r) \to 0$ as
$r \to \infty$, $c$ a hyperbolic zero of $\mathbf{N}$ (Section 1).

The four solutions are independent cases of the checker, declared as the
`cases` `kg0`, `kg1` (the two Klein–Gordon solutions, Theorem 1.1), `sh`
(Theorem 1.2) and `fhn` (Theorem 1.3) of `entry.json`, with `kg0` as
`test_case`; each theorem is proved by the run of its case(s).

Reference: J. B. van den Berg, O. Hénot, J.-P. Lessard, *Constructive proofs
for localised radial solutions of semilinear elliptic systems on
$\mathbb{R}^d$*, Nonlinearity 36 (2023) 6476, doi:10.1088/1361-6544/ad0278;
arXiv:2211.16445. Section, proposition, theorem and lemma numbers below are
those of arXiv **v2** (the original code's messages cite "Theorem 3.2" and
"Corollary 3.3", which are Theorem 3.3 and Corollary 3.5 of v2).

## Setting

Each case is an instance of $\Delta U + \mathbf{N}(U) = 0$, $U \in \mathbb{R}^q$
(KG: $q = 1$, $d = 3$; SH: $q = 2$, $d = 2$, written as a system for
$(U, (\beta_4 + \Delta)U)$; FHN: $q = 3$, $d = 2$), with $\mathbf{N}$ as in
Section 4 and `checker/src/{klein_gordon,swift_hohenberg,fitzhugh_nagumo}.jl`.

- Radial ODE and splitting (Section 3.1): $u(r) = v(r/\ell)$ on $[0, \ell r_*]$
  as a Taylor series $v$; on $[\ell r_*, \ell r_* + L]$ the autonomous system
  $w' = f(w)$, $w = (r^{-1}, u, u') \in \mathbb{R}^{1+2q}$, as Chebyshev series
  $w(s) = w_0 + 2\sum_{n\ge1} w_n T_n(s)$; beyond $r_0 = \ell r_* + L$ the local
  graph $\tilde\alpha$ of the centre-stable manifold of $(0, c, 0)$ in the
  coordinates $w = (0,c,0) + M(x, y, z)$ of Section 2, $M$ built from $\Lambda$, $\Gamma$.
- Space: $X = \mathbb{C}^q \times \mathbb{C}^q \times \mathscr{T}^q \times \mathscr{C}_\nu^{1+2q}$
  with $|\chi|_X = \max(|\eta|_\infty, |\phi|_\infty, \max_i |v_i|_{\mathscr{T}}, \max_i |w_i|_{\mathscr{C}_\nu})$,
  $|a|_{\mathscr{T}} = \sum_n |a_n|$, $|a|_{\mathscr{C}_\nu} = |a_0| + 2\sum_{n\ge1}|a_n|\nu^n$
  (code: `ℓ¹()`, `ℓ¹(GeometricWeight(ν))`, `NormedCartesianSpace(…, ℓ∞())`).
- Map: $F(\chi; \tilde\alpha)$ of Section 3.1 (the zero-finding problem),
  and $F(\chi; 0)$ without the $\tilde\alpha$ terms; code `F`, derivative `DF`
  (`checker/src/proof.jl`, verbatim from the original).
- Constants per case (written in `checker/src/io.jl`, the reals stored in the data):

  | case | $d$ | $q$ | $n_{\mathscr T}$ = num + pad | $n_{\mathscr C}$ = num + pad | $\ell$ | $r_*$ | $L$ | $\nu$ | $\varrho$ | $\mathcal L_x$ | $\mathcal L_y$ |
  |---|---|---|---|---|---|---|---|---|---|---|---|
  | 1 KG $m=0$ | 3 | 1 | 60 + 100 | 85 + 200 | 0.5 | 0.9 | 14 | 1.02 | 1e-7 | 1 | 0.05 |
  | 2 KG $m=1$ | 3 | 1 | 150 + 100 | 200 + 15000 | 0.1 | 0.9 | 15 | 1.02 | 5e-6 | 1 | 0.05 |
  | 3 SH | 2 | 2 | 30 + 100 | 120 + 230 | 1 | 0.9 | 31.1 | 1.05 | 5e-6 | 1 | 0.03 |
  | 4 FHN | 2 | 3 | 60 + 100 | 250 + 10000 | 1 | 0.9 | 42.1 | 1.01 | 1e-7 | 1 | 0.02 |

  $\delta = (\ell r_* + L)^{-1}$, $\mu = |\bar\eta|_\infty + \varrho$,
  $\hat\lambda = \min_i \operatorname{Re}\Lambda_{ii}$; the polynomial order of
  $\mathbf{N}$ and $f$ is $K = 3$.
- Arithmetic: IEEE binary64 interval arithmetic with outward rounding,
  IntervalArithmetic.jl 1.0.12 (CRlibm 1.0.2) through RadiiPolynomial.jl 0.11.13,
  Julia 1.13.3. Interval matrix products use IntervalArithmetic's `:fast` mode
  (Rump's midpoint–radius algorithm, BLAS products under directed rounding in
  the bundled OpenBLAS 0.3.34 built to propagate the rounding mode to all its
  threads); condition 0 tests it. The closed-form enclosures of condition
  $k$.0c are computed in 256-bit BigFloat intervals. The proof was first run
  with the original versions (Julia 1.9.3, RadiiPolynomial 0.7.1,
  IntervalArithmetic 0.20.9); the port keeps the original mathematics: the
  radii-polynomial interval is computed by the algorithm of RadiiPolynomial
  0.7.1 (reimplemented, `checker/src/compat.jl`), and $|\cdot|_{\mathscr B(\mathbb C^q,\mathbb C^q)}$
  is the matrix ∞-norm (maximum row sum) as in the original.

## Inputs (in `data/`)

All inputs are frozen; they are the JLD2 data of the original repository
(`OlivierHnt/EllipticPDELocalizedRadialSol.jl`, commit `c0fda16`), on which the
paper's proofs were run, converted bit for bit to raw binary. The code that
computed them (Newton's method on the truncated $F(\cdot; 0)$, the eigenproblem)
is not in the original repository; `search/convert_jld2.jl` (the `produce`
command) only regenerates `data/*.bin` from the JLD2 files and asserts the
bitwise round trip (the files are exactly the binary64 fields of the JLD2
objects, signed zeros included). Any other data passing the conditions would prove the
theorem as well.

One file per case: `klein_gordon_m_0.bin`, `klein_gordon_m_1.bin`,
`swift_hohenberg.bin`, `fitzhugh_nagumo.bin`. Raw little-endian IEEE binary64,
no header; real interval = (inf, sup); complex interval = (inf re, sup re, inf im, sup im);
complex float = (re, im). In order:

1. 7 scalars: $\ell$, $r_*$, $L$, $\mathcal L_x$, $\mathcal L_y$, $\nu$, $\varrho$.
2. The equation parameters as real intervals, in the order KG $(\beta_1, \beta_2)$;
   SH $(\beta_1, \beta_2, \beta_3, \beta_4)$; FHN $(\epsilon, \beta_1, \beta_2, \beta_3, \beta_4)$.
3. $c$: $q$ real intervals (the equilibrium).
4. $\Lambda$: the $q$ diagonal entries (real intervals; complex intervals for SH).
5. $\Gamma$: $q \times q$, column-major (real; complex for SH).
6. $\bar x = (\bar\eta, \bar\phi, \bar v, \bar w)$ as floats (real; complex for SH):
   $\bar\eta$ ($q$), $\bar\phi$ ($q$), then $\bar v_1, \dots, \bar v_q$ with
   $n_{\mathscr T,\mathrm{num}}+1$ Taylor coefficients each (order 0 first), then
   $\bar w_1, \dots, \bar w_{1+2q}$ with $n_{\mathscr C,\mathrm{num}}+1$ Chebyshev
   coefficients each (order 0 first; $w_1 \approx r^{-1}$, $w_{2..1+q} \approx u$,
   $w_{2+q..1+2q} \approx u'$).

File sizes: 2704, 6184, 11080 and 15944 bytes. SHA-256 digests are in
`data/hashes.json`; the checker verifies them before reading the files.

Recomputed by the checker (not inputs, as in the original code): the
approximate inverse $A = $ `inv(mid(DF(x̄;0)))` (floating point, square part of
$DF$), and the $4\times4$ matrices $A_{\mathrm{eig}}$ of Appendix B. The
conditions hold for whatever $A$ is used ($A$ is only required to be injective,
which $Z_1 < 1$ implies), so their exact values do not matter.

## Conditions to verify

Each is a finite computation in outward-rounded arithmetic on the data above.
The checker prints them in this order, each prefixed `pass` or `fail`, with the
enclosure obtained; conditions marked *(added)* are not in the original code
(they check hypotheses the original takes from the paper).

0. **BLAS rounding** *(added)*: the interval product of two fixed
   $1024 \times 1024$ matrices computed with the configured number of threads
   is bitwise equal to the single-thread one and encloses the exact entries
   sampled (computed in 256-bit BigFloat).

For each case $k = 1$ (KG $m=0$), 2 (KG $m=1$), 3 (SH), 4 (FHN):

- $k$.0a *(added)* Constants satisfy the hypotheses of Theorem 3.3:
  $n_{\mathscr T} \ge 2$, $n_{\mathscr C} \ge 1$, $0 < r_* \le e^{-1/(n_{\mathscr T}+1)}$,
  $\nu > 1$, $\ell, L, \mathcal L_x, \mathcal L_y, \varrho > 0$.
- $k$.0b *(added)* $\bar x = \mathcal S(\bar x)$ (Lemma 3.4): for KG and FHN
  everything is real; for SH $\bar\eta_1 = \overline{\bar\eta_2}$, $\bar\phi, \bar v, \bar w$
  have zero imaginary part, and $\Lambda_{11} = \overline{\Lambda_{22}}$,
  $\Gamma_{\cdot 1} = \overline{\Gamma_{\cdot 2}}$ (as intervals).
- $k$.0c *(added)* The stored intervals enclose the theorem's parameters, an
  exact zero $c$ of $\mathbf N$ and an exact solution of
  $-D\mathbf N(c)\Gamma = \Gamma\Lambda^2$ with $\operatorname{Re}\Lambda_{ii} > 0$
  (Section 2, Section 4): KG $c = 0$, $\Lambda = \Gamma = 1$ exactly; SH $c = 0$,
  $\Lambda_{11} = \sqrt{-\beta_4 - i\sqrt{-\beta_1}}$, $\Gamma = \begin{pmatrix}1&1\\-i\sqrt{-\beta_1}&i\sqrt{-\beta_1}\end{pmatrix}$
  (closed form, Section 4.2, where the paper's $\sqrt{\beta_3}$ reads $\sqrt{\beta_1}$);
  FHN $c_i \ni -(5+\sqrt{145})/20$ (a root of $\upsilon - \upsilon^3 - \epsilon(\beta_1 + (\beta_2+\beta_3)\upsilon)$,
  Section 4.3) and each eigenpair $(\Lambda_{jj}, \Gamma_{\cdot j})$ by the
  Newton–Kantorovich argument of Appendix B ($Z_{\mathrm{eig}} < 1$,
  $\bar\rho_{\mathrm{eig}} = Y_{\mathrm{eig}}/(1-Z_{\mathrm{eig}}) \le 10^{-10}$,
  ball inside the stored intervals), the three $\Lambda_{jj}$ pairwise disjoint.
- $k$.1a–1c **Proposition 2.2**, the three inequalities (first; the
  $\mathcal L_y$ one; the $\mathcal L_x$ one) with $\delta$, $\mu$ above and the
  bound $\hat\psi$ of the original code (`ψ̂` in the equation files; Section 4 formulas).
- $k$.1d *(added)* The same three inequalities with
  $\hat\psi = \tfrac12 B(s)\,|\Lambda^{-1}\Gamma^{-1}|\,|\Gamma|$,
  $s = |\Gamma|(1+\mathcal L_y)\mu$, where $B(s)$ bounds the single nonzero entry
  of $D\mathbf N(c+\zeta) - D\mathbf N(c)$ for $|\zeta|_\infty \le s$ using $|c|$
  (KG: $(2|\beta_1| + 3|\beta_2|(2|c_1|+s))s$; SH: $(2|\beta_2| + 3|\beta_3|(2|c_1|+s))s$;
  FHN: $3\epsilon^{-2}(2|c_1|+s)s$). See the note below.
- $k$.2 **$Y$** (Theorem 3.3, Appendix A): $Y = |A\pi F(\bar x;0)|_X + \max\{\ell^2|\pi^{\infty(n_{\mathscr T}-2)}\mathbf N(\bar v)|/((n_{\mathscr T}+1)(n_{\mathscr T}+d-1)),\ L(\nu+\nu^{-1})|\pi^{\infty(n_{\mathscr C}-1)}f(\bar w)|/(4(n_{\mathscr C}+1))\}$, finite.
- $k$.3 **$Z_1$** (Theorem 3.3, Appendix A), finite; as in the code, the last
  term is $L|Df(\bar w)|/(2\nu^{n_{\mathscr C}+2}((n_{\mathscr C}+2)^2-1))$, half the
  expression printed in the paper (the factor 2 in $|\cdot|_{\mathscr C_\nu}$ allows it).
- $k$.4 **$Z_2$** (Theorem 3.3, Appendix A), finite.
- $k$.5 **Radii polynomial** (Theorem 3.3): with $Q = |A(\Gamma, \Gamma\Lambda, 0, 0)|_{\mathscr B(\mathbb C^q, X)}$,
  `interval_of_existence(Y + Q𝓛_y|η̄|, Z₁ + Q𝓛_y, Z₂, ϱ)` is nonempty, i.e. some
  $\bar\rho \in [0, \varrho]$ satisfies both inequalities of Theorem 3.3
  (RadiiPolynomial convention $Y + (Z_1-1)\rho + Z_2\rho^2/2 \le 0$, $Z_1 + Z_2\rho < 1$).
  Conclusion: a unique zero $\tilde\chi$ of $F(\cdot;\tilde\alpha)$ in
  $\mathrm{cl}(B_{\bar\rho}(\bar x))$, real by Lemma 3.4.
- $k$.6 **$C^0$ error** (Corollary 3.5): with $\bar\rho$ the left end of the
  interval of existence, $\max\{\bar\rho, |\Gamma|_\infty(\bar\rho + (1+\mathcal L_y)\mu)\}$
  (as in the code; it bounds the right-hand side of the estimate of Corollary 3.5) is at
  most the distance stated in the theorem ($2.9\times10^{-7}$, $5.5\times10^{-6}$,
  $2.4\times10^{-5}$, $9.8\times10^{-7}$).

All conditions hold $\Rightarrow$ Theorem (cases 1 and 2 give the two solutions
of Theorem 1.1, distinct because their $\bar\phi = u(0)$ differ by far more
than the error bounds: 2.69 vs 20.25).

**Notes for the reader (and questions to the authors).**
- In the original `ψ̂` for SH and FHN, the argument of $D\mathbf N$ is bounded
  by $(1+\mathcal L_y)\mu$ instead of $|\Gamma|(1+\mathcal L_y)\mu$, and for FHN
  the term $2c_1$ enters with its sign ($c_1 < 0$), making that $\hat\psi$
  negative ($\approx -1.5\times10^{-5}$). Condition $k$.1d recomputes the three
  inequalities with a bound that includes these factors; it passes in all cases.
- In the original `opnorm_Df` for FHN, $|\epsilon^{-2}(1-3w_5^2)|$ appears where
  $Df$ has $\epsilon^{-2}(1-3w_2^2)$; the value used (241.4) exceeds the correct
  one (48.8), so the resulting $Z_1$ is still an upper bound.
