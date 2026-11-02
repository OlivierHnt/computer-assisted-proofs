# Statement and verification conditions

**Theorem (global family of stable periodic orbits; Theorem 1.1 of the paper).**
Consider the system of two predators $X_1, X_2$ competing for one prey $S$,

$$
\dot X_j = \Big(\frac{m_j S}{S + a_j} - d_j\Big) X_j \quad (j = 1, 2), \qquad
\dot S = \Big(\gamma\big(1 - \tfrac{S}{\kappa}\big) - \frac{m_1}{y_1}\frac{X_1}{S + a_1} - \frac{m_2}{y_2}\frac{X_2}{S + a_2}\Big) S,
$$

with $a_1 = 10$, $a_2 = 41$, $d_1 = 0.8$, $d_2 = 0.5$, $m_1 = m_2 = y_1 = y_2 = \gamma = 1$.
There is a global family of positive periodic orbits $\mathcal P(\kappa)$,
analytic in the carrying capacity $\kappa$, and there are
$\hat\kappa_1, \hat\kappa_2 \in [92, 129]$ such that

- (i) for each $\kappa \in (\hat\kappa_1, \hat\kappa_2)$, $\mathcal P(\kappa)$ is positive and stable (its two non-trivial Floquet exponents have negative real part);
- (ii) $\mathcal P(\hat\kappa_1) = \mathcal C_1$, the boundary limit cycle in $Q_1 = \{(X_1, 0, S) : X_1, S > 0\}$;
- (iii) $\mathcal P(\hat\kappa_2) = \mathcal C_2$, the boundary limit cycle in $Q_2 = \{(0, X_2, S) : X_2, S > 0\}$.

The family lies within the radius $r$ of condition 1d (in the norm of $\mathcal X_\nu$, which bounds the $C^0$ distance;
the paper states $r \le 10^{-10}$) of the Fourier–Chebyshev approximation $\bar\chi$ built from the inputs.
Condition 3 gives $\hat\kappa_1 \in (\kappa(h_1^-), \kappa(h_1^+)) \approx (92.66, 94.72)$
and $\hat\kappa_2 \in (\kappa(h_2^-), \kappa(h_2^+)) \approx (124.01, 127.72)$.
This answers the stable connection problem of Butler and Waltman (1981).

Reference: K. E. M. Church, J.-Y. Dai, O. Hénot, P. Lappicy, N. Vassena,
*Global continuation of stable periodic orbits in systems of competing predators*,
arXiv:2504.03058. Theorem, lemma and section numbers below refer to it.

## Setting

- Rescaled system and blow-up (Section 1.2, Section 2): $x_j = \tfrac{1}{\kappa}\tfrac{m_j}{\gamma y_j} X_j$,
  $s = S/\kappa$, time rescaled by $\tau$ so that the period is $2\pi$,
  $\alpha_j = a_j/\kappa$, $\lambda_j = \tfrac{1}{\kappa}\tfrac{a_j d_j}{m_j - d_j}$, $\delta_j = (m_j - d_j)/\gamma$;
  blow-up $x_j = \zeta_j u_j$ ($j = 1,2$), $s = u_3$, normalised by $u_j(0) = 1$.
  The parameter range is $\kappa \in [\kappa_1, \kappa_2] = [92, 129]$, parametrised by
  $\eta \in [-1, 1]$ through $\kappa(\eta) = 2\kappa_1\kappa_2 / (\kappa_1 + \kappa_2 + (\kappa_1 - \kappa_2)\eta)$, so that
  $1/\kappa$, hence $\alpha_j, \lambda_j$, are affine in $\eta$ (a degree-1 Chebyshev series, `vκ⁻¹` in the code).
- Unknown $\chi = (\tau, \zeta_1, \zeta_2, u_1, u_2, u_3)$ with $\tau, \zeta_j$ functions of $\eta$ and $u_j$ functions of $(t, \eta)$.
- Spaces and norms (Section 2.1): $\mathcal P_\nu$ (Chebyshev series in $\eta$, $\|\psi\| = \sum_n |\psi_{|n|}|\nu^{|n|}$),
  $\mathcal W_\nu$ (Fourier in $t$ $\otimes$ Chebyshev in $\eta$, $\ell^1$ in $k$ with weight $\nu^{|n|}$ in $n$;
  `Ell1(GeometricWeight(ν), IdentityWeight())`), $\mathcal X_\nu = \mathcal P_\nu^3 \times \mathcal W_\nu^3$ with the sum norm; $\nu = 1.1$.
  For the Floquet normal form (Section 3, Appendix A): $\mathcal Y_\nu = \mathrm{Mat}_{3\times3}(\mathcal P_\nu) \times \mathrm{Mat}_{3\times3}(\mathcal W_\nu)$, sum norm.
- Truncation: Fourier order $K = 20$ (orbits) and $K_{\text{lin}} = 20$ (Floquet), Chebyshev order $N = 30$.
- Maps: $F(\chi) = (\rho(u), \partial_t u - \tau f(\zeta_1, \zeta_2, u))$ of Section 2.1, with the phase condition
  $\int_0^{2\pi}\langle u, \partial_t\Gamma\rangle\,dt$ against $\partial_t\Gamma = \partial_t \bar u$ and $\ell(u_j) = u_j(0, \cdot)$
  (code: `F_po!`, `DF_po!`, and `construct_W₀_po`, `construct_W₁_po` on the Chebyshev side);
  $G(C, V) = (V(0) - I, \partial_t V + VC - \tilde\tau\,\partial_u f(\tilde\zeta_1, \tilde\zeta_2, \tilde u) V)$ of Section 3, evaluated at the
  validated orbit $\tilde\chi$ of condition 1 (code: `F_lin!`, `DF_lin!`, `construct_W₁_lin`).
- Theorem 2.1 (Newton–Kantorovich): radius threshold $R = 10^{-6}$ for $F$ and $R = \infty$ for $G$ ($G$ is quadratic, $Z_2 = 2\|B\|$).
- Points of Section 2.3 and Section 3: $h_1^- = -0.95$, $h_1^+ = -0.8$, $h_2^- = 0.8$, $h_2^+ = 0.95$; boxes
  $\Omega_1 = [-0.15, -0.05] \times [-0.1, 0.1]\,i$, $\Omega_2 = [-0.1, -0.0025] \times [-0.1, 0.1]\,i$, $\Omega_3 = [-0.04, -0.02] \times [-0.1, 0.1]\,i$.
- Arithmetic: outward-rounded binary64 intervals (RadiiPolynomial.jl on IntervalArithmetic.jl, versions pinned in `checker/Manifest.toml`); the inverses
  $(\bar u_3 + \alpha_j)^{-1}$ of Lemma 2.2 and their error bounds in 256-bit intervals,
  converted outward to binary64. The approximate inverses $A$, $B$ and the Chebyshev
  interpolation of the node values are plain binary64, then used as thin intervals.

## Inputs (in `data/`)

Both inputs are frozen. They were produced by `julia --project=search search/produce.jl`
(the `produce` command of `entry.json`: ODE integration, Newton's method and a continuation in $\kappa$ over 33 Chebyshev nodes),
which is not part of the proof: any other data passing the conditions below would prove the theorem as well.

- $\bar\chi_\ell$, $\ell = 0, \dots, 32$: floating-point approximate zeros of $F$ at the nodes
  $\eta_\ell = -\cos(\ell\pi/32)$ ($\kappa$ from $92$ at $\ell = 0$ to $129$ at $\ell = 32$), truncated to Fourier order $K = 20$.
- $\bar\upsilon_\ell = (\bar C_\ell, \bar V_\ell)$: floating-point approximate Floquet normal forms (zeros of $G$ with $\tilde\chi$ replaced by $\bar\chi_\ell$) at the same nodes, Fourier order 20.

Files (raw little-endian IEEE binary64, no header; every coefficient is a complex
number stored as (re, im)); Fourier coefficients use the convention
$u(t) = \sum_{k=-K}^{K} u_k e^{ikt}$ and are ordered $k = -K, \dots, K$:

- `x_po_grid.bin` ($33 \times 126 \times 16$ bytes): node by node ($\ell = 0$ first), the 126 coefficients
  $\tau, \zeta_1, \zeta_2$, then $(u_1)_{-20..20}$, $(u_2)_{-20..20}$, $(u_3)_{-20..20}$.
- `x_lin_grid.bin` ($33 \times 378 \times 16$ bytes): node by node, the 9 entries of $C$ in column-major order
  ($C_{ij}$ at position $i + 3(j-1)$), then the 9 Fourier series $(V_{ij})_{-20..20}$ in the same order.

SHA-256 digests are in `data/hashes.json`; the checker verifies them before reading the files.

Constructions the checker performs on the inputs (in floating point; the conditions do not depend on their exact values, only on the bounds evaluated on their results):

- polishing (`polish!`: zeroing Fourier coefficients above the geometric-decay envelope) of each node value;
- $\bar\chi$, $\bar\upsilon$: Chebyshev coefficients of order $N = 30$ of the node values (discrete Chebyshev transform of size 64 on the 33 nodes, `grid2cheb`),
  made exactly conjugation-symmetric ($\Sigma(\bar\chi) = \bar\chi$, `conjugacy_symmetry!`, so that Lemma 2.4 gives a real-valued solution) and polished; $\partial_t\Gamma = \partial_t \bar u$;
- $A_{\text{finite}}$ (Section 2.2) and $B_{\text{finite}}$ (Appendix A): Chebyshev interpolants of the floating-point inverses of
  $DF$ and $DG$ at the nodes (LAPACK `inv`, so the exact values may differ across platforms); $A_{\text{tail}}$, $B_{\text{tail}}$ as in the paper.

## Conditions to verify

Each is a finite computation in outward-rounded arithmetic on the data above.
The checker prints them in this order, each prefixed `pass` or `fail`, with the
enclosure obtained.

1. **Family of periodic orbits: Theorem 2.1 for $F$ on $\mathcal X_\nu$, $R = 10^{-6}$.**
   - 1a. **$Y$ (Section 2.2.1, Lemma 2.2).** $Y = \|A\|\,(\|W_0\| + \|\bar\tau f(\bar\zeta_1, \bar\zeta_2, \bar u) - \omega_0\|)$, with $\|A\| = \max(\|A_{\text{finite}}\|, 1/(K+1))$, is a finite bound.
   - 1b. **$Z_1$ (Section 2.2.2).** $Z_1 = \|A W_1 - I\| + \max_j\sum_i\|\omega_1^{(i,j)}\|/(K+1) + \|A\|\,\|DF(\bar\chi) - W_1\|$, the first term computed on the 65 Chebyshev nodes $\cos(\pi j/64)$ (where $A$ and $W_1$ are evaluated, multiplied, and interpolated back to Chebyshev order $2N$), is a finite bound.
   - 1c. **$Z_2$ (Section 2.2.3, Lemma 2.3).** $Z_2 = \|A\|$ times the bound of Section 2.2.3 on $\sup_{B_R(\bar\chi)}\|D^2F\|$, built from the constants $\mathfrak C_j$ of Lemma 2.3, is a finite bound, and the denominators of Lemma 2.3 are positive ($\mathfrak C_j > 0$).
   - 1d. **Contraction.** The interval of existence of Theorem 2.1 determined by $Y, Z_1, Z_2$ and $R$
     (`RadiiPolynomial.interval_of_existence`: $2YZ_2 \le (1 - Z_1)^2$, $Z_1 + Z_2 r < 1$) is non-empty and contained in $[0, R]$,
     and its left end $r$ satisfies $r \le 10^{-10}$, the distance stated in the theorem;
     $r$ is the proved radius: a unique zero $\tilde\chi$ of $F$ in $B_r(\bar\chi)$, real-valued by Lemma 2.4.
2. **Floquet normal form: Theorem 2.1 for $G$ on $\mathcal Y_\nu$, $R = \infty$ (Appendix A).**
   $G$ is evaluated at the validated orbit (the enclosure $\bar\chi + B_r$ of condition 1).
   - 2a. $Y = \|B\|\,\|G(\bar\upsilon)\|$ is a finite bound.
   - 2b. $Z_1 = \|B\hat W_1 - I\| + \max_j\sum_i\|\omega_1^{(i,j)}\|/(K+1) + \|B\|\,\|\tilde\tau\partial_u f - \omega_1\|$ is a finite bound.
   - 2c. $Z_2 = 2\|B\|$ is a finite bound.
   - 2d. The interval of existence is non-empty; its left end is the proved radius of $(\tilde C, \tilde V)$, a zero of $G$.
3. **Positivity and boundary crossings (Section 2.3, Lemmas 2.5 and 2.6).** With $\tilde\tau, \tilde\zeta_j$ enclosed by
   $\bar\tau, \bar\zeta_j \pm r$ (Lemma 2.5), on 100 subintervals of each interval (`prove_crossing`, `prove_interior`):
   - 3a. on $[h_1^-, h_1^+]$: $\tilde\zeta_2 < 0$ at $h_1^-$ or at $h_1^+$; $\tilde\tau > 0$; $\tilde\zeta_1 > 0$; and $\tilde\zeta_2' \ne 0$, i.e. the enclosure of Lemma 2.6 of $\sum_n n (\tilde\zeta_2)_n \sin(n \arccos\eta)$ excludes 0;
   - 3b. on $[h_1^+, h_2^-]$: $\tilde\tau > 0$, $\tilde\zeta_1 > 0$, $\tilde\zeta_2 > 0$;
   - 3c. on $[h_2^-, h_2^+]$: $\tilde\zeta_1 < 0$ at $h_2^-$ or at $h_2^+$; $\tilde\tau > 0$; $\tilde\zeta_2 > 0$; $\tilde\zeta_1' \ne 0$.

   Hence $\tilde\zeta_2$ (resp. $\tilde\zeta_1$) has exactly one zero $\hat\eta_1 \in (h_1^-, h_1^+)$ (resp. $\hat\eta_2 \in (h_2^-, h_2^+)$),
   where the family meets $Q_1$ (resp. $Q_2$), and the orbits are positive for $\eta \in (\hat\eta_1, \hat\eta_2)$; $\hat\kappa_j = \kappa(\hat\eta_j)$.
4. **Floquet exponents (Section 3).** For each of the three intervals $[a, b] = [h_1^-, h_1^+]$, $[h_1^+, h_2^-]$, $[h_2^-, h_2^+]$ with its box
   $\Omega_1, \Omega_2, \Omega_3$ and expected number $1, 2, 1$ of exponents (`prove_inside_box`):
   - 4x.1 the Gershgorin enclosures (column sums after an approximate diagonalisation) of the spectrum of $\tilde C(\eta_0)$ (enclosed by $\bar C(\eta_0) \pm r$) contain exactly
     the expected number of enclosures inside $\Omega$, and these are disjoint from the other enclosures (otherwise Gershgorin's theorem only bounds the total number of eigenvalues, not the number inside $\Omega$). The count is taken at $\eta_0 = 0$ for the middle interval (Section 3, step 2) and at the left end point $\eta_0 = a$ for the two outer intervals; any point of $[a, b]$ is equivalent since 4x.2 shows that the count is constant on $[a, b]$;
   - 4x.2 $(\tilde C(\eta) - zI)^{-1}$ is bounded for all $z \in \partial\Omega$ and all $\eta \in [a, b]$ (adaptive bisection of $\partial\Omega \times [a, b]$ down to width $10^{-6}$),
     so no eigenvalue of $\tilde C(\eta)$ crosses $\partial\Omega$ and the number of Floquet exponents in $\Omega$ stays constant on $[a, b]$.

   With 4b, both non-trivial exponents lie in $\Omega_2 \subset \{\Re < 0\}$ on $[h_1^+, h_2^-]$; with 4a, 4c, one non-trivial exponent stays in
   $\Omega_1$, $\Omega_3 \subset \{\Re < 0\}$ near the boundary planes, and the other one can change sign only through 0, which by the invertibility of $DF(\tilde\chi)$
   happens only on $Q_1$, $Q_2$ (Section 3, step 3). The exponents of the original system are those of $G$ scaled by $\tilde\tau/\gamma > 0$.

All conditions hold $\Rightarrow$ Theorem.
