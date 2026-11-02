# Statement and verification conditions

**Theorem (validated unstable manifolds of three DDEs).** For each of the three
systems below, let $c$ be the equilibrium and $\lambda_1,\dots,\lambda_m$ its
unstable eigenvalues, with eigenvectors $\xi_j$ scaled as stated. The
parameterization of the local unstable manifold of $c$ (Lemmas 1.1, 1.2),

$$P(\sigma,t) = \sum_{\alpha \in \mathbb N^m} A_\alpha\, e^{\langle \alpha,\Lambda\rangle t}\,\sigma^\alpha,
\qquad A_0 = c,\ A_{e_j} = \xi_j,$$

whose coefficients of order $|\alpha| \ge 2$ solve the homological equations
(Theorem 1.3, eq. (10)), converges on the closed unit polydisc; the unstable
eigenvalue satisfies $|\lambda - \bar\lambda| \le \varepsilon$ (and
$\lambda_2 = \bar\lambda_1$ in case I), the coefficients of order
$|\alpha| \le N$ lie in the intervals computed by the checker, and the tail
satisfies $\|\Pi_\infty(A)\|_X = \max_k \sum_{|\alpha| > N} |a^k_\alpha| \le r_0$:

| case | system | $m$ | $N$ | scaling $s$ | $\bar\lambda$ | $\varepsilon$ | $r_0$ |
|---|---|---|---|---|---|---|---|
| I | delayed van der Pol (Section 4.1) | 2 | 99 | 2 | `0.5218379440444361 + 0.8056375092559338i` | `1.6393305131160153e-16` | `1.0350352024895318e-35` |
| II | SEIR (Section 4.2) | 1 | 200 | 1.8 | `0.17847035388447238` | `1.2081032662362984e-16` | `8.073867668146611e-11` |
| III | Mackey–Glass, $\rho = 3$ (Section 4.3) | 1 | 600 | 0.7 | `0.23930016974956494` | `1.2766318821034955e-16` | `1.4805974650552742e-54` |

One hypothesis is taken from the paper and not checked: that the eigenvalue
enclosed in condition 0c (and its conjugate in case I) is the only unstable
eigenvalue of $c$, i.e. the Morse index is $m$ (I: ref. [33]; II: Section 4.2;
III: Section 4.3). Condition 0c proves that there is exactly one eigenvalue in
the disc $|\lambda - \bar\lambda| \le \varepsilon$, not that there are no others
in the right half-plane.

Reference: O. Hénot, J.-P. Lessard, J. D. Mireles James, *Parameterization of
unstable manifolds for DDEs: formal series solutions and validated error
bounds*, J. Dyn. Differ. Equ. 34 (2022), 1285–1324, doi:10.1007/s10884-021-10002-8.
Section, lemma and equation numbers refer to the authors' revised version
(https://www.math.mcgill.ca/jplessard/papers/manifoldValidations_revised.pdf),
which may differ from the journal's.

## Setting

- Spaces and norms: $X = (\ell^1_m)^d$, $d = 2$, with
  $\|A\|_X = \max_{k} \|a^k\|_1$, $\|a\|_1 = \sum_\alpha |a_\alpha|$ (eq. (28),
  with the $\infty$-norm on $\mathbb C^d$, as in the matrix norms of Appendix
  A.1); $\Pi_N$, $\Pi_\infty$ the projections on $|\alpha| \le N$ and
  $|\alpha| > N$ (Section 3.2); $*$ the Cauchy product (Section 2.1). Case I has
  complex coefficients ($m = 2$, $\Lambda = (\lambda, \bar\lambda)$,
  $A_{(\alpha_2,\alpha_1)} = \overline{A_{(\alpha_1,\alpha_2)}}$); II and III
  real ones ($m = 1$). Delayed coefficients
  $B_\alpha = A_\alpha e^{-\langle\alpha,\Lambda\rangle\tau}$ (eq. (19));
  $\|\tau\| = 1$ on $\ell^1_m$.
- Characteristic matrix $\Psi(z) = K_1 + e^{-z\tau}K_2 - zI$, $K_1 = D_1F(c,c)$,
  $K_2 = D_2F(c,c)$ (eq. (2)), and characteristic function
  $f(z) = \det\Psi(z)$; the eigenvalues are the zeros of $f$.
- The systems (polynomial form, $x = x(t)$, $y = x(t-\tau)$):
  - I (eq. (22)): $u' = v$, $v' = -a u - b u(t-\tau) + \mu(1 - u(t-\tau)^2) v$,
    $a = 0.9$, $b = 0.1$, $\mu = 1$, $\tau = 2$; $c = (0,0)$,
    $\xi = s\,(1, \lambda)$ with $s = 2$ (Section 4.1 does not print $s$);
    characteristic matrix (23), $f(z) = z^2 - \mu z + a + b e^{-z\tau}$,
    $f' = 2z - \mu - b\tau e^{-z\tau}$, $f'' = 2 + b\tau^2 e^{-z\tau}$;
    homological equations (25).
  - II (Section 4.2): $S' = \mu - \nu S I - \mu S$, $I' = \nu S(t-\tau) I(t-\tau) e^{-\mu\tau} - (\gamma+\mu) I$,
    $\mu = 0.1$, $\nu = 0.5$, $\gamma = 0.1$, $\tau = 1$; $c = (1, 0)$,
    $\xi = s\,(\nu, -(\lambda+\mu))$;
    $f(z) = (z+\mu)\,g(z)$, $g(z) = z + \gamma + \mu - \kappa e^{-z\tau}$, $\kappa = \nu e^{-\mu\tau}$.
  - III (eq. (12), $v = 1/(1+u^\rho)$): $u' = -\gamma u + \beta u(t-\tau) v(t-\tau)$,
    $v' = -\rho u^{\rho-1} v^2 (-\gamma u + \beta u(t-\tau) v(t-\tau))$,
    $\gamma = 1$, $\beta = 2$, $\tau = 2$, $\rho = 3$; $c = (0, 1)$,
    $\xi = s\,(1, -\rho (u^*)^{\rho-1}) = s\,(1, 0)$ (Section 4.3);
    $f(z) = z\,g(z)$, $g(z) = z + \gamma - \kappa e^{-z\tau}$, $\kappa = \beta$ (eq. (33) times $z$).

  In II and III, $f = q\,g$ with $q = z+\mu$ (II), $q = z$ (III), so
  $f' = g + q g'$, $f'' = 2g' + q g''$, $g' = 1 + \tau\kappa e^{-z\tau}$,
  $g'' = -\tau^2 \kappa e^{-z\tau}$.
- Constants: every decimal parameter is the tightest binary64 interval
  containing it; $r_* = 10^{-3}$ and $r_\lambda = 10^{-10}$ are the binary64
  numbers nearest to these decimals (within one ulp), used exactly. Truncation
  $N$ as in the table.
- The map: $T$ of Section 3.3 on $X^\infty$,
  $T(H) = L\,\Pi_\infty \hat F(\Pi_N(A) + H, \tau(\Pi_N(A) + H))$, with
  $\|L\| \le C_N$, and the bounds $Y$, $Z_1$, $Z_2$ of Lemma 3.3; the radii
  polynomial of Theorem 3.1 is $p(r) = \tfrac{Z_2}{2} r^2 + (Z_1 - 1) r + Y$.
- Arithmetic: outward-rounded binary64 intervals, complex intervals as rectangles
  (RadiiPolynomial.jl on IntervalArithmetic.jl, versions pinned in `checker/Manifest.toml`). Taylor coefficients are sequences on
  `Taylor(N)` (II, III) and `Taylor(N) ⊗ Taylor(N)` (I); products are Cauchy
  products on the full product space, norms the unweighted $\ell^1$ norm.

## Inputs (in `data/`)

- Frozen, one file per case: the approximate unstable eigenvalue $\bar\lambda$,
  a binary64 zero of $f$ found by a floating-point Newton iteration
  (`julia --project=search search/produce.jl`, the `produce` command of
  `entry.json`; not part of the proof: the conditions prove an enclosure around
  the stored value, whatever it is). Byte layout: raw little-endian IEEE
  binary64, no header. Hashes in `data/hashes.json`; the checker verifies them
  before reading the files.

  | file | bytes | content |
  |---|---|---|
  | `vdp_lambda.bin` | 16 | $\operatorname{Re}\bar\lambda_1$, $\operatorname{Im}\bar\lambda_1$ (`0.5218379440444361`, `0.8056375092559338`) |
  | `seir_lambda.bin` | 8 | $\bar\lambda$ (`0.17847035388447238`) |
  | `mg_lambda.bin` | 8 | $\bar\lambda$ (`0.23930016974956494`) |

- Produced: none.
- Constructions the checker performs on the inputs: the eigenvalue enclosure
  (conditions 0a–0c), the eigenvector $\xi$ (0d), the coefficients $\Pi_N(A)$
  and $B = \Pi_N(\tau A)$ (0e), and the bounds $C$, $Y$, $Z_1$, $Z_2$ (1–4).

## Conditions to verify

Each is a finite computation in outward-rounded interval arithmetic on the
inputs above. The checker prints them in this order, each prefixed `pass` or
`fail`, with the enclosure obtained; `<id>` is I, II or III.

0. **Eigenvalue, eigenvector, coefficients.**

   Eigenvalue (0a–0c; `eigenvalue_bounds` in `checker/src/common.jl`):
   Theorem 3.1 in dimension one for the scalar function $f$. Let $A$ be the
   binary64 (complex in I) number nearest to $1/f'(\bar\lambda)$, $A \ne 0$,
   and $T(\lambda) = \lambda - A f(\lambda)$, whose fixed points are the zeros
   of $f$. For $|\lambda - \bar\lambda| \le r \le r_\lambda$, with
   $h = \lambda - \bar\lambda$,
   $$|T(\bar\lambda) - \bar\lambda| = |A f(\bar\lambda)| =: Y,\qquad
   |T'(\lambda)| \le |1 - A f'(\bar\lambda)| + |A|\,|f'(\lambda) - f'(\bar\lambda)| \le Z_1 + Z_2 r,$$
   $$Z_1 := |1 - A f'(\bar\lambda)|,\qquad Z_2 := |A| \sup_{|\lambda - \bar\lambda| \le r_\lambda} |f''(\lambda)|,$$
   since $f'(\lambda) - f'(\bar\lambda) = \int_0^1 f''(\bar\lambda + t h)\, h\, dt$;
   hence $|T(\lambda) - \bar\lambda| \le Y + Z_1 r + \tfrac{Z_2}{2} r^2$. So
   $T$ maps the closed disc of radius $r$ into itself when $p(r) \le 0$ and is
   a contraction there when $Z_1 + Z_2 r < 1$: $f$ has exactly one zero in it.
   $Y$, $Z_1$ are interval evaluations at the thin $\bar\lambda$; the sup of
   $|f''|$ is bounded by the interval evaluation of $f''$ on the square (I) or
   segment (II, III) $\bar\lambda + [-r_\lambda, r_\lambda](+ i[-r_\lambda, r_\lambda])$,
   which contains the disc. In case I the second eigenvalue is $\bar\lambda_1$
   (real coefficients), enclosed by the conjugate square.

   - `<id>.0a` $Y = |A f(\bar\lambda)|$ is a finite bound.
   - `<id>.0b` $Z_1 = |1 - A f'(\bar\lambda)|$ and $Z_2 = |A| \sup_{|\lambda-\bar\lambda|\le r_\lambda}|f''|$,
     $r_\lambda = 10^{-10}$, are finite bounds.
   - `<id>.0c` `interval_of_existence(Y, Z₁, Z₂, r_λ)` succeeds and returns a
     bounded interval $[\varepsilon, \cdot] \subset [0, r_\lambda]$ on which
     $p(r) \le 0$ and $Z_1 + Z_2 r < 1$ (re-evaluated at $\varepsilon$ as in
     condition 5), and $\operatorname{Re}\lambda > 0$ on the enclosure.
     Conclusion: $f$ has exactly one zero $\lambda$ in $|\lambda - \bar\lambda| \le \varepsilon$,
     an unstable eigenvalue; the printed enclosure is the square (I) or
     segment (II, III) of half-side $\varepsilon$, which contains the disc and
     is used in all later conditions.
   - `<id>.0d` Eigenvector (Lemma 1.2): $\Psi(\lambda)\xi \ni 0$ for the
     enclosed $\lambda$ and $\xi$ as in *Setting*.
   - `<id>.0e` Non-resonance and coefficients: $f(\langle\alpha,\Lambda\rangle) \ne 0$
     for every $2 \le |\alpha| \le N$ (the interval excludes zero; its smallest
     modulus is printed), so the homological equations (10) are uniquely
     solvable up to order $N$. The checker computes $\Pi_N(A)$
     (`solve_homological_equations`): $A_0 = c$, $A_{e_j} = \xi_j$, and for
     $2 \le |\alpha| \le N$, order by order,
     $A_\alpha = -\Psi(\langle\alpha,\Lambda\rangle)^{-1}\hat F_\alpha$ with the
     explicit inverse of the $2\times2$ matrix and $\hat F_\alpha$ the hat
     products (eq. (15)) of the enclosures of order $< |\alpha|$ (each factor's
     coefficient of order $\alpha$ set to zero). By induction on $|\alpha|$ and
     inclusion isotonicity, the intervals contain the exact $\Pi_N(A)$ for
     every $\lambda$ in the enclosure of 0c. In case I the equation at
     $(\alpha_2,\alpha_1)$ is the complex conjugate of the one at
     $(\alpha_1,\alpha_2)$, so $A_{(\alpha_2,\alpha_1)} = \overline{A_{(\alpha_1,\alpha_2)}}$
     exactly: the checker computes $\alpha_1 \le \alpha_2$ and sets the mirror
     entry to the conjugate enclosure. $B_\alpha = A_\alpha e^{-\langle\alpha,\Lambda\rangle\tau}$
     is computed alongside. The size of the largest coefficient of order
     $N-2 \le |\alpha| \le N$ is printed (information).

1. **`<id>.1` Resolvent bound (Proposition A.1).** With
   $\zeta = (N+1) \min_j \operatorname{Re}\lambda_j$: for $|\alpha| \ge N+1$,
   $|\langle\alpha,\Lambda\rangle| \ge \operatorname{Re}\langle\alpha,\Lambda\rangle \ge \zeta$,
   and the denominator of (36) increases with $|z|$ and $\operatorname{Re} z$, so
   $C_N = \sup_{|\alpha| \ge N+1} \|\Psi(\langle\alpha,\Lambda\rangle)^{-1}\| \le C := 1/(\zeta - \|K_1\| - e^{-\tau\zeta}\|K_2\|)$
   whenever the denominator is positive (then (34), (35) hold for all these $z$).
   I: $\|K_1\| = \max(1, |a|+|\mu|)$, $\|K_2\| = |b|$; II: $\max(\mu+\nu, \gamma+\mu)$,
   $\nu e^{-\mu\tau}$; III: $\gamma$, $\beta$. The condition: the denominator
   is positive ($0 < C < \infty$), so $\|L\| \le C_N \le C$.
2. **`<id>.2` $Y$ (Lemma 3.3, eq. (30))** is a finite bound, where
   $Y = C \max_j \sum_{|\alpha| > N} |\hat F_j(\Pi_N A, \Pi_N \tau A)_\alpha|$;
   for $|\alpha| > N$ the hat product is the full Cauchy product, computed on
   the full product space (orders up to $3N$, $2N$, $6N$), so the sum runs over
   all $N+1 \le |\alpha| \le \deg(F_j)\,N$.
   I: $\hat F_2 = -\mu\, b^1 * b^1 * a^2$ (eq. (24)); II: $\hat F = (-\nu\, a^1 * a^2,\ \nu e^{-\mu\tau} b^1 * b^2)$;
   III: $\hat F = (\beta\, b^1 * b^2,\ \rho(\gamma (a^1)^{3} * (a^2)^2 - \beta (a^1)^{2} * (a^2)^2 * b^1 * b^2))$.
3. **`<id>.3` $Z_1$ (Lemma 3.3, eq. (31))** is a finite bound, where
   $Z_1 = C \max_j \sum_i \big(\|\partial_{x_i}\hat F_j\|_1 + \|\partial_{y_i}\hat F_j\|_1\big)$,
   the $\ell^1$ norms of the coefficient sequences of the partial derivatives of
   the nonlinear part of $F$ at $(\Pi_N A, \Pi_N \tau A)$ (full products; by
   (27), (29) and $\|\tau\| = 1$ this bounds $\|DT(0)\|$). For III the two
   monomials of $\partial_{x_1}F_2$ and of $\partial_{x_2}F_2$ are summed before
   taking the norm.
4. **`<id>.4` $Z_2$ (Lemma 3.3, eq. (32))** is a finite bound, for
   $r_* = 10^{-3}$, where, for $H \in B_{r_*}(0) \subset X^\infty$,
   $D^2\hat F_j(U)[(V,\tau V),(W,\tau W)] = \sum_{p,q} \partial_p\partial_q F_j(U) * (V,\tau V)_p * (W,\tau W)_q$
   over all ordered pairs $(p,q)$ of the $2d$ variables
   $(x_1,\dots,x_d,y_1,\dots,y_d)$, hence
   $Z_2 = C \max_j \sum_{p,q} \partial_p\partial_q |F_j|(z)$ with
   $z = (\|a^i\|_1 + r_*,\ \|b^i\|_1 + r_*)_i$ (`hessian_sum` in
   `checker/src/common.jl`); each mixed pair $\{x_i, y_j\}$ appears twice (see
   *Notes on the paper*).
   I: $\mu(2 x_2 + 4 y_1)$; II: $\max(2\nu, 2\nu e^{-\mu\tau}) = 2\nu$;
   III: $\max(2\beta,\ \text{Hessian sum of } \rho\gamma x_1^3 x_2^2 + \rho\beta x_1^2 x_2^2 y_1 y_2)$.
   This $Z_2$ bounds the second derivative of $T$ on the whole ball of radius
   $r_*$, as Theorem 3.1 requires, so it enters $p$ with the factor
   $\tfrac12$ and is passed to `interval_of_existence` as is (the library's
   polynomial is $Y + (Z_1 - 1) r + \tfrac{Z_2}{2} r^2$).
5. **`<id>.5` Radii polynomial (Theorem 3.1).** `interval_of_existence(Y, Z₁, Z₂, r*)`
   succeeds and returns a bounded interval $[r_0, r_1] \subset [0, r_*]$ on
   which $p(r) \le 0$ and $Z_1 + Z_2 r < 1$; in addition the checker
   re-evaluates $p(r_0) \le 0$ and $Z_1 + Z_2 r_0 < 1$ at the upper endpoints
   of the enclosures of $Y$, $Z_1$, $Z_2$ (both are increasing in $Y$, $Z_1$,
   $Z_2$ for $r \ge 0$), and $r_0 \le r_*$, exactly in rational arithmetic
   (these are binary64 numbers; an interval evaluation at $r_0$, which is a
   root of $p$ up to a few ulps, would be swamped by its own rounding).
   Conclusion: $T$ has a unique fixed point $\tilde H$ in the ball of radius
   $r_0$ of $X^\infty$, $\Pi_N(A) + \tilde H$ is the coefficient sequence of
   $P$ (Proposition 3.2, with $\Pi_N(A)$ the exact coefficients of 0e), the
   series converges on the closed unit polydisc, and
   $\|\Pi_\infty(A)\|_X \le r_0$.

All conditions hold $\Rightarrow$ Theorem (under the hypothesis taken from the
paper, stated after the table).

Notes on the paper:

- Eq. (32), and the last display of the proof of Lemma 3.3, print each mixed
  second derivative $\partial_{x_i}\partial_{y_j}|F_j|$ once; condition 4 sums
  over all ordered pairs of the $2d$ variables, so each appears twice, because
  the second derivative of $\hat F_j$ in the direction $(V,\tau V)$ is the full
  Hessian in the $2d$ variables.
- Section 4.1 does not print the scaling $s$ of the eigenvector of case I;
  condition 0d and all later conditions use $s = 2$, which gives the bounds
  $Y$, $Z_1$, $Z_2$ printed there.
