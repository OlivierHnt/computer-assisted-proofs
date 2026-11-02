# Statement and verification conditions

**Theorem (period doubling in the Kot–Schaffer model).** The map $f:\mathbb R^{10}\times\mathbb R\to\mathbb R^{10}$
of (6.2), the 10-mode Galerkin projection of the Kot–Schaffer growth-dispersal
model with $\kappa(x) = 0.8 - 0.4\cos x$ and $\beta(x) = 3/(5 - 4\cos x)$,
undergoes a supercritical period-doubling bifurcation at a point
$(a_\star, \lambda_\star)$ such that

$$\lambda_\star \in [2.70423246926594,\ 2.70423246926599],\qquad
c \in [-18.6289496098,\ -18.6289496097],\qquad
\dot\nu_c \in [-0.9888090931,\ -0.9888090930].$$

The bifurcating 2-cycle exists for $\lambda > \lambda_\star$ sufficiently close
to $\lambda_\star$ and is attracting. The corresponding zero of $F_{\mathrm{pd}}$
lies within $\infty$-distance $2.1\cdot10^{-14}$ of the numerical approximation
$\bar X$ stored in `data/X_bar.bin`.

Reference: O. Hénot, J.-P. Lessard, *Rigorous verification of period-doubling
bifurcations in maps*, arXiv:2610.22986 (2026), Theorem 6.1. Section, theorem,
lemma and equation numbers below refer to the arXiv version of that paper,
which defines every object named here.

## Setting

- Space and norm: $\mathbb R^{42} = \mathbb R^n\times\mathbb R^n\times\mathbb R\times\mathbb R^n\times\mathbb R\times\mathbb R^n$,
  $n = 10$, unknown $X = (x, v, \lambda, w, \mu, h)$, with the $\infty$-norm and the
  induced operator norm (Section 5). $\bar B_r(\bar X) = [\bar X - r, \bar X + r]$ componentwise.
- The map: $x = (a_0,\dots,a_9)$ are the cosine coefficients of an even density,
  $a(y) = a_0 + 2\sum_{k\ge1} a_k\cos(ky)$, stored as a `Fourier(9, 1.0)` sequence
  with even symmetry in RadiiPolynomial (`evensym`), one coefficient per $k \ge 0$.
  $\beta_k = 2^{-k}$, $\kappa_0 = 0.8$, $\kappa_1 = -0.2$, $\kappa_k = 0$ ($k\ge 2$).
  $f(x,\lambda) = \lambda\,\beta\odot\bigl(x - M(\kappa,x)\,x\bigr)$ with $M(\kappa,x)$ the
  matrix of multiplication by $\kappa a$ projected on modes $0,\dots,9$; this is (6.2).
  $D_xf = \lambda\,\mathrm{Diag}(\beta)(I - 2M(\kappa,x))$ and
  $B(u,\cdot) = -2\lambda\,\mathrm{Diag}(\beta)M(\kappa,u)$, i.e. $B = \lambda_\star B_1$
  of Section 6; $C \equiv 0$ since $f$ is quadratic in $x$; $f$ is linear in $\lambda$,
  so $f_\lambda = f/\lambda$, $A_\lambda = A/\lambda$, $B_\lambda = B/\lambda$.
- The period-doubling map $F_{\mathrm{pd}}$ of Definition 3.1, (3.1) (code: `F`) and
  its Jacobian at any $X$ (code: `DF`, which reduces to (4.1) at a zero). Code:
  `checker/src/kot_schaffer_map.jl`.
- Constants: $r_\star = 10^{-8}$ (radius of the ball on which $Z$ is computed).
- Arithmetic: outward-rounded binary64 intervals (RadiiPolynomial.jl on IntervalArithmetic.jl, versions pinned in `checker/Manifest.toml`); interval matrix inverses by
  IntervalArithmetic's `inv`, which returns NaI (a failing condition) when it
  cannot certify invertibility.

## Inputs (in `data/`)

- Frozen: `X_bar.bin`: the numerical approximation $\bar X = (\bar x, \bar v, \bar\lambda, \bar w, \bar\mu, \bar h)$
  of a zero of $F_{\mathrm{pd}}$, obtained by Newton's method in floating point
  (the `produce` command of `entry.json`, `julia --project=search search/find_candidate.jl`;
  not part of the proof). Byte layout: 42 IEEE binary64, little-endian, no header,
  in the order
  $\bar x_0..\bar x_9$, $\bar v_0..\bar v_9$, $\bar\lambda$, $\bar w_0..\bar w_9$, $\bar\mu$, $\bar h_0..\bar h_9$
  (336 bytes). SHA-256 in `data/hashes.json`; the checker verifies it before
  reading the file. The model constants $\beta$, $\kappa$ are written in the checker
  ($\beta_k = 2^{-k}$ exactly; $\kappa_0, \kappa_1$ as the tight enclosures `I"0.8"`, `I"-0.2"`).
- Recomputed by the checker (auxiliary, floating point):
  - $Q$ = `inv(DF(X̄))` in binary64 (LAPACK), the approximate inverse of
    Theorem 5.1, then used as a thin interval matrix. Theorem 5.1 holds for any $Q$.
  - $V$ = `eigvecs(mid(A_r))`, approximate eigenvectors of the midpoint of the
    interval matrix $A_r$ below, then used as a thin interval matrix. Theorem 5.2
    holds for any invertible $V$; invertibility is certified by the interval
    inverse $\mathcal W \ni V^{-1}$.

  The conditions do not depend on the exact values of $Q$ and $V$; only the
  printed enclosures do.

## Conditions to verify

Each is a finite computation in outward-rounded interval arithmetic on the data
above. The checker prints them in this order, each prefixed `pass` or `fail`,
with the enclosure obtained.

1. **$Y$ (Theorem 5.1).** $Y \supseteq \|Q\,F_{\mathrm{pd}}(\bar X)\|_\infty$, with
   $F_{\mathrm{pd}}$ evaluated in interval arithmetic at the thin interval $\bar X$,
   is a bounded interval.
2. **$Z$ (Theorem 5.1).** $Z \supseteq \|I - Q\,DF_{\mathrm{pd}}([\bar X - r_\star, \bar X + r_\star])\|_\infty$,
   an enclosure of $\sup_{X\in \bar B_{r_\star}(\bar X)}\|I - Q\,DF_{\mathrm{pd}}(X)\|_\infty$
   (Remark 5.1), is a bounded interval.
3. **Newton–Kantorovich (Theorem 5.1; step (1) of Theorem 5.3).** $\sup Z < 1$ and
   $r_0 := \sup\bigl(Y/(1 - Z)\bigr) < r_\star$ (via `RadiiPolynomial.interval_of_existence(Y, Z, r⋆)`,
   which returns $[r_0, r_\star]$). Conclusion: a unique zero
   $X_\star = (x_\star, v_\star, \lambda_\star, w_\star, -1, h_\star)$ of $F_{\mathrm{pd}}$ in
   $\bar B_{r_0}(\bar X)$, with $DF_{\mathrm{pd}}(X_\star)$ invertible; by Lemma 3.1 and
   Theorem 4.1, (C1), (C2), (C4) and $1\notin\sigma(A)$ hold.

   Below, $(x_r, v_r, \lambda_r, w_r, \mu_r, h_r)$ are the components of
   $[\bar X - r_0, \bar X + r_0] \ni X_\star$, and $A_r = D_xf(x_r, \lambda_r) \ni A$.
4. **Gershgorin (Theorem 5.2; step (2) of Theorem 5.3).** With $\mathcal M = \mathcal W A_r V$
   and $D_k = \mathcal M_{kk} + [-1,1]\sum_{j\ne k}|\mathcal M_{kj}|$, $k_0$ the index
   whose midpoint is closest to $-1$: $V$ is real, so each disk
   $\mathcal D_k = \{z : |z - \operatorname{mid} D_k| \le \operatorname{rad} D_k\}$
   is centred on $\mathbb R$, contains the $k$-th Gershgorin disk of $V^{-1}AV$, and
   has $D_k$ as its trace on $\mathbb R$; disjointness of such disks and
   $\mathcal D_k\subset\{|z|<1\}$ are then equivalent to disjointness of the $D_k$ and
   $\sup|D_k| < 1$.
   - 4a. $-1 \in D_{k_0}$ (and $V$ is real, all $D_k$ are guaranteed enclosures).
   - 4b. $D_{k_0}$ is disjoint from every $D_k$, $k\ne k_0$.
   - 4c. $\sup|D_k| < 1$ for every $k \ne k_0$.

   Conclusion: (C3), and all eigenvalues of $A$ other than $-1$ lie in the open unit
   disk (used for the attraction of the 2-cycle in $\mathbb R^{10}$, Section 2 after
   Proposition 2.2).
5. **Normal-form coefficient (step (3) of Theorem 5.3).** The enclosure of
   $c = -w^TB(v,h)$ ((2.5) with $C\equiv0$; Proposition 2.1, Lemma 3.1) over
   $(w_r, v_r, h_r, \lambda_r)$ satisfies $\sup c < 0$: (C5), supercritical.
6. **Crossing speed (Theorem 4.1, (4.5)).** The enclosure of
   $\dot\nu_c = w^TB(\dot x, v) + w^TA_\lambda v$, with
   $\dot x = -(A - I)^{-1}f_\lambda$ from (2.3) ($(A_r - I)^{-1}$ by interval inverse),
   satisfies $\sup\dot\nu_c < 0$. With 5, Proposition 2.2 places the 2-cycle at
   $\lambda > \lambda_\star$ ($\dot\nu_c(\lambda - \lambda_\star)/c > 0$).
7. **Stated enclosures (Theorem 6.1).** Compared exactly (as rationals) with the
   decimal endpoints of the theorem:
   - 7a. $r_0 \le 2.1\cdot10^{-14}$;
   - 7b. $\lambda_r \subseteq [2.70423246926594, 2.70423246926599]$ (the enclosure
     $a_{\star,0} \in x_r[1]$ is printed alongside);
   - 7c. the enclosure of $c$ is $\subseteq [-18.6289496098, -18.6289496097]$;
   - 7d. the enclosure of $\dot\nu_c$ is $\subseteq [-0.9888090931, -0.9888090930]$.

All conditions hold $\Rightarrow$ Theorem (by Theorem 5.3 and Proposition 2.2).
