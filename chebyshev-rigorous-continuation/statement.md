# Statement and verification conditions

Reference: M. Breden, O. Hénot, *Efficient Rigorous Continuation via Chebyshev
Series Expansion I*, preprint, arXiv:2609.22689 (v1, September 2026). Section,
theorem, remark and equation numbers below are those of arXiv v1; the paper
defines every object named here. The paper states the four results below in its
"Results of the code" paragraphs (they are not numbered theorems): the cube root
(Section 2.2.3) and the square root (Section 3.1.2) are pedagogical
illustrations; the Cahn–Hilliard manifold (Section 4.1) and the SKT arcs
(Section 4.2) are the applications.

This revision of the package corrects the bounds of Sections 4.1 and 4.2 as
computed by the original code (see "Revisions and paper errata" at the end).
With the corrected bounds, the Chebyshev orders of Result 3 ($N_1 = N_2 = 16 \to 20$)
and of the arcs `u_bar2c`, `u_bar3c` of Result 4 ($20 \to 24$, $16 \to 20$) had to be
raised (new inputs, same producer) because with the paper's norm these cases fail
at the paper's orders; the arcs `u_bar2d`, `u_bar3d` were then regenerated from
the new end points of `u_bar2c`, `u_bar3c` so that the arcs still join; and the
truncation order of the approximate inverse of `u_bar3d` was raised from
$K_A = 165$ to $250$, without which that arc is not proved at any Chebyshev order.
Nothing else was tuned. All conditions pass.

All four results are instances of one argument. A continuation problem
$F(x,\lambda) = 0$ (or its pseudo-arclength / reparametrised version) is lifted
to the superposition operator $\mathcal F$ on $\ell^1_\omega(\mathcal X)$, the
space of Chebyshev series in the parameter with coefficients in $\mathcal X$
(Section 2.1, Definition 2.8). Theorem 1.1 (Newton–Kantorovich) applied to
$\mathcal F$ at an approximate zero $\bar x$ with approximate inverse
$\mathcal A = \mathcal M(A)$ gives a unique zero $x^\star$ of $\mathcal F$ with
$\|x^\star - \bar x\|_{\ell^1_\omega(\mathcal X)} \le r$; by Theorem 2.11, for every
parameter value $x^\star(\lambda)$ is then the unique zero of $F(\cdot,\lambda)$ in
the ball of radius $r$ about $\bar x(\lambda)$ in $\mathcal X$, $\lambda \mapsto
x^\star(\lambda)$ is continuous, and no branch bifurcates from $x^\star$
(Remark 3.1 for the pseudo-arclength case: the kernel of $D_{(x,\lambda)}F$
along the curve is at most one-dimensional). Injectivity of $\mathcal A$
follows from $Z_1 < 1$ (Nota Bene 2.4, 2.6, 2.7, 2.12).

**Radii polynomial (convention).** Theorem 1.1, eq. (2), asks for
$r \in (0, r_\star]$ with $Y + r Z_1 + \tfrac{r^2}{2} Z_2 < r$ and
$Z_1 + r Z_2 < 1$, where $Z_2$ is defined by eq. (1c):
$\|A(DF(x) - DF(\bar x))\| \le Z_2 \|x - \bar x\|$ for all $x \in B_{r_\star}(\bar x)$
(a Lipschitz bound of $A\,DF$ at $\bar x$ on the ball; a bound of
$\sup_{B_{r_\star}(\bar x)} \|A\,D^2F\|$ is such a bound by the mean value
inequality, and so is the radial bound used for Example 3). This is exactly
the polynomial certified by RadiiPolynomial's
`interval_of_existence(Y, Z₁, Z₂, r⋆)`, which returns an interval $[r_-, r_+]$
on which $Y + (Z_1 - 1) r + Z_2 r^2/2 \le 0$ and $Z_1 + Z_2 r < 1$. Every $Z_2$
below is such a bound (all nonlinearities are polynomial, and the bounds are
taken over the ball of radius $r_\star$), so the checker passes $Z_2$ directly,
without a factor $2$; the paper uses the same convention. Existence holds in
the ball of radius $r_-$, uniqueness in the ball of radius $r_+$.

## Results

Throughout, $\omega_n \equiv 1$ (Chebyshev weights) and $\chi_k \equiv 1$ (cosine
weights), and $r$ is the left end point $r_-$ of the interval of existence printed
by the checker, rounded up.

**Result 1 (cube root, Section 2.2.3).** Let $\bar x$ be the degree-8 Chebyshev
polynomial built by the checker from `data/cube_root.bin`. There is a unique
$x^\star \in \ell^1(\mathbb R)$ with $x^\star(\lambda)^3 = \lambda + 2$ for all
$\lambda \in [-1,1]$ and $\|x^\star - \bar x\|_{\ell^1} \le r = 6.0362 \times 10^{-7}$;
it is unique in the ball of radius $10^{-6}$.

**Result 2 (square root, Section 3.1.2).** Let $\bar u = (\bar x, \bar\lambda)$ and
$\dot u = (\dot x, \dot\lambda)$ be the degree-32 Chebyshev polynomials in
$s \in [-1,1]$ built by the checker from `data/square_root.bin`. There is a unique
$u^\star = (x^\star, \lambda^\star) \in (\ell^1(\mathbb R))^2$, within
$r = 1.81907 \times 10^{-6}$ of $\bar u$ in the norm
$\|x\|_{\ell^1} + \|\lambda\|_{\ell^1}$, such that for all $s \in [-1,1]$
$x^\star(s)^2 = \lambda^\star(s)$ and
$(x^\star(s) - \bar x(s))\dot x(s) + (\lambda^\star(s) - \bar\lambda(s))\dot\lambda(s) = 0$.
The curve runs from $(1,1)$ through the fold $(0,0)$ towards $(-1,1)$ (Figure 2).

**Result 3 (Cahn–Hilliard, Section 4.1).** Let
$\Lambda = \{ (c,\beta) : 79/2 + 150 c^2 \le \beta \le 41,\ |c| \le 1/10 \}$ and let
$\theta = (c, \beta) : [-1,1]^2 \to \Lambda$ be Coons' map (eq. (coons), Section 3.2)
of the boundary curves $\gamma^{(\text{no pinch})}$ of eq. (parabolacurves1),
resp. $\gamma^{(\text{pinch})}$ of eq. (parabolacurves2). Let $\bar u$ be the
Chebyshev–cosine polynomial of Chebyshev order $(20, 20)$ and cosine order $20$
built by the checker from `data/cahn_hilliard_<case>.bin`, and $\theta$ be
enclosed as explained in the Setting (so the statement is about the exact $\theta$).
For each of the two parameterizations there is a unique
$u^\star \in \ell^{1,\otimes 2}(\ell^1(\mathbb R))$ with
$\Delta u^\star(s) + \beta(s)\,(u^\star(s) - u^\star(s)^3 - c(s)) = 0$ for all
$s \in [-1,1]^2$, within $r^{(\text{no pinch})} = 1.02741 \times 10^{-5}$, resp.
$r^{(\text{pinch})} = 4.63965 \times 10^{-5}$, of $\bar u$; for every $s$,
$u^\star(s)$ is the unique zero of $v \mapsto \Delta v + \beta(s)(v - v^3 - c(s))$
in the ball of radius $r$ about $\bar u(s)$ (and in the ball of radius
$r_+ = 9.7895 \times 10^{-5}$, resp. $3.7706 \times 10^{-4}$), and $s \mapsto u^\star(s)$
is continuous. Since $\theta([-1,1]^2) \supseteq \Lambda$ (Proposition 3.3, see the
Setting), every $(c, \beta) \in \Lambda$ has a Neumann steady state of the
Cahn–Hilliard equation ($\beta = 1/\epsilon^2$, $c = \int_0^1 (v - v^3)$) within $r$ of
$\bar u(s)$ for some $s$ with $\theta(s) = (c, \beta)$.

*On the orders.* The paper uses $N_1 = N_2 = 16$ and reports
$r^{(\text{no pinch})} = 1.28025 \times 10^{-5}$, $r^{(\text{pinch})} = 4.11693 \times 10^{-5}$,
computed with the library norm that weighs the Chebyshev indices of tensor
spaces by $1$ (see the errata). With the paper's norm and the order-16 inputs of
the original scripts the checker gives $Y \le 4.1661 \times 10^{-5}$, $Z_1 \le 0.110983$,
$Z_2 \le 2145.76$, $r = 4.98617 \times 10^{-5}$ for "no pinch", and for "pinch"
$Y \le 1.3472 \times 10^{-4}$, $Z_1 \le 0.37302$, $Z_2 \le 2225.78$, for which the radii
polynomial has no real root ($2 Y Z_2 \approx 0.600 > (1 - Z_1)^2 \approx 0.393$):
"pinch" was not proved at order 16. The inputs were therefore regenerated by the
same producer at order $N_1 = N_2 = 20$ ($K = 20$ unchanged: the cosine
coefficients of the data decay to $10^{-15}$ at $k = 20$, while the Chebyshev
coefficients decay only to $\approx 4 \times 10^{-7}$ at $n = 16$), which gives the
bounds of the following table; the margin of the pinch case is
$2 Y Z_2 = 0.168 \le (1 - Z_1)^2/2 = 0.373$.

| case | $Y \le$ | $Z_1 \le$ (finite part; tail term) | $Z_2 \le$ | $r_\star$ | $r$ | $r_+$ |
|---|---|---|---|---|---|---|
| no pinch | $9.78959 \times 10^{-6}$ | $0.0361355$ ($0.0361355$; $0.0109259$) | $2144.93$ | $10 \sup Y$ | $1.02741 \times 10^{-5}$ | $9.7895 \times 10^{-5}$ |
| pinch | $3.77069 \times 10^{-5}$ | $0.135722$ ($0.135722$; $0.0109594$) | $2223.00$ | $10 \sup Y$ | $4.63965 \times 10^{-5}$ | $3.7706 \times 10^{-4}$ |

**Result 4 (SKT, Section 4.2).** With $d_{11} = d_{21} = d_{22} = 0$, $d_{12} = 3$,
$r_1 = 5$, $r_2 = 2$, $a_1 = a_2 = 3$, $b_1 = b_2 = 1$ and bifurcation parameter
$\lambda = d_1 = d_2$: for each of the fourteen arcs `u_bar1a`, …, `u_bar4c` of
Table 1, let $\bar u = (\bar v_1, \bar v_2, \bar\lambda)$ and $\dot u$ be the
Chebyshev–cosine polynomials in $s \in [-1,1]$ built by the checker from
`data/skt_<arc>.bin`. For every arc marked *proved* below there is a unique
$u^\star = (v^\star, \lambda^\star)$ within the radius $r$ of $\bar u$ in
$\mathcal S_{\text{pa}} = \ell^1(\ell^1(\mathbb R))^2 \times \ell^1(\mathbb R)$ such
that, for all $s$, $v^\star(s)$ is a Neumann steady state of the SKT system at
$\lambda = \lambda^\star(s)$ and $\langle u^\star(s) - \bar u(s), \dot u(s)\rangle = 0$;
for every $s$, $u^\star(s)$ is unique in the ball of radius $r_+$ about $\bar u(s)$
in $\mathcal X_{\text{pa}}$.

| arc | $N$ | arclength | $Y \le$ | $Z_1 \le$ | $Z_2 \le$ | $r$ | verdict |
|---|---|---|---|---|---|---|---|
| `u_bar1a` | 16 | 0.38 | $2.1155 \times 10^{-10}$ | 0.260746 | $3.16014 \times 10^{4}$ | $2.86171 \times 10^{-10}$ | proved |
| `u_bar1b` | 16 | 0.19 | $1.07917 \times 10^{-6}$ | 0.344129 | $2.96781 \times 10^{4}$ | $1.71169 \times 10^{-6}$ | proved |
| `u_bar1c` | 24 | 0.26 | $9.54495 \times 10^{-9}$ | 0.298871 | $1.79899 \times 10^{5}$ | $1.36376 \times 10^{-8}$ | proved |
| `u_bar2a` | 16 | 0.40 | $2.47450 \times 10^{-7}$ | 0.702583 | $9.77966 \times 10^{4}$ | $9.94647 \times 10^{-7}$ | proved |
| `u_bar2b` | 16 | 0.10 | $3.55052 \times 10^{-9}$ | 0.709760 | $1.00582 \times 10^{5}$ | $1.22591 \times 10^{-8}$ | proved |
| `u_bar2c` | **24** | 0.10 | $4.09089 \times 10^{-8}$ | 0.780618 | $1.16947 \times 10^{5}$ | $1.96796 \times 10^{-7}$ | proved |
| `u_bar2d` | 24 | 0.32 | $7.84306 \times 10^{-9}$ | 0.799500 | $5.07479 \times 10^{5}$ | $4.12733 \times 10^{-8}$ | proved |
| `u_bar3a` | 16 | 0.20 | $1.12789 \times 10^{-9}$ | 0.531279 | $1.61172 \times 10^{5}$ | $2.40730 \times 10^{-9}$ | proved |
| `u_bar3b` | 16 | 0.10 | $1.22357 \times 10^{-8}$ | 0.711077 | $1.05098 \times 10^{5}$ | $4.26806 \times 10^{-8}$ | proved |
| `u_bar3c` | **20** | 0.10 | $5.33714 \times 10^{-8}$ | 0.792230 | $1.16668 \times 10^{5}$ | $2.78682 \times 10^{-7}$ | proved |
| `u_bar3d` | 24 | 0.52 | $2.36250 \times 10^{-8}$ | 0.426328 | $1.93883 \times 10^{5}$ | $4.14727 \times 10^{-8}$ | proved ($K_A = 250$) |
| `u_bar4a` | 16 | 0.20 | $8.07141 \times 10^{-9}$ | 0.268635 | $5.33285 \times 10^{4}$ | $1.10406 \times 10^{-8}$ | proved |
| `u_bar4b` | 24 | 0.20 | $2.08562 \times 10^{-7}$ | 0.360556 | $3.30632 \times 10^{4}$ | $3.28958 \times 10^{-7}$ | proved |
| `u_bar4c` | 16 | 0.40 | $2.72124 \times 10^{-10}$ | 0.271833 | $3.84370 \times 10^{4}$ | $3.73717 \times 10^{-10}$ | proved |

The radii are those printed by the checker (condition 4.$i$.4), rounded up; the
uniqueness radii $r_+$ are between $3.9 \times 10^{-7}$ and $2.3 \times 10^{-5}$ (printed
with each condition).

*On the orders.* Table 1 of the paper uses $N = 20$ for `u_bar2c` and $N = 16$ for
`u_bar3c` and reports all fourteen arcs as proved with smaller $Z_1, Z_2, r$; those
values were computed in the wrong norm and with a wrong operator in $Z_2$ (see
the errata) and are not valid bounds. With the paper's norm and the inputs of
the original scripts, `u_bar2c` ($Y \le 7.6164 \times 10^{-7}$, $Z_1 \le 0.780766$,
$Z_2 \le 1.16972 \times 10^5$) and `u_bar3c` ($Y \le 1.11175 \times 10^{-6}$,
$Z_1 \le 0.792428$, $Z_2 \le 1.16700 \times 10^5$) have no real root of the radii
polynomial ($2 Y Z_2 \approx 0.18$, resp. $0.26$, against $(1 - Z_1)^2 \approx 0.048$,
resp. $0.043$); their inputs were regenerated by the same producer at $N = 24$,
resp. $20$, which gives the table above with margins
$2 Y Z_2 = 0.0096 \le (1 - Z_1)^2/2 = 0.024$ and $0.0125 \le 0.0216$. On both arcs
$Z_1$ is its tail term ($0.7806$, resp. $0.7922$; finite part $0.060$, resp. $0.061$),
which does not depend on $N$. For `u_bar3d` ($N = 24$ as in the paper) with the
paper's $K_A = 165$ the tail term is $0.9739$ (finite part $0.053$), so
$(1 - Z_1)^2/2 \approx 3.3 \times 10^{-4}$ against $2 Y Z_2 \approx 9.3 \times 10^{-3}$;
regenerating it at $N = 28$ and $N = 36$ (experiments, not kept) lowers $Y$ to
$5.3 \times 10^{-9}$ and $4.7 \times 10^{-9}$ but leaves $Z_1 = 0.9739$, $0.9735$, and the
arc still fails: raising $N$ cannot prove it. The tail term is
$\|I - \mathcal M(W)\mathcal M(D_v\Phi)\| + \|\mathcal M(W)\|\|\mathcal M(DR)\|/(\pi(K_A+1))^2
= 1.4 \times 10^{-7} + 0.9743$ with $\|\mathcal M(W)\| = 1.20 \times 10^4$ (the
pointwise norms $\|W(s_j)\|$ at the nodes reach $1.09 \times 10^4$) and
$\|\mathcal M(DR)\| = 22.1$; the only lever is $K_A$. The checker therefore uses
$K_A = 250$ for this arc (`SKT.arc_K_A`; $165$ for the thirteen others), which
gives the tail term $0.4263$, finite part $0.119$, $Z_1 = 0.4263$, and the radius of
the table ($2 Y Z_2 = 0.0092 \le (1 - Z_1)^2/2 = 0.164$); the arc takes about $5.5$
minutes instead of $3.5$ with two threads.

*Connection of the arcs.* The arcs of each branch (`1a→1b→1c`, `2a→2b→2c→2d`,
`3a→3b→3c→3d`, `4a→4b→4c`) were computed in sequence, each starting from the
end point of the previous one: at every junction the last node of an arc $a$
($s = -1$) and the first node of the next arc $b$ ($s = 1$) carry the same point
$u^J$ and the same tangent $\dot u^J$, bit for bit in the inputs. For eight
junctions this is inherited from the original data (the producer's Newton step
leaves a converged initial guess unchanged). After the regeneration of `u_bar2c`
and `u_bar3c` at higher orders, with finer pseudo-arclength steps, these arcs
ended $1.06 \times 10^{-4}$, resp. $1.27 \times 10^{-4}$ *before* the stored first node of
`u_bar2d`, resp. `u_bar3d` (the end point of the order-20, resp. order-16, arc of
the original scripts), a gap a thousand times the radii, so the junctions
`2c→2d` and `3c→3d` could not be proved; `u_bar2d` and `u_bar3d` were therefore
regenerated by the producer starting from the last node and tangent of the new
`u_bar2c`, `u_bar3c` (`chained_to` in `search/produce.jl`; no Newton step at
the first node), which restores the bit-identical junction data. Since the
interval interpolants take exactly the node values at the nodes, the maps
$F_{\text{pa}}^{(a)}(\cdot, -1)$ and $F_{\text{pa}}^{(b)}(\cdot, 1)$ on
$\mathcal X_{\text{pa}}$ coincide (same base point, same hyperplane). Hence
$u^\star_b(1)$ is a zero of that map within $r_b$ of $u^J$, and $u^\star_a(-1)$ is
its unique zero within $r_a^+$ of $u^J$: if $r_b \le r_a^+$ (or $r_a \le r_b^+$),
$u^\star_b(1) = u^\star_a(-1)$ and the two validated curves join continuously
(conditions 4.J1–4.J10). Outcome: the junctions `1a→1b`, `1b→1c`, `2a→2b`,
`2b→2c`, `2c→2d`, `3a→3b`, `3b→3c`, `3c→3d`, `4a→4b`, `4b→4c` all hold (in each
case $r_b \le r_a^+$ and $r_a \le r_b^+$), so the four unions
`u_bar1a ∪ u_bar1b ∪ u_bar1c`, `u_bar2a ∪ u_bar2b ∪ u_bar2c ∪ u_bar2d`,
`u_bar3a ∪ u_bar3b ∪ u_bar3c ∪ u_bar3d`, `u_bar4a ∪ u_bar4b ∪ u_bar4c` are connected
validated curves: the four branches of Figure 7. The four branches are separated by the bifurcation points of Figure 7 and are
not connected to each other.

## Setting

- **Chebyshev series** (Section 2.1): $\psi(s) = \sum_{n \in \mathbb Z} \psi_{|n|} T_{|n|}(s)$,
  stored as $(\psi_0, \psi_1, \dots, \psi_N)$; the norm counts every index of
  $\mathbb Z$, $\|\psi\|_{\ell^1} = |\psi_0| + 2\sum_{n \ge 1} |\psi_n|$, and in two
  variables (eq. (ell1_multi)) $\psi(s_1,s_2) = \sum_{n \in \mathbb Z^2} \psi_{|n_1|,|n_2|} T_{|n_1|}(s_1)T_{|n_2|}(s_2)$
  with $\|\psi\| = \sum_{n \in \mathbb Z^2} |\psi_{|n_1|,|n_2|}|$ (weight $2$ per
  nonzero index, $4$ when both are nonzero).
  **Cosine series** (Section 4): $v(y) = \sum_{k \in \mathbb Z} v_{|k|} \cos(\pi k y)$
  on $(0,1)$, stored as $(v_0, \dots, v_K)$ (`evensym(Fourier(K, π))`),
  $\|v\| = |v_0| + 2\sum_{k\ge1}|v_k|$, $\Delta v = (-\pi^2 k^2 v_k)_k$. Products are
  the discrete convolutions of Sections 2.1, 3.2 and 4 (Chebyshev and cosine
  series multiply by the same rule). The norm of
  $\ell^{1,\otimes d}_\omega(\mathcal X)$ is the sum over all Chebyshev indices in
  $\mathbb Z^d$ and all cosine indices in $\mathbb Z$ of the moduli of the
  coefficients (Fubini, Section 2.3.2); for systems the norm is the sum of the
  component norms (Remark 2.14). **In RadiiPolynomial these norms are
  `Ell1(IdentityWeight(), …)` with one weight per factor of the tensor space**
  (`X_CH`, `X_Λ` in `cahn_hilliard.jl`; `X_pa`, `X_cheb_cos` in `skt.jl`): the
  single-weight `norm(·, 1)` / `Ell1()` of RadiiPolynomial 0.11 gives weight $1$
  to every Chebyshev index of a tensor space and is **not** the paper's norm (it
  is what the original code used; see the errata). The checker verified on random
  examples that the per-factor norm equals a BigFloat brute-force evaluation of
  the paper's sums, for sequences and for operators.
- **Operator norms.** The interpolant of a grid of operators $A^{(j)}$ is stored
  by RadiiPolynomial as the operator whose column $l$ is the Chebyshev–cosine
  sequence $(\hat A_{k,l})_{k}$; its induced norm in the per-factor $\ell^1$ norm is
  $\max_l \sum_k \|\hat A_{k,l}\|_{\ell^1_\omega}$, i.e. exactly
  $\|\mathcal M(A)\|_{\mathcal B(\ell^1_\omega(\mathcal X))}$ of eq. (normAgeneral)
  (`opnorm(·, X)`). Norms of multiplication operators by a sequence are the
  sequence norms (eq. (opnormmult)); for $2\times2$ systems the matrix of block
  norms is bounded by its induced 1-norm (eq. (normAgeneral_syst)).
- **Chebyshev nodes and interpolation.** Node values are given on
  $s_j = \cos(\pi j/N)$ in the order expected by RadiiPolynomial's `to_coef`
  (or reversed, as stated per example). Interpolation is by FFT (`to_coef`), and
  values of a Chebyshev polynomial on a finer Chebyshev grid by inverse FFT
  (`to_grid`); products of Chebyshev polynomials needed in the bounds are
  formed on a Chebyshev grid of at least (degree + 1) nodes per variable and
  interpolated back (Section 2.3.2), so every interpolation is exact (the
  node counts are rounded up to $2^k + 1$ for the radix-2 FFT).
- **Arithmetic.** IEEE binary64 interval arithmetic with outward rounding:
  IntervalArithmetic.jl 1.0.12 (correctly rounded elementary functions via
  CRlibm 1.0.2 / CoreMath 0.1.0) under RadiiPolynomial.jl 0.11.13, Julia 1.13
  (`checker/Manifest.toml`). The approximate inverses and the approximate
  interpolants of floating-point data are computed in plain binary64 (LAPACK
  `inv`, FFT) and then enclosed as thin intervals; the interpolants of the
  SKT points and tangents are computed by an interval FFT.

### Example 1: cube root (Section 2.2.3)

- $F(x,\lambda) = x^3 - \lambda - 2$, $\lambda \in [-1,1]$, $\mathcal X = \mathbb R$;
  lifted map $\mathcal F(x) = x^{*3} - e_1 - 2e_0$ on $\ell^1(\mathbb R)$, $e_1 = (0, 1/2)$.
- $N = 8$; node $j$ at $\lambda_j = -\cos(\pi j/N)$, $j = 0,\dots,8$.
- $A_j = 1/(3x_j^2)$ in binary64; $\bar x$, $A$ = interpolants of $(x_j)$, $(A_j)$
  (`to_coef(reverse(·), Chebyshev(8))` in binary64, then thin intervals).
- Bounds of eq. (cubic_bounds_C0), rederived: $Y = \|A * (\bar x^{*3} - e_1 - 2e_0)\|$;
  $Z_1 = \|e_0 - 3A * \bar x^{*2}\|$ (eq. (opnormmult)); $Z_2 = 6\|A\|(\|\bar x\| + r_\star)$
  since $D^2\mathcal F(x)[h,h] = 6 x * h * h$ and $\|x\| \le \|\bar x\| + r_\star$ on the
  ball; $r_\star = 10^{-6}$. The code agrees with the paper line by line.

### Example 2: square root, pseudo-arclength (Section 3.1.2)

- $F(x,\lambda) = x^2 - \lambda$; $F_{\text{pa}}(u, s) = (x^2 - \lambda,\ (x - \bar x(s))\dot x(s) + (\lambda - \bar\lambda(s))\dot\lambda(s))$
  on $\mathcal X_{\text{pa}} = \mathbb R^2$ with the 1-norm.
- $N = 32$; node $j$ ($j = 0,\dots,32$) is the Newton point at pseudo-arclength
  $\tfrac L2(1 - \cos(\pi j/N))$ from $(1,1)$, $L = \sqrt5 + \tfrac12\ln(2+\sqrt5)$,
  i.e. at $s_j = \cos(\pi j/N)$ (interpolated without reversal).
- $A_j = D_uF_{\text{pa}}(\bar u_j, s_j)^{-1}$ (2×2, binary64); $\bar u$, $\dot u$, $A$
  = interpolants (`to_coef(·, Chebyshev(32))` in binary64, then thin intervals).
- Bounds (Section 3.1.2), rederived: $Y = \big\|(\|(\mathcal A\mathcal F_{\text{pa}}(\bar u))_i\|_{\ell^1})_i\big\|_1$,
  $Z_1 = \big\|(\|(I - A * D\mathcal F_{\text{pa}}(\bar u))_{ij}\|_{\ell^1})_{ij}\big\|_{1}$
  (induced 1-norm of the 2×2 matrix of norms, eq. (normAgeneral_syst)),
  $Z_2 = 2(\|A_{11}\|_{\ell^1} + \|A_{21}\|_{\ell^1})$, $r_\star = \infty$.
  $D\mathcal F_{\text{pa}}(u) - D\mathcal F_{\text{pa}}(\bar u)$ has the single
  nonzero block $2\mathcal M(x - \bar x)$ in position $(1,1)$, so only the first
  column of $A$ enters and this $Z_2$ (the paper's) is sharp. The original code
  used the larger $2\max(\|A_{11}\| + \|A_{21}\|, \|A_{12}\| + \|A_{22}\|)$, which
  is what the published $4.98384$ corresponds to.

### Example 3: Cahn–Hilliard, two parameters (Section 4.1, Appendix A)

- $F(v, c, \beta) = \Delta v + \beta(v - v^3 - c)$ on cosine series, $K = 20$;
  lifted map $\tilde{\mathcal F}(u) = \mathcal D^2 u + \beta * (u - u^{*3} - c)$ on
  $\ell^{1,\otimes2}(\mathcal X)$ with $(c, \beta) = \theta$ (Section 4.1, Step 1).
- $N_1 = N_2 = 20$ (the paper: 16, see Result 3); node $(j_1, j_2)$ at
  $s = (-\cos(\pi j_1/20), -\cos(\pi j_2/20))$, $j_1, j_2 = 0, \dots, 20$.
- **Coons' map, enclosed exactly.** The boundary curves of eqs.
  (parabolacurves1)–(parabolacurves2) are polynomials of degree $\le 2$ with
  rational coefficients ($c_{\max} = 1/10$, $g(c) = 79/2 + 150c^2$) and Coons'
  formula (coons) is polynomial in $(s_1, s_2)$ of degree $\le 2$ in each variable.
  The checker evaluates the same Julia functions $\gamma_j$, $\theta$ with the
  coordinate functions $s_1 = T_1(s_1)$, $s_2 = T_1(s_2)$ as interval Chebyshev
  sequences: the result $\bar\theta = (\bar c, \bar\beta)$ is an interval Chebyshev
  polynomial of order $2$ (no interpolation, no rounding other than outward
  rounding) which contains the exact $\theta$ coefficient-wise. All bounds are
  computed with $\bar\theta$, hence hold for the exact $\theta$. Condition 3.$i$.0
  checks the consistency of $\bar\theta$ with the pointwise interval evaluation
  of $\theta$ at the $21 \times 21$ nodes and the corner identities (bdy_rel)
  $\gamma_1(-1) = \gamma_4(1)$, $\gamma_1(1) = \gamma_2(-1)$, $\gamma_2(1) = \gamma_3(-1)$,
  $\gamma_3(1) = \gamma_4(-1)$ (exact rational identities; the check guards the
  implementation).
  **Coverage of $\Lambda$.** With these identities, Coons' map restricted to the
  four edges of $[-1,1]^2$ is $\gamma_1, \gamma_2, \gamma_3(-\cdot), \gamma_4(-\cdot)$
  (a computation from eq. (coons)), which traverse the Jordan curve
  $\partial\Lambda$ (the parabola $\beta = g(c)$, $|c| \le 1/10$, followed by the
  segment $\beta = 41$) exactly once, monotonically in $c$ (for "pinch", two
  edges are constant). Hence $\theta$ is continuous, $\theta(\partial[-1,1]^2) =
  \partial\Lambda$ with winding number $\pm1$ about every interior point, and the
  first part of Proposition 3.3 gives $\theta([-1,1]^2) \supseteq \Lambda$. This
  is a property of the exact formulas, not a numerical check.
- Interpolants (`to_coef(reverse(·; dims=(1,2)), Chebyshev(20) ⊗ Chebyshev(20))`
  in binary64, then thin intervals): $\bar u$ of the node data, $A_{\le K}$ of
  $A^{(j)}_{\le K} = (\Pi_{\le K} D_vF(\bar v_j, \operatorname{mid}\theta(s_j))\Pi_{\le K})^{-1}$
  (21×21, binary64). $A(s) = A_{\le K}(s)\Pi_{\le K} + \Delta_0^{-1}\Pi_{> K}$
  ($\Delta_0^{-1}$ = pseudo-inverse of $\Delta$, zero on $k = 0$).
- **Bounds** (Section 4.1, Step 3, rederived from Appendix A; $N_\theta = 2$ is the
  order of $\bar\theta$, $N = 20$):
  - $Y = \|A\,\tilde{\mathcal F}(\bar u)\|$. $A$ has degree $N$ and
    $\bar\beta\bar u^{*3}$ degree $3N + N_\theta$ in each variable, so the product has
    degree $4N + N_\theta = 82$; it is formed on the $129 \times 129$ grid
    ($129 = 2^7 + 1 \ge 83$) and interpolated in Chebyshev$(128)^{\otimes 2}$, exactly.
    (The cosine degree is $3K$, so $\Pi_{3K}$ of the paper is automatic.)
  - $Z_1 = \max\big(\|\Pi_{\le q} - A\, D\tilde{\mathcal F}(\bar u)\,\Pi_{\le q}\|,\
    \|\bar\beta * (e_{0,0} - 3\bar u^{*2})\| / (\pi(K+1))^2\big)$ with $q = 3K$, as in the paper.
    Derivation: $D\tilde F(\bar u) = \Delta + \mathcal M(\bar\beta(1 - 3\bar u^2))$ is
    banded of bandwidth $2K$ in the cosine index ($\bar u^{*2}$ has $2K$ modes),
    so for $q \ge 3K$: $D\tilde F\,\Pi_{>q} = \Pi_{>K}\, D\tilde F\,\Pi_{>q}$ and
    $A\,\Pi_{>K} = \Delta_0^{-1}\Pi_{>K}$, whence
    $\Pi_{>q} - A D\tilde F \Pi_{>q} = -\Delta_0^{-1}\Pi_{>K}\mathcal M(\bar\beta(1-3\bar u^2))\Pi_{>q}$,
    of norm $\le \|\bar\beta(1 - 3\bar u^2)\| / (\pi^2 (K+1)^2)$ since the diagonal
    $\Delta_0^{-1}\Pi_{>K}$ has entries $-1/(\pi^2k^2)$, $k \ge K+1$ (same weight on
    row and column). The $\ell^1$ operator norm is the maximum over columns, hence
    the `max` of the two column blocks. The finite block has degree
    $3N + N_\theta = 62$ and is formed on the $65 \times 65$ grid; the tail sequence
    has degree $2N + N_\theta$ and is formed on the same grid.
  - $Z_2 = 3\|\bar\beta\|\,\|\mathcal A\|\,(2\|\bar u\| + r_\star)$, $\|\mathcal A\| =
    \max(\|\mathcal A_{\le K}\|, 1/(\pi(K+1))^2)$, $r_\star = 10\sup Y$:
    $\mathcal A(D\tilde{\mathcal F}(u) - D\tilde{\mathcal F}(\bar u)) = -3\mathcal A\mathcal M(\bar\beta * (u + \bar u) * (u - \bar u))$
    and $\|u + \bar u\| \le 2\|\bar u\| + r_\star$ on the ball (a radial Lipschitz bound
    in the sense of eq. (1c)). The columns $> K$ of $A(s)$ are the constant
    $\Delta_0^{-1}\Pi_{>K}$.
- Norms: $\|\cdot\|$ on sequences of $\ell^{1,\otimes2}(\mathcal X)$ is `X_CH`
  (three factors), on $\bar\beta$ it is `X_Λ` (two factors).

### Example 4: SKT, pseudo-arclength, fourteen arcs (Section 4.2)

- $F(v,\lambda) = \Delta\Phi(v,\lambda) + R(v)$, $F_{\text{pa}}(u, s) = (F(v,\lambda),\ \langle u - \bar u(s), \dot u(s)\rangle)$
  with $\langle\cdot,\cdot\rangle$ the Euclidean product of the stored coefficient
  vectors; $u = (v_1, v_2, \lambda)$ on cosine series with $K = 80$
  (`evensym(Fourier(80, π)) × evensym(Fourier(80, π)) × ScalarSpace()`, 163 coefficients).
- Per arc: target arclength $L$ as in Table 1 and Chebyshev order $N$ as in the
  table of Result 4 (Table 1 of the paper, except $N = 24, 20$ for `u_bar2c`,
  `u_bar3c`); node $j$ at pseudo-arclength
  $\tfrac L2(1 - \cos(\pi j/N))$ from the start of the arc, i.e. at $s_j = \cos(\pi j/N)$
  (the arc starts at $s = 1$ and ends at $s = -1$).
- $K_A = 165$ ($250$ for `u_bar3d`, see Result 4; `SKT.arc_K_A`). At each node, in binary64: $A^{(j)}_{\le K_A} = \Pi_{\le K_A}(\Pi_{\le 2K_A} D_uF_{\text{pa}}(\bar u_j, s_j)\Pi_{\le 2K_A})^{-1}\Pi_{\le K_A}$
  and $W^{(j)} \approx D_v\Phi(\bar v_j, \bar\lambda_j)^{-1}$ (pointwise 2×2 inverses on a
  $2K_A+1$ point grid, interpolated in cosine series of order $K_A$).
  $\bar u$, $\dot u$ = interval interpolants of the node data (`to_coef` of thin
  intervals, Chebyshev$(N)$; they take the node values at the nodes); $A_{\le K_A}$, $W$ = binary64 interpolants, then thin
  intervals. $A(s) = A_{\le K_A}(s)\Pi_{\le K_A} + \mathrm{diag}((W(s) - \Pi_{\le K_A}W(s)\Pi_{\le K_A})\Delta_0^{-1}, 0)$;
  its columns $> K_A$ are $\mathrm{diag}(W(s)\Delta_0^{-1}\Pi_{>K_A}, 0)$.
- **Bounds** (Section 4.2, rederived), all norms in $\mathcal S_{\text{pa}}$ (`X_pa`:
  per-factor weights on each component, sum over the components):
  - $Y = \|\mathcal A\,\Pi_{\le 2K}\mathcal F_{\text{pa}}(\bar u)\|$, formed on the
    $3N+1$ grid ($A$ has degree $N$, $\mathcal F_{\text{pa}}(\bar u)$ degree $2N$; its
    cosine degree is $2K$ so $\Pi_{\le 2K}$ is automatic; the last component of
    $\mathcal F_{\text{pa}}(\bar u)$ vanishes).
  - $Z_1 = \max(Z_1^{\text{finite}}, Z_1^{\text{tail}})$ with
    $Z_1^{\text{finite}} = \|\Pi_{\le K+K_A} - \mathcal A D\mathcal F_{\text{pa}}(\bar u)\Pi_{\le K+K_A}\|$
    (degree $2N$, formed on the $2N+1$ grid) and
    $Z_1^{\text{tail}} = \|I - \mathcal M(W)\mathcal M(D_v\Phi(\bar v,\bar\lambda))\| + \|\mathcal M(W)\|\,\|\mathcal M(DR(\bar v))\|/(\pi(K_A+1))^2$.
    Derivation: $D_uF_{\text{pa}}(\bar u)$ is banded of bandwidth $K$ ($D_v\Phi$, $DR$
    are linear in $\bar v$) and $\dot u \in \Pi_{\le K}$, so for $q = K + K_A$ the
    columns $> q$ only meet the columns $> K_A$ of $A$:
    $\Pi_{>q} - A D F_{\text{pa}} \Pi_{>q} = (I - W D_v\Phi)\Pi_{>q} - W\Delta_0^{-1}\Pi_{>K_A} DR\,\Pi_{>q}$
    on the $v$-components (the $\lambda$ row and column vanish on $\Pi_{>q}$). The
    $2\times2$ block norms are bounded by eq. (normAgeneral_syst); $M = I - W D_v\Phi$
    has degree $2N$ and $DR$ degree $N$, both formed on the $2N+1$ grid.
  - $Z_2 = \|\mathcal A\,\mathrm{diag}(\Delta, 0)\|\,C_1 + \|\mathcal A\|\,C_2$, $r_\star = \infty$,
    with $\|\mathcal A\| \le \max(\|\mathcal A\Pi_{\le K_A}\|, \|\mathcal M(W)\|/(\pi(K_A+1))^2)$
    and $\|\mathcal A\,\mathrm{diag}(\Delta,0)\| \le \max(\|\mathcal A\Pi_{\le K_A}\mathrm{diag}(\Delta,0)\|, \|\mathcal M(W)\|)$
    (columns $> K_A$: $W\Delta_0^{-1}\Pi_{>K_A}\Delta = W\Pi_{>K_A}$). Here
    $\mathcal A\Pi_{\le K_A}\mathrm{diag}(\Delta,0)$ is the interpolant of
    $A(s)\,\mathrm{diag}(\Delta,0)\,\Pi_{\le K_A}$ (`AΔ_cheb_trunc`; the original code
    used $\|\mathcal A\Pi_{\le K_A}\|$ in its place, see the errata).
    Constants: $F_{\text{pa}}$ is quadratic, $D^2F_{\text{pa}}[w,h] = (\Delta D^2\Phi[w,h] + D^2R[w,h], 0)$
    with, for $w = (w_1, w_2, w_\lambda)$, $h = (h_1, h_2, h_\lambda)$,
    $D^2\Phi_1[w,h] = w_\lambda h_1 + h_\lambda w_1 + 2d_{11}w_1h_1 + d_{12}(w_1h_2 + w_2h_1)$,
    $D^2\Phi_2[w,h] = w_\lambda h_2 + h_\lambda w_2 + d_{21}(w_1h_2 + w_2h_1) + 2d_{22}w_2h_2$
    ($\Phi$ contains $\lambda v_i$, so the mixed derivative in $(v,\lambda)$ is present),
    $D^2R_1[w,h] = -2a_1w_1h_1 - b_1(w_1h_2 + w_2h_1)$, $D^2R_2[w,h] = -b_2(w_1h_2+w_2h_1) - 2a_2w_2h_2$.
    With $a_i = \|w_i\|$, $b_i = \|h_i\|$ ($\sum a_i = \sum b_i = 1$), the Banach
    algebra property gives $\|D^2\Phi[w,h]\| \le \sum_{ij} M_{ij} a_i b_j$ with
    $M = \begin{pmatrix} 2|d_{11}| & |d_{12}|+|d_{21}| & 1 \\ |d_{12}|+|d_{21}| & 2|d_{22}| & 1 \\ 1 & 1 & 0\end{pmatrix}$,
    whose maximum over the simplices is the largest entry:
    $C_1 = \max(2|d_{11}|, |d_{12}|+|d_{21}|, 2|d_{22}|, 1) = 3$; likewise
    $C_2 = \max(2|a_1|, |b_1|+|b_2|, 2|a_2|) = 6$. These are the constants of the
    original code (the paper's displayed formulas differ, see the errata).

## Inputs (in `data/`)

All inputs are frozen: they are the floating-point outputs of the Newton
continuations of the original scripts (Step 2 of each script), i.e. the
approximate solutions on the Chebyshev nodes and, for the pseudo-arclength
examples, the approximate tangents that define the hyperplanes of
$F_{\text{pa}}$. They were produced by `julia --project=search search/produce.jl`
(the `produce` command of `entry.json`), which is not part of the proof: any
other data passing the conditions below would prove the results as well (with
their own $\bar x$). The SKT continuation starts from the initial guesses of the
original repository (`skt_data.jld2`, converted bit-exactly to
`search/skt_initial_guesses.bin`, 14 × 163 binary64 in the order of the table
below); the Cahn–Hilliard initial guess contains `rand()` (as in the original
script), seeded with `Random.seed!(20260922)`. The cube-root and square-root
files and the ten SKT files that were not regenerated are the node data computed
by the unmodified original scripts (instrumented only to dump them) under the
original environment (Julia 1.12.7, RadiiPolynomial 0.11.5); `produce`
regenerated them bit-exactly under Julia 1.13.2. Under Julia 1.13.3 (the
environment of this revision) `produce` still reproduces the cube-root and
square-root files bit-exactly, but those SKT files only up to $1.6 \times 10^{-10}$
(Newton-tolerance differences of the continuation; they were left as they were).
**New in this revision** are the two Cahn–Hilliard files and the two SKT files
`skt_u_bar2c.bin`, `skt_u_bar3c.bin`, `skt_u_bar2d.bin`, `skt_u_bar3d.bin`: they
are the output of `produce` (`search/produce.jl`, the original continuations
with their binary64 Coons map and initial guesses, orders taken from
`CahnHilliard.N₁, N₂, K` and `SKT.arc_N`) at the raised orders $N_1 = N_2 = 20$,
resp. $N = 24, 20$ (and, for `u_bar2d`, `u_bar3d`, at the paper's order $24$ but
started from the last node of the new `u_bar2c`, `u_bar3c`), under Julia 1.13.3,
RadiiPolynomial 0.11.10, and `produce` regenerated them bit-exactly in that
environment. The first two replace the order-16 (resp. 20, 16) files of the
original scripts because those cases are not proved at the paper's orders in the
paper's norm (Results 3 and 4); the last two replace files that no longer joined
the regenerated arcs (Result 4, connection of the arcs). The checker's Coons map is the exact one (Setting,
Example 3).

**Recomputed by the checker, not stored:** the approximate inverses ($A_j$,
$A^{(j)}_{\le K}$, $A^{(j)}_{\le K_A}$, $W^{(j)}$) at the nodes, the enclosure
$\bar\theta$ of Coons' map, and all Chebyshev interpolants. The approximate
inverses are computed in binary64 as in the original scripts (for Example 3 at
the midpoints of the interval values $\theta(s_j)$); they only need to be some
bounded operators (their injectivity follows from $Z_1 < 1$), so the conditions
do not depend on their exact values, and the interpolants are part of the
definition of $\bar x$, $\bar u$, $\dot u$.

Files (raw little-endian IEEE binary64, no header; "node $j$" as in the Setting):

| file | content | values |
|---|---|---|
| `cube_root.bin` | $x_j$, $j = 0,\dots,8$ | 9 |
| `square_root.bin` | $(x_j, \lambda_j)$ for $j = 0,\dots,32$, then $(\dot x_j, \dot\lambda_j)$ for $j = 0,\dots,32$ | 132 |
| `cahn_hilliard_no_pinch.bin`, `cahn_hilliard_pinch.bin` | cosine coefficients $v_{j_1,j_2,k}$, $k = 0,\dots,20$ fastest, then $j_1 = 0,\dots,20$, then $j_2 = 0,\dots,20$ | 9261 (74 088 bytes) |
| `skt_<arc>.bin` (14 files) | for $j = 0,\dots,N$: $u_j = (v_{1,0..80}, v_{2,0..80}, \lambda)$ (163 values); then for $j = 0,\dots,N$: $\dot u_j$ (163 values) | $2 \cdot 163 (N+1)$ |

with $N$ = 16, 16, 24, 16, 16, **24**, 24, 16, 16, **20**, 24, 16, 24, 16 for
`u_bar1a`, `u_bar1b`, `u_bar1c`, `u_bar2a`, `u_bar2b`, `u_bar2c`, `u_bar2d`,
`u_bar3a`, `u_bar3b`, `u_bar3c`, `u_bar3d`, `u_bar4a`, `u_bar4b`, `u_bar4c`
(`SKT.arc_N`; the paper's Table 1 has 20 and 16 for the two arcs in bold), i.e.
44 336, 54 768 or 65 200 bytes per file.
SHA-256 digests are in `data/hashes.json`; the checker verifies them before
reading the files. Only the two Cahn–Hilliard files and the four SKT files
`skt_u_bar2c.bin`, `skt_u_bar3c.bin`, `skt_u_bar2d.bin`, `skt_u_bar3d.bin` (and
their digests) changed in this revision.

## Conditions to verify

Each is a finite computation in outward-rounded interval arithmetic on the data
above. The checker prints them in this order, each prefixed `pass` or `fail`,
with the enclosure obtained (endpoints at full binary64 precision). For each
application $k$ of Theorem 1.1:

- $k$.1 **$Y$** is enclosed by a bounded nonnegative interval;
- $k$.2 **$Z_1$** is enclosed by a bounded nonnegative interval;
- $k$.3 **$Z_2$** is enclosed by a bounded nonnegative interval;
- $k$.4 **Radii polynomial (Theorem 1.1, conditions (2a)–(2b)).**
  `RadiiPolynomial.interval_of_existence(Y, Z₁, Z₂, r⋆)` returns a non-empty
  interval $[r_-, r_+] \subset [0, r_\star]$ on which $Y + (Z_1 - 1)r + Z_2 r^2/2 \le 0$
  and $Z_1 + Z_2 r < 1$ (it requires $\sup Z_1 < 1$ and a nonnegative discriminant).
  Conclusion: existence in the ball of radius $r_-$, uniqueness in the ball of
  radius $r_+$.

The applications and the additional conditions are, in order:

1. Cube root (Result 1): conditions 1.1–1.4.
2. Square root (Result 2): conditions 2.1–2.4.
3. Cahn–Hilliard (Result 3): 3.$i$.0–3.$i$.4 for $i = 1$ "no pinch", $i = 2$ "pinch",
   where 3.$i$.0 is: **$\bar\theta$ encloses Coons' map**: the corner identities
   (bdy_rel) hold (the interval differences contain $0$) and $\bar\theta$ intersects
   the pointwise interval evaluation of $\theta$ at every node.
4. SKT (Result 4): 4.$i$.1–4.$i$.4 for the $i$-th arc, $i = 1,\dots,14$
   (`u_bar1a`, …, `u_bar4c` in the order of the table above), then
   4.J1–4.J10 **connection** for the ten junctions `1a→1b`, `1b→1c`, `2a→2b`,
   `2b→2c`, `2c→2d`, `3a→3b`, `3b→3c`, `3c→3d`, `4a→4b`, `4b→4c`: both arcs
   proved, the stored last node (point and tangent) of the first arc equals the
   stored first node of the second, and $r_b \le r_a^+$ or $r_a \le r_b^+$.

**Implication.** Conditions 1.1–1.4 $\Rightarrow$ Result 1; 2.1–2.4 $\Rightarrow$
Result 2; 3.1.0–3.1.4 and 3.2.0–3.2.4 $\Rightarrow$ Result 3 (both cases); for
each arc, 4.$i$.1–4.$i$.4 $\Rightarrow$ the statement of Result 4 for that arc,
and 4.J$j$ $\Rightarrow$ the two arcs of junction $j$ join. With the frozen
inputs all 84 conditions (plus the 18 digest checks) pass, and the checker exits
with code 0: Results 1–4 hold in full. (After the last revision, 4.11 and 4.J5,
4.J8, the conditions were re-checked in parts, not in one end-to-end run; the
reported values of this document are the ones printed by those runs.)

## Revisions and paper errata

The original code (and the first version of this package, which reproduced it)
differs from the paper's bounds in the points (a)–(d) below; (e) adds the two
hypotheses that were previously unverified. All changes are implemented in
`checker/src/cahn_hilliard.jl`, `checker/src/skt.jl` and `checker/verify.jl`.

- **(a) Norms.** The original code took all norms of Examples 3 and 4 with
  RadiiPolynomial's `norm(·, 1)` / `Ell1()`, which on a tensor space with a
  Chebyshev factor weighs every Chebyshev index by $1$ instead of $2$ for
  $n \ge 1$ (for instance $\|\bar\beta\| = 40.72 < \sup\bar\beta = 41$ for "no
  pinch"; the coefficients $(1,2,3,4)$ on Chebyshev(1)$^{\otimes 2}$ get norm $10$
  while the paper's norm, eq. (ell1_multi), gives $27 = \psi(1,1)$). The published
  bounds of Sections 4.1–4.2 (the $Y, Z_1, Z_2, r$ of the Cahn–Hilliard
  paragraph and Table 1) are therefore not bounds in the paper's norm and the
  inequality $\sup_s\|\psi(s)\| \le \|\psi\|$ (eq. (C0vsell1)) that Theorem 2.11
  relies on does not hold for them. The checker now uses the per-factor weights.
  Consequence: at the paper's orders the "pinch" case of Result 3 and the arcs
  `u_bar2c`, `u_bar3c`, `u_bar3d` of Result 4 are not proved (hence the raised
  orders of the first three and the raised $K_A$ of `u_bar3d`); the valid radii
  of the other results are larger than the published ones (tables above).
- **(b) $Z_1$ of Example 3.** The original code omitted the tail term
  $\|\bar\beta(1 - 3\bar u^2)\|/(\pi^2(K+1)^2)$ of Section 4.1 (columns $> q$). It is
  smaller than the finite part here ($\approx 1.1 \times 10^{-2}$), so $Z_1$ is unchanged
  numerically, but it is now computed. *Erratum (Appendix A):* the appendix states
  that $DF(\bar v)$ is banded of bandwidth $K$ and splits at $q = 2K$
  ($\Pi_{\le 2K} - \Pi_{\le 3K}A\,DF\,\Pi_{\le 2K}$, $\Pi_{>2K} - A\Pi_{>K}DF\Pi_{>2K}$);
  for the cubic nonlinearity $\bar v^{\odot 2}$ has $2K$ modes, the bandwidth is
  $2K$, and the splitting must be at $q \ge 3K$ (as the main text, Section 4.1,
  correctly uses with $\Pi_{3K}$ and $\Pi_{5K}$). The identities
  $DF(\bar v)\Pi_{\le q} = \Pi_{\le q+K}DF(\bar v)\Pi_{\le q}$ should read $q + 2K$.
  The published numbers are not affected (the original code used $q = 3K+1$,
  the checker uses $q = 3K$; both are $\ge 3K$).
- **(c) Grid of $Y$ in Example 3.** $A\,\tilde{\mathcal F}(\bar u)$ has degree
  $4N + N_\theta$ (with the original binary64 interpolant of $\theta$ of order
  $N$: $5N$) but was sampled on $4N+1$ nodes, so the original $Y$ was the norm of
  an interpolant of the product rather than of the product. The grids are now
  sized by the degree count (Setting). Effect on $Y$: $+0.02\%$.
- **(d) $Z_2$ of Example 4.** The original code computed $\|\mathcal A\,\mathrm{diag}(\Delta,0)\|$
  as $\max(\|\mathcal A\Pi_{\le K_A}\|, \|\mathcal M(W)\|)$, i.e. with the norm of
  $\mathcal A\Pi_{\le K_A}$ instead of $\mathcal A\Pi_{\le K_A}\mathrm{diag}(\Delta,0)$
  (the latter was built but not used). The paper's formula is correct and is now
  used; the published $Z_2$ (and $r$) of Table 1 are too small for this reason as
  well as for (a). *Errata (Section 4.2, $Z_2$):* (i) the derivation bounds
  $\|\mathcal M(D_v^2\Phi)\|$ only, but $\Phi$ depends on $\lambda$ ($\Phi_i \ni \lambda v_i$)
  and $D_u^2 F_{\text{pa}}$ contains the mixed term $w_\lambda h_i + h_\lambda w_i$; the
  constant must be $C_1 = \max(2|d_{11}|, |d_{12}|+|d_{21}|, 2|d_{22}|, 1)$ (the
  original code has the $1$); with the paper's parameters $C_1 = 3$ either way.
  (ii) The displayed $C_1 = \max(\max(|2d_{11}|,|d_{12}|)+|d_{21}|, |d_{12}|+\max(|d_{21}|,|2d_{22}|))$
  and $C_2 = \max(\max(|2a_1|,|b_1|)+|b_2|, |b_1|+\max(|b_2|,|2a_2|))$ are valid
  but not sharp bounds of the bilinear operators; the sharp ones are the largest
  pairing coefficients $\max(2|d_{11}|, |d_{12}|+|d_{21}|, 2|d_{22}|)$ and
  $\max(2|a_1|, |b_1|+|b_2|, 2|a_2|)$. With the paper's parameters the displayed
  $C_2$ equals $7$ while Table 1 was computed with the sharp value $6$ (the
  original code): text and numbers are inconsistent; the checker uses $6$.
- **(e) Previously unverified hypotheses.** (i) Example 3: the proof was about
  the binary64 interpolant of Coons' map; $\bar\theta$ now encloses the exact
  $\theta$ and the coverage $\theta([-1,1]^2) \supseteq \Lambda$ follows from
  Proposition 3.3 (Setting, Example 3; condition 3.$i$.0). *Erratum (Proposition
  3.3, first statement):* continuity and $\theta(\partial D) = \partial\Lambda$ alone do not
  imply $\theta(D) \supseteq \Lambda$ (a map of $D$ into $\partial\Lambda$ whose
  restriction to $\partial D$ runs along $\partial\Lambda$ and back is a
  counterexample); the proof uses, and the statement needs, that $\theta|_{\partial D}$
  has nonzero winding number about the interior points (e.g. traverses the
  Jordan curve $\partial\Lambda$ once), which Coons' map of monotone boundary curves
  satisfies. (ii) Example 4: the connection of consecutive arcs is now a checked
  condition (4.J$j$), using that the stored junction data coincide bit for bit
  (eight of the ten junctions; the two others have a gap, see Result 4).
- **Other remarks on the paper**: Section 3.1.2, $Z_2$: the text's
  $2(\|A_{11}\|+\|A_{21}\|)$ is sharp and is what the checker computes; the
  published $4.98384$ is the coarser $2\|(\|A_{ij}\|)\|_1$ of the original code.
  Section 4.1: $A_{\le K}(s)$ is
  $(K+1)\times(K+1)$, not $K\times K$; in "Results of the code", "(with pinch
  (parabolacurves1) and without (parabolacurves2))" has the labels swapped.
  Section 2.2.3 and 3.1.2 are unaffected by (a) (one-variable Chebyshev norms
  were correct in the library).
