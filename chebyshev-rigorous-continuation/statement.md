# Statement and verification conditions

**Theorem.** The four results below hold. Each one is a continuation problem
$F(x,\lambda) = 0$ (or its pseudo-arclength version) lifted to the map
$\mathcal F$ on $\ell^1_\omega(\mathcal X)$, the space of Chebyshev series in the
parameter with coefficients in $\mathcal X$ (Section 2.1, Definition 2.8), with
weights $\omega_n \equiv 1$ and $\chi_k \equiv 1$. Theorem 1.1 (Newton–Kantorovich)
applied to $\mathcal F$ at an approximate zero $\bar x$ with approximate inverse
$\mathcal A = \mathcal M(A)$ gives a unique zero $x^\star$ of $\mathcal F$ with
$\|x^\star - \bar x\|_{\ell^1_\omega(\mathcal X)} \le r$. By Theorem 2.11, for every
parameter value $x^\star(\lambda)$ is then the unique zero of $F(\cdot,\lambda)$ in
the ball of radius $r$ about $\bar x(\lambda)$ in $\mathcal X$, $\lambda \mapsto
x^\star(\lambda)$ is continuous, and no branch bifurcates from $x^\star$
(Remark 3.1 in the pseudo-arclength case). Injectivity of $\mathcal A$ follows
from $Z_1 < 1$ (Nota Bene 2.4, 2.6, 2.7, 2.12). The radius $r$ is the left end
point $r_-$ of the interval of existence printed by the checker (condition $k$.4),
rounded up; uniqueness also holds in the larger ball of radius $r_+$, its right
end point.

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
built by the checker from `data/cahn_hilliard_<case>.bin`, and let $\theta$ be
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

**Result 4 (SKT, Section 4.2).** With $d_{11} = d_{21} = d_{22} = 0$, $d_{12} = 3$,
$r_1 = 5$, $r_2 = 2$, $a_1 = a_2 = 3$, $b_1 = b_2 = 1$ and bifurcation parameter
$\lambda = d_1 = d_2$: for each of the fourteen arcs `u_bar1a`, …, `u_bar4c`
below, let $\bar u = (\bar v_1, \bar v_2, \bar\lambda)$ and $\dot u$ be the
Chebyshev–cosine polynomials in $s \in [-1,1]$ of Chebyshev order $N$ built by the
checker from `data/skt_<arc>.bin`. For every arc there is a unique
$u^\star = (v^\star, \lambda^\star)$ within the radius $r$ of $\bar u$ in
$\mathcal S_{\text{pa}} = \ell^1(\ell^1(\mathbb R))^2 \times \ell^1(\mathbb R)$ such
that, for all $s$, $v^\star(s)$ is a Neumann steady state of the SKT system at
$\lambda = \lambda^\star(s)$ and $\langle u^\star(s) - \bar u(s), \dot u(s)\rangle = 0$;
for every $s$, $u^\star(s)$ is unique in the ball of radius $r_+$ about $\bar u(s)$
in $\mathcal X_{\text{pa}}$ ($r_+$ is printed with condition 4.$i$.4).

| arc | $N$ | arclength $L$ | $r$ |
|---|---|---|---|
| `u_bar1a` | 16 | 0.38 | $2.86171 \times 10^{-10}$ |
| `u_bar1b` | 16 | 0.19 | $1.71169 \times 10^{-6}$ |
| `u_bar1c` | 24 | 0.26 | $1.36376 \times 10^{-8}$ |
| `u_bar2a` | 16 | 0.40 | $9.94647 \times 10^{-7}$ |
| `u_bar2b` | 16 | 0.10 | $1.22591 \times 10^{-8}$ |
| `u_bar2c` | 24 | 0.10 | $1.96796 \times 10^{-7}$ |
| `u_bar2d` | 24 | 0.32 | $4.12733 \times 10^{-8}$ |
| `u_bar3a` | 16 | 0.20 | $2.40730 \times 10^{-9}$ |
| `u_bar3b` | 16 | 0.10 | $4.26806 \times 10^{-8}$ |
| `u_bar3c` | 20 | 0.10 | $2.78682 \times 10^{-7}$ |
| `u_bar3d` | 24 | 0.52 | $4.14727 \times 10^{-8}$ |
| `u_bar4a` | 16 | 0.20 | $1.10406 \times 10^{-8}$ |
| `u_bar4b` | 24 | 0.20 | $3.28958 \times 10^{-7}$ |
| `u_bar4c` | 16 | 0.40 | $3.73717 \times 10^{-10}$ |

Moreover, consecutive arcs of the same branch join (see "Connection of the
arcs" in the Setting): the four unions `u_bar1a ∪ u_bar1b ∪ u_bar1c`,
`u_bar2a ∪ u_bar2b ∪ u_bar2c ∪ u_bar2d`, `u_bar3a ∪ u_bar3b ∪ u_bar3c ∪ u_bar3d`,
`u_bar4a ∪ u_bar4b ∪ u_bar4c` are connected validated curves: the four branches
of Figure 7. They are separated by the bifurcation points of Figure 7 and are not
connected to each other.

Reference: M. Breden, O. Hénot, *Efficient Rigorous Continuation via Chebyshev
Series Expansion I*, arXiv:2609.22689 (v1). Section, theorem, remark and equation
numbers below refer to it; the paper defines every object named here. The paper
states the four results in its "Results of the code" paragraphs (Sections 2.2.3,
3.1.2, 4.1, 4.2); they are not numbered theorems.

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
  component norms (Remark 2.14). In the checker these norms are
  `Ell1(IdentityWeight(), …)` with one weight per factor of the tensor space
  (`X_CH`, `X_Λ` in `cahn_hilliard.jl`; `X_pa`, `X_cheb_cos` in `skt.jl`).
- **Operator norms.** The interpolant of a grid of operators $A^{(j)}$ is the
  operator whose column $l$ is the Chebyshev–cosine sequence $(\hat A_{k,l})_{k}$;
  its induced norm in the $\ell^1$ norm above is
  $\max_l \sum_k \|\hat A_{k,l}\|_{\ell^1_\omega}$, i.e.
  $\|\mathcal M(A)\|_{\mathcal B(\ell^1_\omega(\mathcal X))}$ of eq. (normAgeneral)
  (`opnorm(·, X)`). The norm of a multiplication operator by a sequence is the
  sequence norm (eq. (opnormmult)); for $2\times2$ systems the matrix of block
  norms is bounded by its induced 1-norm (eq. (normAgeneral_syst)).
- **Chebyshev nodes and interpolation.** Node values are given on
  $s_j = \cos(\pi j/N)$ in the order expected by RadiiPolynomial's `to_coef`
  (or reversed, as stated per example). Interpolation is by FFT (`to_coef`), and
  values of a Chebyshev polynomial on a finer Chebyshev grid by inverse FFT
  (`to_grid`). Products of Chebyshev polynomials needed in the bounds are
  formed on a Chebyshev grid of at least (degree + 1) nodes per variable and
  interpolated back (Section 2.3.2), so every interpolation is exact (the
  node counts are rounded up to $2^k + 1$ for the radix-2 FFT).
- **Radii polynomial.** Theorem 1.1, eq. (2), asks for $r \in (0, r_\star]$ with
  $Y + r Z_1 + \tfrac{r^2}{2} Z_2 < r$ and $Z_1 + r Z_2 < 1$, where $Z_2$ is a
  Lipschitz bound of $A\,DF$ at $\bar x$ on the ball $B_{r_\star}(\bar x)$, eq. (1c).
  Every $Z_2$ below is such a bound (all nonlinearities are polynomial and the
  bounds are taken over the ball of radius $r_\star$). RadiiPolynomial's
  `interval_of_existence(Y, Z₁, Z₂, r⋆)` returns an interval $[r_-, r_+]$ on which
  $Y + (Z_1 - 1) r + Z_2 r^2/2 \le 0$ and $Z_1 + Z_2 r < 1$; existence holds in
  the ball of radius $r_-$, uniqueness in the ball of radius $r_+$.
- **Arithmetic.** Outward-rounded binary64 intervals (RadiiPolynomial.jl on IntervalArithmetic.jl, versions pinned in `checker/Manifest.toml`). The approximate
  inverses and the approximate interpolants are plain binary64, then used as thin
  intervals; the interpolants of the SKT points and tangents are computed by an
  interval FFT.

### Example 1: cube root (Section 2.2.3)

- $F(x,\lambda) = x^3 - \lambda - 2$, $\lambda \in [-1,1]$, $\mathcal X = \mathbb R$;
  lifted map $\mathcal F(x) = x^{*3} - e_1 - 2e_0$ on $\ell^1(\mathbb R)$, $e_1 = (0, 1/2)$.
- $N = 8$; node $j$ at $\lambda_j = -\cos(\pi j/N)$, $j = 0,\dots,8$.
- $A_j = 1/(3x_j^2)$ in binary64; $\bar x$, $A$ = interpolants of $(x_j)$, $(A_j)$
  (`to_coef(reverse(·), Chebyshev(8))` in binary64, then thin intervals).
- Bounds of eq. (cubic_bounds_C0): $Y = \|A * (\bar x^{*3} - e_1 - 2e_0)\|$;
  $Z_1 = \|e_0 - 3A * \bar x^{*2}\|$ (eq. (opnormmult)); $Z_2 = 6\|A\|(\|\bar x\| + r_\star)$
  since $D^2\mathcal F(x)[h,h] = 6 x * h * h$ and $\|x\| \le \|\bar x\| + r_\star$ on the
  ball; $r_\star = 10^{-6}$.

### Example 2: square root, pseudo-arclength (Section 3.1.2)

- $F(x,\lambda) = x^2 - \lambda$; $F_{\text{pa}}(u, s) = (x^2 - \lambda,\ (x - \bar x(s))\dot x(s) + (\lambda - \bar\lambda(s))\dot\lambda(s))$
  on $\mathcal X_{\text{pa}} = \mathbb R^2$ with the 1-norm.
- $N = 32$; node $j$ ($j = 0,\dots,32$) is the Newton point at pseudo-arclength
  $\tfrac L2(1 - \cos(\pi j/N))$ from $(1,1)$, $L = \sqrt5 + \tfrac12\ln(2+\sqrt5)$,
  i.e. at $s_j = \cos(\pi j/N)$ (interpolated without reversal).
- $A_j = D_uF_{\text{pa}}(\bar u_j, s_j)^{-1}$ (2×2, binary64); $\bar u$, $\dot u$, $A$
  = interpolants (`to_coef(·, Chebyshev(32))` in binary64, then thin intervals).
- Bounds (Section 3.1.2): $Y = \big\|(\|(\mathcal A\mathcal F_{\text{pa}}(\bar u))_i\|_{\ell^1})_i\big\|_1$,
  $Z_1 = \big\|(\|(I - A * D\mathcal F_{\text{pa}}(\bar u))_{ij}\|_{\ell^1})_{ij}\big\|_{1}$
  (induced 1-norm of the 2×2 matrix of norms, eq. (normAgeneral_syst)),
  $Z_2 = 2(\|A_{11}\|_{\ell^1} + \|A_{21}\|_{\ell^1})$, $r_\star = \infty$.
  $D\mathcal F_{\text{pa}}(u) - D\mathcal F_{\text{pa}}(\bar u)$ has the single
  nonzero block $2\mathcal M(x - \bar x)$ in position $(1,1)$, so only the first
  column of $A$ enters $Z_2$.

### Example 3: Cahn–Hilliard, two parameters (Section 4.1, Appendix A)

- $F(v, c, \beta) = \Delta v + \beta(v - v^3 - c)$ on cosine series, $K = 20$;
  lifted map $\tilde{\mathcal F}(u) = \mathcal D^2 u + \beta * (u - u^{*3} - c)$ on
  $\ell^{1,\otimes2}(\mathcal X)$ with $(c, \beta) = \theta$ (Section 4.1, Step 1).
- $N_1 = N_2 = 20$; node $(j_1, j_2)$ at
  $s = (-\cos(\pi j_1/20), -\cos(\pi j_2/20))$, $j_1, j_2 = 0, \dots, 20$.
- **Coons' map, enclosed exactly.** The boundary curves of eqs.
  (parabolacurves1)–(parabolacurves2) are polynomials of degree $\le 2$ with
  rational coefficients ($c_{\max} = 1/10$, $g(c) = 79/2 + 150c^2$) and Coons'
  formula (coons) is polynomial in $(s_1, s_2)$ of degree $\le 2$ in each variable.
  The checker evaluates the same Julia functions $\gamma_j$, $\theta$ with the
  coordinate functions $s_1 = T_1(s_1)$, $s_2 = T_1(s_2)$ as interval Chebyshev
  sequences: the result $\bar\theta = (\bar c, \bar\beta)$ is an interval Chebyshev
  polynomial of order $2$ (no interpolation, only outward rounding) which
  contains the exact $\theta$ coefficient-wise. All bounds are computed with
  $\bar\theta$, hence hold for the exact $\theta$. Condition 3.$i$.0 checks that
  $\bar\theta$ agrees with the pointwise interval evaluation of $\theta$ at the
  $21 \times 21$ nodes and that the corner identities (bdy_rel)
  $\gamma_1(-1) = \gamma_4(1)$, $\gamma_1(1) = \gamma_2(-1)$, $\gamma_2(1) = \gamma_3(-1)$,
  $\gamma_3(1) = \gamma_4(-1)$ hold (exact rational identities; the check guards the
  implementation).
  **Coverage of $\Lambda$.** With these identities, Coons' map restricted to the
  four edges of $[-1,1]^2$ is $\gamma_1, \gamma_2, \gamma_3(-\cdot), \gamma_4(-\cdot)$
  (a computation from eq. (coons)), which run along the Jordan curve
  $\partial\Lambda$ (the parabola $\beta = g(c)$, $|c| \le 1/10$, followed by the
  segment $\beta = 41$) exactly once, monotonically in $c$ (for "pinch", two
  edges are constant). Hence $\theta$ is continuous, $\theta(\partial[-1,1]^2) =
  \partial\Lambda$ with winding number $\pm1$ about every interior point, and
  Proposition 3.3 gives $\theta([-1,1]^2) \supseteq \Lambda$. This is a property
  of the exact formulas, not a numerical check.
- Interpolants (`to_coef(reverse(·; dims=(1,2)), Chebyshev(20) ⊗ Chebyshev(20))`
  in binary64, then thin intervals): $\bar u$ of the node data, $A_{\le K}$ of
  $A^{(j)}_{\le K} = (\Pi_{\le K} D_vF(\bar v_j, \operatorname{mid}\theta(s_j))\Pi_{\le K})^{-1}$
  (21×21, binary64). $A(s) = A_{\le K}(s)\Pi_{\le K} + \Delta_0^{-1}\Pi_{> K}$
  ($\Delta_0^{-1}$ = pseudo-inverse of $\Delta$, zero on $k = 0$).
- **Bounds** (Section 4.1, Step 3, and Appendix A; $N_\theta = 2$ is the order of
  $\bar\theta$, $N = 20$):
  - $Y = \|A\,\tilde{\mathcal F}(\bar u)\|$. $A$ has degree $N$ and
    $\bar\beta\bar u^{*3}$ degree $3N + N_\theta$ in each variable, so the product has
    degree $4N + N_\theta = 82$; it is formed on the $129 \times 129$ grid
    ($129 = 2^7 + 1 \ge 83$) and interpolated in Chebyshev$(128)^{\otimes 2}$, exactly.
    (The cosine degree is $3K$, so $\Pi_{3K}$ of the paper is automatic.)
  - $Z_1 = \max\big(\|\Pi_{\le q} - A\, D\tilde{\mathcal F}(\bar u)\,\Pi_{\le q}\|,\
    \|\bar\beta * (e_{0,0} - 3\bar u^{*2})\| / (\pi(K+1))^2\big)$ with $q = 3K$.
    Derivation: $D\tilde F(\bar u) = \Delta + \mathcal M(\bar\beta(1 - 3\bar u^2))$ is
    banded of bandwidth $2K$ in the cosine index ($\bar u^{*2}$ has $2K$ modes),
    so for $q \ge 3K$: $D\tilde F\,\Pi_{>q} = \Pi_{>K}\, D\tilde F\,\Pi_{>q}$ and
    $A\,\Pi_{>K} = \Delta_0^{-1}\Pi_{>K}$, whence
    $\Pi_{>q} - A D\tilde F \Pi_{>q} = -\Delta_0^{-1}\Pi_{>K}\mathcal M(\bar\beta(1-3\bar u^2))\Pi_{>q}$,
    of norm $\le \|\bar\beta(1 - 3\bar u^2)\| / (\pi^2 (K+1)^2)$ since the diagonal
    $\Delta_0^{-1}\Pi_{>K}$ has entries $-1/(\pi^2k^2)$, $k \ge K+1$. The $\ell^1$
    operator norm is the maximum over columns, hence the `max` of the two column
    blocks. The finite block has degree $3N + N_\theta = 62$ and is formed on the
    $65 \times 65$ grid; the tail sequence has degree $2N + N_\theta$ and is formed
    on the same grid.
  - $Z_2 = 3\|\bar\beta\|\,\|\mathcal A\|\,(2\|\bar u\| + r_\star)$, $\|\mathcal A\| =
    \max(\|\mathcal A_{\le K}\|, 1/(\pi(K+1))^2)$, $r_\star = 10\sup Y$:
    $\mathcal A(D\tilde{\mathcal F}(u) - D\tilde{\mathcal F}(\bar u)) = -3\mathcal A\mathcal M(\bar\beta * (u + \bar u) * (u - \bar u))$
    and $\|u + \bar u\| \le 2\|\bar u\| + r_\star$ on the ball. The columns $> K$ of
    $A(s)$ are the constant $\Delta_0^{-1}\Pi_{>K}$.
- Norms: $\|\cdot\|$ on sequences of $\ell^{1,\otimes2}(\mathcal X)$ is `X_CH`
  (three factors), on $\bar\beta$ it is `X_Λ` (two factors).

### Example 4: SKT, pseudo-arclength, fourteen arcs (Section 4.2)

- $F(v,\lambda) = \Delta\Phi(v,\lambda) + R(v)$, $F_{\text{pa}}(u, s) = (F(v,\lambda),\ \langle u - \bar u(s), \dot u(s)\rangle)$
  with $\langle\cdot,\cdot\rangle$ the Euclidean product of the stored coefficient
  vectors; $u = (v_1, v_2, \lambda)$ on cosine series with $K = 80$
  (`evensym(Fourier(80, π)) × evensym(Fourier(80, π)) × ScalarSpace()`, 163 coefficients).
- Per arc: arclength $L$ and Chebyshev order $N$ as in the table of Result 4
  (`SKT.arclength_of`, `SKT.arc_N`); node $j$ at pseudo-arclength
  $\tfrac L2(1 - \cos(\pi j/N))$ from the start of the arc, i.e. at $s_j = \cos(\pi j/N)$
  (the arc starts at $s = 1$ and ends at $s = -1$).
- $K_A = 165$, except $K_A = 250$ for `u_bar3d` (`SKT.arc_K_A`). At each node, in
  binary64: $A^{(j)}_{\le K_A} = \Pi_{\le K_A}(\Pi_{\le 2K_A} D_uF_{\text{pa}}(\bar u_j, s_j)\Pi_{\le 2K_A})^{-1}\Pi_{\le K_A}$
  and $W^{(j)} \approx D_v\Phi(\bar v_j, \bar\lambda_j)^{-1}$ (pointwise 2×2 inverses on a
  $2K_A+1$ point grid, interpolated in cosine series of order $K_A$).
  $\bar u$, $\dot u$ = interval interpolants of the node data (`to_coef` of thin
  intervals, Chebyshev$(N)$; they take the node values at the nodes); $A_{\le K_A}$,
  $W$ = binary64 interpolants, then thin intervals.
  $A(s) = A_{\le K_A}(s)\Pi_{\le K_A} + \mathrm{diag}((W(s) - \Pi_{\le K_A}W(s)\Pi_{\le K_A})\Delta_0^{-1}, 0)$;
  its columns $> K_A$ are $\mathrm{diag}(W(s)\Delta_0^{-1}\Pi_{>K_A}, 0)$.
- **Bounds** (Section 4.2), all norms in $\mathcal S_{\text{pa}}$ (`X_pa`: one
  weight per factor on each component, sum over the components):
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
    $A(s)\,\mathrm{diag}(\Delta,0)\,\Pi_{\le K_A}$ (`AΔ_cheb_trunc`).
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
    $C_2 = \max(2|a_1|, |b_1|+|b_2|, 2|a_2|) = 6$ (see "Notes on the paper").
- **Connection of the arcs.** The arcs of each branch (`1a→1b→1c`,
  `2a→2b→2c→2d`, `3a→3b→3c→3d`, `4a→4b→4c`) are consecutive: at every junction
  the last node of an arc $a$ ($s = -1$) and the first node of the next arc $b$
  ($s = 1$) carry the same point $u^J$ and the same tangent $\dot u^J$, bit for
  bit in the inputs. Since the interval interpolants take exactly the node values
  at the nodes, the maps $F_{\text{pa}}^{(a)}(\cdot, -1)$ and
  $F_{\text{pa}}^{(b)}(\cdot, 1)$ on $\mathcal X_{\text{pa}}$ coincide (same base
  point, same hyperplane). Hence $u^\star_b(1)$ is a zero of that map within $r_b$
  of $u^J$, and $u^\star_a(-1)$ is its unique zero within $r_a^+$ of $u^J$: if
  $r_b \le r_a^+$ (or $r_a \le r_b^+$), $u^\star_b(1) = u^\star_a(-1)$ and the two
  validated curves join continuously (conditions 4.J1–4.J10).

## Inputs (in `data/`)

All inputs are frozen. They are the binary64 outputs of Newton continuations:
the approximate solutions on the Chebyshev nodes and, for the pseudo-arclength
examples, the approximate tangents that define the hyperplanes of
$F_{\text{pa}}$. They were written by `julia --project=search search/produce.jl`
(the `produce` command of `entry.json`), which is not part of the proof: any
other data passing the conditions below would prove the results as well (with
their own $\bar x$). The SKT continuation starts from the initial guesses in
`search/skt_initial_guesses.bin` (14 × 163 binary64, in the order of the table
below); the Cahn–Hilliard initial guess contains `rand()`, seeded with
`Random.seed!(20260922)`.

**Recomputed by the checker, not stored:** the approximate inverses ($A_j$,
$A^{(j)}_{\le K}$, $A^{(j)}_{\le K_A}$, $W^{(j)}$) at the nodes, the enclosure
$\bar\theta$ of Coons' map, and all Chebyshev interpolants. The approximate
inverses are computed in binary64 (for Example 3 at the midpoints of the interval
values $\theta(s_j)$); they only need to be some bounded operators (their
injectivity follows from $Z_1 < 1$), so the conditions do not depend on their
exact values, and the interpolants are part of the definition of $\bar x$,
$\bar u$, $\dot u$.

Files (raw little-endian IEEE binary64, no header; "node $j$" as in the Setting):

| file | content | values |
|---|---|---|
| `cube_root.bin` | $x_j$, $j = 0,\dots,8$ | 9 |
| `square_root.bin` | $(x_j, \lambda_j)$ for $j = 0,\dots,32$, then $(\dot x_j, \dot\lambda_j)$ for $j = 0,\dots,32$ | 132 |
| `cahn_hilliard_no_pinch.bin`, `cahn_hilliard_pinch.bin` | cosine coefficients $v_{j_1,j_2,k}$, $k = 0,\dots,20$ fastest, then $j_1 = 0,\dots,20$, then $j_2 = 0,\dots,20$ | 9261 (74 088 bytes) |
| `skt_<arc>.bin` (14 files) | for $j = 0,\dots,N$: $u_j = (v_{1,0..80}, v_{2,0..80}, \lambda)$ (163 values); then for $j = 0,\dots,N$: $\dot u_j$ (163 values) | $2 \cdot 163 (N+1)$ |

with $N$ = 16, 16, 24, 16, 16, 24, 24, 16, 16, 20, 24, 16, 24, 16 for
`u_bar1a`, `u_bar1b`, `u_bar1c`, `u_bar2a`, `u_bar2b`, `u_bar2c`, `u_bar2d`,
`u_bar3a`, `u_bar3b`, `u_bar3c`, `u_bar3d`, `u_bar4a`, `u_bar4b`, `u_bar4c`
(`SKT.arc_N`), i.e. 44 336, 54 768 or 65 200 bytes per file.
SHA-256 digests are in `data/hashes.json`; the checker verifies them before
reading the files.

## Conditions to verify

Each is a finite computation in outward-rounded interval arithmetic on the data
above. `entry.json` declares the four results as the `cases` `cube_root`,
`square_root`, `cahn_hilliard`, `skt` of the placeholder `example`, with
`cahn_hilliard` as `test_case`; the checker run on one case prints the conditions
of that result, in the order below, each prefixed `pass` or `fail`, with the
enclosure obtained (endpoints at full binary64 precision). The theorem is the
conjunction over the four cases. For each application $k$ of Theorem 1.1:

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

1. Cube root (Result 1, case `cube_root`): conditions 1.1–1.4.
2. Square root (Result 2, case `square_root`): conditions 2.1–2.4.
3. Cahn–Hilliard (Result 3, case `cahn_hilliard`): 3.$i$.0–3.$i$.4 for $i = 1$
   "no pinch", $i = 2$ "pinch", where 3.$i$.0 is: **$\bar\theta$ encloses Coons'
   map**: the corner identities (bdy_rel) hold (the interval differences contain
   $0$) and $\bar\theta$ intersects the pointwise interval evaluation of $\theta$
   at every node.
4. SKT (Result 4, case `skt`): 4.$i$.1–4.$i$.4 for the $i$-th arc, $i = 1,\dots,14$
   (`u_bar1a`, …, `u_bar4c` in the order of the table above), then
   4.J1–4.J10 **connection** for the ten junctions `1a→1b`, `1b→1c`, `2a→2b`,
   `2b→2c`, `2c→2d`, `3a→3b`, `3b→3c`, `3c→3d`, `4a→4b`, `4b→4c`: both arcs
   proved, the stored last node (point and tangent) of the first arc equals the
   stored first node of the second, and $r_b \le r_a^+$ or $r_a \le r_b^+$.

All conditions hold ⇒ Theorem. (Conditions 1.1–1.4 $\Rightarrow$ Result 1;
2.1–2.4 $\Rightarrow$ Result 2; 3.1.0–3.1.4 and 3.2.0–3.2.4 $\Rightarrow$ Result 3;
for each arc, 4.$i$.1–4.$i$.4 $\Rightarrow$ the statement of Result 4 for that arc,
and 4.J$j$ $\Rightarrow$ the two arcs of junction $j$ join.) Cases run with this
checker, all passing: the four examples (84 conditions in all), in one run and
each on its own.

## Notes on the paper

- Sections 4.1–4.2, norms: the numbers printed in the "Results of the code"
  paragraph of Section 4.1 and in Table 1 give weight $1$ to the Chebyshev
  indices of tensor spaces; the conditions use the norm of eq. (ell1_multi)
  (weight $2$ per nonzero index), which eq. (C0vsell1) and Theorem 2.11 need.
- Section 4.1 and Table 1, orders: the paper uses $N_1 = N_2 = 16$ for
  Cahn–Hilliard, $N = 20$ for `u_bar2c`, $N = 16$ for `u_bar3c` and $K_A = 165$ for
  every arc; the inputs here use $N_1 = N_2 = 20$, $N = 24$, $N = 20$ and
  $K_A = 250$ for `u_bar3d`, because at the paper's values the radii polynomial
  of condition $k$.4 has no root in the norm of eq. (ell1_multi).
- Section 3.1.2, $Z_2$: the printed value $4.98384$ does not match the printed
  formula $2(\|A_{11}\| + \|A_{21}\|)$; the condition uses the formula.
- Appendix A: it says $DF(\bar v)$ is banded of bandwidth $K$ and splits at $q = 2K$;
  for the cubic nonlinearity the bandwidth is $2K$ and the split must be at
  $q \ge 3K$ (as Section 4.1 does with $\Pi_{3K}$, $\Pi_{5K}$). The condition uses
  $q = 3K$.
- Section 4.2, $Z_2$: the derivation bounds $\|\mathcal M(D_v^2\Phi)\|$ only, but
  $\Phi$ depends on $\lambda$ ($\Phi_i \ni \lambda v_i$), so $C_1$ must include the
  mixed term: $C_1 = \max(2|d_{11}|, |d_{12}|+|d_{21}|, 2|d_{22}|, 1)$; with the
  paper's parameters $C_1 = 3$ either way.
- Section 4.2, $Z_2$: the printed $C_1 = \max(\max(|2d_{11}|,|d_{12}|)+|d_{21}|, |d_{12}|+\max(|d_{21}|,|2d_{22}|))$
  and $C_2 = \max(\max(|2a_1|,|b_1|)+|b_2|, |b_1|+\max(|b_2|,|2a_2|))$ are valid
  but not sharp (the printed $C_2$ gives $7$); the conditions use the sharp
  values $\max(2|d_{11}|, |d_{12}|+|d_{21}|, 2|d_{22}|, 1) = 3$ and
  $\max(2|a_1|, |b_1|+|b_2|, 2|a_2|) = 6$ (Setting, Example 4).
- Proposition 3.3, first statement: continuity and $\theta(\partial D) = \partial\Lambda$
  alone do not give $\theta(D) \supseteq \Lambda$ (a map of $D$ into $\partial\Lambda$
  whose restriction to $\partial D$ runs along $\partial\Lambda$ and back is a
  counterexample); the proof uses, and the statement needs, that
  $\theta|_{\partial D}$ has nonzero winding number about the interior points,
  which Coons' map of monotone boundary curves satisfies (Setting, Example 3).
- Section 4.1: $A_{\le K}(s)$ is $(K+1)\times(K+1)$, not $K\times K$; in "Results of
  the code", "(with pinch (parabolacurves1) and without (parabolacurves2))" has
  the labels swapped.
