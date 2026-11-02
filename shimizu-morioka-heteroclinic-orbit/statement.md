# Statement and verification conditions

**Theorem (transverse heteroclinic orbit).** For $a = 3/4$, $b = 9/20$, the
Shimizu–Morioka system
$$\dot x = y, \qquad \dot y = x - a y - x z, \qquad \dot z = -b z + x^2$$
has a transverse heteroclinic orbit $w : \mathbb R \to \mathbb R^3$ with
$w(t) \to c_1^\star$ as $t \to -\infty$ and $w(t) \to c_0^\star = (0,0,0)$ as
$t \to +\infty$, where $c_0^\star$ is a saddle with two real stable eigenvalues
and one real unstable eigenvalue, and $c_1^\star \approx (0.67082, 0, 1)$ is a
saddle-focus with one real stable eigenvalue and two complex conjugate unstable
eigenvalues. Precisely: the boundary-value problem of Section 4 on $[0, \tau]$,
$\tau = 13$, has a solution $w(t) = u^\star(2t/\tau - 1)$ with
$w(0) = P^\star(e^{\mathrm i\alpha^\star}, e^{-\mathrm i\alpha^\star}) \in W^u_{loc}(c_1^\star)$ and
$w(\tau) = Q^\star(\theta^\star) \in W^s_{loc}(c_0^\star)$, where
$x^\star = (u^\star, \alpha^\star, \theta_1^\star, \theta_2^\star)$ is the unique zero of the map
$F$ of Section 4, Step 1, in the ball of radius $r$ (condition 6e) around the approximation
$\bar x$ of `data/x_bar.bin` in $\mathcal X_{C,\mu}$.

Reference: O. Hénot, A. Takayasu, *Computer-Assisted Proofs in Dynamical
Systems: A Case Study of a Heteroclinic Orbit in the Shimizu–Morioka System*,
arXiv:2605.07500. Theorem, proposition and lemma numbers below refer to the
arXiv version (v1), which defines every object named here; the paper states no
numbered main theorem, the result is the conclusion of Section 4. The authors'
revised manuscript (v2, in preparation) changes the form of some statements;
where a condition uses the revised form, this is said in "Notes on the paper"
at the end.

## Setting

The proof chains six validations; each uses the enclosures proved by the
previous ones. Each is an application of the Radii Polynomial Theorem 1.1
(Section 1.1): with $Y \ge \|AF(\bar x)\|$ and
$Z(R) \ge \sup_{B(\bar x, R)} \|I - A DF\|$, $Z < 1$ and $Y/(1-Z) \le R$ give a
unique zero in $B(\bar x, r)$, $r = Y/(1-Z)$. Everywhere $R = 10\,\sup Y$, and
the radius proved is $r = \sup(Y/(1-Z))$, the infimum of the interval
returned by `RadiiPolynomial.interval_of_existence(Y, Z, R)` (first-order radii
polynomial $Y + (Z - 1) r \le 0$, $Z < 1$; $R$ is passed as the number
$10\,\sup Y$).

- **Parameters.** $a = 0.75$, $b = 0.45$ enclosed as `I"0.75"`, `I"0.45"`;
  vector field $f$ and Jacobian $Df$ as in Section 2 (`checker/src/system.jl`).
- **Equilibrium $c_1$ (Section 2.1).** $F = f$ on $\mathbb R^3$ with the 1-norm;
  $A = \mathrm{inv}(Df(\bar c_1))$ in floating point;
  $Y = \|A f(\bar c_1)\|_1$, $Z = \|I - A\,Df(\mathrm{Box}(\bar c_1, R))\|_1$.
  $c_0^\star = 0$ is exact and needs no proof.
- **Eigenpairs (Section 2.2; Note 4).** For $M = Df(c_i^\star)$ (an interval matrix, $c_1^\star$
  the box of condition 1c), each eigenpair $(\lambda_j, v_j)$ separately, as zero of
  $F_{\mathrm{eig}}(v, \lambda) = \big(M v - \lambda v,\ v^{(l_*)} - 1\big)$ on
  $X = \mathbb C^3 \times \mathbb C$ ($\mathbb R^3 \times \mathbb R$ for a real
  eigenpair) with the 1-norm $\|(v, \lambda)\|_X = \sum_i |v^{(i)}| + |\lambda|$, $l_*$
  the index of the largest-modulus entry of $\bar v_j$, which is scaled so that
  $\bar v_j^{(l_*)} = 1$ (the enclosed eigenvector is then scaled back by that entry
  of the approximation). Same procedure as the equilibrium: $A = \mathrm{inv}(\mathrm{mid}\,DF_{\mathrm{eig}}(\bar x))$,
  $Y = \|A F_{\mathrm{eig}}(\bar x)\|_1$, $R = 10\sup Y$,
  $Z = \|I - A\,DF_{\mathrm{eig}}(\mathrm{Box}(\bar x, R))\|_1$ (the 1-norm ball lies in
  the box of half-side $R$ in every real and imaginary part), Theorem 1.1 with
  $r = Y/(1-Z)$, a bound on each entry. For the conjugate pair of $c_1^\star$,
  $(\lambda_2, v_2)$ is validated and $(\lambda_3, v_3)$ is defined as its conjugate
  ($M$ is real, and the data are exactly conjugate, condition 3d).
- **Local manifolds (Section 3).** $P$, $Q$ in $\mathcal X_{T,\nu} = (X_{T,\nu}^{\otimes 2})^3$
  (weighted $\ell^1$, weight $\nu^{k_1+k_2}$, sum over the three components), map
  $F(P) = P - \phi - \mathrm{Diag}(\mathcal L_T) f(P)$ of Section 3, Step 1, $\phi(\theta) = c^\star + v_1\theta_1 + v_2\theta_2$,
  $A = A_{\le K}\Pi_{\le K} + \Pi_{>K}$ of Section 3, Step 3.
  - Stable manifold of $c_0^\star$: $Q$, Taylor order $K = 15$ in each variable,
    $\nu = 1.1$, eigenpairs $(\lambda_1, \lambda_2)$ of $c_0^\star$, eigenvectors scaled by $0.5$.
  - Unstable manifold of $c_1^\star$: $P$, Taylor order $K = 22$, $\nu = 1.1$,
    eigenpairs $(\lambda_2, \lambda_3)$ of $c_1^\star$ (the complex pair), eigenvectors scaled by $0.15$.
    The centre $\bar P$ validated is the **symmetrised** data: for every component and
    $m \le n$, $\bar P_{(m,n)} \leftarrow (\bar P_{(m,n)} + \overline{\bar P_{(n,m)}})/2$
    (real part on the diagonal $m = n$) and $\bar P_{(n,m)} \leftarrow \overline{\bar P_{(m,n)}}$,
    so that $\bar P_{(m,n)} = \overline{\bar P_{(n,m)}}$ holds exactly in floating point
    (condition 5e; see "Real-valuedness of $P^\star(\gamma(\alpha))$" below; Note 5).
- **Connecting orbit (Section 4).** $x = (u, \alpha, \theta_1, \theta_2) \in \mathcal X_{C,\mu} = X_{C,\mu}^3 \times \mathbb R^3$,
  $\|u\|_{X_{C,\mu}} = |u_0| + 2\sum_{k\ge1}|u_k|\mu^k$, $\mu = 1.05$, Chebyshev order $K = 125$,
  $\tau = 13$, map $F$ of Section 4, Step 1, $A = A_{\le K}\Pi_{\le K} + \Pi_{>K}$ of Section 4, Step 3.
  $P^\star$, $Q^\star$ and their derivatives are evaluated with the error bounds of Lemma 4.3
  (Notes 2–3): for $\theta \in \mathbb D_\nu^2$, real or complex,
  $P^\star(\theta) \in \bar P(\theta) + [-r, r]$ and
  $\partial_{\theta_j} P^\star(\theta) \in \partial_{\theta_j}\bar P(\theta) + [-r_j, r_j]$,
  $r_j = r\,\beta_\nu(|\theta_j|)$ with $\beta_\nu(\delta) = 1/\nu$ for $\delta \le \nu/e$ and
  $1/(e\delta|\ln(\delta/\nu)|)$ for $\nu/e < \delta < \nu$ (a bound on
  $\sup_{k\ge1} k\delta^{k-1}\nu^{-k}$, attained at $k = 1$ when $\delta \le \nu/e$), the
  radius being added to the real and the imaginary part of the interval evaluation
  of $\bar P$ when $\theta$ is complex. At $\sigma = \gamma(\alpha) = (e^{\mathrm i\alpha}, e^{-\mathrm i\alpha})$,
  $|\sigma_j| = 1 > \nu_P/e$ and $r_j = r_P/(e\ln\nu_P)$, so $DP^\star(\gamma(\alpha))\gamma'(\alpha)$
  carries the radius $2r_P\beta_{\nu_P}(1)$; at $\theta_j = \bar\theta_j \pm R$ (an
  interval, condition 6c) $\beta_\nu$ is bounded by its value at the supremum of the
  interval ($\beta_\nu$ is constant on $[0, \nu/e]$ and increasing on $[\nu/e, \nu)$). The enclosures of $P^\star(\gamma(\alpha))$
  and $DP^\star(\gamma(\alpha))\gamma'(\alpha)$ are complex boxes whose real part is
  taken, which is justified below ("Real-valuedness of $P^\star(\gamma(\alpha))$"); the imaginary
  part is used nowhere.
- **Arithmetic.** Outward-rounded binary64 intervals (RadiiPolynomial.jl on IntervalArithmetic.jl, versions pinned in `checker/Manifest.toml`). The approximate
  inverses $A$ are plain binary64 (LAPACK), then used as thin intervals.

## Inputs (in `data/`)

All inputs are frozen: they are the floating-point outputs of Newton's method,
of the LAPACK eigensolver and of an ODE solve, computed by
`julia --project=search search/produce.jl` (the `produce` command of
`entry.json`). The producer is not part of the proof: any other data passing the
conditions below would prove the theorem as well.

Files: raw little-endian IEEE binary64, no header; complex numbers as (re, im) pairs.

- `c1_bar.bin` (3 numbers): $\bar c_1 = (x, y, z)$, approximate equilibrium.
- `eig_c0.bin` (12 numbers): approximate eigenvalues $\bar\lambda_1, \bar\lambda_2, \bar\lambda_3$ of
  $\mathrm{mid}\,Df(c_0^\star)$, then the $3\times3$ matrix of approximate eigenvectors,
  column-major (column $j$ belongs to $\bar\lambda_j$), real.
- `eig_c1.bin` (12 complex numbers = 24 binary64): the same for $\mathrm{mid}\,Df(c_1^\star)$, complex.
- `Q_bar.bin` (768 numbers): $\bar Q$, real; components $i = 1, 2, 3$ one after the
  other, each the $16^2$ Taylor coefficients $\bar Q^{(i)}_{(k_1,k_2)}$, $0 \le k_1, k_2 \le 15$,
  at position $k_1 + 16 k_2$ ($k_1$ fastest).
- `P_bar.bin` (1587 complex numbers): $\bar P$, the same layout with $23^2$ coefficients per component.
- `x_bar.bin` (381 numbers): $\bar u^{(1)}_0, \dots, \bar u^{(1)}_{125}$, then $\bar u^{(2)}$, $\bar u^{(3)}$
  (Chebyshev coefficients, $u(s) = u_0 + 2\sum_{k\ge1} u_k T_k(s)$), then $\bar\alpha, \bar\theta_1, \bar\theta_2$.

SHA-256 digests are in `data/hashes.json`; the checker verifies them before
reading the files.

Recomputed by the checker (not stored): the approximate inverses $A$ of every
subproblem, in floating point from the frozen data (`inv` of a $3\times3$ matrix, of
$DF_{\mathrm{eig}}$ ($4\times4$, one per validated eigenpair), and of the truncations
$\Pi_{\le K}DF\Pi_{\le K}$ of sizes 768, 1587 (complex) and 381). The theorem
does not depend on their exact value: any $A$ passing the conditions proves it.

## Conditions to verify

Each is a finite computation in outward-rounded interval arithmetic on the
inputs above; the checker prints them in this order with the enclosure obtained.
"Finite bound" means: the enclosure is a non-empty bounded interval.

1. **Equilibrium $c_1^\star$ (Section 2.1, Theorem 1.1).**
   - 1a. $Y = \|A f(\bar c_1)\|_1$ is a finite bound.
   - 1b. $Z = \|I - A\,Df(\mathrm{Box}(\bar c_1, R))\|_1$ is a finite bound.
   - 1c. Theorem 1.1: $Z < 1$ and $Y/(1-Z) \le R$. Then $c_1^\star \in \mathrm{Box}(\bar c_1, r)$.
2. **Eigenpairs of $Df(c_0^\star)$ (Section 2.2)**, each validated eigenpair separately.
   - 2a. $Y = \|A F_{\mathrm{eig}}(\bar x)\|_1$ is a finite bound.
   - 2b. $Z = \|I - A\,DF_{\mathrm{eig}}(\mathrm{Box}(\bar x, R))\|_1$ is a finite bound.
   - 2c. Theorem 1.1: $Z < 1$ and $Y/(1-Z) \le R$; gives the radius $r$ (1-norm on
     $(v, \lambda)$, hence on $\bar\lambda_j$ and on every entry of $\bar v_j$).
   - 2d. $c_0^\star$ is a saddle: the enclosures are real with $\lambda_1, \lambda_2 < 0 < \lambda_3$
     (also the hypothesis $\mathrm{Re}\lambda_1\,\mathrm{Re}\lambda_2 > 0$ of Proposition 3.3 for $Q$).
3. **Eigenpairs of $Df(c_1^\star)$ (Section 2.2)**, on the interval matrix $Df(\mathrm{Box}(\bar c_1, r))$;
   $(\lambda_1, v_1)$ and $(\lambda_2, v_2)$ are validated, $(\lambda_3, v_3)$ is defined as the conjugate of $(\lambda_2, v_2)$.
   - 3a–3c. As 2a–2c.
   - 3d. $c_1^\star$ is a saddle-focus: $\lambda_1 < 0$ real, $\lambda_3 = \overline{\lambda_2}$ (as enclosures),
     $\mathrm{Re}\lambda_2 > 0$, $0 \notin \mathrm{Im}\lambda_2$ (also the hypothesis of Proposition 3.3 for $P$);
     moreover the eigenpairs are conjugate: the approximate data satisfy exactly (bit for bit)
     $\mathrm{Im}\,\bar\lambda_1 = 0$, $\bar v_1 \in \mathbb R^3$, $\bar\lambda_3 = \overline{\bar\lambda_2}$,
     $\bar v_3 = \overline{\bar v_2}$, and the validated enclosures satisfy $v_1$ real and $v_3 = \overline{v_2}$
     (real parts equal, imaginary parts opposite, as intervals).
4. **Local stable manifold $Q^\star$ of $c_0^\star$ (Section 3, Proposition 3.3, Theorem 1.1).**
   - 4a. $Y = \|AF(\bar Q)\|_{\mathcal X_{T,\nu}}$ (Section 3, Step 4; a finite sum) is a finite bound.
   - 4b. $Z_0 = \max\big(\|\Pi_{\le K} - A\,DF(\bar Q)\Pi_{\le K}\|,\ \|Df(\bar Q)\| / ((K+1)\min_i|\mathrm{Re}\lambda_i|)\big)$ of Proposition 3.3 is a finite bound.
   - 4c. $Z_1(R) = 2R\max\big(\|A\,\mathrm{Diag}(\mathcal L_T)\Pi_{\le K}\|,\ 1/((K+1)\min_i|\mathrm{Re}\lambda_i|)\big)$ of Proposition 3.3 is a finite bound.
   - 4d. Theorem 1.1 with $Z = Z_0 + Z_1(R)$: $Z < 1$, $Y/(1-Z) \le R$; gives $r_Q$, $\|Q^\star - \bar Q\|_{\mathcal X_{T,\nu}} \le r_Q$.
5. **Local unstable manifold $P^\star$ of $c_1^\star$** (same section and statements, on the symmetrised $\bar P$): 5a–5d, giving $r_P$.
   - 5e. $\bar P$ is exactly conjugacy-symmetric, $\bar P^{(i)}_{(m,n)} = \overline{\bar P^{(i)}_{(n,m)}}$ for all
     $i, m, n$ (bitwise), and the data of `P_bar.bin` deviate from it by at most $10^{-15}$
     in every coefficient (the deviation found is reported).
6. **Connecting orbit (Section 4, Lemma 4.3, Proposition 4.4, Theorem 1.1).**
   - 6a. The intervals $\bar\theta_j \pm R$, $j = 1, 2$, lie in $(-\nu_Q, \nu_Q)$: Lemma 4.3 is
     applied to $Q^\star$ on the ball (condition 6c) and holds on $\mathbb D_\nu^2$ (Note 6).
   - 6b. $Y = \|AF(\bar x)\|_{\mathcal X_{C,\mu}}$ (Section 4, Step 4), with $P^\star$, $Q^\star$ enclosed via
     $\bar P, r_P$ and $\bar Q, r_Q$, is a finite bound.
   - 6c. $Z_0(R) = \max\big(\sup_{x \in B(\bar x, R)}\|\Pi_{\le 2K+1} - \Pi_{\le 3K+2} A \Pi_{\le 3K+2} B(\alpha, \theta)\Pi_{\le 2K+1}\|,\ M\frac{\tau}{2}\|Df(\bar u)\| + \|A\|\mu^{-(2K+2)}\big)$
     of Proposition 4.4 (Note 1; $B(\alpha, \theta)$ is written $\mathcal J(\alpha, \theta)$ there),
     with $M$ and the operator-norm bounds stated there, is a finite bound. The
     term $\|A\|\mu^{-(2K+2)}$ bounds, on the columns $> 2K+1$, the evaluation functionals
     $\mathcal E_i$ through $A$ ($\|\mathcal E\Pi_{>2K+1}\| = \mu^{-(2K+2)}$); it acts on
     the same columns as the $M$ term, so the two add. $B(\alpha, \theta)$ is $DF(\bar x)$ with its $\alpha$ column
     $-DP^\star(\gamma(\alpha))\gamma'(\alpha)$ and its $\theta$ columns $-\partial_{\theta_j}Q^\star(\theta)$ evaluated at
     the $(\alpha, \theta)$ of the ball (the $u$ block stays at $\bar u$; its variation is the $Z_1$ term).
     Since $\|x - \bar x\| \ge |\alpha - \bar\alpha| + |\theta_1 - \bar\theta_1| + |\theta_2 - \bar\theta_2|$, the supremum
     is enclosed by evaluating these columns, with Lemma 4.3, at the intervals
     $\alpha = \bar\alpha \pm R$, $\theta_j = \bar\theta_j \pm R$. $B(\alpha,\theta)\Pi_{\le 2K+1}$ has range in
     $\Pi_{\le 3K+2}$ ($f$ quadratic, $\bar u$ of order $K$), so the norm is that of a finite matrix.
   - 6d. $Z_1(R) = 2R\frac{\tau}{2}\max(\|A_{\le K}\|, 1)(1 + \mu)$ of Proposition 4.4 is a finite bound.
   - 6e. Theorem 1.1 with $Z = Z_0(R) + Z_1(R)$: $Z < 1$, $Y/(1-Z) \le R$; gives $r$,
     $\|x^\star - \bar x\|_{\mathcal X_{C,\mu}} \le r$. Transversality: Remark 4.2.

### Real-valuedness of $P^\star(\gamma(\alpha))$

The map $F$ of Section 4 uses $P^\star(\gamma(\alpha)) \in \mathbb R^3$ and
$DP^\star(\gamma(\alpha))\gamma'(\alpha) \in \mathbb R^3$, and the checker takes the real
part of their complex enclosures. This is rigorous because $P^\star$ satisfies
$P^{\star(i)}_{(m,n)} = \overline{P^{\star(i)}_{(n,m)}}$, which follows from the
uniqueness part of condition 5d:

Let $S$ be the swap-and-conjugate map $(SP)^{(i)}_{(m,n)} = \overline{P^{(i)}_{(n,m)}}$ on
$\mathcal X_{T,\nu}$, an isometry (the weight $\nu^{m+n}$ is symmetric) and an
involution. The map $F$ of Section 3 satisfies $F(SP) = SF(P)$: $f$ has real
coefficients, so $f(SP) = Sf(P)$; $\phi = c_1^\star + v_2\theta_1 + v_3\theta_2$ with
$c_1^\star \in \mathbb R^3$ and $v_3 = \overline{v_2}$, so $S\phi = \phi$;
$(\mathcal L_T)_{(m,n)} = (m\lambda_2 + n\lambda_3)^{-1}$ with $\lambda_3 = \overline{\lambda_2}$, so
$\mathcal L_T$ commutes with $S$. Here $(\lambda_3, v_3) = (\overline{\lambda_2}, \overline{v_2})$ are the
exact eigenpairs used in $F$: $(\lambda_3, v_3)$ is defined as the conjugate of the
validated $(\lambda_2, v_2)$, which is the zero of $F_{\mathrm{eig}}$ for the conjugate
data (exactly conjugate, condition 3d) in the conjugate ball, since $Df(c_1^\star)$ is
real. By condition 5e, $S\bar P = \bar P$ exactly, hence
$\|SP^\star - \bar P\| = \|S(P^\star - \bar P)\| \le r_P$: $SP^\star$ is a zero of $F$ in the
ball of condition 5d, and by uniqueness $SP^\star = P^\star$. Then
$P^\star(e^{\mathrm i\alpha}, e^{-\mathrm i\alpha}) = \sum_{m,n} P^\star_{(m,n)} e^{\mathrm i(m-n)\alpha}$
is invariant under complex conjugation, hence real, and so is its derivative in
$\alpha$. The complex boxes of Lemma 4.3 contain these real numbers, so their real
parts do too.

All conditions hold $\Rightarrow$ Theorem.

## Notes on the paper

1. Proposition 4.4 in the form of condition 6c, $Z_0 = Z_0(R)$ with
   $\mathcal J(\alpha, \theta)$ over the ball and the tail term $\|A\|\mu^{-(2K+2)}$ (the
   revised manuscript; v1 prints $Z_0 = \max\big(\|\Pi_{\le 2K+1} - \Pi_{\le 3K+2} A \Pi_{\le 3K+2} DF(\bar x)\Pi_{\le 2K+1}\|,\ M\frac{\tau}{2}\|Df(\bar u)\|\big)$,
   without the dependence on $R$ and without the tail term).
2. Lemma 4.3 in the form $r_j = r\,\beta_\nu(|\theta_j|)$ with $\beta_\nu$ as in the Setting
   (the revised manuscript; v1 prints $r_j = r/(e|\theta_j||\ln(|\theta_j|/\nu)|)$ for every
   $|\theta_j| < \nu$, which is not a bound on $k|\theta_j|^{k-1}\nu^{-k}$ at $k = 1$ when $|\theta_j| \le \nu/e$).
3. Lemma 4.3 is applied at the complex point $\gamma(\alpha) \in \mathbb D_{\nu_P}^2$, the
   radius added to the real and the imaginary part (the revised manuscript writes
   this case out explicitly).
4. Section 2.2 in the form of the Setting: each eigenpair as a zero of $F_{\mathrm{eig}}$
   on $\mathbb C^3 \times \mathbb C$ with the 1-norm, a conjugate pair validated once (the
   revised manuscript writes this problem out explicitly).
5. The symmetry $P^\star(\theta_2^*, \theta_1^*) = P^\star(\theta_1, \theta_2)^*$ and the
   symmetrised $\bar P$ of condition 5e are Remark (real_image) of the revised manuscript.
6. Proposition 4.4 does not state the hypothesis $\bar\theta_j \pm R \subset (-\nu_Q, \nu_Q)$
   under which Lemma 4.3 applies on the ball; condition 6a verifies it.
