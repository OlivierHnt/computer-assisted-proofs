# Statement and verification conditions

**Theorem (Marchal's conjecture).** The $P_{12}$ family of the three-body
problem connects the Lagrange equilateral triangle ($\Omega = 1$) to the
figure-eight choreography ($\Omega = 0$): there is a smooth branch of periodic
orbits, parametrised by $\Omega \in [0,1]$, with the symmetries of $P_{12}$,
and the orbit at $\Omega = 0$ has the shape of a figure eight.

Reference: R. Calleja, C. García-Azpeitia, O. Hénot, J.-P. Lessard,
J. D. Mireles James, *From the Lagrange Triangle to the Figure Eight
Choreography: Proof of Marchal's Conjecture*, Trans. Amer. Math. Soc. (2026),
doi:10.1090/tran/9779; arXiv:2406.17564. Section and lemma numbers below refer
to it.

## Setting

- Unknown: $x = (a, \beta, \alpha, u, v, w)$, with $a, \beta, \alpha$ scalars,
  $u = (u_x, u_y, u_z)$ and $v = (v_x, v_y, v_z)$ triples of Fourier series in
  $t$, $w$ a Fourier series in $t$, all depending on $\Omega \in [0,1]$ through
  a Chebyshev expansion (Section 3). A Fourier series is
  $\sum_{|k| \le K} a_k e^{ikt}$. The symmetries of $P_{12}$ (Section 2) put
  each series in one of the subspaces
  EvenCos ($a_{-k} = a_k$, $k$ even), EvenSin ($a_{-k} = -a_k$, $k$ even),
  OddCos, OddSin (same with $k$ odd) and EvenFourier ($k$ even, no other
  relation):
  $u \in$ EvenCos $\times$ EvenSin $\times$ OddCos,
  $v \in$ EvenSin $\times$ EvenCos $\times$ OddSin,
  $w \in$ EvenFourier.
- Symmetric representation: a series in one of these subspaces is stored as
  one complex coefficient per pair $\{k, -k\}$, namely $a_k$ with $k \ge 0$,
  in increasing order of $k$ ($k = 0$ is dropped for the sine subspaces);
  EvenFourier stores $a_k$ for $k = -K, -K+2, \dots, K$. Norms, products,
  shifts, derivatives and projections account for the whole pair.
  For $K = 70$: EvenCos has 36 coefficients, EvenSin, OddCos and OddSin 35
  each, EvenFourier 71, so $x$ has $d = 3 + 106 + 106 + 71 = 286$
  coefficients, in the order $a, \beta, \alpha, u_x, u_y, u_z, v_x, v_y, v_z, w$.
- Truncation: Fourier order $K = 70$, Chebyshev order $N = 20$.
- Banach space: $\ell^1$ with geometric weight $\nu = 1.1$ in the Fourier
  index and no weight in the Chebyshev index, summed over the components
  (Section 3).
- Map: $F(x, \Omega)$ of Section 3 and its derivative $DF$. The shifts by
  $2\zeta$ and $4\zeta$, $\zeta = 2\pi/3$, are $2/3$ and $4/3$ of the period
  and enter as exact roots of unity.
- Radius: $r = 10^{-6}$.
- Arithmetic: IEEE binary64 intervals with outward rounding
  (RadiiPolynomial.jl on IntervalArithmetic.jl, correctly rounded elementary
  functions from CRlibm; versions pinned in `checker/Manifest.toml`). The
  Chebyshev transforms (node values to coefficients, coefficients to a finer
  grid) are done in 256-bit BigFloat intervals.

## Inputs (in `data/`)

- Frozen: `x_nodes.bin`, the values $\bar x_j$ of the approximate solution
  $\bar x$ at the $N + 1 = 21$ Chebyshev–Lobatto nodes
  $\Omega_j = \tfrac12 + \tfrac12\cos(\pi j/N)$, $j = 0, \dots, N$
  (from $\Omega = 1$ down to $\Omega = 0$), each as $d = 286$ complex interval
  coefficients in the symmetric representation above. Byte layout: raw
  little-endian IEEE binary64, no header; node by node, coefficient by
  coefficient, four numbers per coefficient
  $(\inf \operatorname{Re}, \sup \operatorname{Re}, \inf \operatorname{Im}, \sup \operatorname{Im})$;
  file size $21 \cdot 286 \cdot 4 \cdot 8$ bytes. The node $\Omega_0 = 1$ holds
  the Lagrange triangle, the node $\Omega_N = 0$ has $u_x = v_x = 0$, and the
  interior nodes hold thin intervals. Hash in `data/hashes.json`; the checker
  verifies it before reading the file.
- Produced: the frozen file is written by
  `julia --project=search search/find_candidate.jl` (Newton continuation in
  $\Omega$ over the nodes). This step is not part of the proof: any file
  passing the conditions below proves the theorem.
- Constructions the checker performs:
  - $\bar x$ is the Chebyshev interpolant of the node values; its coefficients,
    and its values on the finer grid of $7N + 1$ Chebyshev–Lobatto nodes used
    below, are enclosed by a rigorous DFT. The interpolant is made real-valued
    with $\beta = \alpha = 0$, which condition 0c makes legitimate.
  - $A$ (the $A_{\text{finite}}$ of the paper): on each node, $A_j$ is a
    floating-point inverse of $DF(\operatorname{mid} \bar x_j, \Omega_j)$;
    $A$ is the Chebyshev interpolant of the $A_j$, taken as an exact operator.
    The conditions do not depend on its exact value (which may vary with the
    platform), only on the same $A$ being used in all of them.

## Conditions to verify

Each is a finite computation in outward-rounded arithmetic on the inputs above.
The checker prints them in this order, each prefixed `pass` or `fail`, with the
enclosure obtained.

0. **Input enclosures.**
   - 0a. The node $\Omega = 1$ encloses the Lagrange triangle, whose stored
     coefficients are $a_2(u_x) = 3^{-1/6}/2$, $a_2(u_y) = -i\,3^{-1/6}/2$,
     $a_1(u_z) = 1/2$, $a_2(v_x) = i\,3^{-1/6}$, $a_2(v_y) = 3^{-1/6}$,
     $a_1(v_z) = i/2$, $a_0(w) = 3^{-1/3}$, all others zero.
   - 0b. The node $\Omega = 0$ encloses $u_x = v_x = 0$.
   - 0c. Every node is a real-valued function (cosine coefficients real, sine
     coefficients imaginary, conjugate pairs $a_{-k} = \overline{a_k}$ on $w$)
     with $\beta = \alpha = 0$.
1. **$Y$ (Lemma 3.2).**
   $Y := \|A F(\bar x)\|_\nu + \|\text{tail}_K F(\bar x)\|_\nu / (K+1)$
   is a bounded interval. Here $\text{tail}_K$ keeps the Fourier modes
   $|k| > K$; $F(\bar x)$ is computed without truncation.
2. **$Z_1$ (Lemma 3.3).**
   $Z_1 := \max\big(\|A\,DF(\bar x) - I\|_\nu + \max(\|\text{tail}_K(w^3\Delta u_z \Delta v_z)\|_\nu, M_1)/(K+1),\;
   \|A\|_\nu M_2 + M_1/(K+1)\big)$ is a bounded interval, with $M_1, M_2$ the
   explicit expressions of Lemma 3.3 and $\Delta f = f - f(\cdot + 2\zeta)$.
3. **$Z_2$ (Lemma 3.4).**
   $Z_2 := \max(\|A\|_\nu, 1/(K+1)) \cdot \sup_{\|x - \bar x\|_\nu \le r} \|D^2 F(x)\|_\nu$
   is a bounded interval, with the explicit bound on $\|D^2F\|$ of Lemma 3.4.
4. **Radii polynomial (Section 3).** $r = 10^{-6}$ lies in the interval of
   existence $I \subset [0, r]$ on which $Y + (Z_1 - 1)\,\rho + Z_2\,\rho^2/2 \le 0$
   and $Z_1 + Z_2\,\rho < 1$ for all $\rho \in I$. Conclusion: for every
   $\Omega \in [0,1]$ there is a unique zero of $F$ in the ball of radius $r$
   around $\bar x(\Omega)$.
5. **Figure-eight shape (Lemma 2.2 (iv)).** Let
   $\mu(t) = U_y V_z - U_z V_y$, with $U_y = u_y$, $V_y = v_y$,
   $U_z = \sqrt{a}\,u_z$, $V_z = \sqrt{a}\,v_z$ taken from the zero at
   $\Omega = 0$ (each within $r$ of the stored node value in $\ell^1_\nu$):
   - 5a. $\mu(t) < 0$ for $t \in [0, 3/2]$, verified on 2000 subintervals;
   - 5b. $\mu'''(t) > 0$ for $t \in [3/2, \pi/2]$, verified on 100 subintervals.

   The error of each factor is propagated to $\mu$ by the Banach-algebra
   property of $\ell^1_\nu$, and an $\ell^1_\nu$ error $\varepsilon$ on $\mu$
   gives an unweighted $\ell^1$ error on $\mu'''$ of at most
   $\varepsilon \sup_{k \ge 0} k^3 \nu^{-k} \le \varepsilon\,(3/(e \log \nu))^3$,
   which bounds the pointwise error.

All conditions hold $\Rightarrow$ Theorem.
