# Transverse heteroclinic orbit in the Shimizu–Morioka system — proof package

For $a = 3/4$, $b = 9/20$, the Shimizu–Morioka system has a transverse
heteroclinic orbit from the saddle-focus $c_1 \approx (0.67082, 0, 1)$ to the
saddle at the origin. This package records the computer-assisted proof
(equilibrium, eigenpairs, local invariant manifolds, connecting orbit) as
statement + inputs + checker; see `statement.md` for the theorem and the
conditions.

Check the proof (Julia 1.13):

```
julia --project=checker checker/verify.jl
```

It takes a few minutes (about 1.5 min with 2 BLAS threads). The
checker follows the authors' revised manuscript (version of 28 October 2026) where it differs
from arXiv v1 (the $Z_0(R)$ bound of Proposition 4.4, the error radius of Lemma
4.3 at the complex evaluation point, the conjugacy symmetry of $\bar P$); see the
errata at the end of `statement.md`. Regenerate the inputs (optional, not part
of the proof; the producer runs the original code and does not symmetrise
$\bar P$, which the checker does at load):

```
julia --project=search search/produce.jl
julia --project=checker tools/hash_data.jl
```

O. Hénot, A. Takayasu, *Computer-Assisted Proofs in Dynamical Systems: A Case
Study of a Heteroclinic Orbit in the Shimizu–Morioka System*,
[arXiv:2605.07500](https://arxiv.org/abs/2605.07500).
Original code: ShimizuMoriokaTutorial.jl (O. Hénot).
