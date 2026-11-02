# Validated unstable manifolds of DDEs — proof package

Computer-assisted error bounds for the Taylor series of the parameterization of
the unstable manifold of an equilibrium, for three delay differential equations
(a delayed van der Pol system, an SEIR model, the Mackey–Glass equation with
ρ = 3). This package records the computation as statement + inputs + checker;
see `statement.md` for the theorem, the bounds as rederived from the paper, the
conditions, and the points where the checker departs from the authors'
original notebook and from the printed formulas.

The checker takes a single input per case, the approximate unstable eigenvalue
λ̄ (`data/<case>_lambda.bin`, one binary64 number, two for the complex
eigenvalue of case I). From it, in interval arithmetic, it proves the
eigenvalue enclosure (a radii polynomial argument on the characteristic
function), checks the eigenvector formula, solves the homological equations up
to the paper's orders (N = 99, 200, 600), and evaluates the bounds and the
radii polynomial of the paper.

Check the proof (Julia 1.13, RadiiPolynomial 0.11.13; about 8 minutes with two
threads, `julia -t 2`; the solution of the homological equations of case I at
N = 99 takes most of it):

```
julia --project=checker checker/verify.jl
```

Regenerate the input (optional, not part of the proof: a floating-point Newton
iteration on the characteristic function, started from the paper's values):

```
julia --project=search search/produce.jl
julia --project=checker tools/hash_data.jl
```

The authors' notebook code (interval solution of the homological equations,
IntervalArithmetic 0.17) that earlier versions of this package ran to produce
frozen coefficients is no longer included: the checker recomputes the
coefficients itself. The original notebook remains available at the address
below.

O. Hénot, J.-P. Lessard, J. D. Mireles James, *Parameterization of unstable
manifolds for DDEs: formal series solutions and validated error bounds*,
J. Dyn. Differ. Equ. 34 (2022), 1285–1324,
[doi:10.1007/s10884-021-10002-8](https://doi.org/10.1007/s10884-021-10002-8).
Revised version: <https://www.math.mcgill.ca/jplessard/papers/manifoldValidations_revised.pdf>.
Original code: the notebook `examples_manifolds.ipynb` (April 2020) from
<https://www.math.mcgill.ca/jplessard/code-dde-unstable-manifolds.html>.
