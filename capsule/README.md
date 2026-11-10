# Capsule

A framework for computer-assisted proofs of zeros, on top of [RadiiPolynomial.jl](https://github.com/OlivierHnt/RadiiPolynomial.jl).
A problem is data, the approximate inverse is data, and `validate` is generic:

```julia
cert = validate(problem, x_bar, NewtonKantorovich{1}(inverse))
```

It returns a `Certificate` with one radius per component of the unknown.

## Layout

```
src/
  Capsule.jl               module, exports
  contract.jl              what a problem and an inverse must give; validate and the two bounds, written once
  certificate.jl           Certificate, return codes
  radiipolynomial.jl       additions to RadiiPolynomial.jl, to move there
  finite.jl                numbers: FiniteDimensionalProblem, FiniteInverse
  sequence/                blocks of sequences, with F(u) = L u + B N(u)
    projected_inverse.jl     A = A_K Π + T
    bounds.jl                the cut and the tail of the Z bound
    periodic_orbit.jl, initial_value.jl, steady_state.jl
    quasilinear.jl           F(u) = L Φ(u) + N(u), with the tail M(W) L⁻¹
    pseudo_arclength.jl      a problem with its parameter as an unknown, on a hyperplane
  continuation/            any pointwise problem along a path of parameters
    series.jl                what a block becomes: a series in the continuation variable
    series_operator.jl       what an operator with an error becomes
    continuation.jl          ContinuationProblem: the pointwise F, DF and bounds, run on series
    inverse.jl               ContinuationInverse: a pointwise inverse whose entries are series
test/
  runtests.jl              runs the scripts next to it, each a complete proof, and checks their radii
```

The continuation tier has no bound of its own. A problem along a path is of the kind of its pointwise problem,
and the bound of that kind runs on `Series` blocks and `SeriesOperator` entries.

## Requirements

The unreleased development version of RadiiPolynomial.jl. `Project.toml` points to a checkout next to this
repository, at `../../RadiiPolynomial`.

## Tests

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

Three scripts read data of other packages of this repository: `play.jl` (Shimizu-Morioka), and
`skt_steady_state.jl` and `skt_pseudo_arclength.jl` (Chebyshev continuation). The suite takes about a quarter of an hour.
