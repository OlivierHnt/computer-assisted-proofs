# Marchal's conjecture — proof package

The $P_{12}$ family of the three-body problem connects the Lagrange triangle to
the figure-eight choreography. This package records the computer-assisted proof
as statement + inputs + checker; see `statement.md` for the theorem and the
conditions.

Check the proof:

```
julia --project=checker checker/verify.jl
```

Regenerate the inputs (optional, not part of the proof):

```
julia --project=search search/find_candidate.jl
julia --project=checker tools/hash_data.jl
```

R. Calleja, C. García-Azpeitia, O. Hénot, J.-P. Lessard, J. D. Mireles James,
*From the Lagrange triangle to the figure eight choreography: Proof of Marchal's
conjecture*, Trans. Amer. Math. Soc. (2026), [doi:10.1090/tran/9779](https://doi.org/10.1090/tran/9779);
[arXiv:2406.17564](https://arxiv.org/abs/2406.17564).
Original code: [OlivierHnt/MarchalConjecture.jl](https://github.com/OlivierHnt/MarchalConjecture.jl).
