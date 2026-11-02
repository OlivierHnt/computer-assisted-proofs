module SKT

# Pseudo-arclength continuation of steady states of the Shigesada–Kawasaki–Teramoto
# system (paper, Section 4.2), bifurcation parameter λ = d₁ = d₂, fourteen arcs.
# Adapted from skt.jl of OlivierHnt/ChebyshevContinuation.jl. The Newton
# continuation (Step 2) is removed: its outputs u_grid (points) and u̇_grid
# (tangents) are inputs. The approximate inverses (A_finite, W) are recomputed
# here in floating point from those inputs, exactly as in the original script.
# The bounds are the ones derived in statement.md from Theorem 2.11 and Section
# 4.2; they differ from the original script in two places, marked "revision":
#   (a) the ℓ¹ norms weigh every nonzero Chebyshev index by 2 (sum over n ∈ Z);
#   (d) Z₂ uses the operator norm of 𝒜 Π_{≤K_A} diag(Δ, 0), as the paper's bound
#       requires, instead of that of 𝒜 Π_{≤K_A}.

using RadiiPolynomial, LinearAlgebra

export skt_certificate, ARC_NAMES, arc_N, arc_K_A, K, X_u, X_pa, X_cheb_cos,
    F, DvF, F_pa, DF_pa, tangent, arclength_of, direction_of





# Step 1: Defining the problem

const d₁₁, d₁₂, d₂₁, d₂₂ = exact(0), exact(3), exact(0), exact(0)
const r₁, a₁, b₁ = exact(5), exact(3), exact(1)
const r₂, b₂, a₂ = exact(2), exact(1), exact(3)

Φ(v, λ) = [(λ + d₁₁*v[1] + d₁₂*v[2]) * v[1],
           (λ + d₂₁*v[1] + d₂₂*v[2]) * v[2]]

DvΦ(v, λ) = [exact(2)*d₁₁*v[1]+d₁₂*v[2]+λ                      d₁₂*v[1]
                               d₂₁*v[2]      exact(2)*d₂₂*v[2]+d₂₁*v[1]+λ]

DλΦ(v) = v

R(v) = [(r₁ - a₁*v[1] - b₁*v[2]) * v[1],
        (r₂ - b₂*v[1] - a₂*v[2]) * v[2]]

DR(v) = [r₁-exact(2)*a₁*v[1]-b₁*v[2]                       -b₁*v[1]
                            -b₂*v[2]    r₂-exact(2)*a₂*v[2]-b₂*v[1]]

F(v, λ)   = Laplacian() .* Φ(v, λ) + R(v)

DvF(v, λ) = Laplacian() .* Multiplication.(DvΦ(v, λ)) + Multiplication.(DR(v))

DλF(v)    = LinearOperator.(Laplacian() .* DλΦ(v))

F_pa(u, u̇, w) = [F(u[1:2], u[3][1]) ; adjoint(u̇) * (u - w)]

DF_pa(u, u̇) = [DvF(u[1:2], u[3][1])   DλF(u[1:2])
               adjoint(u̇[1:2])        LinearOperator(u̇[3])]

#--

function tangent(u)
    ub = unpack(u)
    Πseq2 = Projection(space(u)[1:2])
    Π = Projection(space(u))
    D = Πseq2 * [DvF(ub[1:2], ub[3][1]) DλF(ub[1:2])] * Π
    return Sequence(space(u), vec(nullspace(real.(coefficients(D)))))
end

function approx_A(u, u̇, Π_K_A, Π_2K_A)
    DF_ = DF_pa(unpack(u), unpack(u̇))
    return real(Π_K_A * inv(Π_2K_A * DF_ * Π_2K_A) * Π_K_A)
end

function approx_inverse_DvΦ(v, K_A)
    vb = unpack(v)
    U = DvΦ(vb[1:2], vb[3][1])
    m = 2K_A + 1
    U_grid = [real.(to_grid(U[i,j], m)) for i ∈ 1:2, j ∈ 1:2]
    W_pts = [inv([U_grid[i,j][l] for i ∈ 1:2, j ∈ 1:2]) for l ∈ 1:m]
    return [real(to_coef([W[i,j] for W ∈ W_pts], evensym(Fourier(K_A, π)))) for i ∈ 1:2, j ∈ 1:2]
end

#--

const ARC_NAMES = ["u_bar1a", "u_bar1b", "u_bar1c",
                   "u_bar2a", "u_bar2b", "u_bar2c", "u_bar2d",
                   "u_bar3a", "u_bar3b", "u_bar3c", "u_bar3d",
                   "u_bar4a", "u_bar4b", "u_bar4c"]

arclength_of(name) = Dict("u_bar1a" => 0.38, "u_bar1b" => 0.19, "u_bar1c" => 0.26,
                          "u_bar2a" => 0.40, "u_bar2b" => 0.10, "u_bar2c" => 0.10, "u_bar2d" => 0.32,
                          "u_bar3a" => 0.20, "u_bar3b" => 0.10, "u_bar3c" => 0.10, "u_bar3d" => 0.52,
                          "u_bar4a" => 0.20, "u_bar4b" => 0.20, "u_bar4c" => 0.40)[name]

# Chebyshev order per arc. The paper (Table 1) uses 20 and 16 for u_bar2c and u_bar3c;
# with the paper's norm these two arcs fail at those orders, so their inputs were
# regenerated at 24 and 20 (see statement.md, Result 4). u_bar3d fails at the paper's
# order 24 and still at 28 and 36 (Z₁ is dominated by its tail term, independent of
# N), so its input is the original one.
arc_N(name) = Dict("u_bar1a" => 16, "u_bar1b" => 16, "u_bar1c" => 24,
                   "u_bar2a" => 16, "u_bar2b" => 16, "u_bar2c" => 24, "u_bar2d" => 24,
                   "u_bar3a" => 16, "u_bar3b" => 16, "u_bar3c" => 20, "u_bar3d" => 24,
                   "u_bar4a" => 16, "u_bar4b" => 24, "u_bar4c" => 16)[name]

direction_of(name) = Dict("u_bar1a" => -1.0, "u_bar1b" => -1.0, "u_bar1c" => -1.0,
                          "u_bar2a" => -1.0, "u_bar2b" => -1.0, "u_bar2c" => -1.0, "u_bar2d" => 1.0,
                          "u_bar3a" => -1.0, "u_bar3b" => -1.0, "u_bar3c" =>  1.0, "u_bar3d" => 1.0,
                          "u_bar4a" => -1.0, "u_bar4b" => -1.0, "u_bar4c" =>  1.0)[name]

const K = 80

# Truncation order of the approximate inverse, chosen independently from K: 165 as in
# the paper, except for u_bar3d, whose Z₁ tail term ‖𝓜(W)‖‖𝓜(DR)‖/(π(K_A+1))² exceeds
# 0.97 at K_A = 165 whatever the Chebyshev order N (revision: K_A = 250 for that arc).
arc_K_A(name) = name == "u_bar3d" ? 250 : 165

X_u() = evensym(Fourier(K, π)) × evensym(Fourier(K, π)) × ScalarSpace() # space of u = (v₁, v₂, λ)

# ℓ¹ norms (paper, Section 4.2) with ω ≡ 1 and χ ≡ 1: 𝒮_pa = ℓ¹_ω(ℓ¹_χ(R)) × ℓ¹_ω(ℓ¹_χ(R)) × ℓ¹_ω(R)
# with the sum of the component norms; the sums over n ∈ Z and k ∈ Z give weight 2 to
# every nonzero Chebyshev or cosine index. (revision (a): one `IdentityWeight()` per
# factor; the single-weight `Ell1()` of RadiiPolynomial 0.11 gives weight 1 to the
# Chebyshev indices of a tensor space.)

const X_cheb_cos = Ell1(IdentityWeight(), IdentityWeight()) # ℓ¹_ω(ℓ¹_χ(R)): (Chebyshev ⊗ cosine)_sym
const X_cheb     = Ell1(IdentityWeight())                   # ℓ¹_ω(R): Chebyshev
const X_pa       = NormedCartesianSpace((X_cheb_cos, X_cheb_cos, X_cheb), Ell1()) # 𝒮_pa

# utilitary

struct PseudoInverseLaplacian <: AbstractDiagonalOperator end
RadiiPolynomial.getcoefficient(::PseudoInverseLaplacian, (codom, i)::Tuple{SymmetricSpace{<:Fourier},Integer}, (dom, j)::Tuple{SymmetricSpace{<:Fourier},Integer}) =
    (i == j) & !(i == j == 0) ? inv(-(frequency(dom) * exact(i))^2) : zero(frequency(dom))





# Step 2 (the part that does not search): approximate inverses at the nodes and
# Chebyshev interpolants. `u_grid[j]`, `u̇_grid[j]` are the point and tangent at
# arclength 0.5 * (1 - cos(π (j-1)/N)) * arclength along the arc.

function skt_certificate(name, u_grid::Vector{Sequence}, u̇_grid::Vector{Sequence})
    N = arc_N(name)
    K_A = arc_K_A(name)
    length(u_grid) == length(u̇_grid) == N+1 || throw(ArgumentError("expected $(N+1) nodes"))

    A_finite_grid = Vector{LinearOperator}(undef, N+1)

    Π_K_A  = Projection(evensym(Fourier(K_A, π)) × evensym(Fourier(K_A, π)) × ScalarSpace())
    Π_2K_A = Projection(evensym(Fourier(2K_A, π)) × evensym(Fourier(2K_A, π)) × ScalarSpace())
    A_finite_grid[1] = approx_A(u_grid[1], u̇_grid[1], Π_K_A, Π_2K_A)

    for j ∈ 2:N+1
        A_finite_grid[j] = approx_A(u_grid[j], u̇_grid[j], Π_K_A, Π_2K_A)
    end

    W_grid = [approx_inverse_DvΦ(u, K_A) for u ∈ u_grid];

    # retrieve the Chebyshev interpolants

    u_cheb = interval(real(to_coef(interval.(u_grid), Chebyshev(N))));

    u̇_cheb = interval(real(to_coef(interval.(u̇_grid), Chebyshev(N))));

    A_finite_cheb = interval(real(to_coef(A_finite_grid, Chebyshev(N))));

    W_cheb = [interval(real(to_coef(getindex.(W_grid, i, j), Chebyshev(N)))) for i ∈ 1:2, j ∈ 1:2];





    # Step 3: Estimating the bounds

    Δ⁻¹ = PseudoInverseLaplacian()



    #- Y bound

    N_Y = 3N;

    u_grid_Y = real.(to_grid(u_cheb, N_Y+1));

    u̇_grid_Y = real.(to_grid(u̇_cheb, N_Y+1));

    A_finite_grid_Y = real.(to_grid(A_finite_cheb, N_Y+1));

    W_grid_Y_ = [real.(to_grid(W_cheb[i,j], N_Y+1)) for i ∈ 1:2, j ∈ 1:2];
    W_grid_Y = [getindex.(W_grid_Y_, j) for j ∈ eachindex(W_grid_Y_[1])];

    Πseq_K_A = interval(Projection(evensym(Fourier(K_A, π))));
    zero_col = mapreduce(_ -> interval(LinearOperator(ScalarSpace(), evensym(Fourier(0, π)), [0;;])), vcat, 1:2);
    zero_row = mapreduce(_ -> interval(LinearOperator(evensym(Fourier(0, π)), ScalarSpace(), [0;;])), hcat, 1:2);
    zero_corner = interval(LinearOperator(ScalarSpace(), ScalarSpace(), [0;;]));
    A_tail_grid_Y = map(W_grid_Y) do W_
        W = Multiplication.(W_)
        [W .* Δ⁻¹ .- Πseq_K_A .* (W .* Δ⁻¹) .* Πseq_K_A    zero_col
         zero_row                                          zero_corner]
    end;

    A_grid_Y = unpack.(A_finite_grid_Y) .+ A_tail_grid_Y;

    #--

    Π_2K = interval(Projection(evensym(Fourier(2K, π)) × evensym(Fourier(2K, π)) × ScalarSpace()));

    Y_grid = A_grid_Y .* (Π_2K .* F_pa.(unpack.(u_grid_Y), unpack.(u̇_grid_Y), unpack.(u_grid_Y)));
    Y = norm(real(to_coef(Y_grid, Chebyshev(N_Y))), X_pa)



    #- Z₁ bound

    N_Z = 2N;

    u_grid_Z = real.(to_grid(u_cheb, N_Z+1));

    u̇_grid_Z = real.(to_grid(u̇_cheb, N_Z+1));

    A_finite_grid_Z = real.(to_grid(A_finite_cheb, N_Z+1));

    W_grid_Z_ = [real.(to_grid(W_cheb[i,j], N_Z+1)) for i ∈ 1:2, j ∈ 1:2];
    W_grid_Z = [getindex.(W_grid_Z_, j) for j ∈ eachindex(W_grid_Z_[1])];

    A_tail_grid_Z = map(W_grid_Z) do W_
        W = Multiplication.(W_)
        [W .* Δ⁻¹ .- Πseq_K_A .* (W .* Δ⁻¹) .* Πseq_K_A    zero_col
         zero_row                                          zero_corner]
    end;

    A_grid_Z = unpack.(A_finite_grid_Z) .+ A_tail_grid_Z;

    #--

    Π_K_K_A = interval(Projection(evensym(Fourier(K+K_A, π)) × evensym(Fourier(K+K_A, π)) × ScalarSpace()));

    Z₁_finite = opnorm(real(to_coef(Π_K_K_A .- A_grid_Z .* DF_pa.(unpack.(u_grid_Z), unpack.(u̇_grid_Z)) .* Π_K_K_A, Chebyshev(N_Z))), X_pa)

    #-- tail columns (> K + K_A): ‖I - 𝓜(W)𝓜(D_vΦ)‖ + ‖𝓜(W)‖ ‖𝓜(DR)‖ / (π(K_A+1))², each
    #   operator norm bounded by the induced 1-norm of the 2×2 matrix of block norms (eq. (normAgeneral_syst))

    M_grid = Ref(exact(I(2))) .- W_grid_Z .* DvΦ.(getindex.(unpack.(u_grid_Z), Ref(1:2)), getindex.(getindex.(unpack.(u_grid_Z), 3), 1));
    M_cheb = [real(to_coef([W[i,j] for W ∈ M_grid], Chebyshev(N_Z))) for i ∈ 1:2, j ∈ 1:2];

    DR_grid = DR.(getindex.(unpack.(u_grid_Z), Ref(1:2)));
    DR_cheb = [real(to_coef([DR[i,j] for DR ∈ DR_grid], Chebyshev(N_Z))) for i ∈ 1:2, j ∈ 1:2];

    Z₁_tail = opnorm(norm.(M_cheb, X_cheb_cos), 1) +
        opnorm(norm.(W_cheb, X_cheb_cos), 1) * opnorm(norm.(DR_cheb, X_cheb_cos), 1) / (interval(π) * interval(K_A + 1))^2

    #--

    Z₁ = max(Z₁_finite, Z₁_tail)



    #- Z₂ bound

    r_star = Inf # the nonlinearities are quadratic

    A_finite_grid_rig = real.(to_grid(A_finite_cheb, N+1));

    W_grid_rig_ = [real.(to_grid(W_cheb[i,j], N+1)) for i ∈ 1:2, j ∈ 1:2];
    W_grid_rig = [getindex.(W_grid_rig_, j) for j ∈ eachindex(W_grid_rig_[1])];

    Δ = Laplacian();

    opnormW = opnorm(norm.(W_cheb, X_cheb_cos), 1) # ‖𝓜(W)‖ ≤ ‖(‖W_ij‖)_ij‖₁, eq. (normAgeneral_syst)

    A_tail_grid_rig = map(W_grid_rig) do W_
        W = Multiplication.(W_)
        [W .* Δ⁻¹ .- Πseq_K_A .* (W .* Δ⁻¹) .* Πseq_K_A    zero_col
         zero_row                                          zero_corner]
    end;

    A_grid_rig = unpack.(A_finite_grid_rig) .+ A_tail_grid_rig;

    # ‖𝒜‖ ≤ max(‖𝒜 Π_{≤K_A}‖, ‖𝓜(W)‖ / (π(K_A+1))²): the columns > K_A of A(s) are W(s)Δ₀⁻¹Π_{>K_A}
    A_cheb_trunc = real(to_coef(A_grid_rig .* interval(Π_K_A), Chebyshev(N)));
    opnorm_A_trunc = opnorm(A_cheb_trunc, X_pa)
    opnormA = max(opnorm_A_trunc, opnormW / (interval(π) * interval(K_A + 1))^2)

    # ‖𝒜 diag(Δ, 0)‖ ≤ max(‖𝒜 Π_{≤K_A} diag(Δ, 0)‖, ‖𝓜(W)‖): the columns > K_A of A(s) diag(Δ, 0)
    # are W(s)Δ₀⁻¹Π_{>K_A}Δ = W(s)Π_{>K_A} (revision (d): the original script used ‖𝒜 Π_{≤K_A}‖ here)
    AΔ_grid_rig = A_grid_rig .*
                Ref([[Δ exact(0)*I ; exact(0)*I Δ]    zero_col
                     zero_row                         zero_corner]);
    AΔ_cheb_trunc = real(to_coef(AΔ_grid_rig .* interval(Π_K_A), Chebyshev(N)));
    opnorm_AΔ_trunc = opnorm(AΔ_cheb_trunc, X_pa)
    opnormAΔ = max(opnorm_AΔ_trunc, opnormW)

    # Z₂ = ‖𝒜 diag(Δ, 0)‖ C₁ + ‖𝒜‖ C₂ with ‖D²_u (ΔΦ, R)‖: for w, h ∈ 𝒮_pa of unit norm, the
    # bilinear forms D²Φ[w, h] and D²R[w, h] (Φ is quadratic in (v, λ) jointly, R in v) are
    # bounded by the largest coefficient of the pairings (w_i, h_j), see statement.md:
    #   C₁ = max(2|d₁₁|, |d₁₂| + |d₂₁|, 2|d₂₂|, 1)   (the 1 is the cross term λ·v of Φ),
    #   C₂ = max(2|a₁|, |b₁| + |b₂|, 2|a₂|).
    C₁ = interval(max(2*d₁₁, d₁₂ + d₂₁, 2d₂₂, 1))
    C₂ = interval(max(2a₁, b₁ + b₂, 2a₂))

    Z₂ = opnormAΔ * C₁ + opnormA * C₂



    #- radii polynomial (Theorem 1.1: Y + Z₁ r + Z₂ r²/2 < r and Z₁ + Z₂ r < 1, the
    #  convention of `interval_of_existence`)

    ie, proved = interval_of_existence(Y, Z₁, Z₂, r_star)

    return (; Y, Z₁, Z₂, r_star, ie, proved, Z₁_finite, Z₁_tail, opnormW, opnorm_A_trunc, opnorm_AΔ_trunc)
end

end
