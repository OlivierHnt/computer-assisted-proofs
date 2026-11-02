# Newton–Kantorovich bounds of Theorem 2.1: Y, Z₁, Z₂ and the interval of
# existence for the family of periodic orbits (prove_po, Section 2.2) and for
# the Floquet normal form (prove_lin, Appendix A), from the node values of the
# approximations and of the approximate inverses. The Chebyshev variable η is
# the leading factor of the tensor space Chebyshev(N) ⊗ Fourier(K); the
# approximate inverses A, B are arrays of Chebyshev coefficients. The products
# A W₁ are formed on 65 nodes (Chebyshev order 2N) and interpolated back.
# From src/contraction.jl of OlivierHnt/BHWPeriodicCoexistence.jl, commit 6effd7d.

# number of nodes of the discrete Chebyshev transform of order n (33 for n = 30, 65 for n = 60)
_npts(n) = nextpow(2, 2n + 1) ÷ 2 + 1

_vnorm(a::InfiniteSequence) = sequence_norm(a) + total_error(a) # norm of the coefficients plus the tail error

# derivative of a sequence with zero error, kept in the same Banach space
# (the library's `differentiate` would move it to the unweighted ℓ¹)
function _differentiate(a::InfiniteSequence, α)
    iszero(total_error(a)) || throw(ArgumentError("only sequences without error are differentiated here"))
    return InfiniteSequence(differentiate(sequence(a), α), banachspace(a))
end

# operator 1-norm of a matrix of nonnegative entries: largest column sum
function _opnorm1(M::AbstractMatrix)
    nrm = zero(real(eltype(M)))
    for j ∈ axes(M, 2)
        nrmj = zero(real(eltype(M)))
        for i ∈ axes(M, 1)
            nrmj += M[i,j]
        end
        nrm = max(nrm, nrmj)
    end
    return nrm
end

# Chebyshev coefficients of the Fourier mode k of a sequence in Chebyshev(N) ⊗ Fourier(K)
_cheb_slice(a::Sequence, k) = Sequence(Chebyshev(order(a)[1]), a[(:,k)])

# P_ν norms of the Chebyshev series C[:, i, j]
_cheb_norms(C::Array{T,3}, Pν) where {T} = [norm(Sequence(Chebyshev(size(C, 1) - 1), C[:,i,j]), Pν) for i ∈ axes(C, 2), j ∈ axes(C, 3)]

# conjugation symmetry Σ(x̄) = x̄ (Lemma 2.4): ScalarSpace components real,
# Fourier components a_k ← (a_k + conj(a_{-k}))/2, for every Chebyshev coefficient
function _conjugacy_symmetry_cheb!(x)
    for c ∈ eachcomponent(component(x, 1))
        coefficients(c) .= real.(coefficients(c))
    end
    for c ∈ eachcomponent(component(x, 2))
        N, K = order(c)
        A = reshape(coefficients(c), N+1, 2K+1)
        A .= (A .+ conj.(reverse(A; dims = 2))) ./ 2
    end
    return x
end

# `polish!` applied to every Chebyshev series of coefficients
function _polish_cheb!(x)
    for c ∈ eachcomponent(component(x, 1))
        polish!(c)
    end
    for c ∈ eachcomponent(component(x, 2))
        N, K = order(c)
        A = reshape(coefficients(c), N+1, 2K+1)
        for k ∈ 1:2K+1
            s = Sequence(Chebyshev(N), A[:,k])
            polish!(s)
            A[:,k] .= coefficients(s)
        end
    end
    return x
end

#

function construct_approx_inverse(x_po_grid, x_lin_grid, κ₁, κ₂, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂), N)
    npts = length(x_po_grid)
    N_fft = 2(npts - 1)
    _κ⁻¹_grid = [mid(κ₁ + κ₂ - cospi((2j)/N_fft)*(κ₂ - κ₁))/mid(2κ₁*κ₂) for j ∈ npts-1:-1:0];

    A_po_grid = Vector{Matrix{ComplexF64}}(undef, npts);
    A_lin_grid = Vector{Matrix{ComplexF64}}(undef, npts);

    DF_po_cache = zeros(eltype(x_po_grid[1]), space(x_po_grid[1]), space(x_po_grid[1]));
    DF_lin_cache = zeros(eltype(x_lin_grid[1]), space(x_lin_grid[1]), space(x_lin_grid[1]));

    for j ∈ 1:npts
        params = mid.([(m₁-d₁)/γ, (m₂-d₂)/γ, a₁ * _κ⁻¹_grid[j], a₂ * _κ⁻¹_grid[j], a₁ * _κ⁻¹_grid[j] * d₁/(m₁-d₁), a₂ * _κ⁻¹_grid[j] * d₂/(m₂-d₂)])

        # periodic orbit
        approx_tang = differentiate(component(x_po_grid[j], 2))
        DF_po!(DF_po_cache, x_po_grid[j], params, approx_tang)
        A_po_grid[j] = inv(coefficients(DF_po_cache))

        # Floquet normal form
        J = τ∂ᵤf(x_po_grid[j], params)
        DF_lin!(DF_lin_cache, x_lin_grid[j], J)
        A_lin_grid[j] = inv(coefficients(DF_lin_cache))
    end

    return A_po_grid, A_lin_grid
end

#

function construct_W₀_po(x_cheb, approx_tang_cheb, vτf₁, vτf₂, vτf₃, Pν, Wν)
    u₁_cheb = component(component(x_cheb, 2), 1)
    u₂_cheb = component(component(x_cheb, 2), 2)
    u₃_cheb = component(component(x_cheb, 2), 3)

    N, K = order(u₁_cheb)

    #

    ρ₁ = zero(_cheb_slice(u₁_cheb, 0))
    for k ∈ -K:K
        ρ₁ +=
            _cheb_slice(u₁_cheb, k) * conj(_cheb_slice(component(approx_tang_cheb, 1), k)) +
            _cheb_slice(u₂_cheb, k) * conj(_cheb_slice(component(approx_tang_cheb, 2), k)) +
            _cheb_slice(u₃_cheb, k) * conj(_cheb_slice(component(approx_tang_cheb, 3), k))
    end
    ρ₂ = u₁_cheb(nothing, 0) - 1
    ρ₃ = u₂_cheb(nothing, 0) - 1
    s = Chebyshev(N) ⊗ Fourier(K, interval(1))
    W₀_u₁ = project(differentiate(u₁_cheb, (0, 1)) - sequence(vτf₁), s)
    W₀_u₂ = project(differentiate(u₂_cheb, (0, 1)) - sequence(vτf₂), s)
    W₀_u₃ = project(differentiate(u₃_cheb, (0, 1)) - sequence(vτf₃), s)

    # ‖W₀‖ = Σ of the P_ν norms of all its Chebyshev coefficients
    return norm(ρ₁, Pν) + norm(Sequence(Chebyshev(order(ρ₂)[1]), coefficients(ρ₂)), Pν) + norm(Sequence(Chebyshev(order(ρ₃)[1]), coefficients(ρ₃)), Pν) +
        norm(W₀_u₁, Wν) + norm(W₀_u₂, Wν) + norm(W₀_u₃, Wν)
end

# W₁ (Section 2.2.2) at one node η
function construct_W₁_po(approx_tang, f₁, f₂, f₃, τ∂ζ₁f₃, τ∂ζ₂f₃, J)
    K = order(component(approx_tang, 1))

    W₁ = zeros(Complex{Interval{Float64}}, ScalarSpace()^3 × Fourier(2K, interval(1))^3, ScalarSpace()^3 × Fourier(K, interval(1))^3);

    D₂F₁ = component(W₁, 1, 2);

    for j ∈ 1:3, k ∈ -K:K
        component(D₂F₁, 1, j)[1,k] = conj(component(approx_tang, j)[k])
    end

    for k ∈ -K:K
        component(D₂F₁, 2, 1)[1,k] = 1
        component(D₂F₁, 3, 2)[1,k] = 1
    end

    D₁F₂ = component(W₁, 2, 1);
    project!(component(D₁F₂, 1, 1), -f₁);
    project!(component(D₁F₂, 2, 1), -f₂);
    project!(component(D₁F₂, 3, 1), -f₃);
    project!(component(D₁F₂, 3, 2), -τ∂ζ₁f₃);
    project!(component(D₁F₂, 3, 3), -τ∂ζ₂f₃);

    for j ∈ 1:3, i ∈ 1:3
        linop = component(component(W₁, 2, 2), i, j)
        coefficients(linop) .= coefficients(project(Multiplication(-J[i,j]), domain(linop), codomain(linop)))
    end

    linop = component(W₁, 2, 2);
    coefficients(linop) .+= coefficients(project(Derivative(1), domain(linop), codomain(linop)));

    #

    return W₁
end

# Ŵ₁ (Appendix A) at one node η
function construct_W₁_lin(vecC, V, J)
    K_lin = order(V[1,1])

    W₁_lin = zeros(Complex{Interval{Float64}}, ScalarSpace()^9 × Fourier(2K_lin, interval(1))^9, ScalarSpace()^9 × Fourier(K_lin, interval(1))^9);

    k1 = kron(I[1:3,1:3], J);
    k2 = kron(I[1:3,1:3], V);
    transpose_C_lin = [Sequence(Fourier(0, interval(1)), [vecC[j + 3*(i-1)]]) for i ∈ 1:3, j ∈ 1:3]
    k3 = kron(transpose_C_lin, I[1:3,1:3])
    for j ∈ 1:9
        for i ∈ 1:9
            linop = component(component(W₁_lin, 2, 1), i, j)
            project!(linop, k2[i,j])
            linop = component(component(W₁_lin, 2, 2), i, j)
            coefficients(linop) .= coefficients(project(Multiplication(k3[i,j] - k1[i,j]), domain(linop), codomain(linop)))
        end
        linop = component(component(W₁_lin, 1, 2), j, j)
        for k ∈ -K_lin:K_lin
            linop[1,k] += 1
        end
        linop = component(component(W₁_lin, 2, 2), j, j)
        coefficients(linop) .+= coefficients(project(Derivative(1), domain(linop), codomain(linop)))
    end

    #

    return W₁_lin
end

# ‖A W₁ - I‖: A and W₁ at the 65 nodes, product node by node, Chebyshev(2N) interpolation, P_ν norms
function _opnorm_AW₁_I(A_cheb, W₁_grid, N, Pν)
    m = length(W₁_grid)
    A_grid = cheb2grid(A_cheb, m)
    W = first(W₁_grid)
    nrow, ncol = dimension(codomain(W)), dimension(domain(W))
    G = Array{Complex{Interval{Float64}}}(undef, m, nrow, ncol)
    for ℓ ∈ 1:m
        Aℓ = LinearOperator(codomain(W), codomain(W), A_grid[ℓ,:,:])
        G[ℓ,:,:] .= coefficients(project(Aℓ * W₁_grid[ℓ] - I, domain(W), codomain(W), Complex{Interval{Float64}}))
    end
    norms = Matrix{Interval{Float64}}(undef, nrow, ncol)
    for j ∈ 1:ncol, i ∈ 1:nrow
        norms[i,j] = norm(to_coef(G[:,i,j], Chebyshev(2N)), Pν)
    end
    return _opnorm1(norms)
end

# Theorem 2.1 for F (prove_po) and then for G at the validated orbit (prove_lin)

function prove_contraction(x_po_grid_, x_lin_grid_, κ₁, κ₂, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂); K, K_lin, N, ν, R)
    x_po_grid = project.(x_po_grid_, Ref(ScalarSpace()^3 × Fourier(K, 1.0)^3))
    x_lin_grid = project.(x_lin_grid_, Ref(ScalarSpace()^9 × Fourier(K_lin, 1.0)^9))

    for i ∈ 1:length(x_po_grid) # polish Fourier coefficients
        polish!(x_po_grid[i])
        polish!(x_lin_grid[i])
    end

    A_po_grid, A_lin_grid = construct_approx_inverse(x_po_grid, x_lin_grid, κ₁, κ₂, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂), N)

    #

    true_τ, true_ζ₁, true_ζ₂, true_u₁, true_u₂, true_u₃, bounds_po = prove_po(x_po_grid, A_po_grid, κ₁, κ₂, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂), K, N, ν, R)

    true_C, bounds_lin = prove_lin(x_lin_grid, A_lin_grid, true_τ, true_ζ₁, true_ζ₂, true_u₁, true_u₂, true_u₃, κ₁, κ₂, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂), K_lin, N, ν)

    #

    return true_τ, true_ζ₁, true_ζ₂, true_C, bounds_po, bounds_lin
end


function prove_po(x_po_grid, A_po_grid, κ₁, κ₂, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂), K, N, ν, R)
    x_po_cheb_ = grid2cheb(x_po_grid, N)
    _conjugacy_symmetry_cheb!(x_po_cheb_)
    _polish_cheb!(x_po_cheb_) # polish Chebyshev coefficients
    x_cheb = interval(x_po_cheb_)

    approx_tang_cheb = interval(differentiate(component(x_po_cheb_, 2), (0, 1)))

    A_cheb = interval.(grid2cheb(A_po_grid, N))

    # Banach spaces

    Pν = Ell1(GeometricWeight(ν))
    Wν = Ell1(GeometricWeight(ν), IdentityWeight())

    opnorm_A = max(_opnorm1(_cheb_norms(A_cheb, Pν)), inv(interval(K + 1)))

    # data stuff

    τ_ot  = Sequence(Chebyshev(N) ⊗ Fourier(0, interval(1)), coefficients(component(component(x_cheb, 1), 1)))
    ζ₁_ot = Sequence(Chebyshev(N) ⊗ Fourier(0, interval(1)), coefficients(component(component(x_cheb, 1), 2)))
    ζ₂_ot = Sequence(Chebyshev(N) ⊗ Fourier(0, interval(1)), coefficients(component(component(x_cheb, 1), 3)))
    u₁_ot = Sequence(Chebyshev(N) ⊗ Fourier(K, interval(1)), coefficients(component(component(x_cheb, 2), 1)))
    u₂_ot = Sequence(Chebyshev(N) ⊗ Fourier(K, interval(1)), coefficients(component(component(x_cheb, 2), 2)))
    u₃_ot = Sequence(Chebyshev(N) ⊗ Fourier(K, interval(1)), coefficients(component(component(x_cheb, 2), 3)))

    # InfiniteSequence allows to rigorously control the error when inverting a sequence

    vκ⁻¹ = InfiniteSequence(Chebyshev(1) ⊗ Fourier(0, interval(1)), [(inv(κ₁) + inv(κ₂))/2, (inv(κ₂) - inv(κ₁))/4], Wν)
    δ₁, δ₂ = (m₁ - d₁) / γ, (m₂ - d₂) / γ
    vα₁, vα₂ = a₁ * vκ⁻¹, a₂ * vκ⁻¹
    vλ₁, vλ₂ = vα₁ * d₁ / (m₁ - d₁), vα₂ * d₂ / (m₂ - d₂)

    vτ, vζ₁, vζ₂ = InfiniteSequence(τ_ot, Wν), InfiniteSequence(ζ₁_ot, Wν), InfiniteSequence(ζ₂_ot, Wν)
    vu₁, vu₂, vu₃ = InfiniteSequence(u₁_ot, Wν), InfiniteSequence(u₂_ot, Wν), InfiniteSequence(u₃_ot, Wν)

    ###########
    # Y bound #
    ###########

    #---- (ū₃ + αⱼ)⁻¹ and its error bound (Lemma 2.2) in 256-bit intervals, converted outward to binary64
    big_vu₃ = InfiniteSequence(big(u₃_ot), Wν)
    big_vg₁ = inv(big_vu₃ + vα₁)
    big_vg₂ = inv(big_vu₃ + vα₂)
    _f64(z) = complex(interval(Float64, real(z)), interval(Float64, imag(z)))
    vg₁ = InfiniteSequence(Sequence(Chebyshev(N) ⊗ Fourier(K, interval(1)), _f64.(coefficients(big_vg₁))), Wν; total_error = interval(Float64, total_error(big_vg₁)))
    vg₂ = InfiniteSequence(Sequence(Chebyshev(N) ⊗ Fourier(K, interval(1)), _f64.(coefficients(big_vg₂))), Wν; total_error = interval(Float64, total_error(big_vg₂)))
    #----
    vf₁ = δ₁ * (vu₃ - vλ₁) * vg₁ * vu₁
    vf₂ = δ₂ * (vu₃ - vλ₂) * vg₂ * vu₂
    vf₃ = (1 - vu₃ - vζ₁ * vu₁ * vg₁ - vζ₂ * vu₂ * vg₂) * vu₃
    #----
    vτf₁ = vτ * vf₁
    vτf₂ = vτ * vf₂
    vτf₃ = vτ * vf₃
    #----

    error_F = total_error(vτf₁) + total_error(vτf₂) + total_error(vτf₃)
    norm_W₀ = construct_W₀_po(x_cheb, approx_tang_cheb, vτf₁, vτf₂, vτf₃, Pν, Wν)
    Y = opnorm_A * (norm_W₀ + error_F)

    ############
    # Z₁ bound #
    ############

    #----
    v∂ζ₁f₃ = - vu₁ * vu₃ * vg₁
    vτ∂ζ₁f₃ = vτ * v∂ζ₁f₃
    #----
    v∂ζ₂f₃ = - vu₂ * vu₃ * vg₂
    vτ∂ζ₂f₃ = vτ * v∂ζ₂f₃
    #----
    vJ = τ∂ᵤf(vτ, vζ₁, vζ₂, vu₁, vu₂, vu₃, (δ₁, δ₂, vα₁, vα₂, vλ₁, vλ₂); g₁ = vg₁, g₂ = vg₂)
    #----

    error_DF = max(
        total_error(vf₁) + total_error(vf₂) + total_error(vf₃),
        total_error(vτ∂ζ₁f₃),
        total_error(vτ∂ζ₂f₃),
        _opnorm1(total_error.(vJ))
        )

    # the Chebyshev data of W₁, evaluated at the 65 nodes of the products A W₁
    m = _npts(2N) # 65 nodes for N = 30
    tang_grid = [to_grid(project(component(approx_tang_cheb, j), Chebyshev(N) ⊗ Fourier(K, interval(1))), (m,)) for j ∈ 1:3]
    on_grid(v) = to_grid(project(sequence(v), Chebyshev(N) ⊗ Fourier(K, interval(1))), (m,))
    f₁_grid, f₂_grid, f₃_grid = on_grid(vf₁), on_grid(vf₂), on_grid(vf₃)
    τ∂ζ₁f₃_grid, τ∂ζ₂f₃_grid = on_grid(vτ∂ζ₁f₃), on_grid(vτ∂ζ₂f₃)
    J_grid = on_grid.(vJ)
    W₁_grid = [construct_W₁_po(Sequence(Fourier(K, interval(1))^3, mapreduce(j -> coefficients(tang_grid[j][ℓ]), vcat, 1:3)),
                               f₁_grid[ℓ], f₂_grid[ℓ], f₃_grid[ℓ], τ∂ζ₁f₃_grid[ℓ], τ∂ζ₂f₃_grid[ℓ], getindex.(J_grid, ℓ)) for ℓ ∈ 1:m]

    Z₁ = _opnorm_AW₁_I(A_cheb, W₁_grid, N, Pν) + _opnorm1(norm.(sequence.(vJ), Ref(Wν))) / interval(K + 1) + opnorm_A * error_DF

    ############
    # Z₂ bound #
    ############

    α₁_ot, α₂_ot = a₁ * sequence(vκ⁻¹), a₂ * sequence(vκ⁻¹)
    λ₁_ot, λ₂_ot = a₁ * sequence(vκ⁻¹) * d₁ / (m₁ - d₁), a₂ * sequence(vκ⁻¹) * d₂ / (m₂ - d₂)

    rig_g₁ = norm(sequence(vg₁), Wν)/(1 - norm(1 - (sequence(vα₁) + u₃_ot)*sequence(vg₁), Wν) - norm(sequence(vg₁), Wν) * R)
    rig_g₂ = norm(sequence(vg₂), Wν)/(1 - norm(1 - (sequence(vα₂) + u₃_ot)*sequence(vg₂), Wν) - norm(sequence(vg₂), Wν) * R)

    # ||∂τ∇f||
    opnorm_∂τ∇τf₁ = max(
        δ₁ * (norm(u₃_ot - λ₁_ot, Wν) + R) * rig_g₁,
        δ₁ * norm(α₁_ot + λ₁_ot, Wν) * (norm(u₁_ot, Wν) + R) * rig_g₁^2
        )
    opnorm_∂τ∇τf₂ = max(
        δ₂ * (norm(u₃_ot - λ₂_ot, Wν) + R) * rig_g₂,
        δ₂ * norm(α₂_ot + λ₂_ot, Wν) * (norm(u₂_ot, Wν) + R) * rig_g₂^2
        )
    opnorm_∂τ∇τf₃ = max(
        (norm(u₁_ot, Wν) + R) * (norm(u₃_ot, Wν) + R) * rig_g₁,
        (norm(u₂_ot, Wν) + R) * (norm(u₃_ot, Wν) + R) * rig_g₂,
        (norm(ζ₁_ot, Wν) + R) * (norm(u₃_ot, Wν) + R) * rig_g₁,
        (norm(ζ₂_ot, Wν) + R) * (norm(u₃_ot, Wν) + R) * rig_g₂,
        1 + 2 * (norm(u₃_ot, Wν) + R) + norm(α₁_ot, Wν) * (norm(ζ₁_ot, Wν) + R) * (norm(u₁_ot, Wν) + R) * rig_g₁^2 +
                                        norm(α₂_ot, Wν) * (norm(ζ₂_ot, Wν) + R) * (norm(u₂_ot, Wν) + R) * rig_g₂^2
        )
    # ||∂ζ₁∇τf₃||
    opnorm_∂ζ₁∇τf₃ = max(
        (norm(u₁_ot, Wν) + R) * (norm(u₃_ot, Wν) + R) * rig_g₁,
        (norm(τ_ot,  Wν) + R) * (norm(u₃_ot, Wν) + R) * rig_g₁,
        (norm(τ_ot,  Wν) + R) * (norm(u₁_ot, Wν) + R) * norm(α₁_ot, Wν) * rig_g₁^2
        )
    # ||∂ζ₂∇τf₃||
    opnorm_∂ζ₂∇τf₃ = max(
        (norm(u₂_ot, Wν) + R) * (norm(u₃_ot, Wν) + R) * rig_g₂,
        (norm(τ_ot,  Wν) + R) * (norm(u₃_ot, Wν) + R) * rig_g₂,
        (norm(τ_ot,  Wν) + R) * (norm(u₂_ot, Wν) + R) * norm(α₂_ot, Wν) * rig_g₂^2
        )
    # ||∂u₁∇τf₁||
    opnorm_∂u₁∇τf₁ = max(
        δ₁ * (norm(u₃_ot - λ₁_ot, Wν) + R) * rig_g₁,
        δ₁ * (norm(τ_ot, Wν) + R) * norm(α₁_ot + λ₁_ot, Wν) * rig_g₁^2
        )
    # ||∂u₁∇τf₃||
    opnorm_∂u₁∇τf₃ = max(
        (norm(ζ₁_ot, Wν) + R) * (norm(u₃_ot, Wν) + R) * rig_g₁,
        (norm(τ_ot,  Wν) + R) * (norm(u₃_ot, Wν) + R) * rig_g₁,
        (norm(τ_ot,  Wν) + R) * (norm(ζ₁_ot, Wν) + R) * norm(α₁_ot, Wν) * rig_g₁^2
        )
    # ||∂u₂∇τf₂||
    opnorm_∂u₂∇τf₂ = max(
        δ₂ * (norm(u₃_ot - λ₂_ot, Wν) + R) * rig_g₂,
        δ₂ * (norm(τ_ot, Wν) + R) * norm(α₂_ot + λ₂_ot, Wν) * rig_g₂^2
        )
    # ||∂u₂∇τf₃||
    opnorm_∂u₂∇τf₃ = max(
        (norm(ζ₂_ot, Wν) + R) * (norm(u₃_ot, Wν) + R) * rig_g₂,
        (norm(τ_ot,  Wν) + R) * (norm(u₃_ot, Wν) + R) * rig_g₂,
        (norm(τ_ot,  Wν) + R) * (norm(ζ₂_ot, Wν) + R) * norm(α₂_ot, Wν) * rig_g₂^2
        )
    # ||∂u₃∇τf₁||
    opnorm_∂u₃∇τf₁ = max(
        δ₁ * norm(λ₁_ot + α₁_ot, Wν) * norm(u₁_ot, Wν) * rig_g₁^2,
        δ₁ * norm(τ_ot, Wν) * norm(λ₁_ot + α₁_ot, Wν) * rig_g₁^2,
        2δ₁ * norm(τ_ot, Wν) * norm(λ₁_ot + α₁_ot, Wν) * norm(u₁_ot, Wν) * rig_g₁^3
        )
    # ||∂u₃∇τf₂||
    opnorm_∂u₃∇τf₂ = max(
        δ₂ * norm(λ₂_ot + α₂_ot, Wν) * norm(u₂_ot, Wν) * rig_g₂^2,
        δ₂ * norm(τ_ot, Wν) * norm(λ₂_ot + α₂_ot, Wν) * rig_g₂^2,
        2δ₂ * norm(τ_ot, Wν) * norm(λ₂_ot + α₂_ot, Wν) * norm(u₂_ot, Wν) * rig_g₂^3
        )
    # ||∂u₃∇τf₃||
    opnorm_∂u₃∇τf₃ = max(
        1 + 2 * (norm(u₃_ot, Wν) + R) + norm(α₁_ot, Wν) * (norm(ζ₁_ot, Wν) + R) * (norm(u₁_ot, Wν) + R) * rig_g₁^2 +
                                        norm(α₂_ot, Wν) * (norm(ζ₂_ot, Wν) + R) * (norm(u₂_ot, Wν) + R) * rig_g₂^2,
        (norm(τ_ot,  Wν) + R) * norm(α₁_ot, Wν) * (norm(u₁_ot, Wν) + R) * rig_g₁^2,
        (norm(τ_ot,  Wν) + R) * norm(α₂_ot, Wν) * (norm(u₂_ot, Wν) + R) * rig_g₂^2,
        (norm(τ_ot,  Wν) + R) * norm(α₁_ot, Wν) * (norm(ζ₁_ot, Wν) + R) * rig_g₁^2,
        (norm(τ_ot,  Wν) + R) * norm(α₂_ot, Wν) * (norm(ζ₂_ot, Wν) + R) * rig_g₂^2,
        2 * (norm(τ_ot,  Wν) + R) * (1 + norm(α₁_ot, Wν) * (norm(ζ₁_ot, Wν) + R) * (norm(u₁_ot, Wν) + R) * rig_g₁^3 +
                                        norm(α₂_ot, Wν) * (norm(ζ₂_ot, Wν) + R) * (norm(u₂_ot, Wν) + R) * rig_g₂^3)
        )

    Z₂ = opnorm_A * max(
        opnorm_∂τ∇τf₁ + opnorm_∂τ∇τf₂ + opnorm_∂τ∇τf₃,
        opnorm_∂ζ₁∇τf₃,
        opnorm_∂ζ₂∇τf₃,
        opnorm_∂u₁∇τf₁ + opnorm_∂u₁∇τf₃,
        opnorm_∂u₂∇τf₂ + opnorm_∂u₂∇τf₃,
        opnorm_∂u₃∇τf₁ + opnorm_∂u₃∇τf₂ + opnorm_∂u₃∇τf₃
        )

    #

    ie, _ = interval_of_existence(Y, Z₁, Z₂, R)
    r = interval(inf(ie))

    true_τ  = InfiniteSequence(τ_ot,  Wν; total_error = r)
    true_ζ₁ = InfiniteSequence(ζ₁_ot, Wν; total_error = r)
    true_ζ₂ = InfiniteSequence(ζ₂_ot, Wν; total_error = r)

    true_u₁ = InfiniteSequence(u₁_ot, Wν; total_error = r)
    true_u₂ = InfiniteSequence(u₂_ot, Wν; total_error = r)
    true_u₃ = InfiniteSequence(u₃_ot, Wν; total_error = r)

    return true_τ, true_ζ₁, true_ζ₂, true_u₁, true_u₂, true_u₃,
        (; opnorm_A, Y, Z₁, Z₂, rig_g₁, rig_g₂, ie, r)
end

#

function prove_lin(x_lin_grid, A_lin_grid, true_τ, true_ζ₁, true_ζ₂, true_u₁, true_u₂, true_u₃, κ₁, κ₂, (γ, m₁, m₂, y₁, y₂, d₁, d₂, a₁, a₂), K_lin, N, ν)
    x_lin_cheb_ = grid2cheb(x_lin_grid, N)
    _conjugacy_symmetry_cheb!(x_lin_cheb_)
    _polish_cheb!(x_lin_cheb_) # polish Chebyshev coefficients
    x_lin_cheb = interval(x_lin_cheb_)

    A_lin_cheb = interval.(grid2cheb(A_lin_grid, N))

    # Banach spaces

    Pν = Ell1(GeometricWeight(ν))
    Wν = Ell1(GeometricWeight(ν), IdentityWeight())

    opnorm_A_lin = max(_opnorm1(_cheb_norms(A_lin_cheb, Pν)), inv(interval(K_lin + 1)))

    # InfiniteSequence allows to rigorously control the error when inverting a sequence

    vκ⁻¹ = InfiniteSequence(Chebyshev(1) ⊗ Fourier(0, interval(1)), [(inv(κ₁) + inv(κ₂))/2, (inv(κ₂) - inv(κ₁))/4], Wν)
    δ₁, δ₂ = (m₁ - d₁) / γ, (m₂ - d₂) / γ
    vα₁, vα₂ = a₁ * vκ⁻¹, a₂ * vκ⁻¹
    vλ₁, vλ₂ = vα₁ * d₁ / (m₁ - d₁), vα₂ * d₂ / (m₂ - d₂)

    true_J = τ∂ᵤf(true_τ, true_ζ₁, true_ζ₂, true_u₁, true_u₂, true_u₃, (δ₁, δ₂, vα₁, vα₂, vλ₁, vλ₂))

    ###########
    # Y bound #
    ###########

    vecC_lin, vecV_lin = eachcomponent(x_lin_cheb)
    vC_lin = [InfiniteSequence(Sequence(Chebyshev(N) ⊗ Fourier(0, interval(1)),     coefficients(component(vecC_lin, i+3*(j-1)))), Wν) for i ∈ 1:3, j ∈ 1:3]
    vV_lin = [InfiniteSequence(Sequence(Chebyshev(N) ⊗ Fourier(K_lin, interval(1)), coefficients(component(vecV_lin, i+3*(j-1)))), Wν) for i ∈ 1:3, j ∈ 1:3]
    matF1 = evaluate.(vV_lin, Ref((nothing, 0))) - I[1:3,1:3]
    matF2 = _differentiate.(vV_lin, Ref((0, 1))) + vV_lin * vC_lin - true_J * vV_lin

    Y = opnorm_A_lin * (sum(_vnorm, vec(matF1)) + sum(_vnorm, vec(matF2)))

    ############
    # Z₁ bound #
    ############

    m = _npts(2N) # 65 nodes, as for the periodic orbits
    C_grid = [to_grid(project(sequence(vC_lin[i]), Chebyshev(N) ⊗ Fourier(0, interval(1))), (m,)) for i ∈ 1:9] # vec(C) at the nodes
    V_grid = [to_grid(project(sequence(vV_lin[i]), Chebyshev(N) ⊗ Fourier(K_lin, interval(1))), (m,)) for i ∈ 1:9]
    J_grid = [to_grid(project(sequence(true_J[i]), Chebyshev(N) ⊗ Fourier(K_lin, interval(1))), (m,)) for i ∈ 1:9]
    W₁_lin_grid = [construct_W₁_lin([C_grid[i][ℓ][0] for i ∈ 1:9],
                                    reshape([V_grid[i][ℓ] for i ∈ 1:9], 3, 3),
                                    reshape([J_grid[i][ℓ] for i ∈ 1:9], 3, 3)) for ℓ ∈ 1:m]

    Z₁ = _opnorm_AW₁_I(A_lin_cheb, W₁_lin_grid, N, Pν) +
        _opnorm1(norm.(sequence.(true_J), Ref(Wν))) / interval(K_lin + 1) +
        opnorm_A_lin * _opnorm1(total_error.(true_J))

    ############
    # Z₂ bound #
    ############

    Z₂ = 2opnorm_A_lin

    #

    ie, _ = interval_of_existence(Y, Z₁, Z₂, Inf)
    r = interval(inf(ie))

    true_C = [InfiniteSequence(sequence(c), Wν; total_error = r) for c ∈ vC_lin]

    return true_C, (; opnorm_A = opnorm_A_lin, Y, Z₁, Z₂, ie, r)
end
