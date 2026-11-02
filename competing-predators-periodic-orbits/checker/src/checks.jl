# Chebyshev grid <-> coefficient transforms, and the a posteriori checks of
# Sections 2.3 and 3 of the paper: positivity and boundary crossings of the
# family (prove_interior, prove_crossing; sign_derivative is the derivative
# enclosure of Lemma 2.6) and the localisation of the Floquet exponents in a
# box (prove_inside_box). From src/BHWPeriodicCoexistence.jl of
# OlivierHnt/BHWPeriodicCoexistence.jl, commit 6effd7d.

# Chebyshev interpolation, of order N, of node values given at η_ℓ = -cos(ℓπ/(npts-1)),
# ℓ = 0 first; `to_coef` expects η = +1 first, hence the `reverse`
grid2cheb(x_grid::Vector{<:Sequence}, N) = to_coef(reverse(x_grid), Chebyshev(N))

# entrywise Chebyshev interpolation of matrices given on the same nodes; returns
# the array C[n+1, i, j] of the Chebyshev coefficients of order n of entry (i, j)
function grid2cheb(A_grid::Vector{<:AbstractMatrix}, N)
    npts = length(A_grid)
    m, n = size(A_grid[1])
    C = Array{eltype(A_grid[1])}(undef, N+1, m, n)
    vals = Vector{eltype(A_grid[1])}(undef, npts)
    for j ∈ 1:n, i ∈ 1:m
        for ℓ ∈ 1:npts
            vals[ℓ] = A_grid[npts+1-ℓ][i,j]
        end
        C[:,i,j] .= coefficients(to_coef(vals, Chebyshev(N)))
    end
    return C
end

# values of Chebyshev(N) coefficient arrays C[n+1, i, j] at the m nodes cos(πj/(m-1)), j = 0, …, m-1
function cheb2grid(C::Array{T,3}, m) where {T}
    N = size(C, 1) - 1
    G = Array{T}(undef, m, size(C, 2), size(C, 3))
    for j ∈ axes(C, 3), i ∈ axes(C, 2)
        G[:,i,j] .= to_grid(Sequence(Chebyshev(N), C[:,i,j]), m)
    end
    return G
end

# positivity and crossings (Section 2.3)

function prove_interior(τ, ζ₁, ζ₂; a, b, n)
    check = true
    for η ∈ mince(interval(a, b), n)
        check &= strictprecedes(interval(0), real(τ( η, 0 ))) # valid periodic orbit
        check &= strictprecedes(interval(0), real(ζ₁( η, 0 ))) # ζ₁ positive
        check &= strictprecedes(interval(0), real(ζ₂( η, 0 ))) # ζ₂ positive
    end
    return check
end

function sign_derivative(f, η, ν)
    N = order(f)[1]
    val = zero(Interval{Float64})
    err_finite = zero(Interval{Float64})
    for n ∈ 1:N
        val += n * real(f[(n,0)]) * sin(n * acos(η))
        err_finite += n * abs(sin(n * acos(η)))
    end
    return sign(interval(val, total_error(f) * (err_finite + (ν - 1)*N + ν)/((ν - 1)^2 * ν^N); format = :midpoint))
end

function prove_crossing(τ, ζ, ζ_vanish; a, b, n, ν)
    check = strictprecedes(real(ζ_vanish( a, 0 )), interval(0)) | strictprecedes(real(ζ_vanish( b, 0 )), interval(0))
    for η ∈ mince(interval(a, b), n)
        check &= strictprecedes(interval(0), real(τ( η, 0 ))) # valid periodic orbit
        check &= strictprecedes(interval(0), real(ζ( η, 0 ))) # ζ positive on [a, b]
        check &= !in_interval(0, sign_derivative(ζ_vanish, η, ν)) # ζ_vanish strictly monotone on [a, b]
    end
    return check
end

# bisect the rectangle Ix × Iy until f is bounded on every piece (or a piece is thinner than 1e-6)

function evaluate_adaptive(f, Ix, Iy)
    queue = [(Ix, Iy)]

    check = true

    while !isempty(queue) && check
        current_x, current_y = pop!(queue)
        mat = f(current_x, current_y)

        if all(isbounded, mat)
            check = true
        elseif max(diam(current_x), diam(current_y)) > 1e-6
            mid_x = mid(current_x)
            mid_y = mid(current_y)

            push!(queue, (interval(inf(current_x), mid_x), interval(inf(current_y), mid_y)))
            push!(queue, (interval(inf(current_x), mid_x), interval(mid_y, sup(current_y))))
            push!(queue, (interval(mid_x, sup(current_x)), interval(inf(current_y), mid_y)))
            push!(queue, (interval(mid_x, sup(current_x)), interval(mid_y, sup(current_y))))
        else
            check = false
        end
    end

    return check
end

# Floquet exponents in a box (Section 3); `at` ∈ [a, b] is where the eigenvalues are counted
function prove_inside_box(C; a, b, at = a, left_side, right_side, bottom_side, top_side)
    C_at_0 = evaluate.(C, Ref((at, 0)))
    Λ_at_0 = eigvals(C_at_0)

    check = 0
    for j ∈ 1:3
        check += issubset_interval(Λ_at_0[j], complex(interval(left_side, right_side), interval(bottom_side, top_side))) # enclosures inside the box
    end

    #

    evaluate_adaptive((v, w) -> inv(evaluate.( C, Ref((w, 0)) ) - I * complex(left_side, v)), interval(bottom_side, top_side), interval(a, b)) || return check, false
    evaluate_adaptive((v, w) -> inv(evaluate.( C, Ref((w, 0)) ) - I * complex(right_side, v)), interval(bottom_side, top_side), interval(a, b)) || return check, false
    evaluate_adaptive((v, w) -> inv(evaluate.( C, Ref((w, 0)) ) - I * complex(v, bottom_side)), interval(left_side, right_side), interval(a, b)) || return check, false
    evaluate_adaptive((v, w) -> inv(evaluate.( C, Ref((w, 0)) ) - I * complex(v, top_side)), interval(left_side, right_side), interval(a, b)) || return check, false

    #

    return check, true
end
