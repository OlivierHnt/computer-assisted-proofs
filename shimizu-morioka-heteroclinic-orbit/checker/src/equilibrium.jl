# Equilibrium (Section 2.1 of the paper). `validate` takes the approximate zero
# c_bar and returns the Y, Z bounds and the enclosure of Theorem 1.1.
# From ShimizuMoriokaTutorial.jl/src/equilibrium.jl; Step 2 (Newton) is in search/.

struct EquilibriumProblem
    f
    Df
    parameters
end

#-

function validate(pb::EquilibriumProblem, c_bar::Vector{Float64})
    pb_approx = EquilibriumProblem(pb.f, pb.Df, mid.(pb.parameters))

    # Step 2: c_bar is an input (data/c1_bar.bin)

    # Step 3
    A = approximate_inverse(pb_approx, c_bar)

    # Step 4
    #-
    Y = Y_bound(pb, interval(c_bar), A)
    R = interval(10sup(Y))
    Z = Z_bound(pb, interval(c_bar), A, R)
    #-
    ie, proof_success = interval_of_existence(Y, Z, R)

    #
    return (; Y, Z, R, ie, proof_success, c_star = interval.(c_bar, inf(ie); format = :midpoint))
end





# Step 1

F(pb::EquilibriumProblem, u) = pb.f(u, pb.parameters)

DF(pb::EquilibriumProblem, u) = pb.Df(u, pb.parameters)





# Step 3

approximate_inverse(pb_approx::EquilibriumProblem, c_bar) =
    interval(inv(DF(pb_approx, c_bar)))





# Step 4

Y_bound(pb::EquilibriumProblem, c_bar, A) = norm(A * F(pb, c_bar), 1)

function Z_bound(pb::EquilibriumProblem, c_bar, A, R)
    box_c = interval.(c_bar, R; format = :midpoint)
    return opnorm(interval(I) - A * DF(pb, box_c), 1)
end
