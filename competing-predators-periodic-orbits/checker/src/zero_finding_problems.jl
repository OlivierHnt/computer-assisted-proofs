# The zero-finding maps and their derivatives: F_po!, DF_po! for the family of
# periodic orbits (F of Section 2.1), τ∂ᵤf the Jacobian of the vector field,
# and F_lin!, DF_lin! for the Floquet normal form (G of Section 3).
# From src/zero_finding_problems.jl of OlivierHnt/BHWPeriodicCoexistence.jl, commit 6effd7d.

# Periodic orbit

function F_po!(F, x, (δ₁, δ₂, α₁, α₂, λ₁, λ₂), approx_tang)
    τ, ζ₁, ζ₂ = component(x, 1)
    u₁, u₂, u₃ = eachcomponent(component(x, 2))

    #

    component(F, 1)[1] = sum(coefficients(u₁) .* conj.(coefficients(component(approx_tang, 1)))
                          .+ coefficients(u₂) .* conj.(coefficients(component(approx_tang, 2)))
                          .+ coefficients(u₃) .* conj.(coefficients(component(approx_tang, 3)))) # phase condition
    component(F, 1)[2] = u₁(0) - 1
    component(F, 1)[3] = u₂(0) - 1

    #

    g₁ = inv(α₁ + u₃)
    g₂ = inv(α₂ + u₃)
    project!(component(component(F, 2), 1), differentiate(u₁) - τ * (δ₁ * (u₃ - λ₁) * g₁ * u₁))
    project!(component(component(F, 2), 2), differentiate(u₂) - τ * (δ₂ * (u₃ - λ₂) * g₂ * u₂))
    project!(component(component(F, 2), 3), differentiate(u₃) - τ * (1 - u₃ - ζ₁ * u₁ * g₁ - ζ₂ * u₂ * g₂) * u₃)

    #

    return F
end

function DF_po!(DF, x, (δ₁, δ₂, α₁, α₂, λ₁, λ₂), approx_tang)
    τ, ζ₁, ζ₂ = component(x, 1)
    u₁, u₂, u₃ = eachcomponent(component(x, 2))

    coefficients(DF) .= 0

    #

    D₂F₁ = component(DF, 1, 2)

    for j ∈ 1:3, k ∈ -order(component(approx_tang, j)):order(component(approx_tang, j))
        component(D₂F₁, 1, j)[1,k] = conj(component(approx_tang, j)[k])
    end
    project!(component(D₂F₁, 2, 1), Evaluation(0))
    project!(component(D₂F₁, 3, 2), Evaluation(0))

    #

    D₁F₂ = component(DF, 2, 1)
    g₁ = inv(α₁ + u₃)
    g₂ = inv(α₂ + u₃)

    project!(component(D₁F₂, 1, 1), -δ₁ * (u₃ - λ₁) * g₁ * u₁)
    project!(component(D₁F₂, 2, 1), -δ₂ * (u₃ - λ₂) * g₂ * u₂)
    project!(component(D₁F₂, 3, 1), -(1 - u₃ - ζ₁ * u₁ * g₁ - ζ₂ * u₂ * g₂) * u₃)

    project!(component(D₁F₂, 3, 2), τ * u₁ * u₃ * g₁)

    project!(component(D₁F₂, 3, 3), τ * u₂ * u₃ * g₂)

    #

    project!(component(component(DF, 2, 2), 1, 1), Multiplication(-τ * δ₁ * (u₃ - λ₁) * g₁))
    project!(component(component(DF, 2, 2), 1, 3), Multiplication(-τ * δ₁ * (α₁ + λ₁) * g₁^2 * u₁))
    project!(component(component(DF, 2, 2), 2, 2), Multiplication(-τ * δ₂ * (u₃ - λ₂) * g₂))
    project!(component(component(DF, 2, 2), 2, 3), Multiplication(-τ * δ₂ * (α₂ + λ₂) * g₂^2 * u₂))
    project!(component(component(DF, 2, 2), 3, 1), Multiplication(τ * ζ₁ * u₃ * g₁))
    project!(component(component(DF, 2, 2), 3, 2), Multiplication(τ * ζ₂ * u₃ * g₂))
    project!(component(component(DF, 2, 2), 3, 3), Multiplication(-τ * (1 - 2 * u₃ - α₁ * ζ₁ * u₁ * g₁^2 - α₂ * ζ₂ * u₂ * g₂^2)))

    radd!(component(DF, 2, 2), Derivative(1))

    #

    return DF
end



# Jacobian of the vector field

function τ∂ᵤf(x, params; g₁ = inv(params[3] + component(component(x, 2), 3)), g₂ = inv(params[4] + component(component(x, 2), 3)))
    τ, ζ₁, ζ₂ = component(x, 1)
    u₁, u₂, u₃ = eachcomponent(component(x, 2))

    return τ∂ᵤf(τ, ζ₁, ζ₂, u₁, u₂, u₃, params; g₁ = g₁, g₂ = g₂)
end

function τ∂ᵤf(τ, ζ₁, ζ₂, u₁, u₂, u₃, params; g₁ = inv(params[3] + u₃), g₂ = inv(params[4] + u₃))
    δ₁, δ₂, α₁, α₂, λ₁, λ₂ = params

    #
    J₁₁ = τ * δ₁ * (u₃ - λ₁) * g₁
    J₁₃ = τ * δ₁ * (α₁ + λ₁) * g₁^2 * u₁
    #
    J₂₂ = τ * δ₂ * (u₃ - λ₂) * g₂
    J₂₃ = τ * δ₂ * (α₂ + λ₂) * g₂^2 * u₂
    #
    J₃₁ = -τ * ζ₁ * u₃ * g₁
    J₃₂ = -τ * ζ₂ * u₃ * g₂
    J₃₃ = τ * (1 - 2 * u₃ - α₁ * ζ₁ * u₁ * g₁^2 - α₂ * ζ₂ * u₂ * g₂^2)
    #

    return [J₁₁ zero(u₁) J₁₃ ; zero(u₁) J₂₂ J₂₃ ; J₃₁ J₃₂ J₃₃]
end



# Floquet normal form

function F_lin!(F, x, J)
    vecC, vecV = eachcomponent(x)
    #
    C = [vecC[i + 3*(j-1)]            for i ∈ 1:3, j ∈ 1:3]
    V = [component(vecV, i + 3*(j-1)) for i ∈ 1:3, j ∈ 1:3]
    #
    matF1 = evaluate.(V, 0) - I[1:3,1:3]
    matF2 = differentiate.(V) + V * C - J * V
    for i ∈ 1:9
        component(F, 1)[i] = matF1[i]
        project!(component(component(F, 2), i), matF2[i])
    end
    #
    return F
end

function DF_lin!(DF, x, J)
    vecC, vecV = eachcomponent(x)
    #
    C = [vecC[i + 3*(j-1)]            for i ∈ 1:3, j ∈ 1:3]
    V = [component(vecV, i + 3*(j-1)) for i ∈ 1:3, j ∈ 1:3]
    #
    coefficients(DF) .= 0
    #
    k1 = kron(I[1:3,1:3], J)
    k2 = kron(I[1:3,1:3], V)
    k3 = kron(transpose(C), I[1:3,1:3])
    for j ∈ 1:9
        for i ∈ 1:9
            project!(component(component(DF, 2, 1), i, j), k2[i,j])
            project!(component(component(DF, 2, 2), i, j), Multiplication(k3[i,j] - k1[i,j]))
        end
        project!(component(component(DF, 1, 2), j, j), Evaluation(0))
        radd!(component(component(DF, 2, 2), j, j), Derivative(1))
    end
    #
    return DF
end
