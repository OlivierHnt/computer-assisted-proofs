# Newton solver from the authors' code (spiderweb.jl, 2019): the floating-point
# map F (f of the paper, Section 6), its Jacobian DF, and Newton's method. Used
# only by produce.jl; not part of the proof.

module spiderweb

	using LinearAlgebra

	export newton, F_DF

	function newton(x::Vector{T}, m0::T, m::Vector{T}, λ::T, n::Int, l::Int;
				max_iter = 50, atol = 1e-6, rtol = 1e-6) where {T<:Float64}
        # Newton iteration

        F, DF = F_DF(x, m0, m, λ, n, l)

        isapprox(0., norm(F, Inf), atol=atol, rtol=rtol) && return x

        for i in 1:max_iter
            x -= DF\F
			x = sort(abs.(x))
            F, DF = F_DF(x, m0, m, λ, n, l)
            isapprox(0., norm(F, Inf), atol=atol, rtol=rtol) && return x
        end
        error("No convergence after ",max_iter," iterations: |F| = ",norm(F, Inf),". Setup: atol = ",atol,", rtol = ",rtol,".")
	end


	function F_DF(x, m0, m, λ, n, l)
		# compute the map F whose zeros give a spiderweb central configuration
		# and its jacobian DF
		sqrt2 = sqrt(2)
		F = λ*x
		DF = Matrix{Float64}(I, n, n)*λ
		@inbounds for j ∈ eachindex(x)
			xj2 = x[j]^2
			for i ∈ eachindex(x)
				xi2, xixj = x[i]^2, x[i]*x[j]
				xi3 = xi2*x[i]

				if (j==1) F[i] += m0/xi2 end
				if (i==j) DF[i,j] += -2*m0/xi3 end

				for k ∈ 0:l-1
					θ = k*360/l
					cosdθ = cosd(θ)
					cosd2θ = cosd(2*θ)
					if i ≠ j
						d = (xi2 + xj2 - 2*xixj*cosdθ)^(1/2)
						d5 = 2*d^5
						F[i] += m[j]*(x[i] - x[j]*cosdθ)/d^3
						DF[i,i] += -m[j]*(4*xi2 + xj2 - 8*xixj*cosdθ + 3*xj2*cosd2θ)/d5
						DF[i,j] += -m[j]*(-4*(xi2 + xj2)*cosdθ + xixj*(7 + cosd2θ))/d5
					elseif i==j && k>0
						ζ = sqrt2*sqrt(1-cosdθ)
						F[i] += m[i]/(xi2*2*ζ)
						DF[i,j] += -m[i]/(xi3*ζ)
					end
				end
			end
		end
		return F, DF
	end

end
