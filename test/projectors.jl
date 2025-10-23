using LinearAlgebra

"""
converts sparse symmetric matrix to its dense version
"""
toarray(P) = (P + transpose(P) - Diagonal(P)) |> Array

@testset "testing projector properties" begin
	m = 3
	n = 2

	Π = projection_matrices(m, n; load_from_memory = false, save_in_memory = false)

	@testset "projector sum should equal identity" begin
		for _n in 1:2
			d = binomial(_n + m - 1, _n)^2
			@test sum(Π[[m, _n, k]] for k in 0:_n) |> toarray ≈ Matrix(I, d, d) atol = 10e-8
		end
	end

	@testset "trace of projector should equal their dimension" begin
		for _n in 1:2
			for k in 0:_n
				@test Π[[m, _n, k]] |> toarray |> tr ≈ d_λ(m, k) atol = 10e-8
			end
		end
	end

	@testset "projectors should be projectors" begin
		for _n in 1:2
			for k in 0:_n
				P = Π[[m, _n, k]] |> toarray
				@test P^2 ≈ P atol = 10e-6
			end
		end
	end

    @testset "all but projection onto trivial irrep are traceless" begin
        Test.Pass
    end
end