import Base.+, Base.sum, Base./

"""
	DensityMatrixBlock(m::Int, n::Int, ρ::Matrix{ComplexF64})

Stores the density matrix of an ``m``-mode, ``n``-photon Fock state restricted
to the ``n``-photon sector.
"""
struct DensityMatrixBlock
	m::Int
	n::Int
	d::Int
	ρ::Matrix{ComplexF64}

	DensityMatrixBlock(m::Int, n::Int, ρ::Matrix{ComplexF64}) = begin
		@argcheck m >= 0 "m should be greater than 0, got $(m)"
		@argcheck n >= 0 "n should be greater than 0, got $(n)"

		d = binomial(m + n - 1, n)
		@argcheck size(ρ) == (d, d) "incorrect matrix dimension for ρ: expected $((d, d)), got $(size(ρ)). (m = $m, n=$n))"

		return new(m, n, d, ρ)
	end
end

"""
	DensityMatrix(m::Int, n::Int, ρ::Matrix{ComplexF64})
	DensityMatrix(m::Int, n::Int, blocks::Vector{Matrix{ComplexF64}})

Stores the density matrix of an ``m``-mode Fock state with at most ``n``
photons.
"""
struct DensityMatrix
	m::Int
	n::Int
	d::Int
	ρ::Matrix{ComplexF64}

	DensityMatrix(m::Int, n::Int, ρ::Matrix{ComplexF64}) = begin
		@argcheck m >= 0 "m should be greater than 0, got $(m)"
		@argcheck n >= 0 "n should be greater than 0, got $(n)"

		d = binomial(m + n, n)
		@argcheck size(ρ) == (d, d) "incorrect matrix dimension for ρ: expected  $((d, d)), got $(string(size(ρ)))."
		# @argcheck tr(ρ) ≈ 1 "ρ is not normalised"
		spectrum = eigvals(ρ)
		# @argcheck all((real.(spectrum) .>= 0) .|| isapprox.(real.(spectrum), 0.0; atol = 1e-10)) && all(isapprox.(imag.(spectrum), 0.0; atol = 1e-10)) "ρ is not positive semi-definite"
		# @argcheck ishermitian(ρ) "ρ is not self-adjoint" 
		return new(m, n, d, ρ)
	end

	DensityMatrix(m::Int, n::Int, blocks::Vector{Matrix{ComplexF64}}) = begin
		@argcheck m > 0 "m should be greater than 0, got $(m)"
		@argcheck n > 0 "n should be greater than 0, got $(n)"

		for (i, b) in enumerate(blocks)
			d = binomial(m + (i - 1) - 1, (i - 1))
			@argcheck size(b) == (d, d) "incorrect matrix dimension for $(i)-th block: expected $((d, d)), got $(size(b))."
		end

		M = direct_sum(blocks)

		return DensityMatrix(m, n, M / tr(M))
	end

	DensityMatrix(m::Int, n::Int, blocks::Vector{DensityMatrixBlock}) = begin
		return DensityMatrix(m, n, [b.ρ for b in blocks])
	end
end




"""
	DensityVectorBlock(m::Int, n::Int, ρ::Matrix{ComplexF64})

Stores the density matrix of an ``m``-mode, ``n``-photon Fock state restricted
to the ``n``-photon sector in vectorised form.
"""
struct DensityVectorBlock
	m::Int
	n::Int
	d::Int
	ρ::Vector{ComplexF64}

	DensityVectorBlock(m::Int, n::Int, ρ::Vector{ComplexF64}) = begin
		@argcheck m >= 0 "m should be greater than 0, got $(m)"
		@argcheck n >= 0 "n should be greater than 0, got $(n)"

		d = binomial(m + n - 1, n)
		@argcheck size(ρ) == (d^2,) "incorrect matrix dimension for ρ: expected $((d,)), got $(size(ρ))."

		return new(m, n, d, ρ)
	end

end

"""
	extract_blocks(ρ::DensityMatrix)

Extracts the block corresponding to each photon sectors.
"""
function extract_blocks(ρ::DensityMatrix)
	block_sizes = [binomial(ρ.m + l - 1, l) for l in 0:ρ.n]

	block_start_indexes = [
		1 + sum(block_sizes[begin:l]) for l in 0:ρ.n
	]

	block_end_indexes = [
		block_start_indexes[l] + binomial(ρ.m + (l - 1) - 1, (l - 1)) - 1 for l in 1:ρ.n+1
	]

	return [
		DensityMatrixBlock(ρ.m, i - 1, ρ.ρ[start:stop, start:stop])
		for (i, (start, stop)) in enumerate(zip(block_start_indexes, block_end_indexes))
	]
end

"""
	linearise(ρ::DensityMatrixBlock)

Linearise a ``d \\times d`` density matrix block (in ``\\mathcal{H}_m^n``),
where ``d = \\binom{n+m-1}{n}``, to a ``d^2``-dimensional
[`DensityVectorBlock`](@ref).
"""
function linearise(dm::DensityMatrixBlock)
	d = binomial(dm.n + dm.m - 1, dm.n)
	id = 1

	lin_ρ = zeros(ComplexF64, (d^2))

	for j in 1:d
		for i in 1:d
			lin_ρ[id] = dm.ρ[j, i]
			id += 1
		end
	end
	return DensityVectorBlock(dm.m, dm.n, lin_ρ)
end



Base.:+(A::DensityMatrixBlock, B::DensityMatrixBlock) = sum(A, B)
Base.:+(A::DensityMatrix, B::DensityMatrix) = sum(A, B)
Base.:*(A::DensityMatrix, x::Float64) = DensityMatrix(A.m, A.n, A.ρ * x)
Base.:*(C::Matrix{ComplexF32}, A::DensityMatrix) = DensityMatrix(A.m, A.n, C * A.ρ)
Base.:*(C::Matrix{ComplexF64}, A::DensityMatrix) = DensityMatrix(A.m, A.n, C * A.ρ)
Base.:*(C::Matrix{Float64}, A::DensityMatrix) = DensityMatrix(A.m, A.n, C * A.ρ)
Base.:*(C::Matrix{ComplexF64}, A::DensityMatrixBlock) = DensityMatrixBlock(A.m, A.n, C * A.ρ)
Base.:*(C::Matrix{Float64}, A::DensityMatrixBlock) = DensityMatrixBlock(A.m, A.n, C * A.ρ)
Base.:*(A::DensityMatrix, x::Int) = DensityMatrix(A.m, A.n, A.ρ * x)
Base.:*(A::DensityMatrixBlock, x::Int64) = DensityMatrixBlock(A.m, A.n, A.ρ * x)
Base.:*(A::DensityMatrixBlock, x::Float64) = DensityMatrixBlock(A.m, A.n, A.ρ * x)
Base.:/(A::DensityMatrix, x::Int64) = A * (1 / x)
Base.:/(A::DensityMatrix, x::Float64) = A * (1 / x)
Base.:/(A::DensityMatrixBlock, x::Int64) = A * (1 / x)
LinearAlgebra.tr(A::DensityMatrix) = tr(A.ρ)
LinearAlgebra.tr(A::DensityMatrixBlock) = tr(A.ρ)

function Base.sum(A::DensityMatrix, B::DensityMatrix)
	#@TODO: we should be able to sum density matrices with different photon
	#numbers
	return DensityMatrix(A.m, A.n, A.ρ + B.ρ)
end

function Base.sum(A::Vector{DensityMatrixBlock})
	#TODO check all dimensions
	return DensityMatrixBlock(A[1].m, A[1].n, sum(ρ.ρ for ρ in A))
end

function Base.sum(A::DensityMatrixBlock, B::DensityMatrixBlock)
	return DensityMatrixBlock(A.m, A.n, A.ρ + B.ρ)
end

function matricise_subblock(lin_ρ, m, n)
	"""
	Inverse of linearise, i.e., takes a linearised vector back to its matrix
	from.
	"""
	d = binomial(m + n - 1, n)
	ρ = zeros(ComplexF64, (d, d))

	id = 1
	for j in 1:d
		for i in 1:d
			ρ[j, i] = lin_ρ[id]
			id += 1
		end
	end


	return DensityMatrixBlock(m, n, ρ)
end

"""
    apply_channel(ρ::DensityMatrix, Π)

Apply the channel Π to ``\\rho``, i.e., returns ``M(\\rho)`` where Π is the
matrix form of the channel ``M``.
"""
function apply_channel(ρ::DensityMatrix, Π)
	blocks_of_ρ = extract_blocks(ρ)

	m = ρ.m

	ρ_depolarised = [
		matricise_subblock(sparse_sym_mvp(Π[i], linearise(b).ρ), m, i - 1)
		for (i, b) in enumerate(blocks_of_ρ)]


	return DensityMatrix(ρ.m, ρ.n, ρ_depolarised)
end

"""
    apply_nth_channel(ρ::DensityMatrix, Π)

Apply the channel Π to the ``n``-photon sector ``\\rho^{(n)}``,i.e., returns
``M^{(n)}(\\rho^{(n)})`` where Π is the matrix form of the channel ``M``.
"""
function apply_nth_channel(ρ::DensityMatrixBlock, Π)
	m = ρ.m
	n = ρ.n
	return matricise_subblock(sparse_sym_mvp(Π, linearise(ρ)), m, n)
end
