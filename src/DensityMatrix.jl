"""
    DensityMatrix(m::Int, n::Int, ρ::Matrix{ComplexF32})
	DensityMatrix(m::Int, n::Int, blocks::Vector{Matrix{ComplexF32}})

Stores the density matrix of an ``m``-mode Fock state with at most ``n``
photons.
"""
struct DensityMatrix
	m::Int
	n::Int
	d::Int
	ρ::Matrix{ComplexF32}

	DensityMatrix(m::Int, n::Int, ρ::Matrix{ComplexF32}) = begin
		@argcheck m > 0 "m should be greater than 0, got $(m)"
		@argcheck n > 0 "n should be greater than 0, got $(n)"

		d = binomial(m + n, n)
		@argcheck size(ρ) == (d, d) "incorrect matrix dimension for ρ: expected  $((d, d)), got $(string(size(ρ)))."
		@argcheck tr(ρ) ≈ 1 "ρ is not normalised"
		@argcheck all(real.(eigvals(ρ)) .>= 0) && all(imag.(eigvals(ρ)) .≈ 0) "ρ is not positive semi-definite"
		@argcheck ishermitian(ρ) "ρ is not self-adjoint"

		return new(m, n, d, ρ)
	end

	DensityMatrix(m::Int, n::Int, blocks::Vector{Matrix{ComplexF32}}) = begin
		@argcheck m > 0 "m should be greater than 0, got $(m)"
		@argcheck n > 0 "n should be greater than 0, got $(n)"

		for (i, b) in enumerate(blocks)
			d = binomial(m + (i - 1) - 1, (i - 1))
			@argcheck size(b) == (d, d) "incorrect matrix dimension for $(i)-th block: expected $((d, d)), got $(size(b))."
		end
		
		M = direct_sum(blocks)

		return DensityMatrix(m, n, M / tr(M))
	end
end

"""
    DensityMatrixBlock(m::Int, n::Int, ρ::Matrix{ComplexF32})

Stores the density matrix of an ``m``-mode, ``n``-photon Fock state restricted
to the ``n``-photon sector.
"""
struct DensityMatrixBlock
	m::Int
	n::Int
	d::Int
	ρ::Matrix{ComplexF32}

	DensityMatrixBlock(m::Int, n::Int, ρ::Matrix{ComplexF32}) = begin
		@argcheck m > 0 "m should be greater than 0, got $(m)"
		@argcheck n > 0 "n should be greater than 0, got $(n)"

		d = binomial(m + n - 1, n)
		@argcheck size(ρ) == (d, d) "incorrect matrix dimension for ρ: expected $((d, d)), got $(size(ρ))."

		return new(m, n, d, ρ)
	end
end

"""
    DensityVectorBlock(m::Int, n::Int, ρ::Matrix{ComplexF32})

Stores the density matrix of an ``m``-mode, ``n``-photon Fock state restricted
to the ``n``-photon sector in vectorised form.
"""
struct DensityVectorBlock
	m::Int
	n::Int
	d::Int
	ρ::Vector{ComplexF32}

	DensityMatrixBlock(m::Int, n::Int, ρ::Vector{ComplexF32}) = begin
		@argcheck m > 0 "m should be greater than 0, got $(m)"
		@argcheck n > 0 "n should be greater than 0, got $(n)"

		d = binomial(m + n - 1, n)
		@argcheck size(ρ) == (d,) "incorrect matrix dimension for ρ: expected $((d,)), got $(size(ρ))."

		return new(m, n, d, ρ)
	end

end

"""
    extract_blocks(ρ::DensityMatrix)

Extracts the block corresponding to each photon sectors.
"""
function extract_blocks(dm::DensityMatrix)
	block_sizes = [binomial(m + l - 1, l) for l in 0:dm.n]

	block_start_indexes = [
		1 + sum(block_sizes[begin:l]) for l in 0:dm.n
	]

	block_end_indexes = [
		block_start_indexes[l] + binomial(m + (l - 1) - 1, (l - 1)) - 1 for l in 1:dm.n+1
	]

	return [
		DensityMatrixBlock(m, i, dm.ρ[start:stop, start:stop])
		for (start, stop) in zip(block_start_indexes, block_end_indexes)
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

	lin_ρ = zeros(ComplexF32, (d^2))

	for j in 1:d
		for i in 1:d
			lin_ρ[id] = dm.ρ[j, i]
			id += 1
		end
	end
	return DensityVectorBlock(m, n, lin_ρ)
end