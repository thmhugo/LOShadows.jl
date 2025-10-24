"""
    direct_sum(A::Matrix{<:Number}, B::Matrix{<:Number})

Compute the direct sum ``A \\oplus B``
"""
function direct_sum(A::Matrix{<:Number}, B::Matrix{<:Number})
    return [A zeros(ComplexF32, size(A, 1), size(B, 2)); zeros(ComplexF16, size(B, 1), size(A, 2)) B]
end

"""
    direct_sum(M::Vector{Matrix{T}} where T <: Number)

Wrapper of [`direct_sum(A::Matrix{<:Number}, B::Matrix{<:Number})`](@ref) for a
list of matrices, e.g., for M = [A, B, C, ...] returns ``A \\oplus B \\oplus C
\\oplus \\cdots``.
"""
function direct_sum(M::Vector{Matrix{T}} where T <: Number)

    Σ = M[1]

    for m in M[2:end]
        Σ = direct_sum(Σ, m)
    end

    return Σ

end

"""
	sparse_sym_mvp(A::SparseMatrixCSC{Float64, Int64}, x::Vector{<:Number})

Implementation of the sparse-matrix vector product ``Ax`` when ``A`` is
symmetric and given via its upper-triangular form.
"""
function sparse_sym_mvp(A::SparseMatrixCSC{Float64, Int64}, x::Vector{<:Number})
	return A * x + transpose(transpose(x) * A) - Diagonal(A) * x
end

function sparse_sym_mvp(A::SparseMatrixCSC{Float64, Int64}, ρ::DensityVectorBlock)
    return sparse_sym_mvp(A, ρ.ρ) 
end