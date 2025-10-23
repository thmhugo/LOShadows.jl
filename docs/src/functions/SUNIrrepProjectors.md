# SU(m) irreps projectors

```@docs
projection_matrices(m::Int, n::Int; load_from_memory = true, save_in_memory = true, infer_last::Bool = true)

sparse_sym_mvp(A::SparseMatrixCSC{Float64, Int64}, x::Vector{<:Number})

linearised_fock_basis_indices(m::Int, n::Int)
```
