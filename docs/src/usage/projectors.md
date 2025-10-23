# Computing projectors onto irreps of $\phi$

## Compute projectors

Use `projection_matrices` to compute the projection matrices. It will store them
in memory at the end and fetch those you already computed. It will create a
directory `projectors/` to store projectors locally.

```julia
Π = projection_matrices(m, n)
```
It will compute all projectors for $m$ modes and **up to** $n$ photons. Π is a
hashmap `[m, n, k] => SparseMatrixCSC{Float32, Int64}`. For instance for 3 modes and 2 photons:

Start by enabling debugging messages with
`using Logging; global_logger(ConsoleLogger(stdout,Logging.Debug))`.

```julia-repl
julia> Π = projection_matrices(3,2)
 Debug: Infer Last one
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:292
┌ Debug: m=3 n=0 computed projectors 
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:300
┌ Debug: m3-n0 saved at
│   path = "./projectors/P-m3-n0"
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:308
┌ Debug: Compute GC m = 3, n = 1, k = 0
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:75
┌ Debug: Infer Last one
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:292
┌ Debug: m=3 n=1 computed projectors 
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:300
┌ Debug: m3-n1 saved at
│   path = "./projectors/P-m3-n1"
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:308
┌ Debug: Compute GC m = 3, n = 2, k = 0
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:75
┌ Debug: Compute GC m = 3, n = 2, k = 1
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:75
┌ Debug: Infer Last one
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:292
┌ Debug: m=3 n=2 computed projectors 
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:300
┌ Debug: m3-n2 saved at
│   path = "./projectors/P-m3-n2"
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:308
Dict{Any, Any} with 6 entries:
  [3, 2, 1] => sparse([1, 2, 3, 7, 1, 8, 9, 2, 10, 3  …  33, 14, 28, 35, 1, 8, 15, 22, 29, 36], [1, 2, 3, 7, 8, 8, 9…
  [3, 1, 1] => sparse([1, 2, 3, 4, 1, 5, 6, 7, 8, 1, 5, 9], [1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 9, 9], [0.666667, 1.0, 1.…
  [3, 2, 0] => sparse([1, 1, 8, 1, 8, 15, 1, 8, 15, 22  …  8, 15, 22, 29, 1, 8, 15, 22, 29, 36], [1, 8, 8, 15, 15, 1…
  [3, 1, 0] => sparse([1, 1, 5, 1, 5, 9], [1, 5, 5, 9, 9, 9], Float32[0.333333, 0.333333, 0.333333, 0.333333, 0.3333…
  [3, 2, 2] => sparse([1, 2, 3, 4, 5, 6, 7, 1, 8, 9  …  34, 14, 28, 35, 1, 8, 15, 22, 29, 36], [1, 2, 3, 4, 5, 6, 7,…
  [3, 0, 0] => sparse([1], [1], [1.0], 1, 1)
```
All the projectors for less that $2$ photons in $3$ modes are computed, and for
each of them, $k$ ranges from $0$ to $n$. The projectors to the $k$-th irrep of
the $n$-photon state is accessed via `Π[[m, n, k]]` (notice that the hashmap's
key is a **list**). The storage format is detailed below.

If `save_in_memory` is set to `true`, the function call will create the files:
- projectors/P-m3-n0
- projectors/P-m3-n1
- projectors/P-m3-n2
so that another call with now 3 photons displays
```julia-repl
julia> Π = projection_matrices(3, 3)
┌ Debug: m=3 n=0 loaded from memory 
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:281
┌ Debug: m=3 n=1 loaded from memory 
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:281
┌ Debug: m=3 n=2 loaded from memory 
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:281
┌ Debug: Compute GC m = 3, n = 3, k = 0
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:75
┌ Debug: loaded CGC from disk: Irrep[SU₃]("10⁺") ⊗ Irrep[SU₃]("10") → Irrep[SU₃]("1")
└ @ SUNRepresentations .../caching.jl:23
┌ Debug: Compute GC m = 3, n = 3, k = 1
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:75
┌ Debug: loaded CGC from disk: Irrep[SU₃]("10⁺") ⊗ Irrep[SU₃]("10") → Irrep[SU₃]("8")
└ @ SUNRepresentations .../caching.jl:23
┌ Debug: Compute GC m = 3, n = 3, k = 2
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:75
┌ Debug: loaded CGC from disk: Irrep[SU₃]("10⁺") ⊗ Irrep[SU₃]("10") → Irrep[SU₃]("27")
└ @ SUNRepresentations .../caching.jl:23
┌ Debug: Infer Last one
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:292
┌ Debug: m=3 n=3 computed projectors 
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:300
┌ Debug: m3-n3 saved at
│   path = "./projectors/P-m3-n3"
└ @ Main.SUNIrrepProjectors .../SUNIrrepProjectors.jl:308
Dict{Any, Any} with 10 entries:
  [3, 2, 1] => sparse([1, 2, 3, 7, 1, 8, 9, 2, 10, 3  …  33, 14, 28, 35, 1, 8, 15, 22, 29, 36], [1, 2, 3, 7, 8, 8, 9…
  [3, 1, 1] => sparse([1, 2, 3, 4, 1, 5, 6, 7, 8, 1, 5, 9], [1, 2, 3, 4, 5, 5, 6, 7, 8, 9, 9, 9], [0.666667, 1.0, 1.…
  [3, 2, 0] => sparse([1, 1, 8, 1, 8, 15, 1, 8, 15, 22  …  8, 15, 22, 29, 1, 8, 15, 22, 29, 36], [1, 8, 8, 15, 15, 1…
  [3, 1, 0] => sparse([1, 1, 5, 1, 5, 9], [1, 5, 5, 9, 9, 9], Float32[0.333333, 0.333333, 0.333333, 0.333333, 0.3333…
  [3, 3, 3] => sparse([1, 2, 3, 4, 5, 6, 7, 8, 9, 10  …  1, 12, 23, 34, 45, 56, 67, 78, 89, 100], [1, 2, 3, 4, 5, 6,…
  [3, 2, 2] => sparse([1, 2, 3, 4, 5, 6, 7, 1, 8, 9  …  34, 14, 28, 35, 1, 8, 15, 22, 29, 36], [1, 2, 3, 4, 5, 6, 7,…
  [3, 3, 1] => sparse([1, 2, 3, 11, 1, 12, 13, 2, 14, 3  …  77, 88, 99, 1, 12, 34, 56, 67, 89, 100], [1, 2, 3, 11, 1…
  [3, 3, 0] => sparse([1, 1, 12, 1, 12, 23, 1, 12, 23, 34  …  1, 12, 23, 34, 45, 56, 67, 78, 89, 100], [1, 12, 12, 2…
  [3, 3, 2] => sparse([1, 2, 3, 4, 5, 6, 11, 1, 12, 13  …  1, 12, 23, 34, 45, 56, 67, 78, 89, 100], [1, 2, 3, 4, 5, …
  [3, 0, 0] => sparse([1], [1], [1.0], 1, 1)
```

and the following file is created:
- projectors/P-m3-n3


## Usage of the projectors
The projectors are stored in Sparse (CSC) format. As they are symmetric, only
the upper triangular part is stored:

```julia-repl
julia> Π = projection_matrices(2,2)
julia> Π[[2,2,2]]
9×9 SparseMatrixCSC{Float32, Int64} with 14 stored entries:
 0.166667   ⋅    ⋅    ⋅   -0.333333    ⋅    ⋅     ⋅    0.166667
  ⋅        0.5   ⋅    ⋅     ⋅        -0.5   ⋅     ⋅     ⋅ 
  ⋅         ⋅   1.0   ⋅     ⋅          ⋅    ⋅     ⋅     ⋅ 
  ⋅         ⋅    ⋅   0.5    ⋅          ⋅    ⋅   -0.5    ⋅ 
  ⋅         ⋅    ⋅    ⋅    0.666667    ⋅    ⋅     ⋅   -0.333333
  ⋅         ⋅    ⋅    ⋅     ⋅         0.5   ⋅     ⋅     ⋅ 
  ⋅         ⋅    ⋅    ⋅     ⋅          ⋅   1.0    ⋅     ⋅ 
  ⋅         ⋅    ⋅    ⋅     ⋅          ⋅    ⋅    0.5    ⋅ 
  ⋅         ⋅    ⋅    ⋅     ⋅          ⋅    ⋅     ⋅    0.166667
```

We provide a helper function `sparse_sym_mvp(A, x)` that performs the
matrix-vector multiplication $Ax$ by keeping $A$ in sparse and upper-triangular
format.