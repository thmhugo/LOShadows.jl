# export projection_matrices, sparse_sym_mvp, linearised_fock_basis_indices

τ(n::Int, m::Int) = SUNIrrep(Tuple([[n]; zeros(Int, m - 1)]))
λ(k::Int, m::Int) = SUNIrrep(Tuple([[2 * k]; fill(k, m - 2); [0]]))
d_λ(m::Int, k::Int) = binomial(k + m - 2, k)^2 * (2k + m - 1) / (m - 1)

struct FockState
	n::Int
	m::Int
	state::Vector{Int}
	FockState(state) = all(state[:] .>= 0) ? new(sum(state), length(state), state) : error("negative photon counts")
end

struct FockMEl
	ket::FockState
	bra::FockState
	c::Float32
end

@memoize all_fock_states(m::Int, n::Int) = multiexponents(m, n)

""" 
	linearised_fock_basis_indices(m::Int, n::Int)

Indexes the linearised Fock basis used for the projection matrices.
"""
@memoize function linearised_fock_basis_indices(m::Int, n::Int)
	fock_indices = Dict{Vector{Int}, Int}()

	id = 1
	for S in all_fock_states(m, n)
		for T in all_fock_states(m, n)
			fock_indices[[S; T]] = id
			id += 1
		end
	end

	return fock_indices
end



function max_gt_pat(gt::SUNRepresentations.GTPattern)
	d = length(weight(gt))
	max_gt = Base.setindex(gt, gt[1, d], 1, 1)

	for l in 2:d
		for k in 1:l
			max_gt = Base.setindex(max_gt, gt[k, d], k, l)
		end
	end
	return max_gt
end



rowsum(m, l) = l == 0 ? 0 : sum(m[k, l] for k in 1:l)

φ₀(GTpat) = sum([rowsum(GTpat, i) for i in 1:(length(weight(GTpat))-1)])

φ(GTpat) = φ₀(GTpat) - φ₀(max_gt_pat(GTpat))


function projectors(m::Int, n::Int, k::Int)
	"""
		Given (m, n, k), computes the nonzero CG coefficients then express the
		basis state |M⟩ (in the tableau basis, lhs of (T32)) as function of fock
		projectors |s⟩⟨s| 
	"""
	@debug "Compute GC m = $m, n = $n, k = $k"
	#TODO: Should accept other irreps.
	patterns = CGC(τ(n, m), conj(τ(n, m)), λ(k, m))
	all_gt = collect(basis(τ(n, m)))
	all_gt_dual = collect(basis(τ(n, m) |> conj))

	proj = Dict(i => Dict() for i in 1:(d_λ(m, k)|>Int))
	fock_projectors = Dict(i => Vector{FockMEl}() for i in 1:(d_λ(m, k)|>Int))

	for (idx, c) in patterns.data
		proj[idx[3]][[idx[1], idx[2]]] = c
	end

	k = 1
	for i in collect(eachindex(proj))
		p = proj[i]
		k += 1
		fock_projectors[i] = projector_to_fock_state(all_gt, all_gt_dual, p)
	end

	return fock_projectors
end

function projector_to_fock_state(all_gt, all_gt_dual, projector)

	state = Vector{FockMEl}()

	for (state_list, c) in projector
		M1 = all_gt[state_list[1]] |> gt_to_fock
		M2 = all_gt_dual[state_list[2]] |> gt_dual_to_fock
		c *= (-1)^(φ(all_gt_dual[state_list[2]]))

		push!(state, FockMEl(FockState(M1), FockState(M2), c))
	end
	return state
end

function gt_to_fock(gt)
	"""
	Transforms a GT pattern [gt], returns its associated Fock state as per eq.
	(T26)
	"""
	N = gt.N
	first_row = []
	for i in 1:N
		push!(first_row, gt.data[i+(i-1)*N-sum([j for j in 0:i-1])])
	end
	row = first_row |> reverse

	fock_state = [row[1]; [row[i] - row[i-1] for i in 2:length(row)]]

	return fock_state
end

function gt_dual_to_fock(gt)
	"""
	Transforms a dual GT pattern [gt], returns its associated Fock state as per eq.
	(T26)
	"""
	N = gt.N
	n = gt.data[1]

	last_row = []
	for i in 1:N
		push!(last_row, gt.data[(i)*N-sum([j for j in 0:i-1])])
	end

	row = last_row
	fock_state = [row[2]; [row[i] - row[i-1] for i in 3:(length(row))]]
	fock_state = [fock_state; n - sum(fock_state)] |> reverse
	return fock_state
end

function pprint_projector(projectors)
	"""
	Pretty prints a projector in the (linearised) Fock basis
	"""
	for (i, projector) in projectors
		#TODO This could be latexified
		str = ("|M$i\rangle = ")
		for mel in projector
			c = mel.c
			str *= ("$(c > 0 ? "+$c" : c)|$(mel.bra.state...),$(mel.ket.state...)> ")
		end
		print(str)
		println()
	end
end



function projector_to_matrix(projector::Vector{FockMEl}, m::Int, n::Int)
	fock_indices = linearised_fock_basis_indices(m, n)


	entries = Dict{Tuple{Int64, Int64}, Float32}()

	@inbounds for i in 1:length(projector)
		@inbounds for j in i:length(projector)
			mel_1 = projector[i]
			mel_2 = projector[j]

			id1 = fock_indices[[mel_1.ket.state; mel_1.bra.state]]
			id2 = fock_indices[[mel_2.ket.state; mel_2.bra.state]]
			if id2 < id1
				tmp = id1
				id1 = id2
				id2 = tmp
			end

			c = mel_1.c * mel_2.c

			entries[(id1, id2)] = c
		end
	end

	return entries
end


function partial_projector_addition(chunk, projectors::Dict{Int64, Vector{FockMEl}}, m::Int, n::Int)
	Π = Dict{Tuple{Int64, Int64}, Float32}()

	for i in chunk
		p = projector_to_matrix(projectors[i], m, n)
		mergewith!(+, Π, p)
	end

	return Π
end

function projectors_to_matrix(projectors::Dict{Int64, Vector{FockMEl}}, m::Int, n::Int)
	"""
		Returns the matrix form of the [projector] (which is a dict of whose key
		is k and value the projector to |M⟩ associated to λ_k in the Fock basis)
	"""
	"""Best way I found to perform SparseArrays addition is with a dict and
	mergewith (seems better than built-in sparse array addition)"""
	Π = Dict{Tuple{Int64, Int64}, Float64}()

	keys = collect(eachindex(projectors))

	if length(keys) == 0
		return sparse([], [], [])
	end

	if length(keys) < Threads.nthreads()
		for i in keys
			p = projectors[i]
			mergewith!(+, Π, projector_to_matrix(p, m, n))
		end
	else
		chunks = Iterators.partition(keys, length(keys) ÷ Threads.nthreads())

		tasks = map(chunks) do chunk
			Threads.@spawn partial_projector_addition(chunk, projectors, m, n)
		end

		for p in fetch.(tasks)
			mergewith!(+, Π, p)
		end
	end


	"""
	Transform The hashmap version of the matrix 
		(i, j) => Mij
	to a SparseArrays
	"""

	rows = Vector{Int}()
	cols = Vector{Int}()
	values = Vector{Float32}()
	for (k, v) in Π
		push!(rows, k[1])
		push!(cols, k[2])
		push!(values, v)
	end

	return sparse(rows, cols, values)
end

"""
	projection_matrices(m::Int, n::Int; load_from_memory = true, save_in_memory = true, infer_last::Bool = true)

Compute projection matrices for `m` and all `i ≤ n` and `0 ≤ k ≤ i`.

Serialize computed matrix.

Deserialize data on disk when found.

!!! note 
	The implementation supports multithreading with `julia -t <num_threads>`.

!!! warning 
	As of today, the implementation only supports the computation of the
	projectors involving the totally symmetric irrep of ``SU(m)`` and its dual.
"""
function projection_matrices(m::Int, n::Int; load_from_memory = true, save_in_memory = true, infer_last::Bool = true)
	Π = Dict()

	for i in 0:n
		path = "./projectors/P-m$(m)-n$(i)"
		if load_from_memory && isfile(path)
			try
				merge!(Π, deserialize(path))
				@debug "m=$(m) n=$(i) loaded from memory "
			catch e
				@error e
			end
		else
			P = Dict()
			for k in 0:i-1
				P[[m, i, k]] = projectors_to_matrix(projectors(m, i, k), m, i)
			end

			if infer_last
				@debug "Infer Last one"
				PS = sum([P[[m, i, k]] for k in 0:i-1]; init = spzeros(binomial(i + m - 1, i)^2, binomial(i + m - 1, i)^2))
				Id = sparse(1.0I, binomial(i + m - 1, i)^2, binomial(i + m - 1, i)^2)
				P[[m, i, i]] = Id - PS
			else
				P[[m, i, i]] = projectors_to_matrix(projectors(m, i, i), m, i)
			end

			@debug "m=$(m) n=n$(i) computed projectors "

			if save_in_memory
				if !isdir("./projectors/")
					mkdir("./projectors")
				end
				try
					serialize(path, P)
					@debug "m$(m)-n$(i) saved at" path
				catch e
					@error e
				end
			end


			merge!(Π, P)
		end
	end

	return Π
end

"""
	sparse_sym_mvp(A::SparseMatrixCSC{Float64, Int64}, x::Vector{<:Number})

Implementation of the sparse-matrix vector product ``Ax`` when ``A`` is
symmetric and given via its upper-triangular form.
"""
function sparse_sym_mvp(A::SparseMatrixCSC{Float64, Int64}, x::Vector{<:Number})
	return A * x + transpose(transpose(x) * A) - Diagonal(A) * x
end
