"""
    linear_invariant(ρ::DensityMatrixBlock)

Computes the linear optical Lie-algebraic invariant ``I_2`` (Eq.(14) of [this
reference](https://arxiv.org/pdf/2409.12223)).
"""
function linear_invariant(ρ::DensityMatrixBlock)

	m = ρ.m
	n = ρ.n
	return [
		[tr(Oxjk(m, n, i, j) * ρ)^2 for i in 1:m-1 for j in (i+1):m] |> sum,
		[tr(Oyjk(m, n, i, j) * ρ)^2 for i in 1:m-1 for j in (i+1):m] |> sum,
		[tr(Ozj(m, n, i) * ρ)^2 for i in 1:m] |> sum,
	] |> sum
end


function covariance_matrix(ρ::DensityMatrixBlock)
	m = ρ.m
	n = ρ.n

	O = [
		[Oxjk(m, n, i, j) for i in 1:m-1 for j in (i+1):m]...,
		[Oyjk(m, n, i, j) for i in 1:m-1 for j in (i+1):m]...,
		[Ozj(m, n, i) for i in 1:m]...,
	]

	return [tr(oi * ρ) * tr(oj * ρ) - tr((oi * oj + oj * oi) * ρ * 1 / 2) for oj in O, oi in O]
end


"""
    covariance_invariant(ρ::DensityMatrixBlock)

Computes the covariance invariant (Eq.(9) of [this
reference](https://arxiv.org/pdf/2409.12223)).
"""
function covariance_invariant(ρ::DensityMatrixBlock)
	return eigvals(covariance_matrix(ρ))
end

"""
	tangent_invariant(ρ::DensityMatrixBlock)

Computes the tangent invariant (Eq.(20) of [this
reference](https://arxiv.org/pdf/2409.12223)).
"""
function tangent_invariant(ρ::DensityMatrixBlock)
	m = ρ.m
	n = ρ.n

	O = [
		[Ozj(m, n, i) for i in 1:m]...,
		[Oxjk(m, n, i, j) for i in 1:m-1 for j in (i+1):m]...,
		[Oyjk(m, n, i, j) for i in 1:m-1 for j in (i+1):m]...,
	]

	return eigvals(sum([tr(oi * ρ) * oi for oi in O]))
end
