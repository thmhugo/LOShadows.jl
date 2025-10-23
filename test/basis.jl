
@testset "linearised basis indexing" begin
	m = 2
	n = 2

	@test SUNIrrepProjectors.all_fock_states(m, n) |> collect == Any[
		[2, 0], [1, 1], [0, 2],
	]

	@test linearised_fock_basis_indices(m, n) == Dict(
		[2, 0, 2, 0] => 1,
		[2, 0, 1, 1] => 2,
		[2, 0, 0, 2] => 3,
		[1, 1, 2, 0] => 4,
		[1, 1, 1, 1] => 5,
		[1, 1, 0, 2] => 6,
		[0, 2, 2, 0] => 7,
		[0, 2, 1, 1] => 8,
		[0, 2, 0, 2] => 9,
	)
end