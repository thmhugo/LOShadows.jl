@testset "testing matrix direct sum" begin
	A = [
		1 2;
		3 4
	]
	B = [
		5 6;
		7 8
	]
	C = [
		1 2 3;
		4 5 6;
		7 8 9
	]

	M = [A, B, C]

	ApB = [
		1 2 0 0;
		3 4 0 0;
		0 0 5 6;
		0 0 7 8
	]
	ApC = [
		1 2 0 0 0;
		3 4 0 0 0;
		0 0 1 2 3;
		0 0 4 5 6;
		0 0 7 8 9
	]

	@test direct_sum(A, B) == ApB
    @test direct_sum(A, B) != direct_sum(B, A)
	@test direct_sum(A, C) == ApC
	@test direct_sum(M) == direct_sum(direct_sum(A, B), C)
end
