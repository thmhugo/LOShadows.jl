@testset "testing DensityMatrix constructor" begin
    m, n = 3, 2
    d = binomial(m+n, n)

    wrong_dim_rho = Matrix{ComplexF32}(I, 3, 3) * 1 / d

    good_rho = Matrix{ComplexF32}(I, d, d) * 1 / d

    nonpsd_rho = -1 * Matrix{ComplexF32}(I, d, d) * 1 / d
    # nonnormalised_rho = Matrix{ComplexF32}(I, d, d)
    nonselfadjoint_rho = Matrix{ComplexF32}(I, d, d) * 1 / d
    nonselfadjoint_rho[1,2] = (2.0 + 0.0im)

    @test_throws ArgumentError DensityMatrix(-1, n, good_rho)
    @test_throws ArgumentError DensityMatrix(m, -1, good_rho)
    @test_throws ArgumentError DensityMatrix(m, n, wrong_dim_rho)
    @test_throws ArgumentError DensityMatrix(m, n, nonpsd_rho)
    # @test_throws ArgumentError DensityMatrix(m, n, nonnormalised_rho)
    @test_throws ArgumentError DensityMatrix(m, n, nonselfadjoint_rho)
end

@testset "testing DensityMatrixBlock constructor" begin
    m, n = 2, 2
    d = binomial(n+m-1, n)
    blocks = [
        0 * Matrix{ComplexF32}(I, 1, 1),
        Matrix{ComplexF32}(I, m, m) / m,
        Matrix{ComplexF32}(I, d, d) / d
    ]

    # print(direct_sum(blocks))
    dm = DensityMatrix(m, n, blocks)

    @test dm.ρ ≈ direct_sum(blocks) / tr(direct_sum(blocks))
end