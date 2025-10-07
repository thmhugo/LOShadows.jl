using LOShadows
using LinearAlgebra
using ArgCheck
using Test

include("../src/MatrixOperations.jl")
include("../src/DensityMatrix.jl")

@testset "LOShadows.jl" begin
    include("./densitymatrix.jl")
    include("./matrixoperations.jl")
end