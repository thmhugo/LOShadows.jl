struct Snapshot
    m::Int
    n::Int
    s::Vector{Int}
    U::Matrix{ComplexF32}

    Snapshot(m, n, s, U) = begin
        @argcheck length(s) == m
        @argcheck sum(s) == n
        @argcheck all(s .>= 0)
        @argcheck isapprox(U*U', Matrix{ComplexF32}(I, m, m))

        return new(m, n, s, U)
    end

    Snapshot(s, U) = new(length(s), sum(s), s, U)
end