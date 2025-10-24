function χ(η, N, ρ)   
    O = exp(1im*sum([e * n for (e, n) in zip(η, N)]))

    return tr(O*ρ)
end

function N(bins, m, n)
    n_op = [Ozj(m, n, k) for k in 1:m]

    NK = [
        sum([n_op[idx] for idx in bin]) for bin in bins
    ]
    return NK
end

function omega(n, K)
    Ω = collect(0:n)
    return Iterators.product([Ω for _ in 1:K]...) 
end

"""
    binned_probability(k::Vector{Int}, bins::Vector{Vector{Int}}, ρ::DensityMatrixBlock)

Compute the probability of observing the occupation `k` ``=[k_1, \\cdots,
k_{\\#\\text{bins}}]`` when the output bins are given by `bins`.
"""
function binned_probability(k::Vector{Int}, bins::Vector{Vector{Int}}, ρ::DensityMatrixBlock)
    K = length(k)
    NK = N(bins, ρ.m, ρ.n) 
    s = 0
    for l in omega(ρ.n, K)
        nul = (2*π)/(ρ.n+1) * collect(l)
        s += χ(nul, NK, ρ) * exp(-1im * dot(nul, k)) 
    end
    return s / (ρ.n+1)^K |> real
end