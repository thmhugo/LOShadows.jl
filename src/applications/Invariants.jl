using ..FockOperators
using LinearAlgebra
include("../SUN.jl")

function linear_invariant(ρ, m, n)
    """
    Computes the LO invariant of EQ.4 of https://arxiv.org/pdf/2409.12223
    """
    return [
        [tr(Oxjk(m, n, i, j) * ρ)^2 for i in 1:m-1 for j in (i+1):m] |> sum,
        [tr(Oyjk(m, n, i, j) * ρ)^2 for i in 1:m-1 for j in (i+1):m] |> sum,
        [tr(Ozj(m, n, i) * ρ)^2 for i in 1:m] |> sum
    ] |> sum
end


function covariance_matrix(ρ, m, n)
    O = [
        [Oxjk(m, n, i, j) for i in 1:m-1 for j in (i+1):m]...,
        [Oyjk(m, n, i, j) for i in 1:m-1 for j in (i+1):m]...,
        [Ozj(m, n, i) for i in 1:m]...
    ]

    return [tr(oi * ρ) * tr(oj * ρ) - tr((oi * oj + oj * oi) * ρ * 1 / 2) for oj in O, oi in O]
end

function nonlinear_invariant(ρ, m, n)
    # println(covariance_matrix(ρ, m, n))
    return eigvals(covariance_matrix(ρ, m, n))
end

function projection_invariant(ρ, m, n)
    O = [
        [Ozj(m, n, i) for i in 1:m]...,    
        [Oxjk(m, n, i, j) for i in 1:m-1 for j in (i+1):m]...,
        [Oyjk(m, n, i, j) for i in 1:m-1 for j in (i+1):m]...
    ]

    # @info size(O)
    # @info [tr(o*ρ) for o in O]

    M = sum([tr(oi*ρ) * oi for oi in O])|> zchop
    # @info (M)
    return eigvals(
        M
    )
end

function shadow_projection_invariant(dm, m, n)
    Π = inverse_channels(m, n)[end]
    O = [
        [Ozj(m, n, i) for i in 1:m]...,
        [Oxjk(m, n, i, j) for i in 1:m-1 for j in (i+1):m]...,
        [Oyjk(m, n, i, j) for i in 1:m-1 for j in (i+1):m]...
    ]
    # @info tr(apply_nth_channel(O[1], m, n, Π)*dm), tr(apply_nth_channel(O[2], m, n, Π)*dm), tr(apply_nth_channel(O[3], m, n, Π)*dm)
    # @info [tr(apply_nth_channel(dm, m, n, Π)*oi) for oi in O]
    # @info [tr(dm*oi) for oi in O]
    # M = sum([tr(apply_nth_channel(dm, m, n, Π)*oi)*oi for oi in O])
    M = sum([tr(dm*oi)*oi for oi in O])
    # @info M
    return eigvals(M)
end



function mom(X, K, B)
    @assert length(X) % K == 0

    groups = [X[(i*K+1):((i+1)*K)] for i in 0:length(X)÷K-1]
    means = mean.(groups)

    ## Compute median error via empirical bootstraps
    # bootstrap_medians = [
    #     [means[rand(1:length(means))] for _ in 1:length(means)] |> median
    #     for _ in 1:B
    # ]

    # bootstrap_medians_avg = mean(bootstrap_medians)

    # bootstrap_estimate = sum([
    #     (b - bootstrap_medians_avg)^2 for b in bootstrap_medians
    # ]) / (B - 1)

    # @info median(means), mean(X)
    return median(means)
end

# mom([1,2,3,4,5,6,7,8], 4, 1)

function mom_estimate(O, snapshots, K, B)
    """
    Returns a Median-of-Means estimator of O using the snapshots and parameters K and B
    """
    return mom([tr(O * s) |> real for s in snapshots], K, B)
end

@memoize function get_evolved_observables(m, n)
    Π = inverse_channels(m, n)[end]
    return [
        [apply_nth_channel(Oxjk(m, n, i, j), m, n, Π) for i in 1:m-1 for j in (i+1):m]...,
        [apply_nth_channel(Oyjk(m, n, i, j), m, n, Π) for i in 1:m-1 for j in (i+1):m]...,
        [apply_nth_channel(Ozj(m, n, i), m, n, Π) for i in 1:m]...,
    ]
end

function shadow_tom_linear_invariant(ρ, m, n)
    """
    Computes the LO invariant of EQ.4 of https://arxiv.org/pdf/2409.12223
    """
    Π = inverse_channels(m, n)[end]
    # return [
    #     [tr(apply_nth_channel(Oxjk(m, n, i, j), m, n, Π) * ρ)^2 for i in 1:m-1 for j in (i+1):m] |> sum,
    #     [tr(apply_nth_channel(Oyjk(m, n, i, j), m, n, Π) * ρ)^2 for i in 1:m-1 for j in (i+1):m] |> sum,
    #     [tr(apply_nth_channel(Ozj(m, n, i), m, n, Π) * ρ)^2 for i in 1:m] |> sum
    # ] 
    return [tr(o*ρ)^2 for o in get_evolved_observables(m, n)] |> sum
end

function shadow_linear_invariant(snapshots, m, n, K, B)
    """
    Computes the LO invariant of EQ.4 of https://arxiv.org/pdf/2409.12223
    from the classical shadow
    """
    O = get_evolved_observables(m, n)
    return [mean([tr(oi * s) |> real for s in snapshots]) for oi in O ] |> sum
    # return [mom_estimate(o, snapshots, K, B)^2 for o in O]  |> sum
end



function shadow_covariance_matrix(dm, m, n, K, B)
    Π = inverse_channels(m, n)[end]

    O = [
        [Oxjk(m, n, i, j) for i in 1:m-1 for j in (i+1):m]...,
        [Oyjk(m, n, i, j) for i in 1:m-1 for j in (i+1):m]...,
        [Ozj(m, n, i) for i in 1:m]...
    ]
    # 
    return [mom_estimate(apply_nth_channel(oi, m, n, Π), dm, K, B) * mom_estimate(apply_nth_channel(oj, m, n, Π), dm, K, B) - mom_estimate(apply_nth_channel((oi * oj + oj * oi) * 1 / 2, m, n, Π), dm, K, B) for oj in O, oi in O]
end

function shadow_nonlinear_invariant(dm, m, n, K, B)
    # println(covariance_matrix(ρ, m, n))

    # @info covariance_matrix = shadow_covariance_matrix(snapshots, m, n, K, B)
    return eigvals(shadow_covariance_matrix(dm, m, n, K, B))
end
