export Eij, Ozj, Oxjk, Oyjk


# @memoize all_fock_states(m::Int, n::Int) = multiexponents(m, n)

function fock_idx(m, n)
    return Dict(s => i for (i, s) in all_fock_states(m, n) |> enumerate)
end


"""
    Eij(m::Int, n::Int, i::Int, j::Int)

Finite dimensional matrix representation of the operator ``E_{i,j} =
\\hat{a}_i^\\dagger \\hat{a}_j`` acting over ``\\mathcal{H}_m^n``.
"""
function Eij(m::Int, n::Int, i::Int, j::Int)
    d = binomial(n + m - 1, n)
    A = zeros(d, d)
    idx = fock_idx(m, n)
    for s in all_fock_states(m, n)
        # Apply A_i^dagger A_j to s
        t = copy(s)
        t[i] += 1
        t[j] -= 1
        if all(>=(0), t) && all(<=(n), t)
            A[idx[s], idx[t]] = sqrt(t[i]) * sqrt(s[j])
        end
    end
    return A
end


"""
    Ozj(m::Int, n::Int, j::Int)

Finite dimensional matrix representation of the number operator ``O^z_{j} =
\\hat{n}_{j}`` acting over ``\\mathcal{H}_m^n``.
"""
function Ozj(m::Int, n::Int, j::Int)
    return Eij(m, n, j, j)
end

"""
    Oxjk(m::Int, n::Int, j::Int, k::Int)

Finite dimensional matrix representation of the operator ``O^x_{j,k} =
\\frac{1}{\\sqrt{2}}(E_{jk} + E_{k,j})`` acting over ``\\mathcal{H}_m^n``.
"""
function Oxjk(m::Int, n::Int, j::Int, k::Int)
    return (Eij(m, n, j, k) + Eij(m, n, k, j)) * 1 / sqrt(2)
end

"""
    Oyjk(m::Int, n::Int, j::Int, k::Int)

Finite dimensional matrix representation of the operator ``O^y_{j,k} =
\\frac{\\imath}{\\sqrt{2}}(E_{jk} - E_{k,j})`` acting over ``\\mathcal{H}_m^n``.
"""
function Oyjk(m::Int, n::Int, j::Int, k::Int)
    return (Eij(m, n, j, k) - Eij(m, n, k, j)) * im / sqrt(2)
end