function MR_cal(s::S, w::W, p::Int; T1::DataType = ComplexF64, T2::DataType = Float64, digit::Int = 1024, weighted_balanced_truncation = weighted_balanced_truncation) where {S, W}

    @assert iszero(digit % 256)

    if weighted_balanced_truncation
        s, w = Complex{BigFloat}.(s), Complex{BigFloat}.(w)
    end

    s, w, σ = setprecision(digit) do 
        MR(s, w, p) 
    end

    return T1.(s), T1.(w), T2(σ)
end

## WBT MR
function MR(s::Vector{Complex{BigFloat}}, w::Vector{Complex{BigFloat}}, p::Int)

    n = length(s)
    @assert p ≤ n
    A = diagm(- s)
    B = sqrt.(w)
    C = transpose(B)
    P = (B * B') ./ ([s[i] + conj(s[j]) for i in 1:length(s), j in 1:length(s)])
    
    S = cholesky(P).L

    STL = S' * conj(S)
    F = svd(STL)

    H = diagm(inv.(sqrt.(F.S)))
    Tt = S * F.U * H
    St = conj(S) * F.V * H
    At = St' * A * Tt
    Bt = St' * B
    Ct = C * Tt

    σ = 2.0 * sum([F.S[i] for i in p+1:n])

    Ad = At[1:p, 1:p]
    Bd = Bt[1:p]

    Cd = transpose(Ct[1:p])
    
    Ad_eigen = eigen(Ad)
    s = - Ad_eigen.values

    B_mr = inv(Ad_eigen.vectors) * Bd
    C_mr = Cd * Ad_eigen.vectors

    w = [C_mr[i] * B_mr[i] for i in 1:p]

    return s, w, σ
end

# Classical MR
function MR(s::Vector{BigFloat}, w::Vector{BigFloat}, p::Int)

    n = length(s)
    @assert p ≤ n

    A = diagm(- s)
    B = sqrt.(abs.(w))
    C = sign.(w) .* B

    P = lyapc(A, B * B')
    Q = lyapc(A, C * C')

    S = cholesky(P).L
    L = cholesky(Q).L

    STL = S' * L
    F = svd(STL)

    Tt = S * F.U * diagm(inv.(sqrt.(F.S)))

    At = inv(Tt) * A * Tt
    Bt = inv(Tt) * B
    Ct = C' * Tt

    σ = 2.0 * sum([F.S[i] for i in p+1:n])

    Ad = At[1:p, 1:p]
    Bd = Bt[1:p]
    Cd = (Ct[1:p])'

    Ad_eigen = eigen(Ad)
    s = - Ad_eigen.values

    B_mr = inv(Ad_eigen.vectors) * Bd
    C_mr = Cd * Ad_eigen.vectors

    w = [C_mr[i] * B_mr[i] for i in 1:p]

    return s, w, σ
end

function wbt(s::Vector{BigFloat}, w::Vector{BigFloat}, p::Int)

    n = length(s)
    @assert p ≤ n

    A = diagm(- s)
    B = sqrt.(w)
    C = transpose(B)
    #@info size(C), size(B)
    #P = B * B'./(s+s')
    P = (B * B') ./ ([s[i] + conj(s[j]) for i in 1:length(s), j in 1:length(s)])

    #Q = lyapc(A', C' * C)
    #Q = (C' * C) ./ ([conj(s[i]) + s[j] for i in 1:length(s), j in 1:length(s)])
    Q = conj.(P)
    #tmp = C'*C
    #@info P[1], Q[1], tmp[1]
    
    S = cholesky(P).L
#    L = cholesky(Q).L

    STL = S' * conj(S)
    F = svd(STL)

    H = diagm(inv.(sqrt.(F.S)))
    Tt = S * F.U * H
    St = conj(S) * F.V * H
 #   At = inv(Tt) * A * Tt
    At = St' * A * Tt
#    Bt = inv(Tt) * B
    Bt = St' * B
    Ct = C * Tt

    σ = 2.0 * sum([F.S[i] for i in p+1:n])

    Ad = At[1:p, 1:p]
    Bd = Bt[1:p]

    Cd = transpose(Ct[1:p])
    #@info size(Cd), size(Bd)
    
    Ad_eigen = eigen(Ad)
    s = - Ad_eigen.values

    B_mr = inv(Ad_eigen.vectors) * Bd
    C_mr = Cd * Ad_eigen.vectors

    w = [C_mr[i] * B_mr[i] for i in 1:p]

    return s, w, σ
end

function wbt_simplified(s::Vector{T}, w::Vector{T}, b:: T) where{T}
# works for positive weights and nodes
    n = length(s)

    A = diagm(- s)
    B = sqrt.(w)
    C = B'
    #@info size(C), size(B)
    #P = B * B'./(s+s')
    #P = (B * B') .* ([(1-exp.(-(s[i]+s[j])*b))./(s[i] + s[j]) for i in 1:length(s), j in 1:length(s)])
    P = (B * B') ./ ([s[i] + s[j] for i in 1:length(s), j in 1:length(s)])


    F = svd(P)

    S = F.S
    errbnd = T(1.0e-12)
    
    p = sum(S .> errbnd * S[1])
    @info "p: $p"
    
    H = diagm(inv.(sqrt.(F.S)))
    U = F.U
    At = U' * A * U
    Bt = U' * B
    Ct = C * U

    error = 2.0 * sum([S[i] for i in p+1:n])
    @info "error: $error"

    Ad = At[1:p, 1:p]
    Bd = Bt[1:p]

    Cd = transpose(Ct[1:p])
    #@info size(Cd), size(Bd)
    
    Ad_eigen = eigen(Ad)
    s = - Ad_eigen.values

    B_mr = inv(Ad_eigen.vectors) * Bd
    C_mr = Cd * Ad_eigen.vectors

    w = [C_mr[i] * B_mr[i] for i in 1:p]

    return s, w, p, error
end