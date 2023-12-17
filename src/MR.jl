function MR_cal(s::Vector{BigFloat}, w::Vector{BigFloat}, p::Int; T1::DataType = ComplexF64, T2::DataType = Float64, digit::Int = 1024)

    @assert iszero(digit % 256)

    s, w, σ = setprecision(digit) do 
        MR(s, w, p) 
    end

    return T1.(s), T1.(w), T2(σ)
end

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