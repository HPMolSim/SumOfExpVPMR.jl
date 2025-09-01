using LinearAlgebra.LAPACK

# qr decomposition directly using LAPACK
# actually the function LinearAlgebra.qr is also using LAPACK as backend
function lapack_qr(A::Matrix{T}) where{T}
    Q = copy(A)
    m, n = size(Q)
    k = min(m, n)
    tau = zeros(T, k)
    R = zeros(T, m, n)
    LAPACK.geqrf!(Q, tau)
    for i in 1:m
        for j in i:min(n, k)
            R[i, j] = Q[i, j]
        end
    end
    LAPACK.orgqr!(Q, tau, k)
    return Q, R
  end

# Custom hankel constructor for arbitrary precision vectors
function hankel(c::Vector{T}, r::Vector{T}) where{T}
    nrow = length(c)
    ncol = length(r)
    H = Matrix{T}(undef, nrow, ncol)
    for i in 1:nrow
        for j in 1:ncol
            k = i + j - 1
            if k <= length(c)
                H[i, j] = c[k]
            else
                H[i, j] = r[k - length(c) + 1]
            end
        end
    end
    return H
end

function myls2(A::Matrix{T}, b::Vector{T}, eps::T) where{T}

    m, _ = size(A)
    F = qr(A)
    Q = Matrix(F.Q)
    R = Matrix(F.R)

    # the manul qr directly using LAPACK
    # Q, R = lapack_qr(A)

    @info "QR element-wise relative error: $(maximum(abs.((Q * R .- A) ./ A)))"

    s = diag(R)
    r = count(abs.(s) .> eps)
    Qr = Q[:, 1:r]
    Rr = R[1:r, 1:r]

    # Slice b from r+1 to m + r, just like in MATLAB
    b1 = b[r+1 : m+r]  # this is valid if b has length ≥ m + r
    b2 = transpose(Qr) * b1
    x = Rr \ b2

    res = norm(Rr*x-b2) / norm(b2)
    @info "myls2 residual: $res"

    return x
end

# Least-squares using SVD
function myls(A::Matrix{T}, b::Vector{T}, eps::T) where{T}
    _, n = size(A)
    F = svd(A)
    U = F.U
    S = F.S
    V = F.V
    r = count(S ./ S[1] .> eps)

    x = zeros(T, n)
    for i in 1:r
        x += ((U[:, i]' * b) / S[i]) * V[:, i]
    end

    res = norm(A*x-b) / norm(b)
    @info "myls residual: $res"

    return x
end


# Main prony method
function prony(xs::Vector{T}, ws::Vector{T}, errbnd::T) where{T}
    M = length(xs)
    h = Vector{T}(undef, 2*M)
    for j in 1:2*M
        h[j] = sum(xs .^ (j-1) .* ws)
    end

    C = h[1:M]
    R = h[M:2*M-1]
    H = hankel(C, R)


    b = -h
    q = myls2(H, b, errbnd)

    # matlab results
    # q = [2.89538614224549e-09,1.17373935854545e-06,7.83176678450703e-05,0.00201184951196368,0.0259380019882274,0.189145051104628,0.824363460779003,2.17765069246786,3.39749119607388,2.86372834937139]

    r = length(q)
    A = Matrix{T}(undef, 2*M, r)

    Coef = [T(1); reverse(q)]

    # solving the roots of the polynomial using Julia's standard package Polynomials.jl
    xsnew = roots(Polynomial(reverse(Coef)))

    # matlab results
    # xsnew = [-0.7693582923088651, -0.5919121795202147, -0.4553916881934084, -0.35029070128222, -0.2680487783257383, -0.1973409274581753, -0.1304655035382046, -0.07114066567626039, -0.02674408299599483, -0.003035530072314068]

    for j in 1:2*M
        A[j, :] = xsnew .^ (j-1)
    end

    wsnew = myls(A, h, errbnd)

    @assert all(x -> x > 0, wsnew) "Array contains non-positive values"

    return wsnew, xsnew
end