module SumOfExpVPMR
__precompile__()

using LinearAlgebra, SpecialFunctions, GaussQuadrature, SpecialFunctions, MatrixEquations, GenericLinearAlgebra, GenericSchur
using DelimitedFiles, Polynomials

export GaussParameter, Gauss_int
export VP_cal, MR_cal, VPMR_cal, psog_cal
export soe, soe_error, max_error, SoePara
export sog, sog_error, sog_max_error, SOGPara
export prony
export FET1d

include("Gaussian_integral.jl")
include("VP.jl")
include("MR.jl")
include("soe.jl")
include("sog.jl")
include("prony.jl")
include("gsoe_initialguess.jl")
include("psog_initialguess.jl")
include("fast_exp_transform.jl")

function VPMR_cal(f::Function,
    nc::T,
    n::Int,
    N::Int,
    p::Int;
    region::Tuple{TR1, TR2} = (0.0, π),
    T1::DataType = ComplexF64,
    T2::DataType = Float64,
    digit::Int = 512,
    print_info::Bool=false,
    weighted_balanced_truncation::Bool=true) where{T, TR1, TR2}

    @assert iszero(digit % 256)

    s, w = setprecision(digit) do
        VP(f, nc, n, N, region)
    end

    if print_info
        x = big.([0.0:0.001: 10.0...])
        error_VP = max_error(f, s, w, x = x)
        @info "VP error: $error_VP"
    end

    smr, wmr, σ = setprecision(digit) do 
        MR_cal(s, w, p; weighted_balanced_truncation = weighted_balanced_truncation)
	end

    perm = sortperm(wmr,by=abs)
    smr = smr[perm]
    wmr = wmr[perm]
    if print_info
        x = [0.0:0.001: 10.0...]
        error_MR = max_error(f, smr, wmr, x = x)
        @info "MR error: $error_MR"
    end

    Ts, Tw, Tσ = T1.(smr), T1.(wmr), T2(σ)

    x = [0.0:0.001: 10.0...]
    error_TMR = max_error(f, Ts, Tw, x = x)
    @info "Truncated MR error: $error_TMR"

    return SoePara(Ts, Tw), error_TMR
end


function psog_cal(f::Function; T1::DataType = ComplexF64, T2::DataType = Float64, use_sp::Bool = true, digits::Int = 128, initial_tol = 1e-15, print_info::Bool=false) 

    # use_sp for use super precision or not
    T = use_sp ? BigFloat : Float64
    use_sp && setprecision(digits)

    beta0, a0, b0, tol, x = setup_params(T(initial_tol))
    s, w = psog_initialguess(beta0, a0, b0, tol)

    @info "length(s): $(length(s)), maximum(s): $(maximum(s)), minimum(s): $(minimum(s))"
    idx = findall(x -> log(T(2)) <= x <= 10*log(T(10))/a0^2, s)  
    @info "length(idx): $(length(idx))"
    all_idx = collect(1:length(s))  # all valid indices
    cidx = setdiff(all_idx, idx) 

    s1 = s[cidx]
    w1 = w[cidx]

    s0 = s[idx]
    w0 = w[idx]
    
    if print_info
       error0 = sog_max_error(f, s, w, x)
       @info "initial error: $error0"
    end

    smr, wmr, p, σ = wbt_simplified(s, w, b0^2)
    @info "length(smr): $(length(smr))"

    sn = vcat(s1,smr)
    wn = vcat(w1,wmr)
    
    perm = sortperm(wn,by=abs)
    sn = sn[perm]
    wn = wn[perm]
    @info "length(sn): $(length(sn))"
    
    if print_info
        error_MR = sog_max_error(f, sn, wn, x)
        @info "MR error: $error_MR"
    end

    Ts, Tw, Tσ = T2.(smr), T2.(wmr), T2(σ)

    x = T2.(x)
    error_TMR = sog_max_error(f, Ts, Tw, x)
    @info "Truncated MR error: $error_TMR"

    return SOGPara(Ts, Tw), error_TMR
end


function setup_params(initial_tol::T) where T
    beta0 = one(T)
    a0 = T(2)^-10
    b0 = T(10)

    m = 1000
    estart = log10(a0)
    eend = log10(b0)
    rexp = range(estart, eend; length=m)
    x = T.(10.0.^rexp)

    return beta0, a0, b0, initial_tol, x
end

end
