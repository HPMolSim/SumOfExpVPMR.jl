module SumOfExpVPMR

using LinearAlgebra, SpecialFunctions, GaussQuadrature, SpecialFunctions, MatrixEquations, GenericLinearAlgebra, GenericSchur

export GaussParameter, Gauss_int
export VP_cal, MR_cal, VPMR_cal
export soe, soe_error, max_error, SoePara
export FET1d

include("Gaussian_integral.jl")
include("VP.jl")
include("MR.jl")
include("soe.jl")
include("fast_exp_transform.jl")

function VPMR_cal(f::Function, nc::T, n::Int, N::Int, p::Int; region::Tuple{TR1, TR2} = (0.0, π), T1::DataType = ComplexF64, T2::DataType = Float64, digit::Int = 1024, print_info::Bool=false) where{T, TR1, TR2}

    @assert iszero(digit % 256)

    s, w = setprecision(digit) do
        VP(f, nc, n, N, region)
    end

    if print_info
        x = [- nc * log((1.0 + cos(region[1])) / 2.0):0.001: 10.0...]
        error_VP = max_error(f, s, w, x = x)
        @info "VP error: $error_VP"
    end

    smr, wmr, σ = setprecision(digit) do 
        MR(s, w, p)
    end

    if print_info
        x = [- nc * log((1.0 + cos(region[1])) / 2.0):0.001: 10.0...]
        error_MR = max_error(f, smr, wmr, x = x)
        @info "MR error: $error_MR"
    end

    Ts, Tw, Tσ = T1.(smr), T1.(wmr), T2(σ)

    if print_info
        x = [- nc * log((1.0 + cos(region[1])) / 2.0):0.001: 10.0...]
        error_TMR = max_error(f, Ts, Tw, x = x)
        @info "Truncated MR error: $error_TMR"
    end

    return SoePara(Ts, Tw), Tσ
end

end
