using SumOfExpVPMR
using Test
using SpecialFunctions

@testset "SumOfExpVPRM.jl" begin
    include("VP.jl")
    include("VPMR.jl")
    include("fast_exp_transform.jl")
end
