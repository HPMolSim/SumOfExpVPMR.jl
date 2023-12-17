@testset "VPMR, Gaussian" begin
    f = x -> exp(-x^2)

    sw, σ = VPMR_cal(f, 4.0, 40, 100, 8)
    @test max_error(f, sw) < 1e-7

    sw, σ = VPMR_cal(f, 4.0, 40, 100, 12)
    @test max_error(f, sw) < 1e-10

    sw, σ = VPMR_cal(f, 6.0, 60, 200, 16)
    @test max_error(f, sw) < 1e-12
end

@testset "VPMR, error function" begin
    f = x -> erfc(x)

    sw, σ = VPMR_cal(f, 4.0, 40, 100, 8)
    @test max_error(f, sw) < 1e-7

    sw, σ = VPMR_cal(f, 4.0, 40, 100, 12)
    @test max_error(f, sw) < 1e-10

    sw, σ = VPMR_cal(f, 6.0, 60, 200, 16)
    @test max_error(f, sw) < 1e-12
end